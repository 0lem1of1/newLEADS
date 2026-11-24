#!/usr/bin/env python3
"""Check LEADS's time-domain emitted signal (CalcEmittedSignal) against its own
validated spectrum and against an exact multi-particle reference.

Usage (build LEADS first, e.g. `make all` in the repo root, so files/*.o exist):

    python3 -I validation/emitted_signal_check.py [--workdir DIR] [LABEL=FILES_DIR ...]

Each LABEL=FILES_DIR is a LEADS build (a dir holding leads_*.o and *.h); the
default is current=<repo>/files. Case B runs only for the last build given.

Single electron, Model 0, a0 = 0.01, 10 cycles, gamma = 100, detector on axis
(theta = 180 deg, along the electron's motion).

Case A (NP = 1). With f(T) = R*E(T) on the observer-time axis T = t - n.r:
  * Parseval:  int |f|^2 dT  ==  int dP/dOmega dt_ret  ==  (1/pi) int S(w) dw
    (dP/dOmega = |NNU|^2/kappa^5 is RadiationModule's column, S is spectrum.dat).
  * Shape:     |int f e^{iwT} dT|^2 == S(w), since dT = kappa dt and spectrum.dat
    integrates NNU/kappa^2 e^{iwT} dt.  S is already checked against Thomson.

Case B (NP = 20, default bunch shape). In a plane wave every particle's path
is the reference path translated by a in both z and t, a = dc/(1+beta), where
c = z(t0) + beta*t0 is the ballistic intercept. Observer time T = t + z shifts
by 2a, so the exact summed signal is sum_np f_ref(T - 2a_np). Its bursts do not
overlap, so its energy is N * W_ref.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
R_DET = 1.0e8
GAMMA, A0, CYCLES, DTAU = 100.0, 0.01, 10.0, 5e-3
BETA = math.sqrt(1 - 1 / GAMMA ** 2)


def patch(text, params):
    for k, v in params.items():
        pat = re.compile(r"(^[ \t]*%s[ \t]*=[ \t]*)[^;]*;" % re.escape(k), re.M)
        text, n = pat.subn(lambda m, v=v: m.group(1) + v + ";", text)
        assert n == 1, (k, n)
    return text


def run(rundir, objdir, np_, bunch):
    wt = (1 + BETA) / (1 - BETA)
    T = 2 * math.pi * CYCLES
    p = {
        "RadiationReaction": "false", "PulseModel": "0", "Polarization": "0",
        "a0": repr(A0), "Tau_FWHM": "%g * (2*pi)" % CYCLES, "TimeStep": repr(DTAU),
        "NoBeamPart": str(np_), "GamaParticle": repr(GAMMA),
        "PropTheta": "180.0", "PropPhi": "0.0", "z0": "%.6f" % T,
        "BunchLength": "4" if bunch else "0", "BunchRadius": "1" if bunch else "0",
        "CalculateSpectrum": "true", "SignalEmitted": "true", "AngularDistribution": "false",
        "DetecThetaStart": "180.0", "DetecThetaEnd": "180.0",
        "DetecPhiStart": "0.0", "DetecPhiEnd": "0.0", "UseLogScale": "true",
        "SpectrumStart": "%.6f" % (math.log10(wt) - 0.45),
        "SpectrumEnd": "%.6f" % (math.log10(wt) + 0.15), "SpectrumDO": "0.004",
        "StoreTrajectories": "true", "DumpStateBinary": "true",
        "WriteSpatialFieldsProfile": "false",
    }
    subprocess.run(["rm", "-rf", str(rundir)], check=True)
    for d in ("Data", "bench", "files"):
        (rundir / d).mkdir(parents=True)
    subprocess.run("cp %s/*.o %s/*.h %s/files/" % (objdir, objdir, rundir), shell=True, check=True)
    (rundir / "input.cpp").write_text(patch((REPO / "input.cpp").read_text(), p))
    for cmd in ("g++ -std=c++17 -O3 -o leads ./files/*.o input.cpp -larmadillo -lgsl -lgslcblas -lm",
                "./leads 1"):
        r = subprocess.run(cmd, shell=True, cwd=rundir, capture_output=True, text=True)
        if r.returncode:
            sys.exit("%s failed in %s:\n%s" % (cmd, rundir, (r.stdout + r.stderr)[-1500:]))


def read_state(path):
    raw = path.read_bytes()
    NP, nT = np.frombuffer(raw, np.int64, 2, 8)
    d = np.frombuffer(raw, np.float64, offset=72)
    TR = d[:nT]
    posi = d[nT:nT + 3 * NP].reshape(3, NP)          # column-major (x.., y.., z..)
    return int(NP), TR, posi


def load(rundir):
    sig = np.loadtxt(rundir / "Data" / "emitted_signal.dat")
    spec = np.loadtxt(rundir / "Data" / "spectrum.dat")
    return sig, spec


def nudft(T, f, w):
    """|int f(T) e^{iwT} dT|^2 summed over components, trapezoid on the given T."""
    out = np.zeros(len(w))
    for c in range(f.shape[1]):
        out += np.abs(np.trapezoid(f[None, :, c] * np.exp(1j * w[:, None] * T[None, :]), T, axis=1)) ** 2
    return out


def nudft_complex(T, f, w):
    return np.stack([np.trapezoid(f[None, :, c] * np.exp(1j * w[:, None] * T[None, :]), T, axis=1)
                     for c in range(f.shape[1])], axis=1)


def peak(s, y, half=3):
    i = int(np.argmax(y))
    sl = slice(i - half, i + half + 1)
    x = (s[sl] - s[i]) / s[i]
    c = np.polyfit(x, np.log(y[sl]), 2)
    return s[i] * (1 - c[1] / (2 * c[0]))


def case_a(label, rundir):
    sig, spec = load(rundir)
    TR, tlab, E, dP = sig[:, 2], sig[:, 3], sig[:, 4:7], sig[:, 7]
    s, S = spec[:, 2], spec[:, 3]
    f = R_DET * E
    W_ang = np.trapezoid(dP, TR)
    W_sig = np.trapezoid((f ** 2).sum(1), tlab)
    W_spec = np.trapezoid(S, s) / math.pi
    Ssig = nudft(tlab, f, s)
    Rdist = R_DET + (tlab - TR)          # true distance (vT = TR + |R| - r_det)
    used = sig[:, 3] - sig[:, 2]         # what the original code divides by
    print("\n[A: NP=1, %s]" % label)
    print("  min |n.r| (what the pre-fix code divided by): %.3e; true distance %.3e" % (np.abs(used).min(), Rdist.mean()))
    print("  int|R E|^2 dT / int dP/dOmega dt : %.6g" % (W_sig / W_ang))
    print("  (1/pi) int S dw / int dP/dOmega dt : %.6g   (spectrum.dat vs power, code-independent of the bug)" % (W_spec / W_ang))
    print("  peak of |FT signal|^2 / peak of spectrum.dat : %.6f" % (peak(s, Ssig) / peak(s, S)))
    print("  shape: ||FT(signal)/max - S/max|| / ||S/max|| : %.3e" % (
        np.linalg.norm(Ssig / Ssig.max() - S / S.max()) / np.linalg.norm(S / S.max())))
    return dict(tlab=tlab, f=f, s=s, S=S, W=W_sig, W_ang=W_ang)


def intercepts(rundir):
    NP, TR, posi = read_state(rundir / "bench" / "state_init.bin")
    return posi[2] + BETA * TR[0]        # z moves as -beta*t: z = -beta*t + c


def case_b(rundir, ref, refdir):
    sig, spec = load(rundir)
    tlab, E = sig[:, 3], sig[:, 4:7]
    f = R_DET * E
    c = intercepts(rundir)
    c_ref = intercepts(refdir)[0]
    N = len(c)
    a = (c - c_ref) / (1 + BETA)
    W_leads = np.trapezoid((f ** 2).sum(1), tlab)
    # exact energy of sum_i f_ref(T - tau_i): sum_ij C(tau_i - tau_j), C = autocorrelation
    Tr, fr = ref["tlab"], ref["f"]
    fr2 = (fr ** 2).sum(1)
    on = Tr[fr2 > 1e-10 * fr2.max()]
    width = on.max() - on.min()
    tau = 2 * a
    grid = np.linspace(on.min() - width, on.max() + width, 400001)
    fg = np.stack([np.interp(grid, Tr, fr[:, k], left=0, right=0) for k in range(3)], 1)

    def C(d):
        sh = np.stack([np.interp(grid + d, Tr, fr[:, k], left=0, right=0) for k in range(3)], 1)
        return np.trapezoid((fg * sh).sum(1), grid)

    W_true = 0.0
    npairs = 0
    for i in range(N):
        for j in range(N):
            d = tau[j] - tau[i]
            if abs(d) < width:
                W_true += C(d)
                npairs += i != j
    pred = abs(np.exp(1j * (c - c_ref)).sum()) ** 2 / N
    print("\n[B: NP=%d bunch]" % N)
    print("  burst width in T %.3e; overlapping ordered pairs: %d" % (width, npairs))
    print("  exact energy / (N*W_ref)                         : %.4f" % (W_true / (N * ref["W"])))
    print("  LEADS energy / (N*W_ref)                         : %.4f" % (W_leads / (N * ref["W"])))
    print("  predicted for equal-index sum |sum e^{i dc}|^2/N : %.4f" % pred)
    print("  LEADS energy / exact                             : %.4f" % (W_leads / W_true))
    print("  spectrum.dat(NP=%d) / spectrum.dat(NP=1) at peak : %.4f  (incoherent by design: expect N=%d)" % (
        N, spec[np.argmax(ref['S']), 3] / ref["S"].max(), N))
    # exact signal spectrum at a few frequencies near the peak vs LEADS's
    w = ref["s"]
    Fref = nudft_complex(ref["tlab"], ref["f"], w)
    Ftrue = Fref * np.exp(1j * np.outer(w, 2 * a)).sum(1)[:, None]
    Strue = (np.abs(Ftrue) ** 2).sum(1)
    Sleads = nudft(tlab, f, w)
    k = slice(np.argmax(ref["S"]) - 2, np.argmax(ref["S"]) + 3)
    print("  FT^2 near peak, LEADS / exact:", " ".join("%.3g" % x for x in Sleads[k] / Strue[k]))


def default_bunch_ratio():
    NP, TR, posi = read_state(REPO / "TestCase" / "bench" / "state_init.bin")
    b = math.sqrt(1 - 1 / 500.0 ** 2)
    c = posi[2] + b * TR[0]
    return NP, abs(np.exp(1j * (c - c[0])).sum()) ** 2 / NP


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("variants", nargs="*", metavar="LABEL=FILES_DIR")
    ap.add_argument("--workdir", type=Path, default=Path("emitted_signal_runs"))
    a = ap.parse_args()
    variants = dict(v.split("=", 1) for v in (a.variants or ["current=" + str(REPO / "files")]))
    refs = {}
    for label, objdir in variants.items():
        run(a.workdir / ("%s_np1" % label), objdir, 1, False)
        refs[label] = case_a(label, a.workdir / ("%s_np1" % label))
    last = list(variants)[-1]
    run(a.workdir / ("%s_np20" % last), variants[last], 20, True)
    case_b(a.workdir / ("%s_np20" % last), refs[last], a.workdir / ("%s_np1" % last))
    NP, r = default_bunch_ratio()
    print("\nDefault reference bunch (NP=%d, gamma=500): equal-index energy / exact ~ %.0f" % (NP, r))


if __name__ == "__main__":
    main()

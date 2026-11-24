#!/usr/bin/env python3
"""Analytic Thomson / Compton check of the LEADS on-axis radiation spectrum.

Setup (all Model 0 plane wave, linear polarisation, radiation reaction off):
a single on-axis electron of Lorentz factor gamma collides head-on with a weak
(a0 << 1), long Gaussian pulse. LEADS computes the spectrum on the detector
axis theta = 180 deg, i.e. along the electron's direction of motion. In the
a0 -> 0 limit the fundamental peaks at the exact head-on backscatter frequency

    Thomson (classical):  w/w0 = (1+beta)/(1-beta)            ~ 4 gamma^2
    Compton (quantum):    w/w0 = (1+beta)/(1-beta + 2 hbar w0/(gamma m c^2))
                                                         ~ 4 gamma^2 / (1 + 4 gamma hbar w0/(m c^2))

LEADS is purely classical (no recoil), so its ground truth is the THOMSON value.
The Compton edge is printed so the size of the (unmodelled) recoil correction
can be compared with the numerical error: at 800 nm it is 1.2e-3 for gamma=100
and 6.1e-3 for gamma=500.

For a finite pulse the right reference is not (1+beta)/(1-beta) itself. Model 0
builds E = -dA/dt from A = a0 g(eta) cos(eta) with g = exp(-alpha eta^2) and
alpha = 4 ln2 / Tau_FWHM^2, so |E(w)|^2 = w^2 |A(w)|^2 peaks above w0, at

    w_E / w0 = (1 + sqrt(1 + 8 alpha)) / 2      (+0.140% for 10 cycles, +0.559% for 5)

and the backscattered peak sits at Thomson * w_E/w0. The script prints both.
(Before BUGFIXES.md item 2 the envelope term had no sign(eta), the field
spectrum was symmetric about w0, and the peak landed on plain Thomson.) The
measured peak tracks this factor at 5, 10, 20 and 40 cycles, with a residual
of -0.003% to -0.004% each time (2026-11-24, gamma = 100).

Residual systematics even for a correct code (all small for long pulses):
  * a0^2/2 redshift of the fundamental (not modelled; <= 5e-5 for a0 = 0.01,
    ~3e-5 measured after averaging over the pulse envelope).
  * LEADS's spectrum is Jackson 14.65, the acceleration form, which has no w^2
    prefactor. (validation/emitted_signal_check.py checks its normalisation by
    Parseval to 1e-5.) Weighting it by w^2 would count the factor twice.
  * At 40 cycles and longer the default --ds (0.9% per point) is coarser than
    the peak; use e.g. --ds 0.0005 --window -0.05 0.03.
  * Without the grid fix the result is not a fixed 0.5x: the retarded-time step
    error is (2*t_shift/dTau)'s int() truncation plus (1-beta), amplified by
    1/(1-beta), so the peak can land more than a decade below theory (widen
    --window to see it).

Usage (build LEADS first, e.g. `make all` in the repo root, so files/*.o exist):

    python3 -I validation/thomson_compton_check.py \
        --variant before=/path/to/old/files --variant after=/path/to/new/files

    # analyse a spectrum.dat you already have:
    python3 -I validation/thomson_compton_check.py --spectrum Data/spectrum.dat --gamma 100
"""
import argparse
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HC_EV_NM = 1239.841984     # h*c in eV*nm
MC2_EV = 510998.950        # electron rest energy in eV
REPO = Path(__file__).resolve().parent.parent


# ----------------------------------------------------------------- theory --
def beta_of(gamma):
    return math.sqrt(1.0 - 1.0 / gamma ** 2)


def thomson_peak(gamma):
    """Classical head-on backscatter, in units of w0: (1+b)/(1-b) = (1+b)^2 g^2."""
    b = beta_of(gamma)
    return (1.0 + b) ** 2 * gamma ** 2


def pulse_factor(cycles):
    """Peak of the Model 0 E-field spectrum in units of w0 (Gaussian envelope, PulseProfile = 1)."""
    alpha = 4.0 * math.log(2.0) / (2.0 * math.pi * cycles) ** 2
    return (1.0 + math.sqrt(1.0 + 8.0 * alpha)) / 2.0


def compton_edge(gamma, wavelength_nm):
    """Quantum head-on edge with electron recoil, in units of w0 (exact kinematics)."""
    b = beta_of(gamma)
    hw0 = HC_EV_NM / wavelength_nm                       # photon energy, eV
    one_minus_b = 1.0 / (gamma ** 2 * (1.0 + b))         # avoids cancellation
    return (1.0 + b) / (one_minus_b + 2.0 * hw0 / (gamma * MC2_EV))


# ------------------------------------------------------------ peak finder --
def find_peak(s, y, half=3):
    """Peak position of y(s): argmax refined by a parabola in (s, ln y).

    The scattered spectrum of a Gaussian pulse is Gaussian in w, so ln y is a
    parabola in s. (A parabola in ln s is NOT: it is biased low by ~0.03% at
    half=3 and drifts to ~0.1% as the fit window widens. Verified on these runs:
    the (s, ln y) fit gives -0.0030% for every window half-width 1..6.)
    """
    i = int(np.argmax(y))
    if i < half or i > len(s) - 1 - half:
        return float("nan")                              # peak at window edge
    sl = slice(i - half, i + half + 1)
    x = (s[sl] - s[i]) / s[i]                            # centred/scaled for conditioning
    c = np.polyfit(x, np.log(y[sl]), 2)
    return s[i] * (1.0 - c[1] / (2.0 * c[0]))


def load_spectrum(path):
    d = np.loadtxt(path)                                 # th, ph, w/w0, intensity
    return d[:, 2], d[:, 3]


# ------------------------------------------------------------- LEADS runs --
def make_input(template, params):
    text = template
    for name, val in params.items():
        pat = re.compile(r"(^[ \t]*%s[ \t]*=[ \t]*)[^;]*;" % re.escape(name), re.M)
        text, n = pat.subn(lambda m, v=val: m.group(1) + v + ";", text)
        if n != 1:
            raise RuntimeError("input.cpp: expected exactly one '%s = ...;' line, found %d" % (name, n))
    return text


def run_case(label, files_dir, gamma, a, workdir):
    rundir = workdir / ("%s_g%g" % (label, gamma))
    (rundir / "Data").mkdir(parents=True, exist_ok=True)
    subprocess.run(["rm", "-rf", str(rundir / "files")], check=True)
    subprocess.run(["cp", "-r", str(files_dir), str(rundir / "files")], check=True)

    wt = thomson_peak(gamma)
    T = 2.0 * math.pi * a.cycles
    params = {
        "RadiationReaction": "false",
        "PulseModel": "0",
        "Polarization": "0",
        "a0": repr(a.a0),
        "Tau_FWHM": "%g * (2*pi)" % a.cycles,
        "TimeStep": repr(a.dtau),
        "NoBeamPart": "1",
        "GamaParticle": repr(float(gamma)),
        "PropTheta": "180.0",
        "PropPhi": "0.0",
        "z0": "%.6f" % (a.zfac * T),          # window covers the pulse to <~1e-4 of peak
        "CalculateSpectrum": "true",
        "AngularDistribution": "false",
        "SignalEmitted": "false",
        "DetecThetaStart": "180.0", "DetecThetaEnd": "180.0",
        "DetecPhiStart": "0.0", "DetecPhiEnd": "0.0",
        "UseLogScale": "true",
        "SpectrumStart": "%.6f" % (math.log10(wt) + a.window[0]),
        "SpectrumEnd": "%.6f" % (math.log10(wt) + a.window[1]),
        "SpectrumDO": repr(a.ds),
        "StoreTrajectories": "true",
        "DumpStateBinary": "false",
        "WriteSpatialFieldsProfile": "false",
    }
    (rundir / "input.cpp").write_text(make_input((REPO / "input.cpp").read_text(), params))
    # same commands as RunMakefile, then run
    for cmd in ("g++ -std=c++17 -O3 -o leads ./files/*.o input.cpp -larmadillo -lgsl -lgslcblas -lm",
                "./leads 1"):
        r = subprocess.run(cmd, shell=True, cwd=rundir, capture_output=True, text=True)
        (rundir / "run.log").write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            raise RuntimeError("'%s' failed in %s:\n%s" % (cmd, rundir, (r.stdout + r.stderr)[-800:]))
    return rundir / "Data" / "spectrum.dat"


# ----------------------------------------------------------------- report --
def analyse(path, gamma, wavelength_nm, cycles):
    s, y = load_spectrum(path)
    wt, wc = thomson_peak(gamma), compton_edge(gamma, wavelength_nm)
    return wt, wc, wt * pulse_factor(cycles), find_peak(s, y)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", action="append", metavar="LABEL=FILES_DIR",
                    help="a LEADS build to test (dir holding leads_*.o and *.h); repeatable. "
                         "Default: current=<repo>/files (must be freshly built).")
    ap.add_argument("--gammas", type=float, nargs="+", default=[50, 100, 500])
    ap.add_argument("--a0", type=float, default=0.01)
    ap.add_argument("--cycles", type=float, default=10.0, help="Tau_FWHM in laser cycles")
    ap.add_argument("--dtau", type=float, default=5e-3)
    ap.add_argument("--zfac", type=float, default=1.0, help="z0 = zfac * Tau_FWHM (1/omega units)")
    ap.add_argument("--wavelength-nm", type=float, default=800.0)
    ap.add_argument("--window", type=float, nargs=2, default=[-0.45, 0.15],
                    help="log10 spectrum window relative to the Thomson peak")
    ap.add_argument("--ds", type=float, default=0.004, help="log10 spacing of the spectrum grid")
    ap.add_argument("--workdir", type=Path, default=Path("thomson_compton_runs"))
    ap.add_argument("--spectrum", type=Path, help="analyse this spectrum.dat only (needs --gamma)")
    ap.add_argument("--gamma", type=float)
    a = ap.parse_args()

    rows = []
    if a.spectrum:
        if a.gamma is None:
            ap.error("--spectrum needs --gamma")
        rows.append(("file", a.gamma) + analyse(a.spectrum, a.gamma, a.wavelength_nm, a.cycles))
    else:
        variants = [v.split("=", 1) for v in (a.variant or ["current=" + str(REPO / "files")])]
        a.workdir.mkdir(parents=True, exist_ok=True)
        jobs = [(lab, Path(d), g) for lab, d in variants for g in a.gammas]
        with ThreadPoolExecutor(max_workers=min(len(jobs), os.cpu_count() or 1)) as ex:
            futs = [ex.submit(run_case, lab, d, g, a, a.workdir) for lab, d, g in jobs]
            for (lab, d, g), f in zip(jobs, futs):
                rows.append((lab, g) + analyse(f.result(), g, a.wavelength_nm, a.cycles))

    print("\nOn-axis (theta=180) peak, head-on, a0=%g, %g cycles; errors are (peak/theory - 1)" % (a.a0, a.cycles))
    print("%-8s %6s %12s %12s %12s | %12s %9s %9s %9s" %
          ("variant", "gamma", "Thomson", "Thom*pulse", "Compton", "peak", "vs Thom", "vs pulse", "vs Comp"))
    for lab, g, wt, wc, wp, pk in rows:
        print("%-8s %6g %12.5g %12.5g %12.5g | %12.5g %+8.3f%% %+8.3f%% %+8.3f%%" %
              (lab, g, wt, wp, wc, pk, 100 * (pk / wt - 1), 100 * (pk / wp - 1), 100 * (pk / wc - 1)))
    print("pulse factor w_E/w0 for %g cycles: %.6f (see docstring)" % (a.cycles, pulse_factor(a.cycles)))
    print("\nCompton recoil shift (Thomson->Compton edge): " +
          ", ".join("gamma=%g: %+.3f%%" % (g, 100 * (compton_edge(g, a.wavelength_nm) / thomson_peak(g) - 1))
                    for g in a.gammas))
    if any(math.isnan(r[5]) for r in rows):
        print("NaN peak = maximum at the window edge; widen --window.", file=sys.stderr)


if __name__ == "__main__":
    main()

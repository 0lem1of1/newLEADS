# Build the CPU/GPU Boris-pusher benchmark tools.
#
# ARCH selects the CUDA target (default: this dev machine's RTX 3050
# Laptop, compute capability 8.6). Override for other hardware, e.g.:
#   make bench_gpu ARCH=compute_75,code=sm_75

CXX      = g++
CXXFLAGS = -std=c++17 -O3

NVCC      = nvcc
ARCH      = compute_86,code=sm_86
NVCCFLAGS = -std=c++17 -O3 -gencode arch=$(ARCH)

all: bench_cpu bench_gpu convtest statediff

bench_cpu: bench_cpu.cpp leads_boris.cuh ../leads_state.h ../leads_params.h
	$(CXX) $(CXXFLAGS) -o $@ bench_cpu.cpp

bench_gpu: bench_gpu.cu leads_boris.cuh ../leads_state.h ../leads_params.h
	$(NVCC) $(NVCCFLAGS) -o $@ bench_gpu.cu

convtest: convtest.cpp
	$(CXX) $(CXXFLAGS) -o $@ convtest.cpp

statediff: statediff.cpp
	$(CXX) $(CXXFLAGS) -o $@ statediff.cpp

clean:
	rm -f bench_cpu bench_gpu convtest statediff

.PHONY: all clean

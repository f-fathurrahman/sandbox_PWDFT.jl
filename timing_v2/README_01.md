For convenience initialize Hamiltonian first, then run the benchmark for specific functions.


Some results for O2 spinpol
```
julia> bench_PsPotNL(Ham_gth)  # quite fast
Benchmark: 34 samples with 1 evaluation
 min    2.305 ms (13 allocs: 546.585 KiB)
 median 2.502 ms (13 allocs: 546.616 KiB)
 mean   2.969 ms (13 allocs: 546.610 KiB)
 max    17.739 ms (13 allocs: 546.632 KiB)

julia> bench_PsPotNL(Ham_oncv)
Benchmark: 1 sample with 1 evaluation
        334.672 ms (518413 allocs: 20.752 MiB, without a warmup)

julia> bench_PsPotNL(Ham_gbrv)
Benchmark: 1 sample with 1 evaluation
        3.069 s (5283233 allocs: 95.098 MiB, 0.68% gc time, without a warmup)

julia> bench_PsPotNL(Ham_paw_jth)
Benchmark: 1 sample with 1 evaluation
        7.097 s (12761065 allocs: 203.732 MiB, 0.17% gc time, without a warmup)
```

For Pt atom
```
julia> bench_PsPotNL(Ham_gth)
Benchmark: 5 samples with 1 evaluation
 min    22.910 ms (13 allocs: 5.063 MiB)
 median 23.578 ms (13 allocs: 5.063 MiB)
 mean   28.248 ms (13 allocs: 5.063 MiB, 11.95% gc time)
 max    46.772 ms (13 allocs: 5.063 MiB, 51.32% gc time)

julia> bench_PsPotNL(Ham_oncv)
Benchmark: 1 sample with 1 evaluation
        782.042 ms (1032350 allocs: 34.813 MiB, 11.75% gc time, without a warmup)

julia> bench_PsPotNL(Ham_gbrv)
Benchmark: 1 sample with 1 evaluation
        4.152 s (8160394 allocs: 138.200 MiB, 0.14% gc time, without a warmup)

julia> bench_PsPotNL(Ham_paw_jth)
Benchmark: 1 sample with 1 evaluation
        11.173 s (21694959 allocs: 345.127 MiB, 0.13% gc time, without a warmup)
```
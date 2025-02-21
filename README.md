# pdac-efs-bench2024

Restore R libraries required to run the benchmark using `renv`:
```
renv::restore(exclude = c("BiocGenerics", "BiocManager", "BiocVersion", "ComplexHeatmap", "IRanges", "S4Vectors", "MASS", "ipred", "class"))
```


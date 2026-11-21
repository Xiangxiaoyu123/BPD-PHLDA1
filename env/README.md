# Environment

All analyses were run in R 4.3.2 with Seurat 4.4.0.

To record the environment of your own run:

```r
writeLines(capture.output(sessionInfo()), "env/sessionInfo.txt")
write.csv(installed.packages()[, c("Package", "Version")],
          "env/installed_packages.csv", row.names = FALSE)
```

The preprocessing and annotation steps follow the environment of
https://github.com/SucreLab/HumanNeonatal (R 4.3.2, Python 3.10.0).

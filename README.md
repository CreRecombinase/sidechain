# Allotrope HSDS/Shiny Demo
To run this demo, you'll want to make sure you've installed `R`, `Rstudio`, and all the packages in the `Rmd` document.  
Note that `rhdf5client` is a bioconductor package and should be installed as such:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rhdf5")
```

It should also be noted that the python chunk is set to `eval=FALSE`, meaning that it is not evaluated.  You'll want the python library `h5pyd` installed
before trying that out.

Please send questions to Nicholas.Knoblauch@HDFGroup.org

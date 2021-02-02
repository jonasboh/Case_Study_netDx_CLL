# Case_Study_netDx_CLL
This case study illustrates the usage of supervised multi-omics data integration.


Before starting the analysis we generate an output folder in the working directory of the user in order to not overwrite something.

```{r create output files for the anaylsis}
# generate an output file in the home folder for storing the netDx output
dt <- format(Sys.time(),"%y%m%d_%H%M%S")
megaDir <- getwd()
subdir1 <- "netDx_output"
dir.create(file.path(megaDir, subdir1), showWarnings = FALSE)
subdir2 <- "CLL"
outDir <- sprintf("%s/%s/%s_%s",megaDir,subdir1,subdir2,dt)
# check is file exists in order to not overwrite something
if (file.exists(outDir)) unlink(outDir,recursive=TRUE); 
dir.create(outDir)
logFile <- sprintf("%s/log.txt",outDir)
```


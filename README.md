# Case Study on Chronic lymphocytic leukemia (CLL)

This case study illustrates the usage of supervised multi-omics data integration published in our review "A Detailed Catalogue of 
Multi-Omics Methodologies for Identification of Putative Biomarkers and 
Causal Molecular Networks in Translational Cancer Research".
The used data set is included in the MOFAdata package [[1]](#1) and is a part of the review manuskript (Evaluation of multi-omics methodologies and tools for identification of putative biomarkers and causal molecular networks in translational cancer research).  
We aim to illustrate the usage of netDx [[2,3]](#2#3) as an example for supervised multi-omics data integration.

## Start supervised multi-omics analysis with netDx

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

### 1. Preprocessing

Loading and preprocessing of clinical data. Remove samples without known IGHV mutation status.

```{r prepare clinical metadata}
message("Preparing clinical data")
# load clinical data
clinical.data <- read.delim("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt",sep="\t",h=T,as.is=T)
# add variable called ID which includes sample IDs
clinical.data$ID <- clinical.data$sample
  
# generate dataframe with sampleIDs as rownames
clinical.data <- as.data.frame(clinical.data)
rownames(clinical.data) <- clinical.data$sample
  
# filter for NA values in IGHV status
na_idx <- which(is.na(clinical.data$IGHV))
clinical.data <- clinical.data[-na_idx,]
```

```{r include Status variable to clinical metadata for netDx}
# include Status variable for clustering and convert it to string
str <- sprintf("IGHV_%i",clinical.data$IGHV)     ## make a string for status description
clinical.data$STATUS <- str
```

Loading molecular data from MOFAdata package without further preprocessing for missingness.

```{r prepare molekular data}  
# load molecular CLL data from MOFAdata package
message("loading omic data")
utils::data("CLL_data")       
# Sample and feature distribution of CLL omics data  
print("Data")
sapply(names(CLL_data), function(nm){
  message(sprintf("\t%s",nm))
  print(dim(CLL_data[[nm]]))
})
```

Convert Ensembly IDs fo gene symbols for the detection of associated pathways in the gene expression layer.

```{r convert Ensemble Gene IDs to Gene IDs for Pathway prediction}    
message("converting ensembl IDs to hgnc IDs")
# Preprocessing for gene symbols from hg19 (GRCh37)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)         ## connect to Biomart DB
symbs <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), mart=ensembl)               ## list of symbols matching Ensemble IDs
message(sprintf("\t%i genes in Biomart Database",nrow(symbs)))
symbs <- symbs[-which(symbs[,1]==""),]                                                    ## filter for Ensemble IDs without symbols
message(sprintf("\t%i genes with matching gene symbols in Biomart Database",nrow(symbs)))
# Preprocessing for converting Ensemble IDs to gene symbols
idx <- which(!rownames(CLL_data[["mRNA"]]) %in% symbs$ensembl_gene_id)                    ## select Ensemble IDs in mRNA layers
message(sprintf("%i of %i genes from expression layer don't have HGNC symbols", length(idx),nrow(CLL_data$mRNA)))
if (length(idx)>0) {
	message(sprintf("\tremoving RNA features w/o HGNC: %i entries", length(idx)))
	CLL_data$mRNA <- CLL_data$mRNA[-idx,]                                                   ## remove genes without matching gene symbol 
}
midx <- match(rownames(CLL_data[["mRNA"]]),symbs$ensembl_gene_id)                         ## match Ensemble IDs in mRNA layer to gene symbols
if (all.equal(symbs$ensembl_gene_id[midx],rownames(CLL_data$mRNA))!=TRUE) {
	stop("ensembl gene mapping didn't work")
}
rownames(CLL_data$mRNA) <- symbs$hgnc_symbol[midx]                                        ## replace Ensemble IDs by gene symbols
```  

Reorder samples based on ID.

```{r match sample order for all data layer}    
message("checking sample order in all layers")
pID <- clinical.data$ID
for (k in names(CLL_data)) {
	message(k)
	cur <- CLL_data[[k]]
	cur <- cur[,pID]                        ## reorder to match sample order in clinical.data
	CLL_data[[k]] <- cur
}
```

Configure MultiAssayExperiment object as input object for netDx.

```{r creating MultiAssayExperiment from SummarizedExperiment}    
message("Compile experiments")
seList <- list()
for (nm in names(CLL_data)){
	seList[[nm]] <- SummarizedExperiment(CLL_data[[nm]],colData=clinical.data)
}
message("Make MultiAssayExperiment")
MultiAssay <- MultiAssayExperiment(seList,colData=clinical.data)
message("Show Distribution of Samples in MAE")
upsetSamples(MultiAssay)
```
### 2. Configure netDx

Some configuration is needed for specific design of experiments. The gene expression layer will be grouped by pathways where all other layers will have a single feature per layer.

```{r creating grouplist for Patient similarity network construction}    
message("Make group list")
groupList <- list()
message(sprintf("\tRNA layer - pathways"))
# parse GMT file and return the current pathway list
groupList$mRNA <- readPathways(fetchPathwayDefinitions("January",2021,13,verbose=TRUE))
# include grouplist entries for other layers
message("\tother layers - one net per feature")
for (nm in setdiff(names(CLL_data),"mRNA")) {
	blah <- list(tmp=rownames(CLL_data[[nm]]))
	names(blah)[1] <- nm
	groupList[[nm]] <- blah	
}
```

Created custom function based on the design of experiments for this study.

```{r creating function for making Networks}
# creating networks by using RNA expression, Methylation and drug response (excluded Mutations)
makeNets <- function(dataList, groupList, netDir,...) {
  netList <- c() 
	for (nm in c("mRNA","Methylation","Drugs")) {
	  if (!is.null(groupList[[nm]])) { 
	    cur <- makePSN_NamedMatrix(dataList[[nm]],
	                                   rownames(dataList[[nm]]),
	                                   groupList[[nm]],
	                                   netDir,verbose=FALSE, 
	                                   writeProfiles=TRUE,...) 
			message(sprintf("\t%s: %i nets",nm,length(unlist(cur))))
			netList <- c(netList,unlist(cur))
	  }
	}
  return(netList)
}
```

Configure seetings for netDx run and usage of CPUs for faster results.

```{r define parameter for creation of similarity networks suitable for sample size}   
# put in list to save config settings
settings <- list(
	featScoreMax=10L,                       ## maximum score for features (recommended:10)
	featSelCutoff=9L,                       ## threshold to call feature-selected networks for each train/test (recommended:9)
	numSplits=10L,                          ## number of train/test splits (recommended: 100)
	trainProp=0.8,                          ## part of each split for training (default:80/20)
	setSeed=1876,                             ## for reproducibility
	netFunc=makeNets,
	numCores=1                              ## number of used CPUs
)
set.seed(settings$setSeed)                ## make results reproducible
```  

### 3. Run netDx on CLL data

Execute the predictor from netDx and write results to the homedirectory of the user.

```{r configure files for output and run the predictor}
sink(logFile,split=TRUE)
timeStart<- Sys.time()
message("####\tRunning predictor\t####")
out <- buildPredictor(
  dataList=MultiAssay,                    ## patient data (dataList)
  groupList=groupList,                    ## grouping variables into networks (groupList)
  makeNetFunc=makeNets,
  outDir=sprintf("%s/pred",outDir),       ## netDx requires absolute path for output dir
  numSplits=settings$numSplits,
  trainProp=settings$trainProp,
  featScoreMax=settings$featScoreMax,
  featSelCutoff=settings$featSelCutoff,
  numCores=settings$numCores,
  logging="none",
  keepAllData=FALSE,
)
timeEnd<- Sys.time()
difference <- round(difftime(timeEnd, timeStart, units='mins'), digits = 2)
message("####\tPredictor finished\t####\nCalculation time: ",difference," minutes\nUsed CPUs: ", settings$numCores)
```  

### 4. Compile results and evaluate performance

This section is dealing with the evaluation of results. The main part is performance visualization.
```{r compiling results for performance evaluation}   
# The results are stored in the list object returned by the Predictor. 
# This list contains:
# inputNets: all input networks that the model started with.
# Split<i>: a list with results for each train-test split
# predictions: real and predicted labels for test patients
# accuracy: percent accuracy of predictions
# featureScores: feature scores for each label (list with g entries, where g is number of patient labels). Each entry contains the feature selection scores for the corresponding label.
# featureSelected: vector of features that pass feature selection. List of length g, with one entry per label.
message("Compiling results")
summary(out)                              ## get summary of performed splits
# collect results
st <- unique(colData(MultiAssay)$STATUS)
acc <- c()                                ## accuracy
predList <- list()                        ## prediction tables
featScores <- list()                      ## feature scores per class
for (cur in unique(st)) featScores[[cur]] <- list()
for (k in 1:settings$numSplits) { 
	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
	# predictions table
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
	                 sprintf("%s_SCORE",st))]
	predList[[k]] <- tmp 
	# accuracy
	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
	# feature scores
	for (cur in unique(st)) {
	   tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
	}
}
```

Plotting ROC and PR curve

```{r plotting performance evaluation}
# plot ROC and PR curve, compute AUROC, AUPR
predPerf <- plotPerf(predList, predClasses=st)
# save plots as pdf in output folder
pdf(sprintf("%s/perf.pdf",outDir))
tryCatch({
	predPerf <- plotPerf(predList, predClasses=st)
}, error=function(ex){
	print(ex)
},finally={
	dev.off()
})
```

Print Performance values with standard deviation

```{r performance values}  
message("Printing performance values")
auroc <- unlist(lapply(predPerf, function(x) x$auroc))
aupr <- unlist(lapply(predPerf, function(x) x$aupr))
message(sprintf("AUROC = %2.1f%% (%2.1f%%)", mean(auroc*100),sd(auroc*100)))
message(sprintf("AUPR = %2.1f%% (%2.1f%%)", mean(aupr*100),sd(aupr*100)))
``` 

Selecting most important features for further similarity network generation

```{r check for good performing features}   
# get table of feature scores for each split and patient label
featScores2 <- lapply(featScores, getNetConsensus)
summary(featScores2)
# identify features that consistently perform well: features that scored at least 2 out of 2 in 50% or more splits
featSelNet <- lapply(featScores2, function(x) {
    callFeatSel(x, fsCutoff=2, fsPctPass=0.5)
})
```

Save all results to final rda file in the output dir

```{r writing results to file}   
message("saving input and results")
save(MultiAssay,groupList,settings,out,predList,featScores2,featSelNet,
	file=sprintf("%s/results.rda",outDir))
``` 

Prepare the data for plotting the EnrichmentMap in Cytoscape

```{r visualize pathway connected features}   
message("Creating input for EnrichmentMap")
# load recent pathway file 
pathList <- readPathways(fetchPathwayDefinitions("January",2021,13,verbose=TRUE))
Emap_res <- getEMapInput_many(
  featScores2,
  pathList,
  minScore=3,
  maxScore=10,
  pctPass=0.5,
  out$inputNets,
  verbose=FALSE)
gmtFiles <- list()
nodeAttrFiles <- list()
for (g in names(Emap_res)) {
    outFile <- paste(outDir,sprintf("%s_nodeAttrs.txt",g),sep=getFileSep())
    write.table(Emap_res[[g]][["nodeAttrs"]],file=outFile,
        sep="\t",col=TRUE,row=FALSE,quote=FALSE)
    nodeAttrFiles[[g]] <- outFile
    outFile <- paste(outDir,sprintf("%s.gmt",g),sep=getFileSep())
    conn <- suppressWarnings(
         suppressMessages(base::file(outFile,"w")))
    tmp <- Emap_res[[g]][["featureSets"]]
    gmtFiles[[g]] <- outFile
    for (cur in names(tmp)) {
        curr <- sprintf("%s\t%s\t%s", cur,cur,
            paste(tmp[[cur]],collapse="\t"))
        writeLines(curr,con=conn)
    }
close(conn)
}
```

### 5. Plotting EnrichmentMap Similarity Networks

In order to plot the Networks Cytoscape (v.3.8.2) is required.
For installing Cytoscape, take a look at: https://cytoscape.org
The applications EnrichmentMap and AutoAnnotate are also needed to create the similarity networks.
EnrichmentMap and AutoAnnotate are part of the Enrichment Pipeline Collection.
Installation of these applications can be handled inside Cytoscape (Apps/App Manager/Install Apps). 
Follow this guide for further advice: https://cytoscape.org/cytoscape-tutorials/protocols/enrichmentmap-pipeline/#/cytoscape_apps
Cytoscape should run in parallel to the execution of the following code.
Visualisation is handeled to Cytoscape. Further configuration of diagrams can be done in Cytoscape. 

Create networks for clustering connected pathways for each patient
nodes: predictive pathways, edges: similarity in gene-sets. 

```{r creating similarity networks for predictive features for IGHV status 1}
plotEmap(gmtFiles[[1]],nodeAttrFiles[[1]],
         groupClusters=TRUE, hideNodeLabels=TRUE)
```

```{r creating similarity networks for predictive features for IGHV status 0}
plotEmap(gmtFiles[[2]],nodeAttrFiles[],
         groupClusters=TRUE, hideNodeLabels=TRUE)
```

### 5. Plotting Patient Similarity Networks

Setting more strict values for patient similarity network stratification.

```{r filter strict for good performing features}   
# identify features that consistently perform well: features that scored at least 2 out of 2 in all splits
strict_featSelNet <- lapply(featScores2, function(x) {
    callFeatSel(x, fsCutoff=2, fsPctPass=1)
})
print(strict_featSelNet)
```
Generate new list for grouping according to strict filtered values.

```{r creating new groupList for the new strict filtered features} 
topPath <- gsub(".profile","",
        unique(unlist(strict_featSelNet)))
topPath <- gsub("_cont.txt","",topPath)
# create groupList limited to top features
g2 <- list();
for (nm in names(groupList)) {
    cur <- groupList[[nm]]
    idx <- which(names(cur) %in% topPath)
    message(sprintf("%s: %i pathways", nm, length(idx)))
    if (length(idx)>0) g2[[nm]] <- cur[idx]
}
```

Integrating the networks with the mean of the edge weights 

```{r creating patient symilarity networks}
# measure significance with one-sided Wilcoxon-Mann-Whitney test
psn <- suppressMessages(
  plotIntegratedPatientNetwork(
    MultiAssay,
    groupList=g2, 
    makeNetFunc=makeNets,
    aggFun="MEAN",
    prune_pctX=0.30,
    prune_useTop=TRUE,
    numCores=1L,
    calcShortestPath=TRUE,
    showStats=FALSE,
    verbose=FALSE, 
    plotCytoscape=TRUE)
)
```

Visualize the patient similarity network in a tSNE framework

```{r visualize integrated patient symilarity network as a tSNE plot}
tsne <- plot_tSNE(psn$patientSimNetwork_unpruned,colData(MultiAssay))
```

System information for reproducibility

```{r Information about System configuration}   
sessionInfo()
``` 

## References

<a id="1">[1]</a> 
Argelaguet R, Velten B, Arnol D, Buettner F, Huber W, Stegle O (2020). MOFAdata: Data package for Multi-Omics Factor Analysis (MOFA). R package version 1.6.0.

<a id="2">[2]</a> 
Pai, S., Hui, S., Isserlin, R., Shah, M. A., Kaka, H., & Bader, G. D. (2019). netDx: interpretable patient classification using integrated patient similarity networks. Mol Syst Biol, 15(3), e8497. doi:10.15252/msb.20188497

<a id="3">[3]</a> 
Pai, S., Weber, P., Isserlin, R., Kaka, H., Hui, S., Shah, M. A., . . . Bader, G. D. (2020). netDx: Software for building interpretable patient classifiers by multi-'omic data integration using patient similarity networks. F1000Research, 9. doi:10.12688/f1000research.26429.1

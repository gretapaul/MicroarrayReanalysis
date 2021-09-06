# MicroarrayReanalysis
Microarray reanalysis using Bioconductor 
This is a tutorial workflow of the methodology used to perform analysis of data from repositories via Bioconductor. This workflow was influenced by two other workflows (1.) and (2.), as well as drafted from the help pages of Biostars, StackExchange, and reference manuals of packages that are ussed.

## Data Preparation
Set the libraries that will be needed.
```
library(GEOquery)
library(Biobase)
library(dplyr)
library(hgu133a.db)
library(factoextra)
library(limma)
library(readr)
library(topGO)
library(ReactomePA)
```
Download the data from GEO into our workplace. For this analysis, we are using matrix processed data files.
```
gset <- getGEO("GSE68468", GSEMatrix = TRUE, getGPL = FALSE)
````
Next, we see how many platforms are used in the experiment, as GEO queries often are joint. Select the platform for analysis.
```
length(gse)
gset <- gset[[1]]
```
For this session, we are only interested in cancer and healthy samples. The data of interest will be retrieved by filtering out by phenotype/metadata attached to the experiment. To identify which columns present the most apparent differences between the samples.
```
pData(gset.filter)
```
In this case, it is the characteristics_ch1.4 column. Write an aggregate function to filter.
```
filter <- colnames(gset)[gset@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer"|gset@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"]
```
To check if the filter was successful, call the length function. We are looking for 241 samples combined.
```
length(filter)
```
Apply the filter to the expression set.
```
gse.filter <- gset[,filter]
```
Remove the Affymetrix control probes as the GSE entry already came normalised and summarised.
```
gse.filter <- gse.filter[-grep(`^AFFX`, rownames(gse.filter)),]
```
## Log transformation
Log2 transformation is a common method of normalizing the expression data.
```
exprs(gse.filter) <- log2(exprs(gse.filter))
```
## Annotation 
Linking identifiers with annotation Affymetrix probes come with manufacturers names. These identifiers are linked to gene names. Let us connect the probe I.D.s with gene names. This way, the analysis will be more accessible, as if we want to study single gene expression or make sense of the results, we will be able to call out the functions by the annotated gene name. The platform this experiment is based on is the GPL96(or hgu133a). For that, we need GPL96 specific annotation pack. Once we load that, we need to check for key types and columns in the package. The particular naming of these key types and columns are essential to execute successful linking.
```
keytypes(hgu133a.db)
columns(hgu133a.db)
```
Our columns of interest are "GENENAME" and "SYMBOL". In keys, we are linking this GPL96 information via the manufacturers' I.D.s given on a dataset, thus via "PROBEID". On our dataset, PROBEID is listed.
```
featureNames(gse.filter)
```
Make sure that the features names match to our keytype PROBEID.
```
head(keys(hgu133a.db, keytype = "PROBEID"))
```
Write an aggregate function to extract the data we need. For this, we use the AnnotationDbi feature that implements SQLite querying.
```
anno <- AnnotationDbi::select(hgu133a.db,
keys = (featureNames(gse.filter)),
columns = c("SYMBOL", "GENENAME"),
keytype = "PROBEID")
```
Subset the data.
```
anno <- subset(anno, !is.na(SYMBOL))
```
Affymetrix experiments are often designed so that some genes will have multiple dedicated probes that target the same gene. While they are helpful for other purposes, we are removing them for this analysis as they might skew PCA and D.E.
```
annogroup <- group_by(anno, PROBEID)
annosum <- dplyr::summarize(annogroup, no_of_matches = n_distinct(SYMBOL))
```
Separate matches over 1.
```
annofilter <- filter(annosum, no_of_matches = n_distinct(SUMBOL))
nrow(annofilter)
idx <- (featureNames(gse.filter) %in% annofilter&PROBEID)
table(idx)
prep.gse <- subset(gse.filter, !idx)
validObject(prep.gse)
fData(prep.gse)$PROBEID <- rownames(fData(prep.gse))
```
Now only thing is left is to left-join the datasets.
```
fData(prep.gse) <- left_join(fData(prep.gse), anno)
rownames(fData(prep.gse)) <- fData(prep.gse)$PROBEID
```
## Principal Component Analysis (PCA)
Now, we move on to some analysis methods. Reduction of dimensionality is made by performing PCA. It is crucial to transpose the expression set.
```
exps <- log2(Biobase::exprs(annofinal))
PCA <- prcomp(t(exps), scale = FALSE)
```
Calculate the PCA scores in percentages
```
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
```
Calculate standard deviation ratio
```
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
```
Create a data set with the frame of pca results, add metadata
```
dataplot <- data.frame(PC1 = PCA$y[,1], PC2 = PCA$y[,2], Phenotype =  Biobase::pData(annofinal)$characteristics_ch1.1)
Plot the graph.
ggplot(dataplot aes(PC1, PC2)) +
  geom_point(aes(shape = Disease)) +
   ggtitle("PCA plot cancer and normal colon data") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
     ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
       coord_fixed(ratio = sdratio) +
         scale_shape_manual(values = c(3,17)) 
```
Summary of the PCA provides that the dimensionality of the data was reduced to PC241. We started with rotation (n x k) = (21061 x 241)
```
summary(PCA)
biplot(PCA)
library(factoectra)
fviz_eig(PCA, addlabels=TRUE, xlim=c(1,10))
```
Summary the first 6 eigenvectors.
```
eig.val <- get_eigenvalue(PCA)
head(eig.val)
```
Explore the data more. Cosine contribution
```
fviz_cos2(PCA, choice = "var", axes = 1:2, xlim=c(1,10))
Contributions of PC1, PC2. and P.C. 1 and 2.
fviz_contrib(PCA, choice = "var", axes = 1, top = 10)
fviz_contrib(PCA, choice = "var", axes = 2, top = 10)
fviz_contrib(PCA, choice = "var", axes = 1:2, top = 10)
fviz_contrib(PCA, choice = "ind", axes=1:2, top = 10
```
Note that there is more reading available in the factoextra package and other paper too, the cos2 contributions to the first dimensions is a great indication of the variability of genes. The cos2 plot shows which speicifc gene contributes and is represented the best in the principcal components 1 annd 2. It elucidates over expressed or downregulated genes without supervision. 

## Differential Expression
Following the step-by-step provided by the limma documentation. First, to do any testing with limma, a matrix needs to be designed. A matrix will provide the backbone of the linear fitting. We create the matrix for this experiment based on the phenotype we are analysing. Out experiment aim is to compare the differential expression between normal colon mucosa and cancer. Thus, we are performing linear fitting to solve this problem. Depending on the situation itself, the design matrix needs to be altered.
```
metadata <- pData(prep.gse)
design <- model.matrix(~0+metadata$characteristics_ch1.4)
colnames(design) <- c("Cancer", "Normal")
```
It is advised to do filtering based on the median expression. Since we did not apply any of the normalisation techniques in this experiment, I believe it is essential not to skip this recommendation by limma manual.
```
cutoff <- median(exprs(prep.gse))
expressed <- exprs(prep.gse) > cutoff
keep <- rowSums(expressed) > 2
datafit <- prep.gse[keep,]
```
It fits the data. The typical pipeline will fit the data twice. In this experiment, we will also adjust the weights of the data to minimise the false discovery rate.
```
fit <- lmFit(exprs(datafit), design)
contrasts <- makeContrasts(Cancer - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
```
Calculate array weightsa and process outliers.
```
aw <- arrayWeights(exprs(datafit), design)
fit <- lmFit(exprs(datafit), design, weights=aw)
contrasts <- makeContrasts(Cancer - Normal, levles=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
```
Now we can add phenotype information on the fitted tables and export results for further discussion.
```
phenodata <- fData(prep.gse)
phenodata <- select(phenodata, SYMBOL, GENENAME)
fit2$genes <- phenodata
results <- topTable(fit2, number=Inf)
results <- tibble::rownames_to_columns(results, "ID)
filter(results, adj.P.Val < 0.05, abs(logFC) > 1) %>% write_csv(path="DE_Results_Report")
```
## Gene Ontology (G.O.) enrichment analysis
GO-based enrichment analysis gives us a good overview of the down and upregulated pathways between the two datasets.
We will need to perform a new limma sequence, this time without array weights.
```
data_design <- prep.gse[, prep.gse$characteristics_ch1.4 %in% c("histology:
colon cancer", "histology: normal colon mucosa")]

prep.gse$characteristics_ch1.4 <- factor(prep.gse$characteristics_ch1.4)
design <- model.matrix(~ prep.gse$characteristics_ch1.4)
fit0 <- lmFit(prep.gse, design)
fit0 <- eBayes(fit0)
design2 <- model.matrix(~ prep.gse$characteristics_ch1.4 - 1)
colnames(design2) <- c("cancer", "Normal")
fit02 <- lmFit(prep.gse, design2)
contrast.matrix <- makeContrasts("cancer-Normal", levels=design2)
```
Fitting for the second time.

```
fit2c <- contrasts.fit(fit02, contrast.matrix)
fit2c <- eBayes(fit2c)
GOtable <- topTable(fit2c, number = Inf)
back <- subset(GOtable, adj.P.Val < 0.1)$PROBEID
backidx <- genefilter::genefinder(prep.gse, as.character(back), method = "manhattan", scale = "none")
backidx <- sapply(backidx, function(x)x$indices)
backs <- setdiff(backs, back)
intersect(backs, back)
```
Now the data is prepared to run a topGO analysis enrichment.
```
IDs <- rownames(results)
un <- IDs %in% c(backs, back)
sl <- IDs %in% backs
genes <- sl[un]
genes <- factor(as.integer(sl[un]))
names(genes) <- IDs[un]
 GO <- new("topGOdata", ontology="BP", allGenes = genes, nodeSize = 10, annot = annFUN.db, affyLib="hgu133a.db")
 
  result_top_GO_elim <- 
   runTest(GO, algorithm = "elim", statistic = "Fisher")
   
    result_top_GO_classic <- 
     runTest(GO, algorithm = "classic", statistic = "Fisher")
     
  
GO <- GenTable(GO, Fisher.elim = GO,
        Fisher.classic = GO,
        orderBy = "Fisher.elim" , topNodes = 100)
        
top_GO <- printGenes(GO, whichTerms = GO$GO.ID,
                          chip = "hgu133a.db", geneCutOff = 1000)
                          
write_csv(res_top_GO, path ="res_top_GO.csv")
```
## Enrichment analysis via reactomePA
```
entrez <- mapIds(hgu133a.db, 
                 keys = rownames(results),
                 keytype = "PROBEID",
                 colum = "ENTREZID")
reactome <- enrichPathway(gene = entrez[back], 
                                universe = entrez[c(back, 
                                                        backs)],
                                organism = "human",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.9, 
                                readable = TRUE)
        
        
reactome @result$Description <- paste0(str_sub(
                                    reactome @result$Description, 1, 20),
                                    "...")
                                    
head(as.data.frame(reactome_enrich))[1:6]

barplot(reactome_enrich)
````

* possible error messages encountered if using this code and how to solve them. 
Error: The size of the connection buffer (196000) was not large enough
Solve by running: Sys.setenv("VROOM_CONNECTION_SIZE" = 520000 * 2)

 Error: cannot allocate vector of size 1.6 GB
Solve by running: memory.limit(size= 20000) < or whatever memory vector you are using.
 

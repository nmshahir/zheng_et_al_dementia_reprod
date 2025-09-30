################################################################
# Install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("limma")
install.packages("umap")
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

################################################################
# Dataset: GSE122063

ctrl_v_VaD_list <- getGEO("GSE122063") # Importsfile as an ExpressionSet 
ctrl_v_VaD_eData <- ctrl_v_VaD_list[[1]]


###################################
# GEO2R GENERATED CODE
#################################
gset <- getGEO("GSE122063", GSEMatrix =TRUE, AnnotGPL=FALSE)

if (length(gset) > 1) idx <- grep("GPL16699", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))


# group membership for all samples
gsms <- paste0("000000000000000000000000000000000000XXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11111111",
        "111111111111111111111111111111111111")

sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID","GB_ACC","SEQUENCE"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)


# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

#this is the list of patient diagnostics
gset$`patient diagnosis:ch1`

# INSERT COLUMN NAMES 
# SPLIT FIRST COLUMN SO THAT ONLY SAMPLE NAME REMAINS
# REMOVE DUPLICATES TO GET FINAL NUMBERS
EnhancedVolcano(tT2,
  lab = NA,
  x = 'logFC',
  y = 'adj.P.Val')

EnhancedVolcano(tT2,
  x = 'logFC',
  y = 'adj.P.Val')


sig_tT2 <- tT2  |>
  filter(adj.P.Val < 0.05) |>
  filter(abs(logFC) > 1)
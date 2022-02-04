# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)


# load series and platform data from GEO

gset <- getGEO("GSE30639", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL3154", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000001111112222233333344444455555566666"
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) # log2 transform

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("17kL","17kS","40kL","40kS","6kL","6kS","rel606"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c("X6kS-rel606","X6kL-rel606","X17kL-X6kL","X17kS-X6kS","X40kL-X6kL","X40kS-X6kS","X6kS-X6kL","X17kS-rel606","X17kL-rel606","X40kS-rel606","X40kL-rel606")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


write.csv(fit2,'fit2.csv')

cont.matrix <- makeContrasts(contrasts=groups, levels=design)
fit4 <- contrasts.fit(fit, cont.matrix)
write.csv(fit4,'fit4.csv')

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE30639", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE30639", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")



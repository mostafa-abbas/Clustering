#required libraries 
library(EDASeq)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(genefilter)
library(ConsensusClusterPlus)
#read the gene expression matrix
RNASeq.DATA_ <- read.table("Samples_in_one_matrix_cap.tsv", header=TRUE, row.names=1)
#read geneInfo matrix
geneInfo <- read.csv("geneInfo_bimomart_.csv", header=TRUE, row.names=1)
#To be sure the matching between  gene names, we convert all genes to upper case
rownames(RNASeq.DATA_) = toupper(rownames(RNASeq.DATA_))
rownames(geneInfo) = toupper(rownames(geneInfo))
#drop the genes without Information
RNASeq.DATA <- RNASeq.DATA_[rownames(geneInfo),]
# select the genes with very low expression
filter <- apply(RNASeq.DATA,1,function(x) mean(x)>7)
considered_names = rownames(RNASeq.DATA[filter,])
#drop the genes with very low expression
RNASeq.DATA_1 = RNASeq.DATA[considered_names,]
geneInfo_1 = geneInfo [considered_names,]
#convert RNASeq.DATA to matrix
RNASeq.DATA_2 <- data.matrix(RNASeq.DATA_1)
RNASeq.expr.set <- newSeqExpressionSet(RNASeq.DATA_2, featureData = geneInfo_1) # create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "percentage_gc_content"] <- as.numeric(geneInfo_1[, "percentage_gc_content"]) # make sure gcContenet is numeric
fData(RNASeq.expr.set)[, "transcript_length"] <- as.numeric(geneInfo_1[, "transcript_length"]) # make sure gcContenet is numeric
RNASeq.expr.set_withinLanNorm <- withinLaneNormalization(RNASeq.expr.set, "percentage_gc_content", which = "upper", offset = TRUE) #removes lane gene specific effects, for example effects related to  GC content
RNASeq.expr.set_withinLanNorm_1 <- withinLaneNormalization(RNASeq.expr.set_withinLanNorm, "transcript_length", which = "upper", offset = TRUE) #removes lane gene specific effects, for example effects related to gene length 
RNASeq.expr.set_betwLanNorm <- betweenLaneNormalization(RNASeq.expr.set_withinLanNorm_1, which = "upper", offset = TRUE)             #removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.norm_1 = normCounts(RNASeq.expr.set_betwLanNorm) #normalized matrix
write.csv(RNASeq.norm_1, file = "RNA_norm_biomart.csv")
RNASeq.norm_filter = varFilter(RNASeq.norm_1,var.cutoff=0.75) #select the most variant genes
#extracting the heatmap

base_mean = rowMeans(RNASeq.norm_filter)
mat_scaled = t(apply(RNASeq.norm_filter, 1, scale))
type = gsub("s\\d+_", "", colnames(RNASeq.norm_filter))
ha = HeatmapAnnotation(df = data.frame(type = type))
#km means do k-means clustering on rows. If the value is larger than 1, he heatmap will be split by rows according to the k-means clustering. For each row-clusters, hierarchical clustering is applied with parameters above.
pdf("Heatmap.pdf")
Heatmap(mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
		top_annotation = ha, top_annotation_height = unit(4, "mm"), 
		show_row_names = FALSE, show_column_names = FALSE) +Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm")) 
dev.off()
#ConsensusClusterPlus
title="expression"
RNASeq.norm_filter_log = log(RNASeq.norm_filter+1.0) 
results = ConsensusClusterPlus(RNASeq.norm_filter_log,maxK=6,reps=50,pItem=0.8,pFeature=1,
							    title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

 
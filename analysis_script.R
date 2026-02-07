#Cleaning the data file
raw_data <- read.table("GSE53697_RNAseq_AD.txt", header=TRUE, check.names=FALSE, fill=TRUE)

#Keeping only the GeneID and the columns that end in "_raw"
count_matrix <- raw_data[, grepl("_raw", colnames(raw_data))]
rownames(count_matrix) <- raw_data$GeneID

#Loading the metadata
metadata <- read.csv("metadata.csv", row.names=1)
count_matrix <- count_matrix[, rownames(metadata)]

#Rounding the numbers (DESeq2)
count_matrix <- round(count_matrix)


dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "AD", "Control"))

table(res$pvalue < 0.05)

# Sorting by p-value
res_ordered_raw <- res[order(res$pvalue), ]

# Show the top 10 rows using only columns we KNOW exist
# (base DESeq2 columns are baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
head(res_ordered_raw[, c("baseMean", "log2FoldChange", "pvalue")], 10)

#Heatmap---
# Getting the IDs for the top 50 genes based on raw p-value
top50_ids <- rownames(res_ordered_raw)[1:50]

# Transforming the data so colors look balanced (not dominated by outliers)
vsd <- vst(dds, blind=FALSE)
plot_matrix <- assay(vsd)[top50_ids, ]

# Create the Heatmap
#installing the library
install.packages("pheatmap")
library(pheatmap)
#creating the actual plot
pheatmap(plot_matrix, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         show_colnames=TRUE,
         annotation_col=metadata["condition"], 
         main="Top 50 Potential AD Biomarkers",
         scale="row", # This centers the data so 'red' is high and 'blue' is low
         color=colorRampPalette(c("blue", "white", "red"))(100))

png("AD_Heatmap_Final.png", width=1000, height=800, res=120)
pheatmap(plot_matrix, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=metadata["condition"], 
         main="Top 50 Gene Expression: Alzheimer's vs Control",
         scale="row",
         color=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#Preparing the final results 
#Extracting only the genes with p-value less than 0.05
sig_genes_df <- as.data.frame(res[which(res$pvalue < 0.05), ])

#Sorting them so the most significant genes are at the top
sig_genes_df <- sig_genes_df[order(sig_genes_df$pvalue), ]

#Saving to the project folder
write.csv(sig_genes_df, "AD_Significant_Genes_484.csv")

# Invoke by 'Rscript geneCountMunge.R --args path_to_file path_to_destination'
# Writes the voomified output to path_to_destination as a space-seperated table.

library(limma)

filename <- commandArgs(trailingOnly = T)[2]
destination <- commandArgs(trailingOnly = T)[3]
gene_data <- read.table(filename, header=T)
gene_colnames <- colnames(gene_data)

sparse_gene <- function(gene_row) {
	counts <- as.numeric(gene_row[3:length(gene_row)])
	num_zero <- length(which(counts < 1))
	if (num_zero > ncol(gene_data) / 2) {
		return(F)
	} else {
		return(T)
	}
}

keep_row <- by(gene_data, 1:nrow(gene_data), sparse_gene)
gene_data <- gene_data[keep_row,] # remove sparse genes
voom_normalized_counts <- voom(gene_data[,3:ncol(gene_data)]) # voom normalization
gene_data <- cbind(gene_data[,c(1,2)], voom_normalized_counts$weights)
write(gene_colnames, destination, ncolumns = ncol(gene_data), sep="\t")
write.table(gene_data, destination, sep="\t",
	append=T, col.names=F, row.names=F, quote=F)

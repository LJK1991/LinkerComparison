library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(argparse)
library(biomaRt)


parser <- ArgumentParser()

parser$add_argument("-C", "--count_file", help="path to the count file, must be unzipped")
parser$add_argument("-O", "--output_dir", help="path to directory where to store the output")
parser$add_argument("-X", "--ten_X", action="store_true", help="Whether the count files are stored in a 10X format e.g. genes.tsv,matrix.mtx,barcodes.tsv")
parser$add_argument("-L", "--min_cutoff", default = 250, help="miniumum value to use as a cutoff, if lower than this it is ignored")
parser$add_argument("-H", "--max_cutoff", default = 2500, help="maximum value to use as a cutoff, if higher thant this it is ignored")
parser$add_argument("-T", "--cutoff_type", default = 'feature', help="Whether to use feature or counts as a cutoff")
parser$add_argument("-M", "--mito_cutoff", default = 5, help="maximum percentage of mitochondrial reads a cell is allowed to have")

args <- parser$parse_args()

file <- args$count_file
X <- args$ten_X
min <- as.integer(args$min_cutoff)
max <- as.integer(args$max_cutoff)
type <- args$cutoff_type
mito_cutoff <- as.integer(args$mito_cutoff)
output_dir <- args$output_dir


#open the count file

mart_mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://uswest.ensembl.org")

if (X == TRUE){
  counts <- Read10X(file, gene.column = 1)

  #loading the biomaRt
  genes <- data.frame(ensembl_gene_id = rownames(counts), external_gene_name = NA)
  genes_name <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), 
                                  filters = "ensembl_gene_id_version", 
                                  values = genes$ensembl_gene_id, 
                                  mart = mart_mm, 
                                  verbose = FALSE)

  for (id in 1:length(genes$ensembl_gene_id)){
    if ((genes$ensembl_gene_id[id] %in% genes_name$ensembl_gene_id) == TRUE){
      genes$external_gene_name[id] <- genes_name$external_gene_name[which(genes_name$ensembl_gene_id == genes$ensembl_gene_id[id])]
    } 
    else {
      genes$external_gene_name[id] <- genes$ensembl_gene_id[id]
    }
  }
  for (i in 1:nrow(genes)){ if (genes$external_gene_name[i] == ''){genes$external_gene_name[i] <- genes$ensembl_gene_id[i]}}
  
  genes$external_gene_name <- make.names(genes$external_gene_name, unique = TRUE)

  print("This is how many genes are still unnamed")
  print(sum(('' %in% genes$external_gene_name)))
  print("Are all values unique")
  print(length(unique(genes$external_gene_name)) == nrow(genes))

  rownames(counts) <- genes$external_gene_name

  CDS <- CreateSeuratObject(counts, assay = "RNA", min.cells = 3, min.features = 50)

} else{
  counts <- read.table(file=file, sep = "\t", header = TRUE, row.names = "gene")

  #loading the biomaRt

  genes <- data.frame(ensembl_gene_id = rownames(counts), external_gene_name = NA)
  genes_name <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), 
                                  filters = "ensembl_gene_id_version", 
                                  values = genes$ensembl_gene_id, 
                                  mart = mart_mm, 
                                  verbose = FALSE)

  for (id in 1:length(genes$ensembl_gene_id)){
    if ((genes$ensembl_gene_id[id] %in% genes_name$ensembl_gene_id) == TRUE){
      genes$external_gene_name[id] <- genes_name$external_gene_name[which(genes_name$ensembl_gene_id == genes$ensembl_gene_id[id])]
    } 
    else {
      genes$external_gene_name[id] <- genes$ensembl_gene_id[id]
    }
  }

  for (i in 1:nrow(genes)){ if (genes$external_gene_name[i] == ''){genes$external_gene_name[i] <- genes$ensembl_gene_id[i]}}

  genes$external_gene_name <- make.names(genes$external_gene_name, unique = TRUE)

  print("This is how many genes are still unnamed")
  print(sum(('' %in% genes$external_gene_name)))
  print("Are all values unique")
  print(length(unique(genes$external_gene_name)) == nrow(genes))

  rownames(counts) <- genes$external_gene_name

  CDS <- CreateSeuratObject(counts, assay = "RNA", min.cells = 3, min.features = 50)
}

if (colnames(CDS)[1] == 'X.') {
  CDS <- CDS[,-unlist(which(colnames(CDS) == 'X.'))]
} else if (colnames(CDS)[1] == '-'){
  CDS <- CDS[,-unlist(which(colnames(CDS) == '-'))]
}


sprintf("Dimension of CDS: %s %s",dim(CDS)[1], dim(CDS)[2])
print(typeof(CDS$nCount_RNA))
CDS$nCount_RNA <- as.numeric(CDS$nCount_RNA)

print("MT")
#storing %mitcochondrial reads per cell in the metaData
CDS[["percent.mt"]] <- PercentageFeatureSet(CDS, pattern = "^mt.")

print("Violin")
png(paste(output_dir, "_QC_violin.png",sep = ""), width = 960)
print(VlnPlot(CDS, features = c("nFeature_RNA","nCount_RNA", "percent.mt")))
dev.off()

print(min)
print(max)
print(mito_cutoff)

print("subset")

if (type == 'feature'){
  CDS <- subset(CDS, subset = nFeature_RNA > min & nFeature_RNA < max & percent.mt < mito_cutoff)
} else if (type == 'counts'){
  CDS <- subset(CDS, subset = nCount_RNA > min & nCount_RNA < max & percent.mt < mito_cutoff)
}

print("Violin 2")
png(paste(output_dir, "_postQC_violin.png",sep = ""), width = 960)
print(VlnPlot(CDS, features = c("nFeature_RNA","nCount_RNA", "percent.mt")))
dev.off()

#checking some stats
sprintf("Detected a total of %s genes and %s cells.", dim(CDS)[1], dim(CDS)[2])
sprintf("Mean of genes per cell: %s", mean(CDS$nFeature_RNA))
sprintf("Mean of UMI's per cell: %s", mean(CDS$nCount_RNA))
sprintf("Median of genes per cell: %s", median(CDS$nFeature_RNA))
sprintf("Median of UMI's per cell: %s", median(CDS$nCount_RNA))

print("get barcodes")
cells <- colnames(CDS)
writeLines(cells, paste(output_dir, "_cell_barcodes.txt",sep =""))

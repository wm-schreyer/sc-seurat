library(Seurat)
library(ggplot2)

# create pos48hr object
pos48.data <- Read10X(data.dir = "pos_tsv")
pos48.object <- CreateSeuratObject(counts = pos48.data, project = "Pos_48hr")
pos48.object

# create pos72hr object with raw cellranger matrix
pos72.data <- Read10X(data.dir = "pos72_tsv")
pos72.object <- CreateSeuratObject(counts = pos72.data, project = "Pos_72hr")
pos72.object

# merge pos48hr and pos72hr data
pos.times.object <- merge(pos48.object, y = pos72.object, add.cell.ids = c("Pos_48hr", "Pos_72hr"), project = "pos.times")
pos.times.object

# calculate percent MT and RP transcripts, add to object metadata
pos.times.object[["percent.mt"]] <- PercentageFeatureSet(pos.times.object, pattern = "^MT-")
pos.times.object[["percent.rp"]] <- PercentageFeatureSet(pos.times.object, pattern = "^RP")

# visualize QC as vlnplot
VlnPlot(object = pos.times.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0.5)

# filter low quality cells
pos.times.filt <- subset(x = pos.times.object, 
                   subset = (nFeature_RNA >= 100) &
                     (percent.mt < 30))

# comparative QC vlnplot between pos timepoint datasets
VlnPlot(pos.times.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.5) + labs(caption = "30% MT Cutoff")

# filter low quality cells HIGHER THRESH
pos.times.filt2 <- subset(x = pos.times.object, 
                         subset = (nFeature_RNA >= 100) &
                           (percent.mt < 60))

# comparative QC vlnplot between pos timepoint datasets HIGHER THRESH
VlnPlot(pos.times.filt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.5) + labs(caption = "60% MT Cutoff")

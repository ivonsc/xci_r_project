library(rtracklayer)
library(GenomicRanges)
library(GEOquery)
library(dplyr)
library(conflicted)
library(tidyverse)
library(tidyr)
library(tibble)
library(biomaRt)
library(ggplot2)
library(GenomeInfoDb)
conflicts_prefer(dplyr::filter)
#leer metadatos


gse <- getGEO("GSE186458", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))
colnames(meta)

###################################
######## RAW DATA ##################
#################################

files <- list.files("raw_data", pattern = "tsv$", full.names = TRUE)
gsm_ids <- gsub(".*(GSM[0-9]+).*", "\\1", files)


#organise columns

split_index <- which(rownames(meta) == "GSM6810003")
df1 <- meta[1:(split_index - 1), ]       
df2 <- meta[split_index:nrow(meta), ]

#we have two datasets, now rename columns accordingly
df1$characteristics_ch1 <- gsub("tissue: ", "", df1$characteristics_ch1)

meta_files1 <- df1 %>%
  dplyr::mutate(cell_type = paste0(characteristics_ch1.1, "_", characteristics_ch1)) %>%
  dplyr::select(
    geo_accession,
    cell_type,
    age = characteristics_ch1.3,
    sex = characteristics_ch1.4,
    race = characteristics_ch1.5
  )


meta_files2 <- df2 %>%
  dplyr::select(
    geo_accession,
    cell_type = characteristics_ch1,
    sex = characteristics_ch1.2,
    age = characteristics_ch1.3,
    race = characteristics_ch1.6
  )


# join again

metadata <- bind_rows(meta_files1, meta_files2)

# see NANs and missing data

sum(is.na(metadata))
unique(metadata$sex)

# normalize sex variable

metadata$sex <- ifelse(metadata$sex %in% c("Sex: F", "gender: female"), "Female",
                       ifelse(metadata$sex %in% c("Sex: M", "gender: male"), "Male", NA ))

unique(metadata$sex)
metadata <- na.omit(metadata)

metadata$cell_type <- gsub("cell type: ", "", metadata$cell_type)
metadata$age <- gsub("age: ", "", metadata$age)
metadata$race <- gsub("race: ", "", metadata$race)

# DATASET WITH PATIENT ID

metadata_patient_id <- df2 %>%
  dplyr::select(
    geo_accession,
    cell_type = characteristics_ch1,
    patient_id = characteristics_ch1.1,
    sex = characteristics_ch1.2,
    age = characteristics_ch1.3,
    race = characteristics_ch1.6
  )

metadata_patient_id$sex <- ifelse(metadata_patient_id$sex == "gender: female", "Female",
                       ifelse(metadata_patient_id$sex == "gender: male", "Male", NA ))

metadata_patient_id$cell_type <- gsub("cell type: ", "", metadata_patient_id$cell_type)
metadata_patient_id$age <- gsub("age: ", "", metadata_patient_id$age)
metadata_patient_id$race <- gsub("race: ", "", metadata_patient_id$race)
metadata_patient_id$patient_id <- gsub("patient id: ", "", metadata_patient_id$patient_id)

# SUBSETTING WOMEN

metadata_women <- metadata_patient_id[metadata_patient_id$sex == "Female",]
metadata_women_no_patient_id <- metadata[metadata$sex == "Female",]

#GROUPING BY INDIVIDUAL

metadata_women %>%
  group_by(patient_id) %>%
  summarise(n_muestras = n())

# Seems like each patient has two samples (sample id = geo accession)

# Group by cell type

metadata_women %>%
  group_by(cell_type) %>%
  summarise(n_muestras = n())

#METADATA NOT USING PATIENT ID

metadata_women_no_patient_id %>%
  group_by(cell_type) %>%
  summarise(n_muestras = n()) %>%
            print(n = 37)

# EXTRACT GEO ACCESSION INTO A FILE FOR DOWNLOADING SAMPLES
sample_list <- metadata_women_no_patient_id$geo_accession

for ( sample in sample_list) {
  write(sample, file = "output/links.txt", sep = "", append = TRUE)
  
}

# REPRODUCING FIGURE Density - CpG Metilation
listEnsembl()
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# GRCh38.p14 version
filters = listFilters(ensembl)

# STEPS
# for every sample, take chr X, only women
# map every gene of the list to the bigWig coordinates

# BLUE = escapees RED = XCI genes
# have a table stating the genes that are escapees XCI genes
# COLUMNS = SAMPLES, ROWS = BETA-VALUES of EACH GENE


# cols are sample_list
# rows are beta value



# Cargar lista de genes escapees vs XCI
genes <- read.table("raw_data/GenesEscapeXCI.txt", sep="\t", skip = 1, header = TRUE)
genes <- genes[,colnames(genes)[1:14]]
genes <- na.omit(genes)

gr_genes <- GRanges(seqnames = genes$Chr,
                    ranges = IRanges(start = genes$Start.position, end = genes$End.position),
                    gene = genes$Gene.name,
                    status = genes$Combined.XCI.status)


files <- list.files("processed_data", pattern= "tsv$", full.names = TRUE)
chrX_range <- GRanges(seqnames = "chrX", ranges = IRanges(1, 156040895))

bw_list <- lapply(files, function(file) {
  df <- read.table(file, header = FALSE)
  colnames(df) <- c("chr", "start", "end", "meth", "unmeth", "beta_value")
  
  df$beta_value <- (df$meth ) / (df$meth + df$total)
  df[is.na(df)] <- 0
  
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end),
                meth = df$meth,
                unmeth = df$total,
                beta = df$beta_value)
  
  # Filtrar por chrX_range si lo tienes definido
  gr_chrX <- subsetByOverlaps(gr, chrX_range)
  return(gr)
})


names(bw_list) <- gsub("\\.hg38\\tsv$", "", basename(files))

seqlevels(gr_genes)
sample_bw <- bw_list[[1]]
seqlevels(sample_bw)


seqlevelsStyle(gr_genes) <- "UCSC"

# Para cada muestra, calcular beta promedio por gen
overlaps <- findOverlaps(gr_genes, sample_bw)
gene_hits <- split(subjectHits(overlaps), queryHits(overlaps))


# Calcula el promedio de beta por gen
beta_per_gene <- sapply(seq_along(gene_hits), function(i) {
  mean(sample_bw$beta[gene_hits[[i]]], na.rm = TRUE)
})

# Añade los resultados al objeto gr_genes
gr_genes$beta_mean <- NA
gr_genes$beta_mean[as.integer(names(gene_hits))] <- beta_per_gene

# VARIAS MUESTRAS
# Crear matriz: genes en filas, muestras en columnas
# Inicializa matriz con NA
beta_matrix <- matrix(NA, nrow = length(gr_genes), ncol = length(bw_list))
rownames(beta_matrix) <- gr_genes$gene
colnames(beta_matrix) <- names(bw_list)

for (j in seq_along(bw_list)) {
  sample_bw <- bw_list[[j]]
  overlaps <- findOverlaps(gr_genes, sample_bw)
  gene_hits <- split(subjectHits(overlaps), queryHits(overlaps))
  
  for (i in seq_along(gene_hits)) {
    gene_index <- as.integer(names(gene_hits))[i]
    beta_matrix[gene_index, j] <- mean(sample_bw$beta[gene_hits[[i]]], na.rm = TRUE)
  }
}

rownames(beta_matrix) <- gr_genes$gene
df <- data.frame(beta_matrix)
df$status <- gr_genes$status[match(rownames(df), gr_genes$gene)]

df_filtrado <- df[df$status != "variable", ]
df_long <- reshape2::melt(df_filtrado, id.vars = "status")

###############################
#### ONE SAMPLE ###############
###############################

ggplot(df_long, aes(x = value, fill = status)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribución de metilación por gen (chrX)",
       x = "Beta", y = "Densidad") +
  scale_fill_manual(values = c("escape" = "blue", "inactive" = "red")) +
  xlim(0, 1)

##########################################
########## BY CELL TYPE ##################
##########################################


df_beta <- as.data.frame(t(beta_matrix)) %>%
  rownames_to_column(var = "geo_accession") %>%
  pivot_longer(cols = -geo_accession, names_to = "gene", values_to = "beta")


gene_info <- data.frame(gene = gr_genes$gene, status = gr_genes$status )

# Join metadata

df_beta <- df_beta %>%
  left_join(gene_info, by = "gene")

df_beta$geo_accession <- sub("_.*", "", df_beta$geo_accession)
df_beta <- df_beta[!is.na(df_beta$beta), ]

metadata_women_no_patient_id_no_race <- metadata_women_no_patient_id[, !names(metadata_women_no_patient_id) %in% "race"]


df_beta <- df_beta %>%
  left_join(metadata_women_no_patient_id_no_race, by = "geo_accession")

df_beta <- df_beta[!is.na(df_beta$cell_type), ]

# Filter variable status and plot by cell_type

df_beta_filtrado <- df_beta %>%
  filter(!is.na(beta)) %>%
  filter(status %in% c("escape", "inactive")) %>%
  filter(!is.na(beta) & beta >= 0 & beta <= 1)

ggplot(df_beta_filtrado, aes(x = beta, fill = status)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ cell_type) +
  scale_fill_manual(values = c("escape" = "blue", "inactive" = "red")) +
  labs(title = "Distribución de metilación por tipo de gen y tipo celular",
       x = "Beta", y = "Densidad") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 13),
    strip.background = element_rect(fill = "#edf0f7", color = "black", linewidth = 0.5)
  ) +
  xlim(0, 1)

#METADATA NOT USING PATIENT ID

metadata_women_no_patient_id %>%
  group_by(cell_type) %>%
  summarise(n_muestras = n()) %>%
  print(n = 61)


# There are 61 different tissue types


# Free DNA blood samples
df_filtrado <- df[df$status != "variable", ]
df_long <- reshape2::melt(df_filtrado, id.vars = "status")

df_long$geo_accession <- sub("_.*", "", df_long$variable)
df_long <- df_long[!is.na(df_long$value), ]

df_long <- df_long %>%
  left_join(metadata_women_no_patient_id, by = "geo_accession")


# cada muestra

ggplot(df_long, aes(x = value, fill = status)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable) +
  labs(title = "Distribución de metilación por gen (chrX)",
       x = "Beta", y = "Densidad") +
  scale_fill_manual(values = c("escape" = "blue", "inactive" = "red")) +
  xlim(0, 1)



####################################
############ ALL GENOME ############
####################################

library(rtracklayer)
library(GenomicRanges)

# Lista de archivos .tsv con datos de metilación
files <- list.files("processed_data", pattern = "tsv$", full.names = TRUE)

# Crear lista de GRanges para todas las muestras
bw_list_all_chr <- lapply(files, function(file) {
  df <- read.table(file, header = FALSE)
  colnames(df) <- c("chr", "start", "end", "meth", "unmeth", "beta_value")
  
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end),
                meth = df$meth,
                unmeth = df$unmeth,
                beta = df$beta_value)
  
  return(gr)  # No se filtra por cromosoma
})

# Asignar nombres a cada muestra usando el nombre del archivo
names(bw_list_all_chr) <- gsub("\\.hg38\\tsv$", "", basename(files))


genes_all <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
                   filters = "chromosome_name",
                   values = c(1:22, "X", "Y", "MT"),
                   mart = ensembl)

gr_all_genes <- GRanges(seqnames = paste0("chr", genes_all$chromosome_name),
                        ranges = IRanges(start = genes_all$start_position, end = genes_all$end_position),
                        gene = genes_all$external_gene_name)

beta_matrix_all <- matrix(NA, nrow = length(gr_all_genes), ncol = length(bw_list_all_chr))
rownames(beta_matrix_all) <- gr_all_genes$gene
colnames(beta_matrix_all) <- names(bw_list_all_chr)

for (j in seq_along(bw_list_all_chr)) {
  sample_bw <- bw_list_all_chr[[j]]
  overlaps <- findOverlaps(gr_all_genes, sample_bw)
  gene_hits <- split(subjectHits(overlaps), queryHits(overlaps))
  
  for (i in seq_along(gene_hits)) {
    gene_index <- as.integer(names(gene_hits))[i]
    beta_matrix_all[gene_index, j] <- mean(sample_bw$beta[gene_hits[[i]]], na.rm = TRUE)
  }
}

df_all <- data.frame(beta_matrix_all)
df_long_all <- reshape2::melt(df_all)


ggplot(df_long_all, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Distribución de metilación por gen (todo el genoma)",
       x = "Beta", y = "Densidad") +
  xlim(0, 1) +
  theme_minimal()


#########################################################
# Manifest for Illumina's EPIC v2.0 methylation arrays ##
#########################################################
# BiocManager::install("IlluminaHumanMethylationEPICv2manifest")
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)

library(minfi)

# Obtener el objeto de anotación
annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Extraer las posiciones CpG como GRanges
# Convertir a GRanges
epic_gr <- makeGRangesFromDataFrame(annotation,
                                    seqnames.field = "chr",
                                    start.field = "pos",
                                    end.field = "pos",
                                    keep.extra.columns = TRUE)
# Supongamos que tienes tus datos WGBS como un objeto GRanges llamado gr_wgbs:
gr_epic_overlap <- subsetByOverlaps(gr_genes, epic_gr)


par(mar = c(4, 4, 2, 1))  # márgenes: abajo, izquierda, arriba, derecha
hist(gr_epic_overlap$beta_mean, breaks = 50,
     main = "Distribución de beta values en CpGs del array EPIC v2",
     xlab = "Beta value")


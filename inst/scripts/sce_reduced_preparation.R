library(SingleCellExperiment)
library(scuttle)
library(scran)
library(tidyverse)
library(SingleCellExperiment)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
set.seed(1)

# Download for max. 10 minutes
options(timeout=600)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE148842&format=file", "GSE148842_RAW.tar")
untar("GSE148842_RAW.tar")


file_df <- tibble(file = list.files(path = getwd(), full.names = TRUE, pattern= "^GSM"),
                  id = str_match(file, ".+/GSM\\d{7}_(PW[0-9-]*)\\.cts\\.txt\\.gz")[,2])
geo_metadata <- GEOquery::getGEO("GSE148842")
# Parse geo metadata
patient_annotation <- bind_rows(as_tibble(geo_metadata[[1]]),
                                as_tibble(geo_metadata[[2]])) %>%
  unite("description", starts_with("characteristic"), sep = "\n") %>%
  transmute(id = title, description, origin = source_name_ch1) %>%
  separate(id, into = c("patient_id", "treatment_id"), remove = FALSE) %>%
  mutate(age = as.numeric(str_match(description, "age: (\\d{1,3})\\s")[,2]),
         gender = str_match(description, "gender: ([MF])\\s")[,2],
         location = str_match(description, "location:\\s*([^\\n]+)\\s*\\n")[,2],
         diagnosis = str_match(description, "diagnosis:\\s*([^\\n]+)\\s*\\n")[,2],
         treatment = str_match(description, "treatment:\\s*([^\\n]+)\\s*$")[,2]) %>%
  dplyr::select(-description) %>%
  mutate(patient_id = ifelse(str_starts(id, "PW05"), str_sub(id, start = 1L, end = -4L), patient_id),
         treatment_id = ifelse(str_starts(id, "PW05"), str_sub(id, start = -3L, end = -1L), treatment_id)) %>%
  mutate(origin = ifelse(origin == "glioma surgical biopsy", "biopsy", "tissue_slice"),
         recurrent_tumor = str_ends(diagnosis, "recurrent")) %>%
  mutate(condition = case_when(
    treatment == "vehicle (DMSO)" ~ "ctrl",
    treatment == "2.5 uM etoposide" ~ "etoposide",
    treatment == "0.2 uM panobinostat" ~ "panobinostat",
    TRUE ~ "other",
  ))

set.seed(1)
tmp <- patient_annotation %>%
  tidylog::inner_join(file_df) %>%
  filter(condition != "other") %>%
  filter(! recurrent_tumor) %>%
  filter(origin == "tissue_slice")

sces <- map(tmp$file, \(fi){
  id <- str_match(fi, ".+/GSM\\d{7}_(PW[0-9-]*)\\.cts\\.txt\\.gz")[2]
  df <- data.table::fread(file = fi, showProgress = FALSE)
  rowdata <- df[,c("gid", "gene")]
  counts <- df %>%
    dplyr::select(- c(gid, gene)) %>%
    as.matrix() %>%
    as("dgCMatrix")
  rownames(counts) <- rowdata$gid
  rownames(rowdata) <- rowdata$gid
  sce <- SingleCellExperiment(list(counts = counts), colData = data.frame(id = rep(id, ncol(counts))), rowData = as.data.frame(rowdata))
  sce
}, .progress = TRUE)


sce <- do.call(cbind, sces)
colnames <- make.unique(colnames(sce))
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  tidylog::left_join(tmp, by = "id") %>%
  dplyr::select(- c(file, origin, recurrent_tumor)) %>%
  DataFrame()

sce$pat_cond <- paste0(sce$patient_id, "-", sce$condition)


logcounts(sce) <- transformGamPoi::shifted_log_transform(sce)
colnames(sce) <- colnames

# Remove useless genes
sce <- sce[rowSums(counts(sce)) >= 5, ]
# Remove version number of ENSEMBL gene id's
rowData(sce)$gid <- stringr::str_remove(rowData(sce)$gid, "\\.\\d+")
rownames(sce) <- rowData(sce)$gid

ah <- AnnotationHub::AnnotationHub()
annotation <- ah[["AH60085"]] %>%
  as_tibble() %>%
  filter(type == "gene") %>%
  group_by(gene_id) %>%
  summarize(chromosome = unique(seqnames), gene_length = median(width), strand. = strand, source)

table(rowData(sce)$gid %in% annotation$gene_id)

rowData(sce) <- rowData(sce) %>%
  as.data.frame() %>%
  mutate(gene_id = str_remove(gid, "\\.\\d+")) %>%
  left_join(annotation, by = c("gene_id"))

qc_df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = ! is.na(rowData(sce)$chromosome) & rowData(sce)$chromosome == "MT",
                                                       Y_chr = ! is.na(rowData(sce)$chromosome) & rowData(sce)$chromosome == "Y"))

qc_filters <- scuttle::perCellQCFilters(qc_df, sub.fields = c("subsets_Mito_percent"))
rownames(qc_filters) <- rownames(qc_df)

sce <- sce[,!qc_filters$discard]
sce <- sce[,colSums(counts(sce)) > 800 & colSums(counts(sce)) < 12000]

sel_chromosomes <- as.character(1:5)
special_chr <- c("7", "10")
chr_counts <- t(lemur:::aggregate_matrix(t(counts(sce)), group_split = deframe(vctrs::vec_group_loc(rowData(sce)$chromosome))[unique(c(sel_chromosomes, special_chr))], rowMeans2))

ratio7 <- chr_counts["7", ] / colMeans2(chr_counts, rows = 1:5)
ratio10 <- chr_counts["10", ] / colMeans2(chr_counts, rows = 1:5)

simple_label <- ifelse(ratio10 < ratio7, "tumor", "microenvironment")

means <- rowMeans(logcounts(sce)[which(rowData(sce)$chromosome == "7"),])
top_genes <- order(-means)[1:100]
ratio2 <- t(t(logcounts(sce)[which(rowData(sce)$chromosome == "7")[top_genes],]) / means[top_genes])

chr_feature <- cbind("chr10" = colSums(logcounts(sce)[which(rowData(sce)$chromosome == "10"),]),
                     "chr7" = colSums(logcounts(sce)[which(rowData(sce)$chromosome == "7"),]))
chr_clusters <- kmeans(chr_feature, centers = 2)

chr_clusters <- ClusterR::GMM(chr_feature, gaussian_comps = 2)
stopifnot((chr_clusters$centroids[1,1] > chr_clusters$centroids[2,1]) != (chr_clusters$centroids[1,2] > chr_clusters$centroids[2,2]))
cl_1_is_tumor <- chr_clusters$centroids[1,1] < chr_clusters$centroids[2,1]
tumor_label <- case_when(
  colSums2(counts(sce)) < 400 ~ "uncertain",
  predict(chr_clusters, newdata = chr_feature) == 1 & cl_1_is_tumor ~ "tumor",
  predict(chr_clusters, newdata = chr_feature) == 1 & ! cl_1_is_tumor ~ "microenvironment",
  predict(chr_clusters, newdata = chr_feature) == 2 & cl_1_is_tumor ~ "microenvironment",
  predict(chr_clusters, newdata = chr_feature) == 2 & ! cl_1_is_tumor ~ "tumor",
)

set.seed(1)
hvg <- order(-rowVars(logcounts(sce)))
set.seed(1)
subset_pca <- lemur:::pca(logcounts(sce)[hvg[1:2000], ], n = 35)
harm <- harmony::RunHarmony(t(subset_pca$embedding), meta_data = colData(sce), vars_use = c("pat_cond", "treatment_id"), lambda = c(1,1))
set.seed(1)
graph <- bluster::makeKNNGraph(harm, k = 15, BNPARAM = BiocNeighbors::AnnoyParam())
# Warning: this step needs more than 20GB of memory
clustering <- igraph::cluster_walktrap(graph)
set.seed(1)
clusters <- igraph::cut_at(clustering, 4)

cell_type_label <- case_when(
  clusters == 1 ~ "Myeloid cells",
  clusters == 2 ~ "Tumor cells",
  clusters == 3 ~ "Oligodendrocytes",
  clusters == 4 ~ "T cells"
)

colData(sce) <- colData(sce) %>%
  as_tibble() %>%
  mutate(cell_type = cell_type_label, chr_ratio_label = tumor_label,
         chr_10_ratio = ratio10, chr_7_ratio = ratio7) %>%
  DataFrame()

colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
hvg <- order(-rowVars(logcounts(sce)))
sce <- sce[hvg[1:10000],]

sce <- sce[, sce$condition != "etoposide" & sce$patient_id != "PW029"]

# Reduce to two patients only.
sce_2_pat <- sce[,sce$patient_id %in%c("PW030","PW034")]

# Reduce to three cell-types only.
table(sce_2_pat$patient_id,sce_2_pat$cell_type) |> knitr::kable()
sce_2_pat_3ct <-sce_2_pat[,sce_2_pat$cell_type != "T cells"]

# Keep only logcounts assay and needed information for colData.
colData(sce_2_pat_3ct) <- colData(sce_2_pat_3ct)[, c(
  "cell_type",
  "condition",
  "patient_id"
), drop = FALSE]


# Remove useless rowData and reducedDimNames we will anyway compute later.
rowData(sce_2_pat_3ct) <- rowData(sce_2_pat_3ct)[, c(
  "gid",
  "gene"
), drop = FALSE]

reducedDim(sce_2_pat_3ct, "PCA") <- NULL
reducedDim(sce_2_pat_3ct, "UMAP") <- NULL
reducedDim(sce_2_pat_3ct, "harmony_umap") <- NULL
reducedDim(sce_2_pat_3ct, "harmony") <- NULL


# Reduce to top 97 HVGs and 3 genes of interest.
set.seed(1)
pan_affect_genes <- rowData(sce_2_pat_3ct) |>
  as_tibble() |>
  dplyr::select(gid, gene) |>
  dplyr::filter(gene %in% c("HIST3H2A", "PTPRZ1", "HBEGF"))

random_genes <- rowData(sce_2_pat_3ct) |>
  as_tibble() |>
  dplyr::select(gid, gene) |>
  dplyr::filter(!gene %in% pan_affect_genes$gene) |>
  dplyr::slice_sample(n = 97)

selected_genes <- c(pan_affect_genes$gid, random_genes$gid)
sce_reduced <- sce_2_pat_3ct[selected_genes, ]

# Reduce to 3k cells.
set.seed(1)
sce_vignette <- sce_2_pat_3ct[, sample(ncol(sce_2_pat_3ct), 3000)]
saveRDS(sce_vignette, "sce_reduced.rds")

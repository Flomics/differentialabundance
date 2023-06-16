#!/usr/bin/env Rscript

library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

reference <- args[1]
target <- args[2]
samplesheet <- args[3]
counts <- args[4]

count_matrix <- read.table(counts, header = TRUE, sep = "\t", row.names = 1)

metadata <- read.table(samplesheet, header = TRUE, sep = "\t")
# count_matrix <- count_matrix[, names(count_matrix) %in% metadata$sample]
# metadata <- metadata[metadata$sample %in% names(count_matrix), ]

head(count_matrix)
head(metadata)

# group <- metadata$condition
# group <- factor(group)

# y <- DGEList(
#     counts = count_matrix,
#     genes = rownames(count_matrix),
#     group = metadata$condition)

# keep <- filterByExpr(y, y$samples$group)
# table(keep)
# y <- y[keep, , keep.lib.sizes = FALSE]

# y <- calcNormFactors(y)

# design <- model.matrix(~0 + y$samples$group)
# colnames(design) <- levels(y$samples$group)
# design

# y <- estimateDisp(y, design, robust = TRUE)

# fit <- glmQLFit(y, design, robust = TRUE)

# contrast_comp <- makeContrasts(target - reference, levels = design)
# res <- glmQLFTest(fit, contrast = contrast_comp)
# topTags(res)


# ################################################
# ################################################
# ## Generate outputs                           ##
# ################################################
# ################################################



# ################################################
# ################################################
# ## VERSIONS FILE                              ##
# ################################################
# ################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
edger.version <- as.character(packageVersion('edgeR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-edger:', edger.version)
    ),
'versions.yml')

# ################################################
# ################################################
# ################################################
# ################################################

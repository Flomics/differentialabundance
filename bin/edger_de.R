#!/usr/bin/env Rscript

library(edgeR)

args = commandArgs(trailingOnly=TRUE)

reference <- args[1]
target <- args[2]
samplesheet <- args[3]
counts <- args[4]

count_matrix <- read.table(counts, header = T, sep= ",", row.names = 1)

metadata <- read.table(samplesheet, header = T, sep= ",")
# count_matrix= count_matrix[,names(count_matrix)%in%metadata$Sample_name]
# metadata= metadata[metadata$Sample_name %in% names(count_matrix),]

head(count_matrix)
head(meatadata)

# group <- sample.sheet$condition
# group <- factor(group)

# y <- DGEList(
#     counts=count.table,
#     genes=rownames(count.table),
#     group = sample.sheet$condition)

# keep <- filterByExpr(y, y$samples$group)
# table(keep)
# y <- y[keep, , keep.lib.sizes=FALSE]

# y <- calcNormFactors(y)

# design <- model.matrix(~0+y$samples$group)
# colnames(design) <- levels(y$samples$group)
# design

# y <- estimateDisp(y, design, robust=TRUE)

# fit <- glmQLFit(y, design, robust=TRUE)

# contrast_comp <- makeContrasts(target_level-reference_level, levels=design)
# res <- glmQLFTest(fit, contrast=contrast_comp)
# topTags(res)


# if (opt\$control_genes_file != '' && opt\$sizefactors_from_controls){
#     print(paste('Estimating size factors using', length(control_genes), 'control genes'))
#     dds <- estimateSizeFactors(dds, controlGenes=rownames(count.table) %in% control_genes)
# }

# dds <- DESeq(
#     dds,
#     test = opt\$test,
#     fitType = opt\$fit_type,
#     minReplicatesForReplace = opt\$min_replicates_for_replace,
#     useT = opt\$use_t,
#     sfType = opt\$sf_type,
#     parallel=TRUE, BPPARAM=MulticoreParam(opt\$cores)
# )

# comp.results <-
#     results(
#         dds,
#         lfcThreshold = opt\$lfc_threshold,
#         altHypothesis = opt\$alt_hypothesis,
#         independentFiltering = opt\$independent_filtering,
#         alpha = opt\$alpha,
#         pAdjustMethod = opt\$p_adjust_method,
#         minmu = opt\$minmu,
#         contrast = c(
#             contrast_variable,
#             c(opt\$target_level, opt\$reference_level)
#         )
#     )

# if (opt\$shrink_lfc){
#     comp.results <- lfcShrink(dds,
#         type = 'ashr',
#         contrast = c(
#             contrast_variable,
#             c(opt\$target_level, opt\$reference_level)
#         )
#     )
# }

# ################################################
# ################################################
# ## Generate outputs                           ##
# ################################################
# ################################################

# prefix_part_names <- c('contrast_variable', 'reference_level', 'target_level', 'blocking_variables')
# prefix_parts <- unlist(lapply(prefix_part_names, function(x) gsub("[^[:alnum:]]", "_", opt[[x]])))
# output_prefix <- paste(prefix_parts[prefix_parts != ''], collapse = '-')

# contrast.name <-
#     paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
# cat("Saving results for ", contrast.name, " ...\n", sep = "")

# # Differential expression table- note very limited rounding for consistency of
# # results

# write.table(
#     data.frame(
#         gene_id = rownames(comp.results),
#         round_dataframe_columns(data.frame(comp.results, check.names = FALSE)),
#         check.names = FALSE
#     ),
#     file = paste(output_prefix, 'deseq2.results.tsv', sep = '.'),
#     col.names = TRUE,
#     row.names = FALSE,
#     sep = '\t',
#     quote = FALSE
# )

# # Dispersion plot

# png(
#     file = paste(output_prefix, 'deseq2.dispersion.png', sep = '.'),
#     width = 600,
#     height = 600
# )
# plotDispEsts(dds)
# dev.off()

# # R object for other processes to use

# saveRDS(dds, file = paste(output_prefix, 'dds.rld.rds', sep = '.'))

# # Size factors

# sf_df = data.frame(
#     sample = names(sizeFactors(dds)),
#     data.frame(sizeFactors(dds), check.names = FALSE),
#     check.names = FALSE
# )
# colnames(sf_df) <- c('sample', 'sizeFactor')
# write.table(
#     sf_df,
#     file = paste(output_prefix, 'deseq2.sizefactors.tsv', sep = '.'),
#     col.names = TRUE,
#     row.names = FALSE,
#     sep = '\t',
#     quote = FALSE
# )

# # Write specified matrices

# write.table(
#     data.frame(
#         gene_id=rownames(counts(dds)),
#         counts(dds, normalized = TRUE),
#         check.names = FALSE
#     ),
#     file = paste(output_prefix, 'normalised_counts.tsv', sep = '.'),
#     col.names = TRUE,
#     row.names = FALSE,
#     sep = '\t',
#     quote = FALSE
# )

# Write specified matrices

# write.table(
#     normalised_counts,
#     file = 'normalised_counts.tsv',
#     col.names = TRUE,
#     row.names = FALSE,
#     sep = '\t',
#     quote = FALSE
# )

# # Note very limited rounding for consistency of results

# for (vs_method_name in strsplit(opt\$vs_method, ',')){
#     if (vs_method_name == 'vst'){
#         vs_mat <- vst(dds, blind = opt\$vs_blind, nsub = opt\$vst_nsub)
#     }else if (vs_method_name == 'rlog'){
#         vs_mat <- rlog(dds, blind = opt\$vs_blind, fitType = opt\$fit_type)
#     }

#     # Again apply the slight rounding and then restore numeric

#     write.table(
#         data.frame(
#             gene_id=rownames(counts(dds)),
#             round_dataframe_columns(
#                 data.frame(assay(vs_mat), check.names = FALSE)
#             ),
#             check.names = FALSE
#         ),
#         file = paste(output_prefix, vs_method_name,'tsv', sep = '.'),
#         col.names = TRUE,
#         row.names = FALSE,
#         sep = '\t',
#         quote = FALSE
#     )
# }

# ################################################
# ################################################
# ## R SESSION INFO                             ##
# ################################################
# ################################################

# sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
# print(sessionInfo())
# sink()

# ################################################
# ################################################
# ## VERSIONS FILE                              ##
# ################################################
# ################################################

# r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
# edger.version <- as.character(packageVersion('edgeR'))

# writeLines(
#     c(
#         '"${task.process}":',
#         paste('    r-base:', r.version),
#         paste('    bioconductor-edger:', edger.version)
#     ),
# 'versions.yml')

# ################################################
# ################################################
# ################################################
# ################################################

process EDGER_DIFFERENTIAL {
    tag "$meta"
    label 'process_medium'

    conda "bioconda::bioconductor-edger=3.40.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-edger:3.40.0--r42hf17093f_1' :
        'quay.io/biocontainers/bioconductor-edger:3.40.0--r42hf17093f_1' }"

    input:
    path(edger)
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    // tuple val(meta), path("*.edger.results.tsv")              , emit: results
    // tuple val(meta), path("*.edger.dispersion.png")           , emit: dispersion_plot
    // tuple val(meta), path("*.dds.rld.rds")                     , emit: rdata
    // tuple val(meta), path("*.edger.sizefactors.tsv")          , emit: size_factors
    tuple val(meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    // tuple val(meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    // tuple val(meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    // tuple val(meta), path("*.R_sessionInfo.log")               , emit: session_info
    //path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp $edger edger.R
    Rscript edger.R $reference $target $samplesheet $counts
    """
}

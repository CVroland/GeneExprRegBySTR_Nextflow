process GET_FASTA_BEDTOOLS{
    conda params.condaEnv

    input:
    tuple val(strClass), val(moduleId), path(bedFile)
    path refFasta
    path refFastaFaidx

    output:
    tuple val(strClass), val(moduleId), path("${bedFile.baseName}.fasta")

    script:
    """
    bedtools getfasta -fi ${refFasta} -bed ${bedFile} > ${bedFile.baseName}.fasta
    """
}
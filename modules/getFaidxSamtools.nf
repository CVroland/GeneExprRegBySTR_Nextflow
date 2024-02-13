process GET_FAIDX_SAMTOOLS{
    conda params.condaEnv

    input :
    path(fastaFile)

    output:
    path("${fastaFile}.fai")

    script:
    """
    samtools faidx ${fastaFile}
    """
}
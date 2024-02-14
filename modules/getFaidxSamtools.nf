process GET_FAIDX_SAMTOOLS{

    input :
    path(fastaFile)

    output:
    path("${fastaFile}.fai")

    script:
    """
    samtools faidx ${fastaFile}
    """
}
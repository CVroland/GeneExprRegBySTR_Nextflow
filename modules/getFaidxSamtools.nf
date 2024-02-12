process GET_FAIDX_SAMTOOLS{
    //TODO : modify conda env and add to this process

    input :
    path(fastaFile)

    output:
    path("${fastaFile}.fai")

    script:
    """
    samtools faidx ${fastaFile}
    """
}
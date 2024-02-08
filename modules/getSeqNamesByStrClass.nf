process GET_SEQ_NAMES_BY_STR_CLASS{
    conda params.condaEnv

    input:
    val strClass
    path seqNamesFile

    output:
    tuple val(strClass), path("${strClass}_seqNames.npy")

    script:
    """
    filterSeqNames.py ${seqNamesFile} ${strClass} ${strClass}_seqNames.npy
    """
}
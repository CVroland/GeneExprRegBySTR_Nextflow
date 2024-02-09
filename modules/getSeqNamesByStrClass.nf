process GET_SEQ_NAMES_BY_STR_CLASS{
    conda params.condaEnv

    input:
    val strClass
    path mergedResultsFile

    output:
    tuple val(strClass), path("${strClass}_seqNames.npy")

    script:
    """
    filterSeqNames.py ${mergedResultsFile} ${strClass} ${strClass}_seqNames.npy
    """
}
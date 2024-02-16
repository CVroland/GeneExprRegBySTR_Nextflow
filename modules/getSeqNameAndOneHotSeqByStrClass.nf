process GET_SEQ_NAMES_AND_ONE_HOT_BY_STR_CLASS{
    publishDir "$params.resultsDir/$strClass", mode: 'copy'

    input:
        val strClass
        path seqNameFile
        path oneHotSeqFile
        path mergedResultsFile

    output:
        tuple val(strClass), path("${strClass}_seqNames.npy")
        tuple val(strClass), path("${strClass}_oneHotSeqs.npy")

    script:
    """
    filterSeqNameAndOneHotSeq.py ${seqNameFile} ${oneHotSeqFile} ${mergedResultsFile} ${strClass} ${strClass}_seqNames.npy ${strClass}_oneHotSeqs.npy
    """
}
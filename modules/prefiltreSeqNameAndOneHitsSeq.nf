process PREFILTRE_SEQ_NAMES_AND_ONE_HOT{

    input:
        path seqNameFile
        path oneHotSeqFile
        path mergedResultsFile

    output:
       path "prefiltered_seqNames.npy"
       path "prefiltered_oneHotSeqs.npy"

    script:
    """
    filterSeqNameAndOneHotSeq.py ${seqNameFile} ${oneHotSeqFile} ${mergedResultsFile} prefiltered_seqNames.npy prefiltered_oneHotSeqs.npy
    """
}
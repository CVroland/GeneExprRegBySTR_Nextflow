process COMPUTE_MNN_RESULTS{
    conda params.condaEnv

    input:
    val strClass
    path oneHotSeqFile
    path seqNameFile
    path mnnModelHParams
    path mnnModelParams

    output:
    path "mnnResultsArray.npy"

    script:
    """
    getMnnResults.py ${oneHotSeqFile} ${seqNameFile} ${mnnModelHParams} ${mnnModelParams} --output mnnResultsArray.npy
    """
}

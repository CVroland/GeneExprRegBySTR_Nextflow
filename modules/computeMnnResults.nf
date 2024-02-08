process COMPUTE_MNN_RESULTS{
    conda params.condaEnv

    input:
    tuple val(strClass), path(mnnModelHParams), path(mnnModelParams), path(strSeqNameFile)
    path oneHotSeqFile
    path seqNameFile

    output:
    tuple val(strClass), path("mnnResultsArray.npy")

    script:
    """
    getMnnResults.py ${oneHotSeqFile} ${seqNameFile} ${mnnModelHParams} ${mnnModelParams} --seqNameList ${strSeqNameFile} --output mnnResultsArray.npy
    """
}

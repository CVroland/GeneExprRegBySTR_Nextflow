process COMPUTE_MNN_RESULTS{
    conda params.condaEnv
    publishDir "$params.resultsDir/$strClass", mode: 'copy'

    input:
    tuple val(strClass), path(mnnModelHParams), path(mnnModelParams), path(strSeqNameFile), path(strOneHotSeqFile)

    output:
    tuple val(strClass), path("mnnResultsArray.npy")

    script:
    """
    getMnnResults.py ${strOneHotSeqFile} ${strSeqNameFile} ${mnnModelHParams} ${mnnModelParams} --output mnnResultsArray.npy
    """
}

process GET_STR_CLASS_BED_FILES{
    conda params.condaEnv

    input:
    tuple val(strClass), path(strClassSeqNames), path(mnnResultsArray), path(modelHParams), path(modelParams)

    output:
    tuple val(strClass), path("positiveMnnHits.bed")
    tuple val(strClass), path("negativeMnnHits.bed")

    script:
    """
    mnnResultBedFilsGenerator.py --outputDir . ${mnnResultsArray} ${strClassSeqNames} ${modelHParams} ${modelParams} --margin 0 --offset 450 --allNegHits
    """
}
process GET_MNN_RESULTS_BEDFILE{
    input: 
    val strClass
    path mnnResultsArray
    path mnnSeqNames
    path modelHParams
    path modelParams
    path strResultsDir

    output:
    path "${strResultsDir}/positiveMnnHits.bed"
    path "${strResultsDir}/negativeMnnHits.bed"

    script:
    """
    mnnResultBedFilsGenerator.py --outputDir ${strResultsDir} ${mnnResultsArray} ${mnnSeqNames} ${modelHParams} ${modelParams} --margin 0 --offset 450 --allNegHits
    """
}
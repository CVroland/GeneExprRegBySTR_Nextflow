process PLOT_MNN_SCORE{
    errorStrategy 'ignore'
    publishDir "$params.resultsDir/$strClass/$moduleId", mode: 'copy'

    input:
    tuple val(strClass), val(moduleId), path(mnnResultsArray), path(modelHParams), path(modelParams)
    val poolFunction

    output:
    tuple val(strClass), val(moduleId), path("moduleActivation_${poolFunction}.svg")

    script:
    """
    plotMnnScore.py --mnnResultsArray ${mnnResultsArray} --moduleId ${moduleId} --mnnHParams ${modelHParams} --mnnParams ${modelParams} --fig moduleActivation_${poolFunction}.svg --poolFunction ${poolFunction}
    """
}
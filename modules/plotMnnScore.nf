process PLOT_MNN_SCORE{
    errorStrategy 'ignore'
    publishDir "$params.resultsDir/$strClass/$moduleId", mode: 'copy'

    input:
    tuple val(strClass), val(moduleId), path(mnnResultsArray), path(modelHParams), path(modelParams)

    output:
    tuple val(strClass), val(moduleId), path("moduleActivation.svg")

    script:
    """
    plotMnnScore.py --mnnResultsArray ${mnnResultsArray} --moduleId ${moduleId} --mnnHParams ${modelHParams} --mnnParams ${modelParams} --fig moduleActivation.svg
    """
}
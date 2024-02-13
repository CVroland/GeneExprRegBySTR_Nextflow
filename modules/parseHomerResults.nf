process PARSE_HOMER_RESULTS{
    conda params.condaEnv

    input:
    tuple val(strClass), val(moduleId), path(homerResultsDir)
    val subName

    output:
    tuple val(strClass), val(moduleId), path("${subName}_homerResults.csv")

    script:
    """
    homerResultsToCsv.py ${homerResultsDir} -o "${subName}_homerResults.csv" --strClass ${strClass} --moduleId ${moduleId}
    """
}
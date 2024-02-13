process CONCATE_HOMER_RESULTS{

    input:
    path(strModuleHomerCsv, stageAs: "?/*")
    val subName

    output:
    path("${subName}_homerResults.csv")

    shell:
    def homerCsvList= strModuleHomerCsv instanceof List ? strModuleHomerCsv.join(" ") : strModuleHomerCsv
    '''
    awk 'NR==FNR||FNR>1' !{homerCsvList} | awk '!/^[[:space:]]*$/' > ${subName}_homerResults.csv
    '''
}
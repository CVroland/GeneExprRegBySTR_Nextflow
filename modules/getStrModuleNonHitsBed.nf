process GET_STR_MODULE_NON_HITS_BED{
    input: 
    tuple val(strClass), val(moduleId), path(strNegativeHits)
    
    output:
    tuple val(strClass), val(moduleId), path("negativeHits.bed")

    shell:
    '''
    cat !{strNegativeHits} | awk '$4==!{moduleId}' >  "negativeHits.bed"
    '''

}
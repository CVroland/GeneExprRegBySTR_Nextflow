process GET_STR_MODULE_HITS_BED{
    input: 
    tuple val(strClass), val(moduleId), path(strPositiveHits)
    
    output:
    tuple val(strClass), val(moduleId), path("positiveHits.bed")

    shell:
    '''
    cat !{strPositiveHits} | awk '$4==!{moduleId}' >  "positiveHits.bed"
    '''

}
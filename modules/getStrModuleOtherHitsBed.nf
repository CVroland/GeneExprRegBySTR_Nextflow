process GET_STR_MODULE_OTHER_HITS_BED{
    conda params.condaEnv

    input: 
    tuple val(strClass), val(moduleId), path(strPositiveHits)
    
    output:
    tuple val(strClass), val(moduleId), path("otherHits.bed")

    shell:
    '''
    cat !{strPositiveHits} | sort -k1,1 -k2,2 | bedtools groupby -g 1,2 -c 3,4,5,6 -o max,collapse,max,first | awk '$4 !~ /^!{moduleId}$|^!{moduleId},|,{moduleId},|{moduleId}$/ {print $0}' > otherHits.bed 
    '''
    // explanation : 
    // use bedtools groupby to group the hits by the first two columns (chr, start) and keep the max end, collapse the moduleIds (separate with comma), keep the max score, keep any strand (they are all positive).
    // use awk to filter out the hits that does not contain the moduleId
    // * cases : alone (^!{moduleId}$), start (^!{moduleId},), middle (,!{moduleId},), end (,!{moduleId}$) -> regex : ^!{moduleId}$|^!{moduleId},|,!{moduleId},|,!{moduleId}$
}

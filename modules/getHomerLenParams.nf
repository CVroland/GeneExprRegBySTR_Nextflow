process GET_HOMER_LEN_PARAM{
    input:
    tuple val(strClass), val(moduleId), path(strModuleHitsFasta)

    output:
    tuple val(strClass), val(moduleId), env(mlens)

    shell:
    '''
    # get the parameter for the option len as : len=[5,filterSize] if filterSize>5 else filterSize. Smallest motif in jaspar custom : 5
    filterSize=$(awk ' NR == 2 {print; exit}' !{strModuleHitsFasta} | tr -d '\n\r' | wc -m)
    #mlens=5 6 7 ... filterSize
    mlens=$(seq 5 $filterSize)
    #if filterSize<5, set mlens to filterSize
    mlens="${mlens:-$filterSize}"
    mlens=$(echo $mlens | tr -s '[:blank:]' ',' )
    '''
}
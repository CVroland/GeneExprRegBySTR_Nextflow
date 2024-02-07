process LIST_STR_CLASS{
    input:
    path mnnModelsDir // Channel.fromPath(${mnnModelsDir})

    output:
    path "X_StrClassList.txt"

    shell:
    '''
    ls !{mnnModelsDir}/*.pt | cut -d"_" -f3 | sort > X_StrClassList.txt
    '''

    
}
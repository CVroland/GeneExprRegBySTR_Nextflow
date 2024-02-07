process LIST_STR_MODULE{
    input:
    path strClassList
    path mnnModelsDir
    
    output:
    path "X_strBlockPairs.txt"

    shell:
    '''
    > X_strBlockPairs.txt
    # read all hyper parameters files and add pair (strClass, block) to the output file
    cat !{strClassList} | while read strClass
    do
        nbBlock=$(cat "!{mnnModelsDir}/MNN_ranks_${strClass}_params.npy" | head -1 | awk -F'[:(,]' '{ print $7 }')
        for ((blockId=0; blockId<nbBlock; blockId++))
        do
            echo -e "${strClass}\t${blockId}" >> X_strBlockPairs.txt
        done
    done
    '''
}
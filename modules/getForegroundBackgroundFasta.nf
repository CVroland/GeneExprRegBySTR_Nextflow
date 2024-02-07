process GET_FOREGROUND_BACKGROUND_FASTA{
    input:
    val strClass 
    val moduleId
    path positiveMnnHits
    path negativeMnnHits
    path strFasta


    output:
    path "positiveMnnHits.bed"
    path "positiveMnnHits.fasta"
    path "negativeMnnHits.bed"
    path "negativeMnnHits.fasta"
    path "otherPositiveMnnHits.bed"
    path "otherPositiveMnnHits.fasta"

    script:
    """
    cat ${positiveMnnHits} | awk "\$4 == ${moduleId}" > positiveMnnHits.bed
    bedtools getfasta -fi ${strFasta} -bed positiveMnnHits.bed > positiveMnnHits.fasta
    cat ${negativeMnnHits} | awk "\$4 == ${moduleId}" > negativeMnnHits.bed
    bedtools getfasta -fi ${strFasta} -bed negativeMnnHits.bed > negativeMnnHits.fasta
    cat ${positiveMnnHits} | awk "\$4 != ${moduleId}" | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3,4,5,6 -o max,collapse,max,first >  otherPositiveMnnHits.bed
    if [ -f otherPositiveMnnHits.bed ]  && [ -s otherPositiveMnnHits.bed ] && [ -f positiveMnnHits.bed ] && [ -s positiveMnnHits.bed ]; then
        removeSameChrStart.py otherPositiveMnnHits.bed positiveMnnHits.bed tmpOtherPositiveMnnHits.bed
        mv tmpOtherPositiveMnnHits.bed otherPositiveMnnHits.bed
        bedtools getfasta -fi ${strFasta} -bed otherPositiveMnnHits.bed > otherPositiveMnnHits.fasta
    else :
        echo "class:$strClass,block:$blockId : no positive hits or other hits found"
    fi
    """
}
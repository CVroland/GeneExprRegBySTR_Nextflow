process RENAME_SEQ_IN_FASTA{
    input:
    path inputFasta

    output:
    path "hg38.hipstr_reference.cage.500bp.around3end.2.fa"

    shell:
    '''
    cat !{inputFasta} | sed '/^>/s/|.*$//' > hg38.hipstr_reference.cage.500bp.around3end.2.fa
    '''
}
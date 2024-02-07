process MEME_TO_HOMER_FORMAT{
    conda params.condaEnv

    input:
    path memeFile

    output:
    path "${memeFile.baseName}.motif"

    script:
    """
    echo ${params.condaEnv}
    conda list
    pwm2homer.py -i ${memeFile} -m fpr --pValue 0.0001 -o ${memeFile.baseName}.motif
    """
}
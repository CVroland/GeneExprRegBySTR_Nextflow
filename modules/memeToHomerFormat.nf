process MEME_TO_HOMER_FORMAT{

    input:
    path memeFile

    output:
    path "${memeFile.baseName}.motif"

    script:
    """
    pwm2homer.py -i ${memeFile} -m fpr --pValue 0.0001 -o ${memeFile.baseName}.motif
    """
}
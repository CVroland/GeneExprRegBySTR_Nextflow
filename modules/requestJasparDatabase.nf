process REQUEST_JASPAR_DATABASE{

    input:

    output:
    path "jasparMotif_custom_r2022_cCore_gVertebrates_fMeme.txt" 

    script:
    """
    requestJasparDatabase.py -r 2022 -c CORE -g Vertebrates -V latest -f meme > jasparMotif_r2022_cCore_gVertebrates_fMeme.txt
    requestJasparDatabase.py -r 2020 -c POLII  -V latest -u "https://jaspar2020.genereg.net/api/v1/" -f meme -a >jasparMotif_r2020_cPolII_fMeme_a.txt
    cat jasparMotif_r2022_cCore_gVertebrates_fMeme.txt jasparMotif_r2020_cPolII_fMeme_a.txt > jasparMotif_custom_r2022_cCore_gVertebrates_fMeme.txt
    """
}
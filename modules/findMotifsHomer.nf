process FIND_MOTIFS_HOMER{
    publishDir "$params.resultsDir/$strClass/$moduleId/$subName", mode: 'copy'

    input:
    tuple val(strClass), val(moduleId), path(foregroundFasta), path(backgroundFasta), val(lenParam)
    path(jasparDatabaseHomer)
    val subName 

    output:
    tuple val(strClass), val(moduleId), path("homer")

    script:
    def outDir = "homer"

    """
    mkdir p "${outDir}"
    findMotifs.pl ${foregroundFasta} fasta ${outDir} -fasta ${backgroundFasta} -len ${lenParam} -norevopp -mcheck ${jasparDatabaseHomer} -mknown ${jasparDatabaseHomer} 
    """
}
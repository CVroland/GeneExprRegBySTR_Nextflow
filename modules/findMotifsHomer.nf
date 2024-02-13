process FIND_MOTIFS_HOMER{
    conda params.condaEnv

    input:
    tuple val(strClass), val(moduleId), path(foregroundFasta), path(backgroundFasta), val(lenParam)
    path(jasparDatabaseHomer)
    val subName 

    output:
    tuple val(strClass), val(moduleId), path("${subName}")

    script:
    def outDir = "${subName}"

    """
    mkdir p "${outDir}"
    findMotifs.pl ${foregroundFasta} fasta ${outDir} -fasta ${backgroundFasta} -len ${lenParam} -norevopp -mcheck ${jasparDatabaseHomer} -mknown ${jasparDatabaseHomer} 
    """
}
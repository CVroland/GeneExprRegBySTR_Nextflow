include {LIST_STR_CLASS} from './modules/listStrClass.nf'
include {LIST_STR_MODULE} from './modules/listModule.nf'
include {REQUEST_JASPAR_DATABASE} from './modules/requestJasparDatabase.nf'
include {MEME_TO_HOMER_FORMAT} from './modules/memeToHomerFormat.nf'
include {RENAME_SEQ_IN_FASTA} from './modules/renameSeqIn1001ncFasta.nf'
include{GET_SEQ_NAMES_AND_ONE_HOT_BY_STR_CLASS} from './modules/getSeqNameAndOneHotSeqByStrClass.nf'
include {COMPUTE_MNN_RESULTS} from './modules/computeMnnResults.nf'
include {GET_STR_CLASS_BED_FILES} from './modules/getStrClassBedFiles.nf'
include {GET_STR_MODULE_HITS_BED} from './modules/getStrModuleHitsBed.nf'
include {GET_FASTA_BEDTOOLS} from './modules/getFastaBedtools.nf'
include {GET_FAIDX_SAMTOOLS} from './modules/getFaidxSamtools.nf'

workflow{
    /* 
    ## Get the list of classes and moduleIds
    */
    mnnModelsDir=Channel.fromPath(params.mnnModelsDir, type:'dir')
    strClassListPath=LIST_STR_CLASS(mnnModelsDir)
    strClass=strClassListPath.splitText().map{it -> it.trim()}
    strModuleListPath=LIST_STR_MODULE(strClassListPath, mnnModelsDir)
    //strClassModule=strModuleListPath.splitCsv(sep:'\t', header:['class', 'moduleId'])
    strClassModule=strModuleListPath.splitCsv(sep:'\t')

    /* 
    ## create results directory
    */
    // TODO : Nothing?

    /*
    ## get the jaspar database
    */
    jasparDatabase=REQUEST_JASPAR_DATABASE()
    jasparDatabaseHomer=MEME_TO_HOMER_FORMAT(jasparDatabase)

    /*
    ## Prepare the fasta files of 1001 bp sequences
    */
    originalHipStr1001bpFasta=Channel.fromPath(params.originalHipStr1001bpFasta)
    hipStr1001bpFasta=RENAME_SEQ_IN_FASTA(originalHipStr1001bpFasta).first()
    hipStr1001bpFaidx=GET_FAIDX_SAMTOOLS(hipStr1001bpFasta).first()
    
    /*
    ## get MNN scores for each STR class (npy files of 3D arrays (module, seq, position))
    */
    //get the input files
    oneHotSeqFile=Channel.fromPath(params.oneHotSeqFile).first()
    seqNameFile=Channel.fromPath(params.seqNameFile).first()
    mergedResultsFile=Channel.fromPath(params.mergedResultsFile).first()
    // XXX: Nexflow does not ensure the order of the (output) channels, so we need to keep all the channels indexed by strClass
    mnnModelParams=strClass.map(strClass -> [strClass, file("${params.mnnModelsDir}/MNN_ranks_${strClass}_.pt")])
    mnnModelHParams=strClass.map(strClass -> [strClass, file("${params.mnnModelsDir}/MNN_ranks_${strClass}_params.npy")])
    // get seqNameFile and oneHotSeqFile grouped by strClass
    (strSeqNameFile, strOneHotSeqFile)=GET_SEQ_NAMES_AND_ONE_HOT_BY_STR_CLASS(strClass, seqNameFile, oneHotSeqFile, mergedResultsFile)
    //join input channel by strClass 
    // computeMnnResultsJoinedParameters : [strClass, mnnModelHParams, mnnModelParams, strSeqNameFile]
    computeMnnResultsJoinedParameters = strClass.join(mnnModelHParams).join(mnnModelParams).join(strSeqNameFile).join(strOneHotSeqFile)
    mnnResultsArray=COMPUTE_MNN_RESULTS(computeMnnResultsJoinedParameters)
    
    /*
    ## Get the fasta files of background sequences and foreground sequences
    */
    // get the "positive" and the "negative" bed files for each STR class. For sorting purpose, the blockId is store in the name column of the bed file.
    getStrClassBedFilesJoinedParameters = strClass.join(strSeqNameFile).join(mnnResultsArray).join(mnnModelHParams).join(mnnModelParams)
    (strPositiveHits, strNegativeHits) = GET_STR_CLASS_BED_FILES(getStrClassBedFilesJoinedParameters)
    // now we have general foreground and background bed files, we can make foreground and background for each module
    // Foreground : positive hits for each (strClass,module)
    getStrModuleHitsBedParameters=strPositiveHits.cross(strClassModule).map(it -> [it[1][0], it[1][1], it[0][1]]) //join and remap to get tuples (strClass, ModuleId, strPositiveHits)
    strModuleHitsBed=GET_STR_MODULE_HITS_BED(getStrModuleHitsBedParameters)
    strModuleHitsFasta=GET_FASTA_BEDTOOLS(strModuleHitsBed, hipStr1001bpFasta, hipStr1001bpFaidx)
    // Background : non hits for each (strClass,module)
    /*
    strModuleNonHitsBed
    strModuleNonHitsFasta
    // Background : for each (strClass, module) : non hits for the module but hits in other modules of the same STR class
    strModuleOtherHitsBed
    strModuleOtherHitsFasta
    */
}
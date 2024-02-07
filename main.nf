include {LIST_STR_CLASS} from './modules/listStrClass.nf'
include {LIST_STR_MODULE} from './modules/listModule.nf'
include {REQUEST_JASPAR_DATABASE} from './modules/requestJasparDatabase.nf'
include {MEME_TO_HOMER_FORMAT} from './modules/memeToHomerFormat.nf'
include {RENAME_SEQ_IN_FASTA} from './modules/renameSeqIn1001ncFasta.nf'
include {COMPUTE_MNN_RESULTS} from './modules/computeMnnResults.nf'

workflow {
    mnnModelsDir=Channel.fromPath(params.mnnModelsDir, type:'dir')
    // get the list of classes and moduleIds
    strClassListPath=LIST_STR_CLASS(mnnModelsDir)
    strClass=strClassListPath.splitText().map{it -> it.trim()}
    strModuleListPath=LIST_STR_MODULE(strClassListPath, mnnModelsDir)
    strClassModule=strModuleListPath.splitCsv(sep:'\t', header:['class', 'moduleId'])
    // create results directory
    // get the jaspar database
    jasparDatabase=REQUEST_JASPAR_DATABASE()
    jasparDatabaseHomer=MEME_TO_HOMER_FORMAT(jasparDatabase)
    // prepare the fasta files of 1001 bp
    originalHipStr1001bpFasta=Channel.fromPath(params.originalHipStr1001bpFasta)
    hipStr1001bpFasta=RENAME_SEQ_IN_FASTA(originalHipStr1001bpFasta)
    // get MNN scores for each STR class (npy files of 3D arrays (module, seq, position))
    oneHotSeqFile=Channel.fromPath(params.oneHotSeqFile)
    seqNameFile=Channel.fromPath(params.seqNameFile)
    mnnModelParams=strClass.map(strClass -> file("${params.mnnModelsDir}/MNN_ranks_${strClass}_.pt"))
    mnnModelHParams=strClass.map(strClass -> file("${params.mnnModelsDir}/MNN_ranks_${strClass}_params.npy"))
    mnnRsultsArray=COMPUTE_MNN_RESULTS(strClass.first(), oneHotSeqFile.first(), seqNameFile.first(), mnnModelHParams.first(), mnnModelParams.first())
}
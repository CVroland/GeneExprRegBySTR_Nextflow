# Identify candidate motifs responsible of STR-initiating RNA

## Introduction

This is a part of the study "Impact of transcription initiation at microsatellites on gene expression" by [M. Grapotte et al.](https://plmlatex.math.cnrs.fr/project/6352614c3e03ed92d5ea13ea). In this part, we aim to identify candidate motifs responsible of STR-initiating RNA with HOMER findMotifs.pl.

## Requirements

[Nextflow](https://www.nextflow.io/) need to be installed on your system. If it is not the case,  you can install it with the following the instructions [here](https://www.nextflow.io/docs/latest/getstarted.html#installation).

[Conda](https://docs.conda.io/en/latest/) need to be installed on your system. If it is not the case, you can install it with the following the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Download the data

The data used in this study is available [here](liketotheData). The data is composed of 4 files and 1 directory:

- `mnnModels/`: a directory containing the hyperparameters and the parameters of the different MNN models.
- `hg38.hipstr_reference.cage.500bp.around3end.fa`: a fasta file containing the sequences of the 500bp around the 3' end of all STRs.
- `hg38all_names_raw.npy` : a numpy file containing the names of STRs with a correct prediction of the CAGE signal.
- `hg38all_seqs_raw.npy` : a numpy file containing the one-hot encoded sequences of STRs with a correct prediction of the CAGE signal.

Those files need to be downloaded and placed in the `data` directory.

## Launch the pipeline

To launch the pipeline, you need to execute the following command in the repository where the main.nf file is located:

```bash
nextflow main.nf
```

A `local` profile is used to run the pipeline on your local machine (option `-profile local`). It launches the pipeline with a single process. Some processes need a large amount of memory and can crash if you run the pipeline with too much parallel executions or on a machine with limited memory. The `-resume` option allows you to resume the pipeline from where it stopped if it was stopped for any reason.

## Results

The results of the pipeline are located in the `results` directory. Pregenerated results are available [here](linktotheResults).
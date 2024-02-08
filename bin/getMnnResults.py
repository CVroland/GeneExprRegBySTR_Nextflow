#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""


Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 07/04/2023
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '07/04/2023'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import os
import argparse
import mnnProcess
import mnnPseudoModel

import numpy.typing as npt
from typing import Union
import numpy as np
import pandas as pd

def loadData(
    oneHotSeqFilePath:os.PathLike,
    namesFilePath:os.PathLike,
    seqNameList:Union[None,npt.ArrayLike]=None,
)->tuple[npt.NDArray[np.str_], npt.NDArray[np.integer]]:
    """
    Load the one-hot encoded sequences and their corresponding names from files.

    Parameters
    ----------
    oneHotSeqFilePath : PathLike
        Path to the file containing the one-hot encoded sequences.
    namesFilePath : PathLike
        Path to the file containing the sequence names.
    seqNameList : Optional[NDArray[np.str_]], optional
        Optional array of sequence names to filter the data, by default None.

    Returns
    -------
    Tuple[NDArray[np.str_], NDArray[np.integer]]
        A tuple containing the sequence names and the one-hot encoded sequences.
    """
    #load data
    seqNames, oneHotSeqs=mnnProcess.loadMnnOneHotSequences(oneHotSeqFilePath,namesFilePath)
    if seqNameList is not None:
        seqNames, oneHotSeqs=mnnProcess.filterMnnOneHotSequencesBySeqNames(seqNames, oneHotSeqs, seqNameList)
    return seqNames, oneHotSeqs
    
def loadModel(
    hParamsPath:os.PathLike, 
    paramsPath:os.PathLike
)->mnnPseudoModel.Net:
    """
    Load the MNN model.

    Parameters
    ----------
    hParamsPath : PathLike
        Path to the file containing the hyperparameters of the MNN model.
    paramsPath : PathLike
        Path to the file containing the parameters of the MNN model.

    Returns
    -------
    mnnPseudoModel.Net
        The loaded MNN model.
    """
    mnnModel:mnnPseudoModel.Net=mnnPseudoModel.load_model(hParamsPath, paramsPath)
    return mnnModel

def getMnnResults(
    oneHotSeqs:npt.NDArray[np.integer],
    mnnModel:mnnPseudoModel.Net
)->tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]:
    """
    Get the MNN results for the given one-hot encoded sequences using the MNN model.

    Parameters
    ----------
    oneHotSeqs : NDArray[np.integer]
        Array of one-hot encoded sequences.
    mnnModel : mnnPseudoModel.Net
        MNN model.

    Returns
    -------
    Tuple[NDArray[np.floating], NDArray[np.floating]]
        A tuple containing the MNN results array and the MNN max results array.
    """
    # load model
    blockList=mnnPseudoModel.getBlockList(mnnModel)
    filterLengthList=mnnPseudoModel.getFilterLengthList(blockList)
    # compute results
    mnnResultsArray, mnnMaxResultsArray=mnnProcess.getBlocksResultsArray(blockList, oneHotSeqs, filterLengthList)
    return mnnResultsArray, mnnMaxResultsArray



def parseArgs() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process MNN results.")
    parser.add_argument("oneHotSeqFilePath", type=str, help="Path to the file containing the one-hot encoded sequences.")
    parser.add_argument("namesFilePath", type=str, help="Path to the file containing the sequence names.")
    parser.add_argument("hParamsPath", type=str, help="Path to the file containing the hyperparameters of the MNN model.")
    parser.add_argument("paramsPath", type=str, help="Path to the file containing the parameters of the MNN model.")
    parser.add_argument("-l","--seqNameList", type=str, help="Path to a file containing a list of sequence names to filter the data.")
    parser.add_argument("-o","--output", type=argparse.FileType('wb'), default="-", help="Path to the output file. Use '-' for stdout. Default: stdout")
    return parser.parse_args()

def main():
    args = parseArgs()
    oneHotSeqFilePath=args.oneHotSeqFilePath
    namesFilePath=args.namesFilePath
    hParamsPath=args.hParamsPath
    paramsPath=args.paramsPath
    seqNameListPath=args.seqNameList

    if seqNameListPath is not None:
        if seqNameListPath.endswith((".npy", ".npz")):
            seqNameList = np.load(seqNameListPath)
        else: # assume it's a plain text file
            seqNameList = pd.read_csv(seqNameList, sep="\t", header=None)[0].to_numpy()
    else :
        seqNameList=None
    seqNames, oneHotSeqs = loadData(oneHotSeqFilePath, namesFilePath, seqNameList=seqNameList)
    mnnModel = loadModel(hParamsPath, paramsPath)
    mnnResultsArray, _ = getMnnResults(oneHotSeqs, mnnModel)
    np.save(args.output, mnnResultsArray)

if __name__ == "__main__":
    main()
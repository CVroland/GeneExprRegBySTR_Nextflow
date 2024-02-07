#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
get results data from MNN model

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 06/13/2023
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '06/13/2023'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import os
import io
import argparse
import pathlib

import typing
import numpy.typing as npt
from typing import Union, Tuple, TypeVar, Callable, Any, Sequence
PathLikeOrBuffer=typing.Union[os.PathLike, io.TextIOWrapper]


import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F

import mnnPseudoModel 

# model:Net=load_model(paramsPath, keysPath)
# blockList:nn.ModuleList=getBlockList(model)
# filterLengthList:list[int]=getFilterLengthList(blockList)

def loadMnnOneHotSequences(
    oneHotSeqFilePath:PathLikeOrBuffer,
    namesFilePath:PathLikeOrBuffer
)->tuple[npt.NDArray[np.str_], npt.NDArray[np.integer]]:
    """
    Load MNN one-hot encoded sequences and their corresponding names from file.

    Parameters
    ----------
    oneHotSeqFilePath : PathLikeOrBuffer
        The path or buffer to the one-hot encoded sequences file.
    namesFilePath : PathLikeOrBuffer
        The path or buffer to the sequence names file.

    Returns
    -------
    tuple
        A tuple containing the sequence names and the one-hot encoded sequences arrays.

    Raises
    ------
    FileNotFoundError
        If the specified file path does not exist.
    """
    oneHotSeqs = np.load(oneHotSeqFilePath)
    seqNames = np.load(namesFilePath)
    return seqNames, oneHotSeqs

def filterMnnOneHotSequencesBySeqNames(
    seqNames:npt.ArrayLike, 
    oneHotSeqs:npt.NDArray[np.str_],
    seqNameList:npt.NDArray[np.str_]
)->tuple[npt.NDArray[np.str_], npt.NDArray[np.str_]]:
    """
    Filter MNN one-hot encoded sequences and their corresponding names based on a list of sequence names.

    Parameters
    ----------
    seqNames : ArrayLike
        Array of sequence names.
    oneHotSeqs : NDArray[str]
        Array of one-hot encoded sequences.
    seqNameList : NDArray[str]
        List of sequence names to filter.

    Returns
    -------
    tuple
        A tuple containing the filtered sequence names and the filtered one-hot encoded sequences arrays.
    """
    uniqSeqNameList=np.unique(seqNameList)
    # the following solution is too slow !
    # mask=np.isin(seqNames, uniqSeqNameList, assume_unique=True)
    # so i will use the pythonic version with a set (search=O(1))
    uniqSeqNameSet=set(uniqSeqNameList)
    mask=np.array([seqName in uniqSeqNameSet for seqName in seqNames])
    return seqNames[mask], oneHotSeqs[mask]

def getBlocksResultsArray(
    blockList:nn.ModuleList,
    oneHotSeqs:npt.NDArray[np.integer],
    filterLengthList:list[np.integer]=None,
) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]:
    """process the convolution for each block of `blockList` for each seq in `oneHotSeqs` and return the results and the max of the result for each seq.
    
    Parameters
    ----------
    blockList : nn.ModuleList
        list of BlockNet from the pseudo model
    oneHotSeqs : ArrayLike[int]
        a array containing the sequences One-Hot encoded. The shape should be (nbSeq, seqSize, alphabetSize)
    filterLengthList : list[int], optional
        list of kernel size for each convolution in `blockList`, default None.
    
    Returns
    -------
    np.array
        the results of the convolution on the sequences. The shape should be `(nbBlock, nbSeq, seqSize)`
    np.array
        the max of the results of the convolution on the sequences by sequence. The shape should be `(nbBlock,nbSeq)`
    """
    if filterLengthList is None :
        filterLengthList=mnnPseudoModel.getFilterLengthList(blockList)
    # Adding 1 dim as a single Channel because BlockNet use a 2D convolution with a kernel of size (filterLength, alphabet) and a 0 padding. 
    # oneHotSeqs dim (seq#, letter, char)
    # seqs dim (seq#, 1, letter, char)
    seqs = torch.from_numpy(np.expand_dims(oneHotSeqs, axis=1)).float()
    # get cnn result for all block : 
    # block(seqs) return a tensor of dim (seq#, letter-(filter_size-1))
    # add Padding with 0 on the right according to the filterLength
    mnnResultsListTorchList=[
        F.pad(block(seqs), (0,filterLength-1),"constant", 0) 
            for block, filterLength in zip(blockList, filterLengthList)
    ]
    # stack to form a single tensor of dim (block,seq,seqSize)
    mnnResultsTorch=torch.stack(mnnResultsListTorchList)
    # compute the max value foreach sequence along the seqSize axis
    # the dimensions of mnnResultsMaxTorch should be (block, seq)
    mnnResultsMaxTorch=torch.amax(mnnResultsTorch, dim=-1)
    # Now transform the torch tensors into numpy arrays. The graph should be compute here
    mnnResultsArray=mnnResultsTorch.detach().cpu().numpy()
    mnnMaxResultsArray=mnnResultsMaxTorch.detach().cpu().numpy()
    return mnnResultsArray, mnnMaxResultsArray

def getMnnIntervals(
    mnnResultsArray:npt.NDArray,
    filterLengthList:Sequence[int],
    sequenceNames:Sequence[str]=None,
    blockNames:Sequence[str]=None
):
    """
    Compute the hit intervals of MNN results based on the maximum results array.

    Parameters
    ----------
    mnnResultsArray : npt.NDArray
        The MNN results array of shape (nbBlock, nbSeq, seqSize).
    filterLengthList : ArrayLike[int]
        List of filter lengths for each block.
    sequenceNames : ArrayLike[str], optional
        Array of sequence names, by default None.
    blockNames : ArrayLike[str], optional
        Array of block names, by default None.
    """
    if sequenceNames is None :
        sequenceNames = [str(i) for i in range(np.shape(mnnResultsArray)[1])]
    sequenceNames=np.asarray(sequenceNames)
    if blockNames is None :
        blockNames=[str(i) for i in range (len(mnnResultsArray))]
    blockNames=np.asarray(blockNames)
    blockIdx, seqIdx, matchIdx=np.nonzero(mnnResultsArray>0)
    sequences=sequenceNames[seqIdx]
    blocks=blockNames[blockIdx]
    filterLengthList=np.asarray(filterLengthList)
    filterLengths=filterLengthList[blockIdx]
    begin=matchIdx
    end=matchIdx+filterLengths
    score=mnnResultsArray[(blockIdx, seqIdx, matchIdx)].reshape(-1)
    return pd.DataFrame(
        data={
            "sequence_name":sequences,
            "start":begin,
            "stop":end,
            "block":blocks,
            "score":score
        }
    )


def getMnnScoreMatrix(
    mnnMaxResultsArray:npt.NDArray,
    sequenceNames:Sequence[str]=None,
    blockNames:Sequence[str]=None
):
    if sequenceNames is None :
        sequenceNames = [str(i) for i in range(np.shape(mnnMaxResultsArray)[1])]
    if blockNames is None :
        blockNames=[str(i) for i in range (len(mnnMaxResultsArray))]
    return pd.DataFrame(
        data=mnnMaxResultsArray.T,
        index=sequenceNames,
        columns=blockNames
    )


def getMnnHitPfm(
    mnnIntervalsDf:pd.DataFrame,
    seqNames:npt.ArrayLike, 
    oneHotSeqs:npt.NDArray[np.str_],
    weightByScore=False
):
    sorter=np.argsort(seqNames)
    mnnIntervalsSeqIndex=sorter[np.searchsorted(seqNames, mnnIntervalsDf["sequence_name"], sorter=sorter)]
    mnnIntervalsDf["seqIdx"]=mnnIntervalsSeqIndex
    def getMnnHitPfmPerBlock_(grp:pd.DataFrame):
        indices=np.asarray([np.arange(start,stop) for start,stop in zip(grp["start"], grp["stop"])])
        weight=grp["score"].to_numpy()  if weightByScore else np.asarray([1])
        return np.sum(
            np.take_along_axis(
                oneHotSeqs[grp["seqIdx"]]*weight[..., np.newaxis, np.newaxis], 
                indices=indices[...,np.newaxis], 
                axis=1
            ),
            axis=0
        )
    blockIdArray=mnnIntervalsDf["block"].unique()
    pfmList={}
    for blockId in blockIdArray :
        pfmList[blockId]=getMnnHitPfmPerBlock_(mnnIntervalsDf.loc[mnnIntervalsDf["block"]==blockId])
    return pfmList
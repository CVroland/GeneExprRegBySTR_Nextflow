#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute a correlation between MNN Score and max Fimo -log(p-value) by including non MNN hit as score of 0 and non Fimo hit as score of 0.
Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 07/13/2023
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '07/19/2023'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import pathlib
import argparse
import warnings
import numpy as np
#setup a random generator with a fix seed
npRandomGen=np.random.default_rng(seed=42)
import pandas as pd

from typing import Any, Sequence, Union
import numpy.typing as npt

import mnnPseudoModel

def getScore(mnnResultsArray: np.ndarray, blockIdx: np.ndarray, seqIdx: np.ndarray, matchIdx: np.ndarray) -> np.ndarray:
    """
    Get the scores from the MNN results array based on the given block, sequence, and match indices.

    Parameters
    ----------
    mnnResultsArray : numpy.ndarray
        The 3D NumPy array representing the MNN results.
    blockIdx : numpy.ndarray
        An array of integers representing block indices for each match.
    seqIdx : numpy.ndarray
        An array of integers representing sequence indices for each match.
    matchIdx : numpy.ndarray
        An array of integers representing match indices.

    Returns
    -------
    numpy.ndarray
        An array containing the scores corresponding to the provided indices.
    """
    return mnnResultsArray[blockIdx, seqIdx, matchIdx]

def getMnnHitPos(mnnResultsArray: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get the indices of positive hits from the MNN results.

    Parameters
    ----------
    mnnResultsArray : numpy.ndarray
        The 3D NumPy array representing the MNN results.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple containing three arrays representing block indices, sequence indices, and match indices of positive hits.
    """
    blockIdx, seqIdx, matchIdx=np.nonzero(mnnResultsArray>0)
    return blockIdx, seqIdx, matchIdx

def getMnnNonHitPos(mnnResultsArray: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get the indices of non-hit positions from the MNN results.

    Parameters
    ----------
    mnnResultsArray : numpy.ndarray
        The 3D NumPy array representing the MNN results.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple containing three arrays representing block indices, sequence indices, and match indices of non-hit positions.
    """
    blockIdx, seqIdx, matchIdx=np.nonzero(mnnResultsArray<=0)
    return blockIdx, seqIdx, matchIdx

def _randomDrawByBlock(
        
):
    pass
    

def randomDraw(
    blockIdx: np.ndarray, 
    seqIdx: np.ndarray, 
    matchIdx: np.ndarray,
    filterLengthList: np.ndarray,
    sequenceLength: int,
    blockIdxLike: np.ndarray, 
    seqIdxLike: np.ndarray, 
    warning:bool=True #XXX
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Randomly draws matches from given block and sequence indices

    Parameters
    ----------
    blockIdx : numpy.ndarray
        An array of integers representing block indices for each match.
    seqIdx : numpy.ndarray
        An array of integers representing sequence indices for each match.
    matchIdx : numpy.ndarray
        An array of integers representing match indices.
    filterLengthList : numpy.ndarray
        An array of integers representing the filter lengths for each block.
    sequenceLength : int
        The length of the sequences.
    blockIdxLike : numpy.ndarray
        An array of integers representing block indices used for drawing matches.
    seqIdxLike : numpy.ndarray
        An array of integers representing sequence indices used for drawing matches.
    warning : bool, optional
        Whether to issue a warning and draw additional matches with replacement if there are not enough negative hits 
        (default is True).

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        drawBlockIdx : numpy.ndarray
            An array of integers representing the drawn block indices.
        drawSeqIdx : numpy.ndarray
            An array of integers representing the drawn sequence indices.
        drawMatchIdx : numpy.ndarray
            An array of integers representing the drawn match indices.

    Notes
    -----
    - The function draws matches based on the specified block and sequence indices provided in `blockIdxLike` and `seqIdxLike`.
    - The number of drawn matches from each block and sequence is determined by the occurrence count of the corresponding block and sequence indices.
    - The function performs sampling without replacement (i.e., `replace=False`), and it may draw additional matches with replacement if there are not enough negative hits, based on the `warning` parameter.
    - If `warning` is True, a warning will be issued if there are not enough negative hits to fulfill the drawing requirements, and the function will attempt to draw additional matches with replacement to meet the specified number of draws.
    - If `warning` is False, a `ValueError` will be raised if there are not enough negative hits to fulfill the drawing requirements.

    Raises
    ------
    ValueError
        If `warning` is False and there are not enough negative hits to fulfill the drawing requirements.
    """
    # count the number of hits in a (block, seq) groups
    blockSeqLikeStack, count=np.unique(
        np.stack([blockIdxLike, seqIdxLike], axis=-1),
        return_counts=True,
        axis=0
    )
    blockIds=blockSeqLikeStack[:,0]
    seqIds=blockSeqLikeStack[:,1]
    drawnMatchIdxList=[]
    for i in range(0, len(blockSeqLikeStack)):
        blockId=blockIds[i]
        #get subMatchIdx that correspond to all matchIdx with the same blockIds and the same seqIds
        subMatchIdx=matchIdx[(blockIdx==blockId) & (seqIdx==seqIds[i])]
        #remove also matches where end is beyond the length of the sequence: their's dummy values (0) on the right side due to the padding after the convolution.
        subMatchIdx=subMatchIdx[subMatchIdx<=sequenceLength-filterLengthList[blockId] ]
        # pass if their's no subMatch : impossible to draw
        if len(subMatchIdx) == 0 : 
            msg="No value for block {}, sequence {} : {} needed, {} available".format(blockIds[i], seqIds[i], count[i], len(subMatchIdx))
            if warning :
                warnings.warn(msg)
                pass
            else :
                raise ValueError(msg)
        try :
            drawnMatchIdxList.append(
                npRandomGen.choice(
                    subMatchIdx, # get the same number of draw by block and by sequences
                    size=count[i], 
                    replace=False,
                    shuffle=False
                )
            )
        except ValueError as e :
            msg="Not enough value for block {}, sequence {} : {} needed, {} available".format(blockIds[i], seqIds[i], count[i], len(subMatchIdx))
            if warning :
                warnings.warn(msg)
                #add the available and draw with replacement the rest
                drawnMatchIdxList.append(subMatchIdx)
                additionalCount=count[i]-len(subMatchIdx)
                additionalDrawnMatch=npRandomGen.choice(
                    subMatchIdx,
                    size=additionalCount,
                    replace=True
                )
                drawnMatchIdxList.append(additionalDrawnMatch)
            else :
                raise ValueError(msg)

    drawMatchIdx=np.concatenate(drawnMatchIdxList, axis=None) #axis=0 should work too
    #now recreate drawBlockIdx and drawSeqIdx
    drawBlockSeqStack=np.repeat(blockSeqLikeStack, count, axis=0)
    drawBlockIdx=drawBlockSeqStack[:,0]
    drawSeqIdx=drawBlockSeqStack[:,1]
    return drawBlockIdx, drawSeqIdx, drawMatchIdx


def processOffset(offset: Union[int, tuple[int, int]]) -> tuple[int, int]:
    """
    Process the offset parameter.

    Parameters
    ----------
    offset : int or tuple[int, int]
        Offset to process.

    Returns
    -------
    tuple[int, int]
        Processed left and right offsets.
    """
    if hasattr(offset, '__getitem__') and len(offset) >=1:
        leftOffset=offset[0]
        rightOffset=offset[1] if len(offset)>=2 else offset[0]
    else :
        leftOffset = rightOffset = offset
    return leftOffset, rightOffset

def processMargin(margin: Union[int, tuple[int, int]]) -> tuple[int, int]:
    """
    Process the margin parameter.

    Parameters
    ----------
    margin : int or tuple[int, int]
        Margin to process.

    Returns
    -------
    tuple[int, int]
        Processed left and right margins.
    """
    if hasattr(margin, '__getitem__') and len(margin) >=1:
        leftMargin=margin[0]
        rightMargin=margin[1] if len(margin)>=2 else margin[0]
    else :
        leftMargin = rightMargin = margin
    return leftMargin, rightMargin

def getBed(
    blockIdx: np.ndarray,
    seqIdx: np.ndarray,
    matchIdx: np.ndarray,
    filterLengths: Sequence[int],
    sequenceLength: int,
    margin: Union[int, tuple[int, int]] = 0,
    offset: Union[int, tuple[int, int]] = 0,
    scores: Any = 0,
    blockNames: Any = None,
    seqNames: Any = None
) -> pd.DataFrame:
    """
    Generate a DataFrame in BED format containing the positions of matches. Return an empty DataFrame if there's no entry.

    Parameters
    ----------
    blockIdx : numpy.ndarray
        An array of integers representing block indices for each match.
    seqIdx : numpy.ndarray
        An array of integers representing sequence indices for each match.
    matchIdx : numpy.ndarray
        An array of integers representing match indices.
    filterLengths : Sequence[int]
        List of integers representing the filter lengths for each match.
    sequenceLength : int
        The length of the sequences.
    margin : int or tuple[int, int], optional
        Margin to add on both sides of the match positions (default is 0).
    offset : int or tuple[int, int], optional
        Offset to add to the match positions (default is 0).
    scores : Any, optional
        Scores corresponding to each match (default is 0).
    blockNames : Any, optional
        Names for each block (default is None).
    seqNames : Any, optional
        Names for each sequence (default is None).

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the BED data for the matches.
    """
    # case where there's no entry : create an empty dataframe
    if len(blockIdx) == 0 :
        return pd.DataFrame(columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand"]).astype(
            {
                "chrom":str,
                "chromStart":int,
                "chromEnd":int,
                "name":str,
                "score":float,
                "strand":str
            }
        )
    if blockNames is None :
        blockNames=np.arange(np.max(blockIdx)+1)
    else :
        blockNames=np.asarray(blockNames)
    if seqNames is None :
        seqNames=np.arange(np.max(seqIdx))
    else :
        seqNames=np.asarray(seqNames)
    leftMargin, rightMargin = processMargin(margin)
    leftOffset, rightOffset = processOffset(offset)
    minStart=0
    maxEnd=sequenceLength+rightOffset
    filterLengths=np.asarray(filterLengths)
    chrom=seqNames[seqIdx]
    # compute start and stop in the bedfile
    # add margin and change offset to correspond to new sequences
    # verify if not overflow
    chromStart=matchIdx-leftMargin+leftOffset 
    chromStart=np.maximum(chromStart, minStart)
    chromEnd=matchIdx+filterLengths[blockIdx]+rightMargin+leftOffset
    chromEnd=np.minimum(chromEnd, maxEnd)
    names=blockNames[blockIdx]
    bedDf=pd.DataFrame(
        {
            "chrom":chrom,
            "chromStart":chromStart,
            "chromEnd":chromEnd,
            "name":names,
            "score":scores,
            "strand":["+"]*len(chrom)
        }
    ).astype(
        {
            "chrom":str,
            "chromStart":int,
            "chromEnd":int,
            "name":str,
            "score":float,
            "strand":str
        }
    )
    return bedDf

def generateMnnResultBedFiles(
    mnnResultsArray: npt.NDArray,
    filterLengthList: Sequence[int],
    seqNames: Sequence[str] = None,
    margin: Union[int, tuple[int, int]] = 0,
    offset: Union[int, tuple[int, int]] = 0,
    allNegHits=False
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate BED files for positive and randomly selected negative hits from the MNN results.

    Parameters
    ----------
    mnnResultsArray : numpy.typing.NDArray
        3D NumPy array representing the MNN results.
    filterLengthList : Sequence[int]
        List of integers representing the filter lengths for each match.
    seqNames : Sequence[str], optional
        List of sequence names (default is None).
    margin : int or tuple[int, int], optional
        Margin to add on both sides of the match positions (default is 0).
    offset : int or tuple[int, int], optional
        Offset to add to the match positions (default is 0).

    Returns
    -------    
    Tuple[pandas.DataFrame, pandas.DataFrame]
        A tuple containing two DataFrames:
        - The first DataFrame contains the BED data for positive hits.
        - The second DataFrame contains the BED data for randomly selected negative hits.
    """
    # get positive mnnResult hits
    posBlockIdx, posSeqIdx, posMatchIdx=getMnnHitPos(mnnResultsArray)
    posScores=getScore(mnnResultsArray, posBlockIdx, posSeqIdx, posMatchIdx)
    # get all negative mnnResult hits
    allNegBlockIdx, allNegSeqIdx, allNegMatchIdx=getMnnNonHitPos(mnnResultsArray)
    # draw negative mnnResults hits
    sequenceLength=mnnResultsArray.shape[2]
    if allNegHits :
        #remove matches where end is beyond the length of the sequence: their's dummy values (0) on the right side due to the padding after the convolution.
        filterLengthRepeated=np.asarray(filterLengthList)[allNegBlockIdx]
        outOfBoundMask=(allNegMatchIdx<=sequenceLength-filterLengthRepeated)
        negBlockIdx, negSeqIdx, negMatchIdx=allNegBlockIdx[outOfBoundMask], allNegSeqIdx[outOfBoundMask], allNegMatchIdx[outOfBoundMask]
    else :
        negBlockIdx, negSeqIdx, negMatchIdx=randomDraw(allNegBlockIdx, allNegSeqIdx, allNegMatchIdx, filterLengthList, sequenceLength, posBlockIdx, posSeqIdx)
    # generate bed files
    posBedDf=getBed(
        posBlockIdx,
        posSeqIdx,
        posMatchIdx,
        filterLengthList,
        sequenceLength,
        margin=margin,
        offset=offset,
        scores=posScores,
        seqNames=seqNames
    )
    negBedDf=getBed(
        negBlockIdx,
        negSeqIdx,
        negMatchIdx,
        filterLengthList,
        sequenceLength,
        margin=margin,
        offset=offset,
        scores=0,
        seqNames=seqNames
    )
    return posBedDf,negBedDf


def parseArgs():
    parser = argparse.ArgumentParser(description="Generate BED files for positive and randomly selected negative hits from MNN results.")
    parser.add_argument("mnnResultsArray", type=pathlib.Path, help="Path to the .npy file containing the MNN results array.")
    parser.add_argument("seqNames", type=pathlib.Path, help="Path to the .npy file containing the sequence names.")
    parser.add_argument("modelHParam", type=pathlib.Path, help="path to hyper-parameters of the MNN model")
    parser.add_argument("modelParam", type=pathlib.Path, help="path to parameters of the MNN model")
    parser.add_argument("--outputDir", type=pathlib.Path, default=pathlib.Path.cwd(), help="Path to the output directory for saving the BED files.")
    parser.add_argument("--margin", type=int, nargs="+", default=[0], help="Margin to add on both sides of the match positions (default is 0).")
    parser.add_argument("--offset", type=int, nargs="+", default=[0], help="Offset to add to the match positions (default is 0).")
    parser.add_argument("--allNegHits", action="store_true", help="Return all negative hits instead of a subset.")
    return parser.parse_args()

def main():
    args = parseArgs()

    # Load data from files
    mnnResultsArray = np.load(args.mnnResultsArray)
    seqNames = np.load(args.seqNames)

    # Load model and get filter length list
    model = mnnPseudoModel.load_model(args.modelHParam, args.modelParam)
    blockList = mnnPseudoModel.getBlockList(model)
    filterLengthList = mnnPseudoModel.getFilterLengthList(blockList)
    del blockList
    del model

    # Generate BED files
    margin = args.margin
    offset = args.offset
    allNegHits=args.allNegHits
    posBedDf, negBedDf = generateMnnResultBedFiles(mnnResultsArray, filterLengthList, margin=margin, offset=offset, seqNames=seqNames, allNegHits=allNegHits)

    # Save the BED files
    outputDir=args.outputDir
    outputDir.mkdir(parents=True, exist_ok=True)
    posBedDf.to_csv(outputDir / "positiveMnnHits.bed", sep="\t", index=False, header=False)
    del posBedDf
    negBedDf.to_csv(outputDir / "negativeMnnHits.bed", sep="\t", index=False, header=False)
    del negBedDf

if __name__ == "__main__":
    main()
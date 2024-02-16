#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter sequenceNames file (hg38all_names_raw.npy) and oneHot sequences file (hg38all_seqs_raw.npy) to keep only the sequences sequences for a given STR class.

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 08/02/2024
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '08/02/2024'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import argparse
import numpy as np
import pandas as pd
from miscFct import splitEStrHeader, rcDnaSeq

def filterSeqNames(
    seqNamesArray:np.ndarray,
    strClass:str
)->np.ndarray:
    """
    Filter sequenceNames file to keep only the sequences for a given STR class.

    Parameters
    ----------
    seqNamesArray : np.ndarray
        The sequence names array.
    strClass : str
        The STR class sequence of the STR class to keep.

    Returns
    -------
    np.ndarray
        The filtered sequence names array.

    """
    rcStrClass=rcDnaSeq(strClass)
    seqHeaderDf=splitEStrHeader(seqNamesArray)
    strClassMask=np.asarray(
        ((seqHeaderDf["seq"]==strClass) & (seqHeaderDf["strand"]=='+')) |
        ((seqHeaderDf["seq"]==rcStrClass) & (seqHeaderDf["strand"]=='-'))
    )
    return np.unique(seqNamesArray[strClassMask])

def getSeqNamesArray(mergedResultsFilePath:str)->np.ndarray:
    """
    Get the sequence names array from the merged results file.

    Parameters
    ----------
    mergedResultsFilePath : str
        Path to the merged results file.

    Returns
    -------
    np.ndarray
        The sequence names array.

    """
    # merged_results.txt is a TSV file where the first column is the STR-gene association separated by a ':'
    strSeries=pd.read_csv(mergedResultsFilePath, sep=':', header=None, index_col=False, usecols=[0])[0]
    return strSeries.to_numpy(dtype=str)


def getStrMask(allSeqNames:np.ndarray, strSeqNames:np.ndarray)->np.ndarray:
    """
    Get the mask for the sequence names of a given STR class.

    Parameters
    ----------
    allSeq
    Names : np.ndarray
        The sequence names array for all STR classes.
    strSeqNames : np.ndarray
        The sequence names array for the STR class.

    Returns
    -------
    np.ndarray
        The mask for the sequence names of the given STR class.

    """
    # using a set and python is way more efficient than using the np.isin function
    strSeqNamesSet=set(strSeqNames)
    mask = np.array([seqName in strSeqNamesSet for seqName in allSeqNames])
    return mask

def filterSeqNamesAndOneHotSeqArray(allSeqNamesArray:np.ndarray, allOneHotSeqArray:np.ndarray, strSeqNames:np.ndarray)->tuple[np.ndarray, np.ndarray]:
    """
    Filter sequenceNames file and oneHotSeq file to keep only the sequences sequences for a given STR class.

    Parameters
    ----------
    allSeqNamesArray : np.ndarray
        The sequence names array for all STR classes.
    allOneHotSeqArray : np.ndarray
        The oneHotSeq array for all STR classes.
    strSeqNames : np.ndarray
        The sequence names array for the STR class.

    Returns
    -------
    (np.ndarray, np.ndarray)
        The filtered sequence names array and the filtered oneHotSeq array.

    """
    mask = getStrMask(allSeqNamesArray, strSeqNames)
    return allSeqNamesArray[mask], allOneHotSeqArray[mask]

def main():
    #parse arguments
    parser=argparse.ArgumentParser(description="Filter sequenceNames file to keep only the sequences for a given STR class.")
    parser.add_argument("allSeqNamesFilePath", type=str, help="Path to the `hg38all_names_raw.npy` file.")
    parser.add_argument("allOneHotSeqFilePath", type=str, help="Path to the `hg38all_seqs_raw.npy` file.")
    parser.add_argument("mergedResultsFilePath", type=str, help="Path to the merged_results.txt file.")
    parser.add_argument("strClass", type=str, help="The STR class sequence to keep.")
    parser.add_argument("outputSeqNamesFilePath", type=str, help="Path to the output sequenceNames file.")
    parser.add_argument("outputOneHotSeqFilePath", type=str, help="Path to the output oneHotSeq file.")
    args=parser.parse_args()

    # get the sequence names array for all STR classes
    #load data
    mergedResultsSeqNamesArray=getSeqNamesArray(args.mergedResultsFilePath)
    #filter data
    unSortedUniqStrSeqNames=filterSeqNames(mergedResultsSeqNamesArray, args.strClass)
    del mergedResultsSeqNamesArray
    # get the seqNamesArray and oneHotSeqArray for the STR class
    allSeqNamesArray=np.load(args.allSeqNamesFilePath)
    allOneHotSeqArray=np.load(args.allOneHotSeqFilePath)
    #filter data
    seqNamesArray, oneHotSeqArray = filterSeqNamesAndOneHotSeqArray(allSeqNamesArray, allOneHotSeqArray, unSortedUniqStrSeqNames)
    #save data
    np.save(args.outputSeqNamesFilePath, seqNamesArray)
    np.save(args.outputOneHotSeqFilePath, oneHotSeqArray)


if __name__ == "__main__":
    main()
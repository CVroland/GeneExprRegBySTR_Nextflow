#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter sequenceNames file (hg38all_names_raw.npy) to keep only the sequences for a given STR class.

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 08/02/2024
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '08/02/2024'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import os
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


def main():
    #parse arguments
    parser=argparse.ArgumentParser(description="Filter sequenceNames file to keep only the sequences for a given STR class.")
    parser.add_argument("mergedResultsFilePath", type=str, help="Path to the sequenceNames file.")
    parser.add_argument("strClass", type=str, help="The STR class sequence of the STR class to keep.")
    parser.add_argument("outputFilePath", type=str, help="Path to the output file.")
    args=parser.parse_args()
    #load data
    seqNamesArray=getSeqNamesArray(args.mergedResultsFilePath)
    #filter data
    filteredSeqNamesArray=filterSeqNames(seqNamesArray, args.strClass)
    #save data
    np.save(args.outputFilePath, filteredSeqNamesArray)

if __name__=="__main__":
    main()

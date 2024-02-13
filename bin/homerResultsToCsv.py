#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scrap HOMER results directory and compile results into a CSV file

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 08/02/2024
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '08/02/2024'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import argparse
import io
import os
import sys
import pathlib
import re

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from typing import List, Tuple, Dict, Union

PathLike=Union[str, pathlib.Path, os.PathLike]

def parseMotifInfoHtml(htmlFilePath: PathLike) -> Tuple[dict, List[dict]]:
    """
    Parse motif information from an HTML file.

    Parameters
    ----------
    htmlFilePath : PathLike
        Path to the HTML file containing motif information.

    Returns
    -------
    Tuple[dict, List[dict]]
        A tuple containing a dictionary with motif information and a list of dictionaries with match information.
    """
    soup=BeautifulSoup(open(htmlFilePath), features="lxml")
    #get info table
    motifTitle=soup.select("html > body > h2")[0].text
    motifNameRegex=re.compile(r"Information for (.*) \(.*\)")
    motifName=motifNameRegex.match(motifTitle).group(1)
    motifInfoTable=soup.select("html > body > table")[0]
    motifInfoRawDf=pd.read_html(io.StringIO(str(motifInfoTable)))[0]
    motifInfo={
        "name":motifName,
        "p":np.float128(motifInfoRawDf.loc[0,1]),
        "logP":float(motifInfoRawDf.loc[1,1]),
        "infoContentPerBp":float(motifInfoRawDf.loc[2,1]),
        "t":float(motifInfoRawDf.loc[3,1]),
        "pT":float(motifInfoRawDf.loc[4,1][:-1]), #remove '%' char at the end
        "b":float(motifInfoRawDf.loc[5,1]),
        "pB":float(motifInfoRawDf.loc[6,1][:-1]), #remove '%' char at the end
    }

    #get match Name list
    matchNameList=[e.text for e in soup.select("html > body > table > tr > td > h4")]
    #get match information tables list
    matchInfoTableList=soup.select("html > body > table > tr > td > table > tr > td > table")
    matchInfoDfList=pd.read_html(io.StringIO(str(matchInfoTableList)))
    matchInfoList=[]
    for i,matchInfoDf in enumerate(matchInfoDfList) :
        matchInfo={
            "matchName":matchNameList[i],
            "matchRank":int(matchInfoDf.loc[0,1]),
            "matchScore":float(matchInfoDf.loc[1,1]),
            "matchStrand":'+' if matchInfoDf.loc[3,1]=='forward strand' else '-'
        }
        matchInfoList.append(matchInfo)
    return motifInfo, matchInfoList

def motifInfoHtmlToMatchDf(motifInfo: dict, matchInfoList: List[dict]) -> pd.DataFrame:
    """
    Join parsed motif and match information into a DataFrame.

    Parameters
    ----------
    motifInfo : dict
        Dictionary containing motif information.
    matchInfoList : List[dict]
        List of dictionaries containing match information.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing match information with added motif information as columns.

    """
    matchInfoDf=pd.DataFrame(matchInfoList)
    return matchInfoDf.assign(**motifInfo)

def readMatchDfInMotifInfoHtml(htmlFilePath: PathLike) -> pd.DataFrame:
    """
    Read match information from an HTML file and convert it to a DataFrame.

    Parameters
    ----------
    htmlFilePath : PathLike
        Path to the HTML file containing match information.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing match information with added motif information as columns.

    """
    motifInfo, matchInfoList=parseMotifInfoHtml(htmlFilePath)
    matchInfoDf=motifInfoHtmlToMatchDf(motifInfo, matchInfoList)
    return matchInfoDf

def getMotifInfoHtmlFilePathList(homerResultsDirPath: PathLike) -> List[PathLike]:
    """
    Get a list of motif information HTML file paths in a HOMER results directory.

    Parameters
    ----------
    homerResultsDirPath : PathLike
        Path to the HOMER results directory.

    Returns
    -------
    List[PathLike]
        A list of file paths to motif information HTML files.
    """
    homerResultsDirPath=pathlib.Path(homerResultsDirPath)
    motifInfoHtmlFileNameRegex=re.compile(r"motif\d+.info.html")
    motifInfoHtmlFilePathList=[
        path for path in homerResultsDirPath.iterdir()
        if path.is_file() and motifInfoHtmlFileNameRegex.match(path.name) is not None
    ]
    return motifInfoHtmlFilePathList

def readMatchDfInHomerResultsDir(homerResultsDirPath: PathLike) -> pd.DataFrame:
    """
    Read match information from HTML files in a HOMER results directory and convert it to a DataFrame.

    Parameters
    ----------
    homerResultsDirPath : str
        Path to the HOMER results directory.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing match information with added motif information as columns.

    Example
    -------
    match_df = readMatchDfInHomerResultsDir('path/to/homer/results')
    """
    homerResultsDirPath = pathlib.Path(homerResultsDirPath)
    homerResultsSubDirPath = homerResultsDirPath / "homerResults"
    # Check if the homerResultsSubDirPath exists, if not return an empty dataframe for Nextflow compatibility
    if not homerResultsSubDirPath.exists():
        return pd.DataFrame(columns=["matchName", "matchRank", "matchScore", "matchStrand", "name", "p", "logP", "infoContentPerBp", "t", "pT", "b", "pB"])
    motifInfoHtmlFilePathList=getMotifInfoHtmlFilePathList(homerResultsSubDirPath)
    matchDfList=[
        readMatchDfInMotifInfoHtml(motifInfoHtmlFilePath) for motifInfoHtmlFilePath in motifInfoHtmlFilePathList
    ]
    return pd.concat(matchDfList).reset_index(drop=True)

def main():
    parser = argparse.ArgumentParser(description="Scrap HOMER results directory and compile results into a CSV file")
    parser.add_argument("homerResultsDirPath", type=str, help="Path to the HOMER results directory")
    parser.add_argument("-o", "--output", type=str, default="-", help="Path to the output CSV file. Use '-' for stdout. Default: stdout")
    parser.add_argument("--strClass", type=str, default=None, help="STR class of the results")
    parser.add_argument("--moduleId", type=str, default=None, help="Module ID of the results")
    args = parser.parse_args()
    matchDf=readMatchDfInHomerResultsDir(args.homerResultsDirPath)
    if args.moduleId is not None:
        matchDf.insert(0, "moduleId", args.moduleId)
    if args.strClass is not None:
        matchDf.insert(0, "strClass", args.strClass)
    output=args.output if args.output != "-" else sys.stdout
    matchDf.to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    main()
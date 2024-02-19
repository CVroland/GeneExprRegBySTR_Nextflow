#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
return a histogram score for a specific motif.

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 30/01/2024
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '19/02/2024'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"
 
# python plotMnnScore.py --mnnResultsArray mnnResultsArray.npy --moduleId moduleId --mnnHParams mnnHParams  --mnnParams mnnParams --fig OUTPUT_FIGURE_PATH --values OUTPUT_VALUES_PATH --bias


import functools
import argparse
import pathlib

import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

from typing import Generator, NewType, Union, Tuple

FString=NewType("FString", str)
PdQuery=NewType("PdQuery", str)
PathLike=Union[str, pathlib.Path]

import mnnPseudoModel


def getMnnModuleResultsIdx(
    mnnModuleResultsArray:np.ndarray, 
    threshold:float=0
)->Tuple[np.ndarray, np.ndarray]: ##tuple of 1D array of idx
    """
    Get the mnn results idx for a module. 
    Return the idx of the score above the threshold.

    Parameters
    ----------
    mnnModuleResultsArray : np.ndarray
        The mnn results array. Score of the convolution in a 2D array (seq, pos)
    threshold : float, optional
        The threshold, by default 0

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        The idx of the score above the threshold for the module.
    """
    seqIdx,PosIdx=np.nonzero(mnnModuleResultsArray>threshold)
    return seqIdx, PosIdx

@functools.lru_cache(maxsize=4)
def loadMnnModel(
    paramsPath:PathLike,
    keysPath:PathLike
)->mnnPseudoModel.Net:
    """
    Load the mnn model from the given paths.

    Parameters
    ----------
    paramsPath : PathLike
        Path to the mnn model parameters.
    keysPath : PathLike
        Path to the mnn model keys.

    Returns
    -------
    mnnPseudoModel.Net
        The mnn model.

    """
    return mnnPseudoModel.load_model(
        paramsPath=paramsPath, 
        keysPath=keysPath
    )

def getPosCoefArray(mnn:mnnPseudoModel.Net, moduleId, seqSize=101)->np.ndarray: #1D array of weight (pos)
    """
    Get the position coefficient array for a given module.

    Parameters
    ----------
    mnn : mnnPseudoModel.Net
        The mnn model.
    moduleId : int
        The module id.
    seqSize : int, optional
        The sequence size, by default 101

    Returns
    -------
    np.ndarray
        The position coefficient array for the module.

    """
    moduleList=mnnPseudoModel.getBlockList(mnn)
    posCoefList=[module.dense.weight.detach().cpu().numpy().flatten() for module in moduleList]
    posCoef=posCoefList[moduleId]
    posCoefArray=np.zeros((seqSize))
    posCoefArray[0:len(posCoef)]=posCoef
    return posCoefArray

def getPosBias(mnn:mnnPseudoModel.Net, moduleId)->np.ScalarType: #A scalar (float) of bias
    """
    Get the position bias for a given module.

    Parameters
    ----------
    mnn : mnnPseudoModel.Net
        The mnn model.
    moduleId : int
        The module id.

    Returns
    -------
    np.ScalarType
        The position bias for the module.

    """
    moduleList=mnnPseudoModel.getBlockList(mnn)
    posBiasList=[module.dense.bias.detach().cpu().numpy() for module in moduleList]
    posBiasArray=np.asarray(posBiasList).flatten()
    return posBiasArray[moduleId]

def getModuleWeight(
    mnn:mnnPseudoModel.Net,
    moduleId:int
)->np.ScalarType: #a scalar (float) of weight (module)
    """
    Get the module weight for a given module.

    Parameters
    ----------
    mnn : mnnPseudoModel.Net
        The mnn model.
    moduleId : int
        The module id.

    Returns
    -------
    np.ScalarType
        The module weight for the module.

    """
    return mnn.linear.weight.detach().cpu().numpy().flatten()[moduleId]

def getMeanPosActivationScore(
    mnnResultsArray:np.ndarray,
    mnn:mnnPseudoModel.Net,
    moduleId:int,
    threshold:float=0,
    bias:bool=False,
    seqSize:int=101
)->np.ndarray: #1D array of mean score (pos)
    """
    Get the mean position activation score for a given module.

    Parameters
    ----------
    mnnResultsArray : np.ndarray
        The mnn results array.
    mnn : mnnPseudoModel.Net
        The mnn model.
    moduleId : int
        The module id.
    threshold : float, optional
        The threshold (ReLU), by default 0
    bias : bool, optional
        If True, add the bias, by default False

    Returns
    -------
    np.ndarray
        The mean position activation score for the module.

    """
    x=mnnResultsArray[moduleId]
    #apply ReLU
    x[x<threshold]=0
    #mean over the sequence
    x=np.mean(x, axis=0)
    #apply position coefficient
    posCoefArray=getPosCoefArray(mnn, moduleId)
    x*=posCoefArray
    #apply bias if needed 
    if bias:
        posBias=getPosBias(mnn, moduleId)
        #add the bias, normalize by the sequence size to spread the bias all over the position
        x+=posBias/seqSize
    # apply the module weight
    moduleWeight=getModuleWeight(mnn, moduleId)
    x*=moduleWeight
    return x

def drawAxe(
    ax:plt.Axes,
    resultsArray:np.ndarray,
    title:str=None,
    xlabel:str="Distance to STR 3' end",
    ylabel:str=None
)-> plt.Axes:
    """
    Draw an histogram with seaborn.

    Parameters
    ----------
    ax : plt.Axes
        The axe to draw.
    resultsArray : np.ndarray
        The results array. 1D array of score (pos).
    title : str, optional
        The title, by default None
    xlabel : str, optional
        The xlabel, by default "Distance to STR 3' end"
    ylabel : str, optional
        The ylabel, by default None
    """
    nbPos=len(resultsArray)
    middlePos=nbPos//2
    xValues=np.arange(-middlePos, middlePos+1) #[-50,51[
    sns.histplot(x=xValues, weights=resultsArray, discrete=True, ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is not None :
        ax.set_title(title)
    return ax

def main():
    parser = argparse.ArgumentParser(description='Plot the MNN score for a given module.')
    parser.add_argument('--mnnResultsArray', type=pathlib.Path, required=True, help='Path to the mnn results array.')
    parser.add_argument('--moduleId', type=int, required=True, help='The module ID.')
    parser.add_argument('--mnnHParams', type=pathlib.Path, required=True, help='Path to the mnn hyperparameters.')
    parser.add_argument('--mnnParams', type=pathlib.Path, required=True, help='Path to the mnn parameters.')
    parser.add_argument('--bias', action='store_true', help='Add the bias to the score.')
    parser.add_argument('--fig', type=pathlib.Path, default=None, help='Output figure path.')
    parser.add_argument('--values', type=pathlib.Path, default=None, help='Output values path.')
    args = parser.parse_args()

    ylabel="Mean of the positional activation score"
    mnnResultsArray=np.load(args.mnnResultsArray)
    mnn=loadMnnModel(args.mnnHParams, args.mnnParams)
    resultsArray=getMeanPosActivationScore(mnnResultsArray, mnn, args.moduleId, bias=args.bias)
    if args.fig is not None:
        fig, ax = plt.subplots()
        drawAxe(ax, resultsArray, ylabel=ylabel)
        plt.savefig(args.fig)
    if args.values is not None:
        np.save(args.values, resultsArray)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pseudo Model for MNN

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
from typing import Union, Tuple, TypeVar, Callable, Any
PathLikeOrBuffer=typing.Union[os.PathLike, io.TextIOWrapper]

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

class BlockNet(nn.Module):
    def __init__(self, filter_size, size=101):
        super(BlockNet, self).__init__()
        self.conv = nn.Conv2d(1, 1, (filter_size, 4), bias=False)
        self.dense = nn.Linear(size - (filter_size - 1), 1)

    def forward(self, x):
        # x  dim (seq#, channel=1, letter, alphabet)
        x = self.conv(x)
        # -> dim (seq#, channel=1, letter-(filter_size-1),  alphabet-3=1)
        # x = F.relu(x)
        # -> dim (seq#, channel=1, letter-(filter_size-1),  alphabet-3=1)
        x = x.reshape(x.shape[0], x.shape[2])
        # -> dim (seq#,letter-(filter_size-1))
        # apply the factor of the dens layer and bias without the sum.
        # x= x * self.dense.weight + self.bias
        output = self.dense(x)
        # output dim (seq#, 1)
        # XXX : shunt the dense layer for result analysis purpose
        return x


class Net(nn.Module):
    def __init__(self, filter_size, block_type):
        super(Net, self).__init__()
        if block_type == "max":
            self.last_block = BlockMax(filter_size)
        else:
            self.last_block = BlockNet(filter_size)
        self.blocks = nn.ModuleList()
        self.len = 0
        self.linear = nn.Linear(self.len + 1, 1)

    def add_block(self, filter_size, block_type):
        self.last_block.require_grad = False
        self.blocks.append(self.last_block)
        if block_type == "max":
            self.last_block = BlockMax(filter_size)
        else:
            self.last_block = BlockNet(filter_size)
        self.len = len(self.blocks)
        self.linear = nn.Linear(self.len + 1, 1)

    def forward(self, x, block="last"):
        if block == "last":
            return self.last_block(x)
        else: 
            for i, l in enumerate(self.blocks):
                if i == block:
                    return self.blocks[i](x)
        raise KeyError(str(block) + " not in block list or out of range")



def build_modular(
        params:npt.ArrayLike
    )->Net:
    netmodel = Net(int(params[0][0]), 'net')
    for i in range(1, len(params)):
        netmodel.add_block(int(params[i][0]), 'net')
    return netmodel


def load_model(
    paramsPath:os.PathLike, 
    keysPath:os.PathLike
)->Net:
    params = np.load(paramsPath)
    model = build_modular(params)
    model.load_state_dict(torch.load(keysPath, map_location=torch.device('cpu')))
    return model

def getBlockList(model:Net) -> nn.ModuleList:
    """
    Get a list of blocks from the given model.

    Parameters
    ----------
    model : Net
        The neural network model.

    Returns
    -------
    nn.ModuleList
        The list of blocks in the model.
    """
    blockList=nn.ModuleList()
    blockList.extend(model.blocks)
    blockList.append(model.last_block)
    # i am ensuring that gradient is not computed. 
    for block in blockList :
        block.require_grad = False
    return blockList

def getFilterLengthList(blockList:nn.ModuleList)->list[int]:
    """
    Get a list of filter lengths from the given block list.

    Parameters
    ----------
    blockList : nn.ModuleList
        The list of blocks.

    Returns
    -------
    List[int]
        The list of filter lengths.
    """
    filterLengthList=[block.conv.kernel_size[0] for block in blockList]
    return filterLengthList
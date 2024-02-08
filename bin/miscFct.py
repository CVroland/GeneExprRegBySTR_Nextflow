#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Different useful tools.

Auteur : Mathys Grapotte, Christophe Vroland and Charles Lecellier
Date : 06/13/2023
"""

__authors__ = ("Mathys Grapotte", "Christophe Vroland", "Charles Lecellier")
__date__ = '06/13/2023'
__email__ = ("mathysgrapotte@gmail.com", 'christophe@cvroland.com', "charles.lecellier@igmm.cnrs.fr")
__status__ = 'Prototype'
__version__ = "0.0.1"

import numpy as np
import pandas as pd

import numpy.typing as npt
from typing import Sequence


def splitEStrHeader(seqHeaderArray:Sequence[str])->pd.DataFrame:
    seqHeaderDf=pd.Series(seqHeaderArray).str.split(pat=";", expand=True)
    seqHeaderDf=seqHeaderDf.rename(columns={0:"motifId",1:"seq",2:"strand"})
    return seqHeaderDf

COMPLEMENT_DNA_TABLE:bytes=bytes.maketrans("ACGTacgt".encode("ASCII"),"TGCAtgca".encode("ASCII"))
"""
COMPLEMENT_DNA_TABLE: bytes
    Translate table for complementing DNA sequences.
    
    This table is used to complement the DNA sequences, converting each nucleotide 
    to its complementary nucleotide (A to T, C to G, G to C, and T to A). The characters 
    'N' and 'n' are not changed.

    Notes
    -----
    The table is created using the `maketrans()` method from the `bytes` class.
"""

def complementDnaSeq(seq:str)->str:
    """
    Complements a DNA sequence.

    Parameters
    ----------
    seq : str
        The input DNA sequence.

    Returns
    -------
    str
        The complemented DNA sequence.

    """
    # encode the str into bytes with ASCII code. 
    # Use the translate table define below to transform a character by its complement. 
    # Decode the bytes into str with ASCII code. 
    return seq.encode("ASCII").translate(COMPLEMENT_DNA_TABLE).decode("ASCII")

def rcDnaSeq(seq:str)->str:
    """
    Generates the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        The input DNA sequence.

    Returns
    -------
    str
        The reverse complement of the DNA sequence.

    """
    # Reverse the complement sequence.
    return complementDnaSeq(seq)[::-1]
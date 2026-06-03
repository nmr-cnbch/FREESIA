"""
This module can be used for amino acid residues of a peptide or protein.

Functions:
    amino_acid_symbol(amino_acid)
    create_residue_list_from_sequence_file(sequence_file)

Classes:
    Residue: for residues of the peptide/protein.
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)


amino_acid_symbols = {
    "A": ("Ala", "alanine"),
    "C": ("Cys", "cysteine"),
    "D": ("Asp", "aspartate"),
    "E": ("Glu", "glutamate"),
    "F": ("Phe", "phenylalanine"),
    "G": ("Gly", "glycine"),
    "H": ("His", "histidine"),
    "I": ("Ile", "isoleucine"),
    "L": ("Leu", "leucine"),
    "K": ("Lys", "lysine"),
    "M": ("Met", "methionine"),
    "N": ("Asn", "asparagine"),
    "P": ("Pro", "proline"),
    "Q": ("Gln", "glutamine"),
    "R": ("Arg", "arginine"),
    "S": ("Ser", "serine"),
    "T": ("Thr", "threonine"),
    "V": ("Val", "valine"),
    "W": ("Trp", "tryptophan"),
    "Y": ("Tyr", "tyrosine"),
}


def amino_acid_symbol(amino_acid: str) -> tuple[str, str, str]:
    """
    Return one letter symbol, three letter symbol and name of an amino acid.

    Required arguments:
        amino_acid (str): Name or symbol of an amino acid.
    
    Returns:
        A tuple (one_letter, three_letters, full_name)

    ToDo:
        *add logging and raising errors when name not found
        *change the output to only return one value,
            based on the input parameter (0, 1 or 3; or default None
            when all three would be the output - the tuple)
    """
    found = False
    if len(amino_acid) == 1:
        if amino_acid in amino_acid_symbols:
            found = True
            one_letter = amino_acid
            three_letters, full_name = amino_acid_symbols[amino_acid]
    else:
        for key, value in amino_acid_symbols.items():
            if amino_acid in value:
                found = True
                one_letter = key
                three_letters, full_name = amino_acid_symbols[key]
    if found == False:
        #TODO Error!
        pass
    return one_letter, three_letters, full_name


def create_residue_list_from_sequence_file(
    sequence_file: os.PathLike
    ) -> list[Residue]:
    """
    Create a list of 'Residue' class objects from the protein sequence file.

    The protein sequence file has to contain only an ordered sequence
    of one letter symbols of subsequent amino acid residues.
    It can include any number of line breaks and whitespaces.

    Required arguments:
        sequence_file (os.PathLike): Path of the file with the sequence.

    Returns:
        list[Residue]: List of 'Residue' class objects with their number
            and one letter amino acid symbol.

    Raises:
        TypeError: description of the error

    ToDo:
        *add raising errors?
    """
    
    sequence = ""
    residue_list = []
    with open(sequence_file) as file:
        lines = file.readlines()
    logger.info(f"Finished reading the sequence file: {sequence_file}")
    for line in lines:
        line = line.strip().replace(" ","")
        sequence += line
    for index, char in enumerate(sequence):
        residue_list.append(Residue(index + 1, char))
    return residue_list



class Residue:
    """
    The amino acid residue of a peptide or a protein.

    Attributes:
        number (int): The number of the residue within the sequence
            of the peptide.
        type (str): The type of the residue, preferrably full name(?)
        resonances (dict[str:float]): Dictionary that assigns a chemical shift
            to a specific nucleus of the residue (in ppm).

    Methods:
        __init__(number, amino_acid_type)
        create_residue_list_from_sequence_file(sequence_file)
    """
    

    def __init__(self,
        number: int,
        amino_acid_type: str,
        resonances: dict[str:float] = None,
        ) -> None:
        """
        Assign the number and type of the residue.

        Required arguments:
            number (int): Number of the residue within the peptide sequence.
            amino_acid_type (str): Type of the amino acid, its name or symbol.
        
        ToDo:
            *?
        """
        self.number = number # for example 56
        self.type = amino_acid_type # "G" or "A" etc.
        if resonances != None:
            self.resonances = resonances
        else:
            self.resonances = {} # "CO: 175.6; HN: 114.3"
        logger.debug(f"Created Residue object with number: {self.number} and amino acid type: {self.type} ({self.number}{self.type})")

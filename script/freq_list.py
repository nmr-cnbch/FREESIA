"""
This module can be used to generate or read the list of frequencies.

It can be use on its own to generate a specific list of frequencies.
"""
#TODO: add documentation to the functions; add raising errors instead of exit
import logging
import os
import argparse
from typing import Any

logger = logging.getLogger(__name__)



def read_freq_list_file(file_path: os.PathLike):
    freqs = []
    with open(file_path) as file:
        lines = file.readlines()
    if lines[0].strip() != "bf ppm":
        exit("Unknown format of the frequency list file")
    for line in lines:
        line = line.strip()
        if line == "bf ppm":
            continue
        freqs.append(float(line))
    return freqs



def generate_freq_list(
    count: int,
    file_path: os.PathLike = None,
    ppm_range_start: float = None,
    ppm_range_end: float = None,
    protein: str = None,
    ) -> list[float]:

    proteins_dict = {
        "ubiq": (168.5, 183.5),
        "MBP": (169.5, 183),
        "Mpro": (169.5, 181.7),
    }

    if ppm_range_start == None and ppm_range_end == None and protein == None:
        exit("Neither protein nor the range have been specified!")
    elif protein in proteins_dict:
        ppm_range_start, ppm_range_end = proteins_dict[protein]
        print(f"Ppm range has been selected for the protein: {protein}")
        print(f"Ppm range is: {ppm_range_start} - {ppm_range_end}")
    elif protein not in proteins_dict and protein != None:
        print(f"There is no preset range for the protein: {protein}")
    elif ppm_range_start == None or ppm_range_end == None:
        print(f"The range is not fully specified")

    freq_increment = (ppm_range_end - ppm_range_start) / (count - 1)
    #freq_increment = (end - begin) / (increments-1)
    file_contents = ["bf ppm"]
    freq_list = []
    for i in range(0, count):
        freq = ppm_range_start + freq_increment * i
        file_contents.append(str(freq))
        freq_list.append(freq)

    if file_path != None:
        with open(file_path, "w") as file:
            for j in range(0, len(file_contents)):
                file.write(file_contents[j])
                file.write("\n")
    return freq_list #file_contents[1:]



if __name__ == "__main__":


    # ***** COMMAND LINE ARGUMENTS *****
    
    parser = argparse.ArgumentParser(
        description="The program %(prog)s extracts chemical shift information \
            from accordion NMR experiments with decoupling frequency \
            incremented along the evolution time in one of the dimensions \
            (called the accordion dimension).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file_path',
        help="Path of the output file")
    parser.add_argument('count',
        type = int,
        help="Number of frequencies in the output")
    parser.add_argument('ppm_range_start',
        type = float, nargs = '?', default=None,
        help="The low limit of the frequencies range (in ppm)")
    parser.add_argument('ppm_range_end',
        type = float, nargs = '?', default=None,
        help="The high limit of the frequencies range (in ppm)")
    parser.add_argument('protein', nargs = '?', default=None,
        help="Name of the protein for which the range is to be established. \
            If the protein was included in the presets, the frequencies list \
            will correspond to the range of chemical shift of CO nuclei \
            of the protein.")
    args = parser.parse_args() #after this line the script demands the positional arguments to exist

    generate_freq_list(args.count, args.file_path, args.ppm_range_start,
        args.ppm_range_end, args.protein)
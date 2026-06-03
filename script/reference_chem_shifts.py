"""
This module can be used for extracting reference chemical shifts 
and comparing them with measured chemical shifts.

Functions:
    extract_reference_chemical_shifts(reference_file_path, peaklist, residue_list, coupling_CO)

Classes:
    -
"""

import logging
import os

import peaks
import residues

logger = logging.getLogger(__name__)

def extract_reference_chemical_shifts(
    reference_file_path: os.PathLike,
    peaklist: peaks.PeakList,
    residue_list: list[residues.Residue] = None,
    coupling_CO: str = None,
    ) -> dict[str:float]:
    """
    Extract reference chemical shifts from the file.
    
    Adds 'acc_ref_shift' attribute to the peak objects from the peaklist.
    
    Required arguments:
        reference_file_path (os.PathLike): Path of the reference file.
            It should be a two-column text file with the name of the peak
            in the first column and its reference chemical shift in the second.
            Alternatively, for CO detection in proteins it can be
            a HNCO peaklist file in the Sparky format, the file type is detected
            based on the name (the name has to end with 'HNCO.list').
            If HNCO peaklist file is used, then optional arguments
            'residue_list' and 'coupling_CO' must be provided!
        peaklist (peaks.PeakList): Object of the 'PeakList' class with the list
            of peaks for which the profiles and chemical shifts were obtained.
    
    Optional arguments:
        residue_list (list[residues.Residue]): List of objects of 'Residue'
            class, containing all residues of the protein. Defaults to None.
            Required only if 'HNCO.list' is passed as the reference file.
            Can be created using function residues.create_residue_list_from_sequence_file.
        coupling_CO (str): Type of J-coupling which was used in the protein
            experiment to create the coherence with CO magnetization before
            decoupling. Can be either 'CACO' or 'NCO'. Defaults to None.
            Required only if HNCO peaklist is passed as the reference file.

    Returns:
        reference_shifts (dict[str:float]): Dictionary with reference chemical
            shifts assigned to the names of the peaks.

    Raises:
        TypeError: If HNCO peaklist file was passed but the residue_list
            or coupling_CO arguments are missing.
        ValueError: If the value of the coupling_CO argument is incorrect.
    """
    reference_shifts = {}
    allowed_couplings = ["CACO", "NCO"]
    
    if os.path.basename(reference_file_path).endswith("HNCO.list"):
        if (residue_list == None or coupling_CO == None):
            raise TypeError("HNCO peaklist file was provided for reference "
                "chemical shifts but the residue_list or coupling_CO "
                "arguments are missing!")
        
        elif coupling_CO not in allowed_couplings:
            raise ValueError("Unrecognized type of coupling was passed: "
                f"'{coupling_CO}'. Must be 'CACO' or 'NCO'.")
    
    if os.path.basename(reference_file_path).endswith("HNCO.list") == False:
        with open(reference_file_path) as file:
            lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line == "":
                continue
            words = line.split()
            peak_name = words[0]
            if words[1] == "None" or words[1] == "unknown":
                reference_chemical_shift = None
            else:
                reference_chemical_shift = float(words[1])
            reference_shifts[peak_name] = reference_chemical_shift
        
        for i, peak in enumerate(peaklist.peaklist):
            ref_shift = reference_shifts.get(peak.name)
            if ref_shift != None:
                peak.acc_ref_shift = ref_shift
                log_msg = (f"Reference chem. shift of {peak.acc_ref_shift} "
                    f"has been assigned to the peak {peak.name}.")
                logger.info(log_msg)
            else:
                peak.acc_ref_shift = "unknown"
                log_msg2 = (f"Reference chem. shift of peak {peak.name} "
                    "was not found.")
                logger.info(log_msg2)
            
            
    else:
        HNCO_peaklist = peaks.PeakList()
        HNCO_peaklist.create_from_peaklist_sparky(reference_file_path)
        HNCO_resonance_list = peaks.ResonanceList()
        HNCO_resonance_list.create_from_peaklist(HNCO_peaklist)
        
        for resonance in HNCO_resonance_list.resonances:
            residue_number = resonance.residue_number
            if resonance.nucleus_type not in residue_list[residue_number - 1].resonances:
                residue_list[residue_number - 1].resonances[resonance.nucleus_type] = resonance.chem_shift
                #print("CHECKIF: ",residue_number, residue_list[residue_number - 1].resonances, resonance.nucleus_type, resonance.chem_shift)
            #print("CHECK: ",residue_number, residue_list[residue_number - 1].resonances, resonance.nucleus_type, resonance.chem_shift)
    
        for peak in peaklist.peaklist:
            peak_resonances = peaks.ResonanceList()
            peak_resonances.create_from_peak(peak)
            
            ref_shift = None
            if coupling_CO == "NCO":
                for resonance in peak_resonances.resonances:
                    if resonance.nucleus_type == "N":
                        N_resonance = resonance
                for resonance in HNCO_resonance_list.resonances:
                    if resonance.residue_number == N_resonance.residue_number - 1 and resonance.nucleus_type == "CO":
                        CO_resonance = resonance
                        ref_shift = CO_resonance.chem_shift
            
            elif coupling_CO == "CACO":
                for resonance in peak_resonances.resonances:
                    if resonance.nucleus_type == "CA":
                        CA_resonance = resonance
                for resonance in HNCO_resonance_list.resonances:
                    if resonance.residue_number == CA_resonance.residue_number and resonance.nucleus_type == "CO":
                        CO_resonance = resonance
                        ref_shift = CO_resonance.chem_shift
            
            reference_shifts[peak.name] = ref_shift
            if ref_shift != None:
                peak.acc_ref_shift = ref_shift
                log_msg = (f"Reference chem. shift of {peak.acc_ref_shift} "
                    f"has been assigned to the peak {peak.name}.")
                logger.info(log_msg)
            else:
                peak.acc_ref_shift = "unknown"
                log_msg2 = (f"Reference chem. shift of peak {peak.name} "
                    "was not found.")
                logger.info(log_msg2)
            
    return reference_shifts



if __name__ == "__main__":
    HNCO_peaklist_path = "/net/synology/volume1/nhome/home/bartek/Documents/Doktorat/Python_scripts/classes/testinput_110/MBP_deut_2025_05_exp41_HNCO.list"
    sequence_file_path = "/net/synology/volume1/nhome/home/bartek/Documents/Doktorat/Python_scripts/classes/testinput_110/Protein_sequence"
    peaklist_path = "/net/synology/volume1/nhome/home/bartek/Documents/Doktorat/Python_scripts/classes/testinput_110/peak.list"
    ref_peaklist_path = "/net/synology/volume1/nhome/home/bartek/Documents/Doktorat/Python_scripts/classes/testinput_110/ref_peaklist.list"

    residue_list = residues.create_residue_list_from_sequence_file(sequence_file_path)

    peaklist = peaks.PeakList()
    peaklist.create_from_peaklist_sparky(peaklist_path)
    #ref_shifts = extract_reference_chemical_shifts(ref_peaklist_path, peaklist)
    ref_shifts = extract_reference_chemical_shifts(HNCO_peaklist_path, peaklist, residue_list, "CACO")

    output_file = "/net/synology/volume1/nhome/home/bartek/Documents/Doktorat/Python_scripts/classes/testinput_110/new_ref_peaklist.list"
    with open(output_file, "w") as file:
        #file.write("Assignment   w1 \n")
        for key, item in ref_shifts.items():
            file.write(f"{key}   {item}\n")

    #ref_shifts = extract_reference_chemical_shifts(HNCO_peaklist_path, peaklist, residue_list,"CCO")
    #ref_shifts = extract_reference_chemical_shifts(HNCO_peaklist_path, peaklist, residue_list, "CCO")
    print(ref_shifts)
    for peak in peaklist.peaklist:
        print(peak.name, peak.acc_ref_shift)
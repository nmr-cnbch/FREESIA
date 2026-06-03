"""
The program freesia.py extracts chemical shift information from FREESIA 
NMR experiments with irradiation pulse frequency co-incremented with 
the evolution time in one of the indirect dimensions 
(called the accordion dimension).
"""
import argparse
import logging
import os
import sys
import time
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

import cross_sections
import figures
import freq_list
#import parameters
import peaks
import reference_chem_shifts
import residues
import spectra
import weighting



logger = logging.getLogger(__name__)
#TODO: There must be a better way to do this
# to suppress debug logs from imported library modules:
logging.getLogger('matplotlib').setLevel(logging.WARNING)



def check_input_file_path(
    input_file_path: os.PathLike,
    file_type_name: str = None,
    ) -> os.PathLike:
    """
    Check if a path leads to a file and return its absolute path.

    Required arguments:
        input_file_path (os.PathLike): Path of the input file.
    
    Optional arguments:
        file_type_name (str): Name of the type of the input file.
            Affects only the text of the potential error message.
            Defaults to None.

    Returns:
        abs_input_file_path (os.PathLike): Absolute path of the input
            file.

    Raises:
        FileNotFoundError: If the 'input_file_path' does not lead to
            any file.
    """
    if input_file_path is None:
        return None
    else:
        abs_input_file_path = os.path.abspath(input_file_path)
    
    if not os.path.isfile(abs_input_file_path):
        if file_type_name is not None:
            error_msg = (f"This path was searched for the {file_type_name} "
                f"file but no file was found: {abs_input_file_path}")
        else:
            error_msg = (f"This path was searched for an input file but "
                f"no file was found: {abs_input_file_path}")

        logger.critical(error_msg)
        raise FileNotFoundError(error_msg)
    return abs_input_file_path



def process_cross_section(
    cross_section: cross_sections.CrossSection,
    processing_method: tuple[str, Any],
    weighting_values: list[float] = None,
    #weighting_values_dict: dict[tuple[str, Any]: list[float]] = None,
    processing_methods_dict: dict[tuple[str, Any]: str] = None,
    ) -> cross_sections.CrossSection:
    """
    Process the cross section according to the provided method.

    Required arguments:
        cross_section (cross_sections.CrossSection): 
            Cross section to process.
        processing_method (tuple[str, Any]): 
            Method in which the cross section will be processed.
            First element of the tuple indicates the method type,
            the second element is its appropriate parameter
            (such as the name assigned to a weighting function).
    
    Optional arguments:
        weighting_values (list[float]): List of values by which every 
            value of the cross section gets multiplied. Has to have
            the same length as the values of the cross section.
        processing_method_str (str): Text that could be used in output
            or error messages.
        
        weighting_values_dict: dict[tuple[str, Any]: list[float]]:
            Dictionary with the lists of weighting values assigned to
            the processing method. Required if the cross section will be
            multiplied by a weighting function. Defaults to None.
        processing_methods_dict: dict[tuple[str, Any]: str]:
            Dictionary which assigns a text to every processing method,
            the text may be used in the output or error messages.
            
    Returns:
        cross_section (cross_sections.CrossSection): Processed cross section.

    Raises:
        TypeError: If no list with weighting values was provided
            despite weighting being the selected processing method.
        ValueError: If the provided processing method is unknown.
        
        TypeError: If no dictionary with weighting values was provided
            despite weighting being the selected processing method.
        ValueError: If the selected processing method
    """
    if processing_method[0] == "s": #TODO: add somewhere a way of processing mode=r
        corr_peak_index, _ = cross_section.correct_peak_pos()
        cross_section.center_on_peak(peak, corr_peak_index)
        #cross_section.center_on_peak(peak)
    
    elif processing_method[0] == "trim":
        corr_peak_index, _ = cross_section.correct_peak_pos()
        cross_section.center_on_peak(peak, corr_peak_index)
        #cross_section.center_on_peak(peak)
        trimming_value = int(processing_method[1])
        cross_section.trim(trimming_value)
    
    elif processing_method[0] == "weighting":
        if weighting_values is not None:
            corr_peak_index, _ = cross_section.correct_peak_pos()
            cross_section.center_on_peak(peak, corr_peak_index)
            #cross_section.center_on_peak(peak)
            #weighting_values = weighting_values_dict[processing_method]
            cross_section.new_weighting(weighting_values)
        else:
            raise TypeError("No list with weighting values was provided")
    
    else:
        raise ValueError("Unknown processing method")
        
    return cross_section



if __name__ == "__main__":
    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    #default_input_dir = os.path.join(parent_dir, "input")
    #default_output_dir = os.path.join(parent_dir, "output")
    default_weighted_spectrum = "weighted_spectrum.ucsf"
    
    

    # ***** COMMAND LINE ARGUMENTS *****
    
    parser = argparse.ArgumentParser(
        description=("The program %(prog)s extracts chemical shift "
            "information from accordion NMR experiments with "
            "irradiation pulse frequency co-incremented with "
            "the evolution time in one of the indirect dimensions "
            "(called the accordion dimension)."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # *** Positional arguments: ***
    parser.add_argument('accordion_dim', type=int, #default="1",
        help=("Number of the indirect dimension whose evolution time was "
            "co-incremented with the irradiation frequency.")
        )
    parser.add_argument('input_dir', #nargs = '?', #default=default_input_dir,
        help="Path of the directory which contains input files."
        )
    parser.add_argument('output_dir', #nargs = '?', default=default_output_dir,
        help=("Path of the directory where output files will be saved. "
            "The directory will be created if it does not exist. If it "
            "already exists, the files in the directory may be overwritten!")
        )
    
    # *** Optional arguments: ***
    parser.add_argument('-w', "--weighting", nargs = '*', default="g[5]",
        help=("Weighting function by which the cross sections get "
            "multiplied, written as function[parameter], for example "
            "'c[4]'. Four types of functions are supported: c, s, e, g. "
            "They are defined in Section 3 of the Supporting Information of "
            "the publication (doi: ). Multiple functions can be used, "
            "each function will be used on the data separately "
            "and more subfolders will be created with output "
            "for every used weighting function."),
        )
    parser.add_argument('-spectrum', "--spectrum_files", nargs = '*',
        default="ucsf", 
        help=("Name or path of the spectrum file. Use this argument "
            "if you wish to use a different file name or directory "
            "than default (or want to use the sum of multiple spectra). "
            "Spectrum used should have been created from data which was "
            "submitted to apodization in every dimension except the accordion "
            "dimension. Apodization in accordion dimension will lead to worse "
            "or incorrect results. Multiple spectra can be provided, "
            "they will be added together, but they all must have exactly "
            "the same parameters and numbers of points in each dimension - "
            "otherwise the data in the output might be incorrect!"),
        metavar='spectrum_file_name'
        )
    parser.add_argument('-fws', "--fully_weighted_spectrum_files",
        nargs = '*', default=default_weighted_spectrum,#"fully_weighted_spectrum_real", 
        help=("Name or path of the fully weighted spectrum file. Use this "
            "argument if you wish to use a different file name or directory "
            "than default (or want to use the sum of multiple spectra). "
            "This spectrum must have been created from data which was "
            "submitted to apodization in every dimension, including "
            "the accordion dimension. Multiple spectra can be provided, "
            "they will be added together, but they all must have exactly "
            "the same parameters and numbers of points in each dimension."),
            #"If not provided, the main spectrum will be used to search for "
            #"the best cross section and to validate obtained results, but it "
            #"may lead to worse results."),
        metavar='fully_weighted_spectrum_real_file_name'
        )
    parser.add_argument('-plf', "--peak_list_file", default="peak.list",
        help=("Name or path of the peak list file. Use this argument "
            "if you wish to use a different file name or directory "
            "than default."), 
        metavar='peak_list_file_name'
        )
    parser.add_argument('-flf', "--freq_list_file", 
        default="fq1list", 
        help=("Name or path of the frequency list file. Use this argument "
            "if you wish to use a different file name or directory "
            "than default."), 
        metavar='frequency_list_file_name'
        )
    parser.add_argument('-tol', "--tolerance", type=float, 
        default="0.25",
        help=("Acceptable tolerance (chemical shift difference in ppm) "
            "between the detected chemical shift and the reference "
            "chemical shift. Requires the --reference argument with "
            "the appropriate file.")
        )
    parser.add_argument('-ref', "--reference_file", 
        help=("Name or path of the file with reference chemical shifts "
            "which will be used to compare the detected chemical "
            "shifts. Detected shifts are interpreted as correct "
            "if they match the reference values within the tolerance "
            "passed in the --tolerance argument (0.25 ppm by default). "
            "It should be a two-column text file with the name of "
            "the peak in the first column and its reference chemical "
            "shift in the second (columns separated by white spaces). "
            #"Alternatively, for CO detection "
            #"in proteins it can be a HNCO peaklist file in the Sparky "
            #"format, the file type is detected based on the name "
            #"(the name has to end with 'HNCO.list'). If HNCO peaklist "
            #"is used, then optional arguments --sequence_file "
            #"and --coupling_CO must be provided!"),
            ),
        metavar="reference_file_name"
        )
    parser.add_argument("--latex", action='store_true',
        help=("Create a text file that will include a latex code for "
            "a table containing main results of the script: peak names, "
            "their positions in the spectrum (chemical shifts), "
            "chemical shifts determined from the amplitude profiles for "
            "the accordion dimension, reference chemical shifts "
            "(if -ref argument was provided) and the assigned profile group")
        )
    parser.add_argument("--figures",
        help=("Available options: individual - a plot of every cross section "
            "and profile is saved into a separate png file; "
            "combined - plots for profiles are combined into "
            "fewer figures, one figure having 35 profile plots; "
            "both - both individual and combined options are used"), 
        choices = ['individual', 'combined', 'both']
        )
    parser.add_argument("--text",
        help=("Available options: cs - values for cross sections are saved "
            "as text files; profiles - values for profiles are saved as text "
            "files; both - values for both cross sections and profiles are "
            "saved as text files"),
        choices = ['cs', 'profiles', 'both']
        )
    
    parser.add_argument('-v', '--verbosity', action='count', default=0,
        help='Increase verbosity level for the logging'
        )
    #TODO: if there are many more arguments, 
    #       might be wise to use parser.add_argument_group() ?

    #HIDDEN: mode argument is currently unused, 
    # TODO: allow for s and r?;; default None?
    parser.add_argument('-m', "--mode", default='s',
        help=argparse.SUPPRESS,#("Determines how the cross sections of the spectrum get "
            #"processed before the inverse Fourier transformation: "
            #"r - raw cross sections, "
            #"t - trimmed, "
            #"w - weighted, "
            #"s - shifted"),
        choices=['r', 't', 'w', 's']
        )
    
    #HIDDEN: trimming is not described in paper, it is not recommended
    parser.add_argument('-t', "--trimming", type = int, nargs = '*',
        default=None, help=argparse.SUPPRESS
        # help = ("Number of points remaining within the cross sections "
            #"after trimming them around the peak position")
        )
    
    #HIDDEN: use of parameters file is no longer needed as we use ucsf
    # and mp_acc is taken as the number of frequencies in freq_list
    parser.add_argument('-pf', "--parameters_file", 
        default=None,#"parameters.txt",
        help=argparse.SUPPRESS,#("Name or path of the parameters file. Use this argument "
            #"if you wish to use a different file name or directory "
            #"than default."),
        metavar='parameters_file_name'
        )
    
    #HIDDEN: sequence file is not needed as it is recommended to use reference
    # shifts file instead of extracting values from HNCO
    parser.add_argument('-seq', "--sequence_file",
        help=argparse.SUPPRESS,#"Name or path of the protein sequence file",
        metavar='sequence_file_name'
        )
    
    #HIDDEN: coupling info not needed, recommended to use reference shifts file
    # instead of extracting values from experiment peaklist
    parser.add_argument("--coupling_CO", 
        help=argparse.SUPPRESS,
            #("Type of J-coupling which was used in the protein "
            #"experiment to create the coherence with CO magnetization "
            #"before decoupling. Can be either 'CACO' or 'NCO'. "
            #"Required only if HNCO peaklist is passed as "
            #"the reference file."),
        choices=['CACO', 'NCO']
        )



    # All parser arguments will exist now:
    args = parser.parse_args()
    
    
        
    # ***** MAIN OUTPUT DIRECTORY *****
    output_dir = os.path.abspath(args.output_dir)
    
    main_output_textfile_name = "FREESIA_output_summary.txt"
    main_output_textfile_path = os.path.join(
        output_dir, main_output_textfile_name)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    '''else: #code that asks the user if they want to save data in existing directory
        user_input_msg = (f"WARNING: The output directory {output_dir} "
            "already exists and it could contain data from previous uses "
            "of the program. If you continue, the data may be overwritten. "
            "Do you want to continue? (Y/N) ")
        overwrite_decision = input(user_input_msg)
        if overwrite_decision == "Y" or overwrite_decision == "y":
            print("\nYou chose to continue...\n")
        elif overwrite_decision == "N" or overwrite_decision == "n":
            exit("Stopped the program!")
        else:
            exit("Stopped the program!")'''

        
        
    # ***** LOGGING *****
    logging.addLevelName(999, "BASIC")
    log_file_path = os.path.join(output_dir, "log")

    if args.verbosity > 0:
        log_level = max(logging.DEBUG, logging.WARNING - args.verbosity * 10)
    else:
        log_level = logging.WARNING
    
    if log_level == logging.DEBUG:
        log_format = (
            "[%(asctime)s]   %(levelname)s: %(message)s"
            "   (file: %(filename)s, function: %(funcName)s, "
            "line: %(lineno)d, time: %(relativeCreated)d ms)")
    else:
        log_format = (
            "[%(asctime)s]   %(levelname)s: %(message)s")

    """Logging levels from lowest to highest:
        debug, info, warning, error, critical, (basic)"""
    logging.basicConfig(
        filename = log_file_path,
        encoding = "utf-8",
        level = log_level,
        format = log_format
        )

    with open(log_file_path, "a") as log_file:
        log_file.write("-"*50+"\n")
        
    script_path = os.path.abspath(os.path.realpath(__file__))
    used_command = f"{' '.join(sys.argv)}"
    log_script_use = ("The program was started using the following command:\n"
                      f"{used_command}\n"
                      f"The program's path is: {script_path}")
    logger.log(999, log_script_use)    



    # ***** INPUT FILES *****
    input_paths = {
        "NMR spectrum": [],
        "Weighted NMR spectrum": [],
        "peak list": os.path.join(args.input_dir, args.peak_list_file),
        #"parameters": os.path.join(args.input_dir, args.parameters_file),
        "freq list": os.path.join(args.input_dir, args.freq_list_file),
        }
    optional_files = {
        "protein sequence": args.sequence_file,
        "reference shifts": args.reference_file,
        "parameters": args.parameters_file,
        }
    
    for file_type, file_name in optional_files.items():
        if file_name is not None:
            input_paths[file_type] = os.path.join(args.input_dir, file_name)
    
    if isinstance(args.spectrum_files, list):
        for index, file_name in enumerate(args.spectrum_files):
            input_paths["NMR spectrum"].append(
                os.path.join(args.input_dir, file_name))
    else:
        input_paths["NMR spectrum"].append(
            os.path.join(args.input_dir, args.spectrum_files))
    
    if isinstance(args.fully_weighted_spectrum_files, list):
        for index, file_name in enumerate(args.fully_weighted_spectrum_files):
            input_paths["Weighted NMR spectrum"].append(
                os.path.join(args.input_dir, file_name))
    else:
        input_paths["Weighted NMR spectrum"].append(
            os.path.join(args.input_dir, args.fully_weighted_spectrum_files))
    
    #TODO: how to implement other spectra formats, change in the future if
    # this option will be desired
    '''if input_paths["NMR spectrum"][0].endswith("ucsf"):
        spectrum_file_type = "ucsf"
        create_cs_func = cross_sections.CrossSection.create_from_ucsf
    else:
        spectrum_file_type = "spectrum_real"
        create_cs_func = cross_sections.CrossSection.create_from_spectrum_real'''
    
    # for now only ucsf format will be supported
    spectrum_file_type = "ucsf"
    create_cs_func = cross_sections.CrossSection.create_from_ucsf
    

    # *** Verification of the input paths: ***
    input_files_msg = "\nLOCATED INPUT FILES:\n"
    for file_type, file_path in input_paths.items():
        if isinstance(file_path, list):
            input_files_msg += f"{file_type}:\n"
            for i, path in enumerate(file_path):
                input_paths[file_type][i] = check_input_file_path(
                    path, file_type)
                input_files_msg += f" - {input_paths[file_type][i]}\n"
        else:
            input_paths[file_type] = check_input_file_path(
                file_path, file_type)
            input_files_msg += f"{file_type}:   {input_paths[file_type]}\n"
    
    print(input_files_msg)
    time.sleep(1) #TODO: remove? change time?
    logger.log(999, input_files_msg)


    weighted_spectra = input_paths["Weighted NMR spectrum"]


    # ***** PROCESSING METHODS *****
    #processing_methods_dict = {(args.mode, ""): f"{args.mode}"}
    processing_methods_dict = {}

    if args.trimming != None:
        for trim_value in args.trimming:
            processing_methods_dict[("trim", trim_value)] = \
                f"trim_{trim_value}"
    '''if args.weighting != "gauss[5]": #None:
        for weighting_function in args.weighting:
            if weighting_function == 'none':
                processing_methods_dict[("s", "")] = f"noweighting"
            else:
                func_str = weighting_function.replace("[","_")
                func_str = func_str.replace("]","_")
                processing_methods_dict[("weighting", weighting_function)] = \
                    f"weighting_{func_str}"#f"weighting_{weighting_function}"
    else:
        processing_methods_dict[("weighting", "gauss[5]")] = "weighting_gauss_5_"'''
    #if args.weighting != "gauss[5]": #None:
    weighting_functions_list = []
    if not isinstance(args.weighting, list):
        weighting_functions_list.append(args.weighting)
    else:
        for function in args.weighting:
            weighting_functions_list.append(function)
    
    for weighting_function in weighting_functions_list:
        if weighting_function == 'none':
            processing_methods_dict[("s", "")] = f"noweighting"
        else:
            func_str = weighting_function.replace("[","")
            func_str = func_str.replace("]","")
            processing_methods_dict[("weighting", weighting_function)] = \
                f"weighting_{func_str}"#f"weighting_{weighting_function}"
    #else:
    #    processing_methods_dict[("weighting", "gauss[5]")] = "weighting_gauss_5_"
    processing_methods = processing_methods_dict.keys()



    # ***** OUTPUT SUBDIRECTORIES *****
    #TODO: change this section so that subdirectories are created 
    #       only when there will be some output inside 
    #       (based on text/figures arguments); 
    #       first have to change the final output files, 
    #       to stop them from being created inside the subdirectories
    output_dir_cs = os.path.join(output_dir, "cross_sections_before_processing")
    if not os.path.exists(output_dir_cs):
        os.makedirs(output_dir_cs)

    output_dirs_methods = {}
    output_dirs_methods_cs = {}
    output_dirs_methods_profiles = {}

    for method in processing_methods:
        output_dir_method = os.path.join(
            output_dir, f"{processing_methods_dict[method]}")#f"profiles_{processing_methods_dict[method]}")
        output_dir_method_cs = os.path.join(
            output_dir_method, "cross_sections")#f"cross_sections_{processing_methods_dict[method]}")
        output_dir_method_profiles = os.path.join(
            output_dir_method, "profiles")
        
        output_dirs_methods[method] = output_dir_method
        output_dirs_methods_cs[method] = output_dir_method_cs
        output_dirs_methods_profiles[method] = output_dir_method_profiles
        if not os.path.exists(output_dir_method):
            os.makedirs(output_dir_method)
        if not os.path.exists(output_dir_method_cs):
            os.makedirs(output_dir_method_cs)
        if not os.path.exists(output_dir_method_profiles):
            os.makedirs(output_dir_method_profiles)



    # ***** PREPARATIONS *****
    start_time = time.time()
    accordion_dim = args.accordion_dim
    tolerance = args.tolerance


    # *** List of irradiation frequencies: ***
    freqs = freq_list.read_freq_list_file(input_paths["freq list"])


    # *** Parameters of the spectrum: ***
    #spectrum_pars = spectra.spectrum_params_types
    #params_types = parameters.ParamsTypes(spectrum_pars)
    #spectrum_params = spectra.SpectrumParams(params_types = params_types)

    if spectrum_file_type == "ucsf":
        spectrum = spectra.UCSF_Spectrum.read_ucsf_file(
            input_paths["NMR spectrum"][0])
        spectrum_params = spectrum.params
        dim = spectrum_params["dim"]
        mp_acc = len(freqs) #TODO: might have to be adjusted if different sampling is desired in the future
    elif spectrum_file_type == "spectrum_real":
        spectrum_params = spectra.SpectrumParams()
        spectrum_params.read_parameters_file(input_paths["parameters"])
        dim = spectrum_params.establish_dimensionality()
        spectrum_params.check_reffrq()
        spectrum_params.calc_center()
        spectrum_params.calc_measured_points()
        mp_acc = spectrum_params.get_value("mp", accordion_dim)
        spectrum = spectra.Spectrum(spectrum_params)


    # *** Peak list for the cross sections: ***
    peaklist_cs = peaks.PeakList()
    peaklist_cs.create_from_peaklist_sparky(input_paths["peak list"])


    # *** Preparing a list of correct/expected CO chemical shifts: ***
    if "reference shifts" not in input_paths:
        ref_comparison_flag = False
    else:
        reference_file_path = input_paths["reference shifts"]
        ref_comparison_flag = True
        if "protein sequence" in input_paths:
            residue_list = residues.create_residue_list_from_sequence_file(
                input_paths["protein sequence"])
        else:
            residue_list = None
        coupling_CO = args.coupling_CO 
        ref_shifts = reference_chem_shifts.extract_reference_chemical_shifts(
            reference_file_path, peaklist_cs, residue_list, coupling_CO)
        
        
    # *** Preparing the weighting functions: ***
    points_dir = []
    if spectrum_file_type == "ucsf":
        cs_len = spectrum_params.get_value("fnz", accordion_dim)
        pass
    else:
        for i in range(0, dim):
            if f"points_dir{i}" in spectrum_params.params:
                points_dir.append(spectrum_params.get_value("points_dir", i))
            else:
                points_dir.append(spectrum_params.get_value("fnz", i))
        cs_len = points_dir[accordion_dim]
    
    #if num == 0:
    weighting_values_dict = {}
    for method in processing_methods:
        if method[0] == "weighting":
            weighting_function = method[1]
            weighting_func, _ = weighting.select_weighting_function_from_str(
                weighting_function, cs_len
            )
            weighting_values = weighting.weight_func_vals(weighting_func, cs_len)
            weighting_values_dict[method] = weighting_values
    
    
    # *** Preparing the pandas table for all the data: ***
    table_columns = ["peak_name"]
    dims = [i for i in range(1, dim)]
    dims.append(0)
    for i in range(0, dim):
        dim_nucleus = spectrum_params.get_value("name", i)
        table_columns.append(f"chem_shift_{i}__{dim_nucleus}")
        
    if ref_comparison_flag:
        table_columns.append("ref_chem_shift")
    
    table_columns.append("cross_section")
    
    for method in processing_methods:
        method_str = processing_methods_dict[method]
        table_columns.append(f"{method_str}_cs")
        table_columns.append(f"{method_str}_profile")
        table_columns.append(f"{method_str}_chem_shift")
        if ref_comparison_flag:
            table_columns.append(f"{method_str}_correctness")
        #table_columns.append(f"{method_str}_group")
        table_columns.append(f"{method_str}_validity")

    results_table = pd.DataFrame(columns=table_columns)
    
    for peak in peaklist_cs.peaklist:
        coords = peak.calc_peak_coords(spectrum_params)
        #peak.establish_widths(spectrum, spectrum_real_file = input_paths["NMR spectrum"][0])
        #peak.establish_height(spectrum, spectrum_real_file = input_paths["NMR spectrum"][0])
        



    # ***** PROCESSING DATA *****
    print("Processing spectrum...")
    for num, peak in enumerate(peaklist_cs.peaklist):
        peaklist_len = len(peaklist_cs.peaklist)
        cur_progress = num / peaklist_len * 100
        cur_progress = round(cur_progress, 2)
        #print("", end='\r')
        #sys.stdout.write("\033[K")
        print("\033[K", end="") #ANSI escape code, clears the line
        print((f"Progress: {cur_progress:.2f}% - "
            f"processing peak {peak.name} (number {num+1} out of "
            f"{len(peaklist_cs.peaklist)})"), end='\r')
        results_table.at[num, "peak_name"] = peak.name
        if ref_comparison_flag:
            results_table.at[num, "ref_chem_shift"] = peak.acc_ref_shift

        for i in range(0, dim):
            dim_nucleus = spectrum_params.get_value("name", i)
            results_table.at[num, f"chem_shift_{i}__{dim_nucleus}"] = \
                peak.chem_shifts[i]


        # ***** EXTRACTING CROSS SECTION *****
        coords = peak.calc_peak_coords(spectrum_params)
        
        # *** Establishing the least overlapped cross section: ***
        #fully_weighted_spectrum_cs = cross_sections.CrossSection()
        #fully_weighted_spectrum_cs = cross_sections.CrossSection.create_from_spectrum_real(
        fully_weighted_spectrum_cs = cross_sections.CrossSection.create_from_any_spectrum(
            weighted_spectra[0], spectrum, accordion_dim, coords, peak.name)
        additional_weighted_spectrum_cs = []
        for spectrum_file in weighted_spectra[1:]: #TODO: change weighted_spectra to appropriate thing later
            #additional_cs = cross_sections.CrossSection()
            #additional_cs = cross_sections.CrossSection.create_from_spectrum_real(spectrum_file,
            additional_cs = cross_sections.CrossSection.create_from_any_spectrum(spectrum_file,
                spectrum, accordion_dim, coords, peak.name)
            additional_weighted_spectrum_cs.append(additional_cs)
        fully_weighted_spectrum_cs.add_cross_sections(
            additional_weighted_spectrum_cs)
        
        
        for method in processing_methods:
            method_str = processing_methods_dict[method]
            fws_processed_cs = cross_sections.CrossSection(fully_weighted_spectrum_cs)
            
            if method[0] == "weighting":
                #fws_processed_cs.center_on_peak(peak)
                weighting_function = method[1]
                weighting_values = weighting_values_dict[method]
            else:
                weighting_values = None
            fws_processed_cs = process_cross_section(fws_processed_cs, 
                method, weighting_values, method_str)
            
            # *** Handling all processing methods: ***
            """
            if method[0] == "s": #TODO: add somewhere a way of processing mode=r
                fws_processed_cs.center_on_peak(peak)
            
            elif method[0] == "trim":
                fws_processed_cs.center_on_peak(peak)
                trimming_value = int(method[1])
                fws_processed_cs.trim(trimming_value)
            
            elif method[0] == "weighting":
                fws_processed_cs.center_on_peak(peak)
                weighting_function = method[1]
                weighting_values = weighting_values_dict[method]
                fws_processed_cs.new_weighting(weighting_values)
            
            else:
                exit("Unknown processing method") #TODO: improve"""
            
            """overlap_results = fws_processed_cs.check_weighted_spectrum_overlap_v3()
            interfering_peaks = overlap_results[0]
            main_peak_index = overlap_results[1]
            overlap_factor = overlap_results[2]
            main_overlap_factor = overlap_results[3]
            # correct coords so that it matches the shifted cs
            if overlap_factor == 0:
                final_cs = fws_processed_cs
                print(f"No overlap on the initial peak: {main_overlap_factor}.")
            else:
                centered_coords = coords
                print(coords)
                #centered_coords[accordion_dim] = int(len(fws_processed_cs.values)/2)
                cs_to_compare_list = fws_processed_cs.find_neighbouring_cross_sections(
                    coords, weighted_spectra, spectrum)
                for compared_cs in cs_to_compare_list:
                    compared_cs = process_cross_section(compared_cs, 
                        method, weighting_values, method_str)
                    """
                    
            """
                    if method[0] == "s": #TODO: add somewhere a way of processing mode=r
                        compared_cs.center_on_peak(peak)
                    
                    elif method[0] == "trim":
                        compared_cs.center_on_peak(peak)
                        trimming_value = int(method[1])
                        compared_cs.trim(trimming_value)
                    
                    elif method[0] == "weighting":
                        compared_cs.center_on_peak(peak)
                        weighting_function = method[1]
                        weighting_values = weighting_values_dict[method]
                        compared_cs.new_weighting(weighting_values)
                    
                    else:
                        exit("Unknown processing method") #TODO: improve
                        """
            r"""
                    
                    
                main_factors = []
                for j, cs in enumerate(cs_to_compare_list):
                    interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = cs.check_weighted_spectrum_overlap_v3()
                    main_factors.append(abs(main_overlap_factor))
                    
                    output_dir_cs_overlap_check = os.path.join(
                        output_dirs_methods_cs[method], "overlap_check")
                    
                    if not os.path.exists(output_dir_cs_overlap_check):
                        os.makedirs(output_dir_cs_overlap_check)
                        
                    x_values = [i for i in range(0, len(cs.values))]
                    y_values = cs.values
                    #x_values = [i for i in range(0, len(fully_weighted_spectrum_cs.values))]
                    #y_values = fully_weighted_spectrum_cs.values
                    #x_values = [i for i in range(0, len(fws_processed_cs.values))]
                    #y_values = fws_processed_cs.values
                    
                    xlim = (min(x_values), max(x_values))
                    xticks = None #[len(x_values)/2]
                    xticks_labels = None #["0"]
                    yticks = [0]
                    cs_file_name = (f"CS_p{num}_{peak.name}__{j}")
                    figures.create_2Dplots(
                        list_of_x_values = [x_values],
                        list_of_y_values = [y_values],
                        fig_file_dir = output_dir_cs_overlap_check,
                        fig_file_basename = cs_file_name,
                        transparent = True,
                        xticks = xticks,
                        xticks_labels = xticks_labels,
                        yticks = yticks,
                        xlim = xlim,
                        xlabel = "$\omega_1$",
                        ylabel = "intensity",
                        list_of_scatter_x_values = [interfering_peaks],#[0],#[value_index],
                        list_of_scatter_y_values = [[cs.values[k] for k in interfering_peaks]],#[noise],#[value, noise],
                        plot_settings = figures.plot_settings_cross_section_for_paper
                    )
                    
                    
                    
                    
                lowest_overlap = min(main_factors)
                ind = main_factors.index(lowest_overlap)
                #cross_section_for_overlap = cs_to_compare_list[ind]
                final_cs = cs_to_compare_list[ind]
                interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = final_cs.check_weighted_spectrum_overlap_v3()
                print(f"Final peak, {method_str}: {main_overlap_factor}, {ind}, {lowest_overlap}, {main_factors}")
            """
            cs_to_compare_list = fws_processed_cs.find_neighbouring_cross_sections(
                coords, weighted_spectra, spectrum)
            for compared_cs in cs_to_compare_list:
                compared_cs = process_cross_section(compared_cs, 
                    method, weighting_values, method_str)
            final_cs = cross_sections.select_best_cross_section(cs_to_compare_list, peak)
            #final_cross_section = cross_sections_to_compare[ind]
            #if len(interfering_peaks) == 0:
            #    final_cross_section = cross_section_for_overlap
            del cs_to_compare_list
            final_coords = final_cs.slicing_point_coords
            #print(method_str, peak.name, "best_coords: ", final_coords)
            if num == 20:
                pass
                #temp_txt_path = os.path.join(output_dir, f"{method_str}_p20_cs_ift.txt")
                #cs_ift.create_text_file(temp_txt_path)
                #temp_txt_path = os.path.join(output_dir, f"{method_str}_p20_cross_section.txt")
                #cross_section.create_text_file(temp_txt_path)
                #temp_txt_path = os.path.join(output_dir, f"{method_str}_p20_final_cs.txt")
                #final_cs.create_text_file(temp_txt_path)
        
        
        
        
        
            #cross_section = cross_sections.CrossSection()
            #cross_section = cross_sections.CrossSection.create_from_spectrum_real(input_paths["NMR spectrum"][0],
            cross_section = cross_sections.CrossSection.create_from_any_spectrum(input_paths["NMR spectrum"][0],
                spectrum, accordion_dim, final_coords, peak.name)


            additional_cross_sections = []
            for spectrum_file in input_paths["NMR spectrum"][1:]:
                #additional_cs = cross_sections.CrossSection()
                #additional_cs = cross_sections.CrossSection.create_from_spectrum_real(spectrum_file,
                additional_cs = cross_sections.CrossSection.create_from_any_spectrum(spectrum_file,
                    spectrum, accordion_dim, final_coords, peak.name)
                additional_cross_sections.append(additional_cs)
            cross_section.add_cross_sections(additional_cross_sections)
            
            
            
            #interfering_peaks, interfering_peak_heights = cross_section.check_overlap()
            #print(f"{peak.name}: {interfering_peaks}, {interfering_peak_heights}")
            
            # *** Cross sections output: ***
            results_table.at[num, "cross_section"] = cross_section
            cs_file_name = (f"CS_p{num}_{peak.name}")
            
            # text output
            if args.text == 'cs' or args.text == 'both':
                cs_txt_name = cs_file_name + ".txt"
                output_dir_cs_text = os.path.join(output_dir_cs, "txt_files")
                if not os.path.exists(output_dir_cs_text):
                    os.makedirs(output_dir_cs_text)
                cs_txt_path = os.path.join(output_dir_cs_text, cs_txt_name)
                cross_section.create_text_file(cs_txt_path)
        
            # figures
            if args.figures == 'individual' or args.figures == 'both':
                output_dir_cs_figs = os.path.join(output_dir_cs, "png_files")
                if not os.path.exists(output_dir_cs_figs):
                    os.makedirs(output_dir_cs_figs)
                x_values = [i for i in range(0, len(cross_section.values))]
                y_values = cross_section.values
                xlim = (min(x_values), max(x_values))
                xticks = [len(x_values)/2]
                xticks_labels = ["0"]
                yticks = [0]
            
                figures.create_2Dplots(
                    list_of_x_values = [x_values],
                    list_of_y_values = [y_values],
                    fig_file_dir = output_dir_cs_figs,#output_dir_cs,
                    fig_file_basename = cs_file_name,
                    transparent = True,
                    xticks = xticks,
                    xticks_labels = xticks_labels,
                    yticks = yticks,
                    xlim = xlim,
                    xlabel = r"$\omega_1$",
                    ylabel = "intensity",
                    plot_settings = figures.plot_settings_cross_section_for_paper
                )
    
        
        
        
        
        # ***** PROCESSING CROSS SECTION *****
        #for method in processing_methods:
            method_str = processing_methods_dict[method]
            cs_ift = cross_sections.CrossSection(cross_section)
            
            
            if method[0] == "weighting":
                #fws_processed_cs.center_on_peak(peak)
                weighting_function = method[1]
                weighting_values = weighting_values_dict[method]
            else:
                weighting_values = None
            cs_ift = process_cross_section(cs_ift, 
                method, weighting_values, method_str)
            if num == 20:
                temp_txt_path = os.path.join(output_dir, f"{method_str}_p20_cs_ift.txt")
                #cs_ift.create_text_file(temp_txt_path)
            # *** Handling all processing methods: ***
            """if method[0] == "s": #TODO: add somewhere a way of processing mode=r
                cs_ift.center_on_peak(peak)
            
            elif method[0] == "trim":
                cs_ift.center_on_peak(peak)
                trimming_value = int(method[1])
                cs_ift.trim(trimming_value)
            
            elif method[0] == "weighting":
                cs_ift.center_on_peak(peak)
                weighting_function = method[1]
                weighting_values = weighting_values_dict[method]
                cs_ift.new_weighting(weighting_values)
            
            else:
                exit("Unknown processing method") #TODO: improve
            """
                
            #WIP
            #cross_sections.check_cs_overlap(spectrum, accordion_dim, coords, peak, num, output_dir_cs)
            #if not weighting_values:
            #    weighting_values = None
            
            #selected_cross_section = cross_sections.check_cs_overlap_v3(spectrum, accordion_dim, coords, peak, num, output_dir_cs, weighting_values)
            #new_coords = selected_cross_section.peak_coords
            """
            cross_section = cross_sections.CrossSection()
            cross_section.create_from_spectrum_real(input_paths["NMR spectrum"][0],
                spectrum, accordion_dim, new_coords, peak.name)
            additional_cross_sections = []
            for spectrum_file in input_paths["NMR spectrum"][1:]:
                additional_cs = cross_sections.CrossSection()
                additional_cs.create_from_spectrum_real(spectrum_file,
                    spectrum, accordion_dim, new_coords, peak.name)
                additional_cross_sections.append(additional_cs)
            cross_section.add_cross_sections(additional_cross_sections)
            """
            # *** Processed cross sections output: ***
            #results_table.at[num, f"{method_str}_cs"] = cs_ift
            
            output_dir_method = output_dirs_methods[method]
            output_dir_method_cs = output_dirs_methods_cs[method]
            
            cs_file_name = f"CS_p{num}_{peak.name}__{method_str}"

            # text output
            if args.text == 'cs' or args.text == 'both':
                cs_txt_name = cs_file_name + ".txt"
                output_dir_method_cs_text = os.path.join(output_dir_method_cs, "txt_files")
                if not os.path.exists(output_dir_method_cs_text):
                    os.makedirs(output_dir_method_cs_text)
                cs_txt_path = os.path.join(output_dir_method_cs_text, cs_txt_name)
                cs_ift.create_text_file(cs_txt_path)

            # figures
            if args.figures == 'individual' or args.figures == 'both':
                output_dir_method_cs_figs = os.path.join(output_dir_method_cs, "png_files")
                if not os.path.exists(output_dir_method_cs_figs):
                    os.makedirs(output_dir_method_cs_figs)
                cs_x_values = [i for i in range(0, len(cs_ift.values))]
                cs_y_values = cs_ift.values
                xlim = (min(cs_x_values), max(cs_x_values))
                xticks = [len(cs_x_values)/2]
                xticks_labels = ["0"]
                yticks = [0]
            
                figures.create_2Dplots(
                    list_of_x_values = [cs_x_values],
                    list_of_y_values = [cs_y_values],
                    fig_file_dir = output_dir_method_cs_figs,
                    fig_file_basename = cs_file_name,
                    transparent = True,
                    xticks = xticks,
                    xticks_labels = xticks_labels,
                    yticks = yticks,
                    xlim = xlim,
                    xlabel = r"$\omega_1$",
                    ylabel = "intensity",
                    plot_settings = figures.plot_settings_cross_section_for_paper
                )
            
            # *** Calculating profiles: ***
            cs_processed = cross_sections.CrossSection(cs_ift)
            results_table.at[num, f"{method_str}_cs"] = cs_processed
            
            zerofilling_factor = cs_ift.zerofilling(2)#TODO add zerofilling argument to argparse
            points_count_new = int(zerofilling_factor * mp_acc)  #-2 is best; TODO: log previous (mp_acc) and new point counts?
            
            cs_length = int(len(cs_ift.values)) #TODO: unused, remove?
            accdim_shifts_amount = int(len(cs_ift.values)/2) #TODO: unused, remove?
            
            ift = cs_ift.inverse_Fourier_transformation()
            
            accdim_frequencies = freq_list.generate_freq_list(points_count_new,
                ppm_range_start = freqs[0], ppm_range_end = freqs[-1],)

            profile = cross_sections.Profile(cs_ift, accdim_frequencies, points_count_new)
            
            #profile.rescale_percentage()

            
            # *** Profiles output: ***
            results_table.at[num, f"{method_str}_profile"] = profile
            profile_file_name = f"Profile_p{num}_{peak.name}__{method_str}"
            
            # text output
            if args.text == 'profiles' or args.text == 'both':
                profile_txt_name = profile_file_name + ".txt"
                output_dir_method_profiles_text = os.path.join(output_dir_method_profiles, "txt_files")
                if not os.path.exists(output_dir_method_profiles_text):
                    os.makedirs(output_dir_method_profiles_text)
                profile_txt_path = os.path.join(output_dir_method_profiles_text, profile_txt_name)
                profile.create_text_file(profile_txt_path, reverse=True)
            
            # figures
            if args.figures == 'individual' or args.figures == 'both':
                output_dir_method_profiles_figs = os.path.join(output_dir_method_profiles, "png_files")
                if not os.path.exists(output_dir_method_profiles_figs):
                    os.makedirs(output_dir_method_profiles_figs)
                x_values = profile.x_values
                y_values = profile.y_values
                xlim = (min(x_values), max(x_values))
                xticks = list(np.linspace(min(x_values), max(x_values), num = 4))
                yticks = None#[0, 50, 100]
            
                figures.create_2Dplots(
                    list_of_x_values = [x_values],
                    list_of_y_values = [y_values],
                    fig_file_dir = output_dir_method_profiles_figs,
                    fig_file_basename = profile_file_name,
                    transparent = True,
                    xticks = xticks,
                    yticks = yticks,
                    xlim = xlim,
                    xlabel = "chemical shift [ppm]",
                    ylabel = "intensity [arb. units]",#"relative intensity [%]",
                    plot_settings = figures.plot_settings_profile_for_paper
                )
                        
            
            min_index, min_value = profile.establish_minima()
            
            #noise, verifying_minimum = profile.calc_noise()
            #quality = profile.verify_quality() #TODO: change, remove, use it?
            results_table.at[num, f"{method_str}_chem_shift"] = round(min_index, 3)
            
            
            #results_table_file = os.path.join(
            #    output_dir, "results_table_file.txt")
            #results_table.to_csv(results_table_file, index=False)
            

            #if profile.correct_chemical_shift == None:
            #TODO: Verification of the profile quality goes here:!!!
            
            # verifying overlap
            peak_overlap = peaks.PeakOverlap(peak)
            #peak.widths = [0.05, 1, 0.5] #ppm
            
            #print(peak.name, peak.chem_shifts, peak.widths)
            #for other_peak in peaklist_cs.peaklist:
                #other_peak.widths = [0.05, 1, 0.5]
                #peak_overlap.evaluate_overlap(other_peak, accordion_dim, 40)
                #peak_overlap.evaluate_overlap_cross_section(other_peak, accordion_dim, 40)
            #print(peak.name, peak_overlap.overlap, peak_overlap.overlap_strong, peak_overlap.overlap_strong_acc, peak_overlap.overlap_weak, peak_overlap.overlap_weak_acc)
            
            # *** Determination of correct shifts: ***
            if ref_comparison_flag:
                if (peak.acc_ref_shift == None or 
                    peak.acc_ref_shift == 'unknown'
                    ):
                    profile.correctness = "unknown"
                    #unknown_counter += 1
                elif abs(peak.acc_ref_shift - profile.chemical_shift) <= tolerance:
                    profile.correctness = True
                    #correct_counter += 1
                else:
                    profile.correctness = False
                    #wrong_counter += 1
    
                results_table.at[num, f"{method_str}_correctness"] = profile.correctness
            #cs_processed.values = list(cs_processed.values)
            cs_processed.check_weighted_spectrum_overlap_v3()
            cross_sections.calc_sino(cs_processed)
            #group = cross_sections.assign_to_group(cs_processed, profile, peak, None)
            validity = cross_sections.determine_result_validity(cs_processed, profile, peak, None)
            
            #results_table.at[num, f"{method_str}_group"] = group
            results_table.at[num, f"{method_str}_validity"] = validity
            
            results_table.at[num, f"{method_str}_profile"] = profile
            
            
    

    print("\033[K", end="") #ANSI escape code, clears the line
    print("Progress: 100% - processed all peaks. Finalizing output...\n")#, end='\033[1F')
    
    #TODO: add creating figures with a processed cross section and below the profile obtained from it; multiple such plots on one page; 6 rows and 5 columns, 1 row is cross sections, the 2nd is their profiles etc. with arrow in between; maybe cross section and profile as subplots within one plot, but one subplot underneath the other
    # ***** FINALIZING OUTPUT *****
    final_results = []
    for method in processing_methods:
        method_str = processing_methods_dict[method]
        profiles = results_table.get(f"{method_str}_profile").tolist()
        processed_cs_list = results_table.get(f"{method_str}_cs").tolist()
        
        
        
        # *** Figures with multiple cross sections in one image: ***
        if args.figures == 'acombined' or args.figures == 'aboth': #TODO: this process is way too slow for cross sections
            list_of_x_values = [[i for i in range(0, len(cs.values))] for cs in processed_cs_list]
            list_of_y_values = [cs.values for cs in processed_cs_list]
            fig_basename = f"Cross_sections_combined_{method_str}"
            #fig_file_path = os.path.join(output_dirs_methods[method], fig_name)
            fig_file_dir = output_dirs_methods_cs[method]
            names = [cs.name for cs in processed_cs_list]
            xticks = [len(cs.values)/2 for cs in processed_cs_list]
            xticks_labels = [0 for cs in processed_cs_list]
            yticks = [0 for cs in processed_cs_list]
            
            figures.create_2Dplots( #TODO: take list of 2D/cross_section/profile objects?
                list_of_x_values,
                list_of_y_values,
                fig_file_dir,
                fig_file_basename = fig_basename,
                #correct_points = None,
                names = None,
                #titles = None,
                titles = names,
                number_of_rows = 6,
                number_of_columns = 5,
                figure_width = 17,
                figure_height = 24,
                transparent = True,
                #xticks: dict[float: float|str] = None, #tick_value: label
                xticks = xticks,
                xticks_labels = xticks_labels,
                yticks = yticks,
                yticks_labels = None,
                #yticks: dict[float: float|str] = None, #tick_value: label
                xlim = None,
                ylim = None,
                xlabel = None,
                ylabel = None,
                list_of_scatter_x_values = None,
                #list_of_scatter_y_values = ,
                plot_settings = figures.plot_settings_combined_plots_for_paper
            )
            del list_of_x_values
            del list_of_y_values
            
        
        # *** Figures with multiple profiles in one image: ***
        if args.figures == 'combined' or args.figures == 'both':
            list_of_x_values = [profile.x_values for profile in profiles]
            list_of_y_values = [profile.y_values for profile in profiles]
            fig_basename = f"Profiles_combined_{method_str}"
            #fig_file_path = os.path.join(output_dirs_methods[method], fig_name)
            fig_file_dir = output_dirs_methods_profiles[method]
            names = [profile.name for profile in profiles]
            if ref_comparison_flag:
                correct_shifts = results_table.get("ref_chem_shift").tolist()
                refs = []
                for shift in correct_shifts:
                    if shift == 'unknown':
                        refs.append([None])
                    else:
                        refs.append([shift])
            else:
                refs = None
            #refs = [[shift] for shift in correct_shifts]

            #works to make y axis scale the same for all, but it does not look good because on many profiles the minima are not visible:
            """maxima, minima = [], []
            for i in list_of_y_values:
                maximum = max(i)
                minimum = min(i)
                maxima.append(maximum)
                minima.append(minimum)
            highest = max(maxima)
            lowest = min(minima)
            ylim = (lowest, highest)"""
            ylim=None
            x_values = list_of_x_values[0]
            plot_colors = []
            for num in range(0, len(profiles)):
                validity = results_table.at[num, f"{method_str}_validity"]
                if validity:
                    plot_color = "b"
                else:
                    plot_color = "r"
                plot_colors.append(plot_color)
            figures.create_2Dplots_v2( #TODO: take list of 2D/cross_section/profile objects?
                list_of_x_values,
                list_of_y_values,
                fig_file_dir,
                fig_file_basename = fig_basename,
                #correct_points = None,
                names = None,
                #titles = None,
                titles = names,
                number_of_rows = 6,
                number_of_columns = 5,
                figure_width = 17,
                figure_height = 24,
                transparent = True,
                #xticks: dict[float: float|str] = None, #tick_value: label
                xticks = list(np.linspace(min(x_values), max(x_values), num = 4)),#None,
                xticks_labels = None,
                yticks = None,
                yticks_labels = None,
                #yticks: dict[float: float|str] = None, #tick_value: label
                xlim = None,
                ylim = ylim,#None,
                xlabel = None,
                ylabel = None,
                list_of_scatter_x_values = refs,
                #list_of_scatter_y_values = ,
                plot_settings = figures.plot_settings_combined_plots_for_paper,
                xlabel_global = "chemical shift [ppm]",
                ylabel_global = "intensity [arb. units]",
                plot_colors = plot_colors,
                first_figure_plots_number = 25,
            )
            del list_of_x_values
            del list_of_y_values
        
        detected_shifts = []
        correct_shifts = []
        correct_profiles_detected = []
        correct_profiles_correct = []
        correct_counter = 0
        wrong_counter = 0
        unknown_counter = 0
        if ref_comparison_flag:
            for profile in profiles:
                if profile.correctness == True:
                    correct_counter += 1
                elif profile.correctness == False:
                    wrong_counter += 1
                else:
                    unknown_counter += 1
            all_counter = correct_counter + wrong_counter + unknown_counter
            correct_percentage = correct_counter * 100 / all_counter
            correct_percentage = round(correct_percentage, 4)
            wrong_percentage = wrong_counter * 100 / all_counter
            wrong_percentage = round(wrong_percentage, 4)
            unknown_percentage = unknown_counter * 100 / all_counter
            unknown_percentage = round(unknown_percentage, 4)
            correctness_msg = (f"For processing method {method_str}:""\n"
                f"Correct: {correct_percentage}% ({correct_counter})."
                f" Incorrect: {wrong_percentage}% ({wrong_counter})."
                f" Unknown: {unknown_percentage}% ({unknown_counter}).""\n")
            #print(f"For processing method {method}:\n{correctness_msg}")
            print(correctness_msg)
            final_results.append(correctness_msg)

            # INSERTING NEW PART
            cor_and_valid_sum = 0
            cor_and_invalid_sum = 0
            incor_and_valid_sum = 0
            incor_and_invalid_sum = 0
            unknown_and_valid_sum = 0
            unknown_and_invalid_sum = 0
            for num, profile in enumerate(profiles):
                correctness = profile.correctness
                validity = results_table.at[num, f"{method_str}_validity"]
                if correctness == True and validity:
                    cor_and_valid = True
                    cor_and_valid_sum += 1
                elif correctness == True and not validity:
                    cor_and_invalid = True
                    cor_and_invalid_sum += 1
                elif correctness == False and validity:
                    incor_and_valid = True
                    incor_and_valid_sum += 1
                elif correctness == False and not validity:
                    incor_and_invalid = True
                    incor_and_invalid_sum += 1
                elif correctness == 'unknown' and validity:
                    unknown_and_valid_sum += 1
                elif correctness == 'unknown' and not validity:
                    unknown_and_invalid_sum += 1


            new_cor_message = (f"For method {method_str}:\n"
                f"Correct and valid: {cor_and_valid_sum}\n"
                f"Correct but invalid: {cor_and_invalid_sum}\n"
                f"Incorrect but valid: {incor_and_valid_sum}\n"
                f"Incorrect and invalid: {incor_and_invalid_sum}\n"
                f"Unknown and valid: {unknown_and_valid_sum}\n"
                f"Unknown and invalid: {unknown_and_invalid_sum}\n")
            print(new_cor_message)
            final_results.append(new_cor_message)
            
            # END OF NEW PART

        output_dir_method = output_dirs_methods[method]
        output_textfile_path = os.path.join(output_dir_method, "output.txt")
        for i, profile in enumerate(profiles):
            detected_shifts.append(profile.chemical_shift)
            if ref_comparison_flag:
                correctness = profile.correctness
                correct_shift = results_table.at[i, "ref_chem_shift"]
            else:
                correct_shift = 'unknown'
            
            if ref_comparison_flag:
                if correct_shift != None or correct_shift != 'unknown':
                    #detected_shifts.append(profile.chemical_shift)
                    correct_shifts.append(correct_shift)
                if correctness == True:
                    correct_profiles_detected.append(profile.chemical_shift)
                    correct_profiles_correct.append(correct_shift)

        '''with open(output_textfile_path, "w") as file: #TODO: Should this be a function that takes a list of profiles as the argument?
            #file.write("Assignment"+" "*20 + "Global_min_shift" + " "*14 + "Relative intensity [%] \n \n")
            file.write("Assignment"+" "*20 + "Global_min_shift" + " "*14 + "Intensity [arb. units] \n \n")
            for i, profile in enumerate(profiles):
                
                #print(record.profiles)
                #profile = record.profiles[processing_method]
            #for profile in profiles:
                name = profile.name
                chemical_shift = profile.chemical_shift
                cs = results_table.at[i, f"{method_str}_cs"]
                overlap_factor = cs.overlap_factor
                cs_sino = cs.sino
                profile_sino = profile.sino
                #group = results_table.at[i, f"{method_str}_group"]
                validity = results_table.at[i, f"{method_str}_validity"]
                profile_minimum_intensity = profile.global_minimum
                if ref_comparison_flag:
                    correctness = profile.correctness
                #correct_shift = profile.correct_chemical_shift
                #correct_shift = profile_obj.peak.correct_CO_shift
                if ref_comparison_flag:
                    #correct_shift = profile_obj.peak.acc_ref_shift
                    #correct_shift = results_table.at[i, table_columns["reference_chem_shift"]]
                    correct_shift = results_table.at[i, "ref_chem_shift"]
                else:
                    correct_shift = 'unknown'

                #if correct_shift != None:
                if ref_comparison_flag:
                    if correct_shift != None or correct_shift != 'unknown':
                        detected_shifts.append(chemical_shift)
                        correct_shifts.append(correct_shift)
                    if correctness == True:
                        correct_profiles_detected.append(chemical_shift)
                        correct_profiles_correct.append(correct_shift)
                #file.write(f"{name}"+" "*20+f"{chemical_shift}"+" "*14+f"{profile_minimum_intensity}"+" "*5+f"{correctness}" + " "*5 + f"{correct_shift}" + f" {profile_obj.profile.quality}")
                #file.write(f"{name}"+" "*20+f"{chemical_shift}"+" "*14+f"{profile_minimum_intensity}"+" "*5+f"{correctness}" + " "*5 + f"{correct_shift}")# + f" {profile.quality}")
                    file.write(f"{name} {chemical_shift} "
                        f"{profile_minimum_intensity} {correctness} " 
                        f"{correct_shift} {overlap_factor} {cs_sino} "
                        f"{profile_sino} {validity}")# + f" {profile.quality}")
                    file.write("\n")
            if ref_comparison_flag:
                file.write(correctness_msg)'''
        if ref_comparison_flag:
            logger.log(999, correctness_msg)


        
        # Creating a graph of detected vs correct CO chemical shifts:
        if ref_comparison_flag:
            #detected_vs_correct_CO_shifts = (peak, correct_shift, chemical_shift)
            #x_values = correct_shifts
            correct_shifts = results_table.get(f"ref_chem_shift").tolist()
            #groups = results_table.get(f"{method_str}_group").tolist()
            validities = results_table.get(f"{method_str}_validity").tolist()
            indexes_to_remove = []
            for i, correct_shift in enumerate(correct_shifts):
                if (correct_shift == 'unknown' or 
                    (validities[i] !=True)):
                    #(groups[i] !="A" and groups[i] !="B")):
                    #correct_shifts.pop(i)
                    #detected_shifts.pop(i)
                    #groups.pop(i)
                    indexes_to_remove.append(i)
            
            indexes_to_remove.reverse()
            for ind in indexes_to_remove:
                correct_shifts.pop(ind)
                detected_shifts.pop(ind)
                #groups.pop(ind)
                validities.pop(ind)
            #i=0
            """while i < len(correct_shifts):
                correct_shift = correct_shifts[i]
                if correct_shift == 'unknown':
                    correct_shifts.pop(i)
                    detected_shifts.pop(i)
                    i -= 1
                i += 1"""
            #print(correct_shifts)
            x_values = [float(correct_shift) for correct_shift in correct_shifts]
            y_values = detected_shifts
            linregress_results = scipy.stats.linregress(correct_profiles_correct, correct_profiles_detected)
            slope = linregress_results.slope
            intercept = linregress_results.intercept
            r_squared = (linregress_results.rvalue)**2
            linregress_equation = f"y = {slope}*x + {intercept}; R2 = {r_squared}"
            #print(linregress_equation)
            logger.log(999, linregress_equation)
            with open(output_textfile_path, "a") as file:
                file.write(linregress_equation)
            #plt.text(
            fig_name = f"Ref_shifts_comparison_{method_str}.png"
            fig_path = os.path.join(output_dir_method, fig_name)
            #fig = plt.figure(layout = "tight")
            fig = plt.figure()
            fig_plot = fig.add_subplot()
            fig_plot.set_title(method)
            fig_plot.set_xlabel("Correct chemical shift [ppm]")
            fig_plot.set_ylabel("Detected chemical shift [ppm]")
            ax = fig.gca()
            #fig_plot.text(0, 1, linregress_equation, horizontalalignment='center',
            #    verticalalignment='center',
            #    transform=ax.transAxes)
            #x = f.gca()
            #ax = plt.gca()
            
            fig_plot.set_aspect("equal")
            fig_plot.plot([170,183],[170,183], color = "g")
            #y_regress = []
            #for x_value in correct_profiles_correct:
            #    y_regress.append(slope*x_value+intercept)
            #y_regress = [(float(slope) * x_value + float(intercept)) for x_value in x_values] 
            fig_plot.plot(x_values, y_values, color = "b", marker = 'o', ms = 1, ls = '') #!!!
            fig_plot.set_xticks(ticks=[170, 183])
            fig_plot.plot([170,183],[170*slope+intercept,183*slope+intercept], color = "r")
            ax.set_xlim(170, 183)
            fig_plot.set_xlim((170, 183))
            fig_plot.set_ylim((170, 183))
            fig_plot.set_xticks(ticks=[170, 183])
            #fig_plot.plot(correct_profiles_correct, y_regress, color = "r", ls='dashed')
            fig.savefig(fig_path)
            #print("saved")
            #plt.close(fig)
            plt.cla()
            #plt.clear()
            fig.clear()
            plt.close('all')
        
        #if ref_comparison_flag:
            cm = 1/2.54
            fig_width, fig_height = 8, 6
            xlim = (167, 183)
            ylim = (167, 183)
            fig = plt.figure(figsize = (fig_width * cm, fig_height * cm))
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(x_values, y_values, color = "b", marker = 'o', ms = 1, ls = '')
            ax.plot([170, 183], [170*slope+intercept, 183*slope+intercept], color = 'r', linewidth = 1)
            ax.plot([170,183],[170,183], color = '#00ff0080', linewidth = 0.5)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.xaxis.set_inverted(False)
            ax.yaxis.set_inverted(False)
            fig.savefig(fig_path, dpi = 300, transparent = False)
            plt.cla()
            fig.clear()
            plt.close()

            shift_axis_limits = None#(min(freqs), max(freqs))
            ticks = list(np.linspace(min(freqs), max(freqs), num = 4))
            figures.create_2Dplots(
                figure_width = 8.4,
                figure_height = 8.4,#6.3,
                list_of_x_values = [[min(freqs), max(freqs)]],
                list_of_y_values = [[min(freqs), max(freqs)]],
                fig_file_dir = output_dir_method,
                fig_file_basename = fig_name,
                xlabel = "reference chemical shift [ppm]",
                ylabel = "detected chemical shift [ppm]",
                xlim = shift_axis_limits,
                ylim = shift_axis_limits,
                list_of_scatter_x_values=[x_values],
                list_of_scatter_y_values=[y_values],
                plot_settings = figures.plot_settings_reference_comparison,
                xticks=ticks,
                yticks=ticks,
            )

        method_output_columns = ["peak_name"]

        for i in dims:
            dim_nucleus = spectrum_params.get_value("name", i)

            method_output_columns.append(f"chem_shift_{i}__{dim_nucleus}")
        if ref_comparison_flag: method_output_columns.append("ref_chem_shift")


        method_str = processing_methods_dict[method]
        method_output_columns.append(f"{method_str}_chem_shift")
        #latex_wanted_columns.append(f"{method_str}_group")
        method_output_columns.append(f"{method_str}_validity")
        #method_str = processing_methods_dict[("weighting", "gauss[5]")]
        #latex_wanted_columns.append(f"{method_str}_chem_shift")
        #latex_wanted_columns.append(f"{method_str}_validity")
        
        method_output_table = results_table[method_output_columns]
        csv_file = os.path.join(output_dir_method, "output_table.csv")
        method_output_table.to_csv(csv_file,sep=";", index=False)

    with open(main_output_textfile_path, "w") as file:
        for results in final_results:
            file.write(results)
            file.write("\n")

    # *** Creating latex code for a table with results: ***
    if args.latex:
        latex_wanted_columns = ["peak_name"]

        for i in dims:
            dim_nucleus = spectrum_params.get_value("name", i)

            latex_wanted_columns.append(f"chem_shift_{i}__{dim_nucleus}")
        if ref_comparison_flag: latex_wanted_columns.append("ref_chem_shift")

        for method in processing_methods:
            method_str = processing_methods_dict[method]
            latex_wanted_columns.append(f"{method_str}_chem_shift")
            #latex_wanted_columns.append(f"{method_str}_group")
            latex_wanted_columns.append(f"{method_str}_validity")
        #method_str = processing_methods_dict[("weighting", "gauss[5]")]
        #latex_wanted_columns.append(f"{method_str}_chem_shift")
        #latex_wanted_columns.append(f"{method_str}_validity")
        
        results_table_for_latex = results_table[latex_wanted_columns]
        csv_file = os.path.join(output_dir, "output_table.csv")
        results_table_for_latex.to_csv(csv_file,sep=";", index=False)
        #IMPORTANT
        """
            results_table_for_latex.rename(columns={old_name: new_name,
                                          f"{method_str}_cs": f"Cross section after processing with {method_str}",
                                          f"{method_str}_profile": f"Profile from processing method: {method_str}",
                                          f"{method_str}_chem_shift": f"Chem. shift determined from processing method: {method_str} [ppm]",
                                          f"{method_str}_group": f"Group assigned to the profile from processing method: {method_str}",
                                          f"{method_str}_correctness": f"Chem. shift correctness from processing method: {method_str}"})
        """
        
        latex_code=results_table_for_latex.to_latex(
            float_format="%.3f", index=False)
        latex_file = os.path.join(output_dir, "latex_table.txt")
        with open(latex_file, "w") as file:
            file.write(latex_code)

# TODO add documentation to all newly created functions and methods
# TODO move around some methods in cross_section.py so that they fit the class


    print(f"\nOutput saved in the following directory: {output_dir}")

    end_time = time.time()
    runtime = round(end_time - start_time, 3)
    logger.info("Program finished, it was running for: " + str(runtime) +
        " seconds.")
    print(f"Finished! Running for {round(runtime/60.0, 3)} minutes.")
    
    logger.info("Execution finished!")
    with open(log_file_path, "a") as log_file:
        log_file.write("-"*50+"\n")
        
        
        
        
        
        
# TODO: better noise calculation (skipping the big peaks)
# TODO: when overlap is the same check the signal to noise ratio; or the intensity of the signal itself

"""
1) Create cross sections from weighted spectrum
1.5) Find the main peak
2) Shift them
3) Do the weighting
4) Evaluate overlap, if there is some, find neighbouring cross sections.
5) Shift and do the weighting on the found neighbouring cross sections.
6) Evaluate overlap on the processed neighbouring cross sections and choose the best.
7) Based on the coords of the best cross section create the proper cross section for the actual spectrum.
8) Process it as normal.
    
Evaluation should be based on the noise calculated from the initial cross section before processing


From the center_on_peak method obtain the shifting factor, and then can shift all the peaks found in ssig.find_peaks? 
"""

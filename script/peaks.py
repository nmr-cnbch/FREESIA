"""
This module can be used for resonances and peaks of an NMR spectrum.

Functions:
    chemical_shift_from_ppm_to_point

Classes:
    Resonance: for a single resonance - single nucleus.
    ResonanceList: a list of resonances.
    Peak: for peaks of the spectrum.
    PeakList: for a list of peaks of the spectrum.
    PeakOverlap: for the detection of the peak overlap.
"""

from __future__ import annotations
import logging
import os

import parameters
import spectra

logger = logging.getLogger(__name__)



def chemical_shift_from_ppm_to_point(
    spectrum_params: spectra.SpectrumParams,
    dim: int,
    chem_shift_ppm: float,
    fnz: int = None,
    reffrq: float = None,
    sw: float = None,
    center: float = None,
    point_dir_begin: int = None,
    ) -> int:
    """
    Calculate the point number corresponding to the chemical shift.

    Required arguments:
        spectrum_params (SpectrumParams): Object of the 'SpectrumParams' 
            class, containing information about the parameters 
            of the spectrum.
        dim (int): Number of the spectrum's dimension.
        chem_shift_ppm (float): Value of the chemical shift in ppm.
    
    Optional arguments:
        fnz (int): Number of frequency points in the same dimension
            of the spectrum. If not specified, it is obtained from
            'spectrum_params'.
        reffrq (float): Reference frequency of the dimension in MHz.
            If not specified, it is obtained from 'spectrum_params'.
        sw (float): Spectral width of the dimension in Hz. 
            If not specified, it is obtained from 'spectrum_params'.
        center (float): Chemical shift of the center of the dimension 
            in ppm. If not specified, it is obtained from 
            'spectrum_params'.
        point_dir_begin (int): Number of the point from which 
            the desired range of the spectrum starts. It is non-zero 
            when only a part of the spectrum is being investigated, 
            it indicates the location of the first point of this part 
            within the entire dimension.
            If not specified, it is obtained from 'spectrum_params'. 
            If not present there, it defaults to 0.

    Returns:
        chem_shift_point (float): Value of the chemical shift 
            in point number.
    """
    # verifying input:
    if fnz == None:
        fnz = spectrum_params.get_value("fnz", dim)
    if reffrq == None:
        reffrq = spectrum_params.get_value("reffrq", dim)
    if sw == None:
        sw = spectrum_params.get_value("sw", dim)
    if center == None:
        center = spectrum_params.get_value("center", dim)
    if point_dir_begin == None:
        try:
            point_dir_begin = spectrum_params.get_value("point_dir_begin", dim)
        except KeyError:
            point_dir_begin = 0
    
    # calculation:
    chem_shift_point = fnz/2 - (chem_shift_ppm - center) * reffrq * fnz / sw
    chem_shift_point = round(chem_shift_point - point_dir_begin)
    #chem_shift_point = int(chem_shift_point - point_dir_begin)
    
    log_msg = (f"Calculated the point number ({chem_shift_point}) for"
        f" chemical shift equal to: {chem_shift_ppm} ppm.")
    logger.debug(log_msg)
    log_msg2 = (f"Dim {dim}, chemical shift in ppm: {chem_shift_ppm}; "
        f"point_number = {chem_shift_point} = {fnz}/2 - "
        f"({chem_shift_ppm} - {center}) * {reffrq} * {fnz} / {sw} - "
        f"{point_dir_begin}")
    logger.debug(log_msg2)
    
    return chem_shift_point



def new_chemical_shift_from_ppm_to_point(
    spectrum_params: spectra.SpectrumParams,
    dim: int,
    chem_shift_ppm: float,
    fnz: int = None,
    reffrq: float = None,
    sw: float = None,
    center: float = None,
    point_dir_begin: int = None,
    ) -> int:
    """
    Calculate the point number corresponding to the chemical shift.

    Required arguments:
        spectrum_params (SpectrumParams): Object of the 'SpectrumParams' 
            class, containing information about the parameters 
            of the spectrum.
        dim (int): Number of the spectrum's dimension.
        chem_shift_ppm (float): Value of the chemical shift in ppm.
    
    Optional arguments:
        fnz (int): Number of frequency points in the same dimension
            of the spectrum. If not specified, it is obtained from
            'spectrum_params'.
        reffrq (float): Reference frequency of the dimension in MHz.
            If not specified, it is obtained from 'spectrum_params'.
        sw (float): Spectral width of the dimension in Hz. 
            If not specified, it is obtained from 'spectrum_params'.
        center (float): Chemical shift of the center of the dimension 
            in ppm. If not specified, it is obtained from 
            'spectrum_params'.
        point_dir_begin (int): Number of the point from which 
            the desired range of the spectrum starts. It is non-zero 
            when only a part of the spectrum is being investigated, 
            it indicates the location of the first point of this part 
            within the entire dimension.
            If not specified, it is obtained from 'spectrum_params'. 
            If not present there, it defaults to 0.

    Returns:
        chem_shift_point (float): Value of the chemical shift 
            in point number.
    """
    # verifying input:
    if fnz == None:
        fnz = spectrum_params.get_value("fnz", dim)
    if reffrq == None:
        reffrq = spectrum_params.get_value("reffrq", dim)
    if sw == None:
        sw = spectrum_params.get_value("sw", dim)
    if center == None:
        center = spectrum_params.get_value("center", dim)
    if point_dir_begin == None:
        try:
            point_dir_begin = spectrum_params.get_value("point_dir_begin", dim)
        except KeyError:
            point_dir_begin = 0
    
    # calculation:
    sw_sparky = sw - (sw / fnz)
    sw_sparky_ppm = sw_sparky / reffrq
    sw_per_bin = sw_sparky_ppm / (fnz - 1)
    downfield = center + fnz/2 * sw_per_bin
    chem_shift_point = round((downfield - chem_shift_ppm) / sw_per_bin)
    chem_shift_point = chem_shift_point - point_dir_begin

    
    log_msg = (f"Calculated the point number ({chem_shift_point}) for"
        f" chemical shift equal to: {chem_shift_ppm} ppm.")
    logger.debug(log_msg)
    log_msg2 = (f"Dim {dim}, chemical shift in ppm: {chem_shift_ppm}; "
        f"point_number = {chem_shift_point} = {fnz}/2 - "
        f"({chem_shift_ppm} - {center}) * {reffrq} * {fnz} / {sw} - "
        f"{point_dir_begin}")
    logger.debug(log_msg2)
    
    return chem_shift_point



def ucsf_chemical_shift_from_ppm_to_point(
    spectrum_params: spectra.SpectrumParams,
    dim: int,
    chem_shift_ppm: float,
    fnz: int = None,
    reffrq: float = None,
    sw: float = None,
    center: float = None,
    point_dir_begin: int = None,
    ) -> int:
    """
    Calculate the point number corresponding to the chemical shift.

    Required arguments:
        spectrum_params (SpectrumParams): Object of the 'SpectrumParams' 
            class, containing information about the parameters 
            of the spectrum.
        dim (int): Number of the spectrum's dimension.
        chem_shift_ppm (float): Value of the chemical shift in ppm.
    
    Optional arguments:
        fnz (int): Number of frequency points in the same dimension
            of the spectrum. If not specified, it is obtained from
            'spectrum_params'.
        reffrq (float): Reference frequency of the dimension in MHz.
            If not specified, it is obtained from 'spectrum_params'.
        sw (float): Spectral width of the dimension in Hz. 
            If not specified, it is obtained from 'spectrum_params'.
        center (float): Chemical shift of the center of the dimension 
            in ppm. If not specified, it is obtained from 
            'spectrum_params'.
        point_dir_begin (int): Number of the point from which 
            the desired range of the spectrum starts. It is non-zero 
            when only a part of the spectrum is being investigated, 
            it indicates the location of the first point of this part 
            within the entire dimension.
            If not specified, it is obtained from 'spectrum_params'. 
            If not present there, it defaults to 0.

    Returns:
        chem_shift_point (float): Value of the chemical shift 
            in point number.
    """
    # verifying input:
    if fnz == None:
        fnz = spectrum_params.get_value("fnz", dim)
    if reffrq == None:
        reffrq = spectrum_params.get_value("reffrq", dim)
    if sw == None:
        sw = spectrum_params.get_value("sw", dim)
    if center == None:
        center = spectrum_params.get_value("center", dim)
    '''if point_dir_begin == None:
        try:
            point_dir_begin = spectrum_params.get_value("point_dir_begin", dim)
        except KeyError:
            point_dir_begin = 0'''
    
    sw_sparky = sw - (sw / fnz)
    sw_sparky_ppm = sw_sparky / reffrq
    sw_per_bin = sw_sparky_ppm / (fnz - 1)
    downfield = center + fnz/2 * sw_per_bin
    chem_shift_point = round((downfield - chem_shift_ppm) / sw_per_bin)
    
    '''log_msg = (f"Calculated the point number ({chem_shift_point}) for"
        f" chemical shift equal to: {chem_shift_ppm} ppm.")
    logger.debug(log_msg)
    log_msg2 = (f"Dim {dim}, chemical shift in ppm: {chem_shift_ppm}; "
        f"point_number = {chem_shift_point} = {fnz}/2 - "
        f"({chem_shift_ppm} - {center}) * {reffrq} * {fnz} / {sw} - "
        f"{point_dir_begin}")
    logger.debug(log_msg2)'''
    
    return chem_shift_point


class Resonance:
    """
    The resonance in an NMR spectrum - chemical shift of a nucleus.
     
    It holds information about the nucleus type, the residue 
    it belongs to and the chemical shift.

    Attributes:
        residue_number (int): Number of the residue within the peptide 
            sequence.
        residue_name (str): Name of the amino acid residue. For example:
            "G", "Gly" or "glycine".
        nucleus_type (str): Type of the nucleus. For example: "CO" 
            or "HN".
        chem_shift (float): Chemical shift value in ppm.

    Methods:
        __init__
        __str__
    """
    def __init__(self,
        residue_number: int,
        residue_name: str,
        nucleus_type: str,
        chem_shift: float,
        ) -> Resonance:
        """
        Initialize a new instance of the class.

        Required arguments:
            residue_number (int): Number of the residue within the peptide
                sequence.
            residue_name (str): Name of the amino acid residue. For example:
                "G", "Gly" or "glycine".
            nucleus_type (str): Type of the nucleus. For example: "CO" or "HN".
            chem_shift (float): Chemical shift value in ppm.
        """
        self.residue_number = residue_number
        self.residue_name = residue_name
        self.nucleus_type = nucleus_type
        self.chem_shift = chem_shift

    def __str__(self) -> str:
        resonance_info = f"Resonance {self.residue_name}{self.residue_number}{self.nucleus_type}: {self.chem_shift}"
        return resonance_info



class ResonanceList:
    """
    The list of resonances in an NMR spectrum.

    Attributes:
        resonances (list[Resonance]): List of resonances.

    Methods:
        create_from_peaklist(peaklist)

    ToDo:
        *?
    """
    def __init__(self,
        resonances: list[Resonance] = None,
        ) -> None:
        """
        Class constructor.

        Required arguments:
            resonances (list[Resonance]): List of resonances.

        Returns:
            None

        ToDo:
            *
        """
        self.resonances = resonances
    
    def create_from_peak(self,
        peak: Peak,
        ) -> ResonanceList:
        """
        Create and return a list of resonances based on the peak.  

        Required arguments:
            peak (Peak): Peak of the spectrum.

        Returns:
            self (ResonanceList): List of created resonances.

        Raises:
            TypeError: description of the error

        ToDo:
            *
        """
        resonances = []
        #for peak in peaklist.peaklist:
        #for peak in peaklist.peaks:
        name = peak.name
        chem_shifts = peak.chem_shifts
        #print(len(chem_shifts))
        peak_resonances = name.split("-")
        for i in range(1, len(chem_shifts) + 1):
            resonance = peak_resonances[i - 1]
        #for index, resonance in enumerate(peak_resonances):
            final = 0
            #print(index, resonance, chem_shifts)
            #print(i, resonance, chem_shifts)
            for j in range(1, len(resonance)):
                if resonance[1:j].isdigit():
                    residue_number = int(resonance[1:j])
                    residue_name = resonance[0]
                    final = j
            nucleus_type = resonance[final:].split("i")[0]
            if i == len(chem_shifts):
                chem_shift = chem_shifts[0]
            else:
                chem_shift = chem_shifts[i]
            resonance = Resonance(residue_number, residue_name,
                nucleus_type, chem_shift)
            resonances.append(resonance)
            #print(peak_name, residue_number, residue_name,nucleus_type, resonance_chem_shift)
            #self.resonances_list.append([residue_number, residue_name, nucleus_type, resonance_chem_shift])
            #for character in resonance[1:]:
            #    number = "".join(character )
        self.resonances = resonances
        return self

    def create_from_peaklist(self,
        peaklist: PeakList,
        ) -> ResonanceList:
        """
        Create and return a list of resonances based on the list of peaks.  

        Required arguments:
            peaklist (PeakList): List of peaks in the spectrum.

        Returns:
            self (ResonanceList): List of created resonances.

        Raises:
            TypeError: description of the error

        ToDo:
            *
        """
        if self.resonances == None:
            resonances = []
        else:
            resonances = self.resonances
        for peak in peaklist.peaklist:
        #for peak in peaklist.peaks:
            name = peak.name
            chem_shifts = peak.chem_shifts
            #print(len(chem_shifts))
            peak_resonances = name.split("-")
            for i in range(1, len(chem_shifts) + 1):
                resonance = peak_resonances[i - 1]
            #for index, resonance in enumerate(peak_resonances):
                final = 0
                #print(index, resonance, chem_shifts)
                #print(i, resonance, chem_shifts)
                for j in range(1, len(resonance)):
                    if resonance[1:j].isdigit():
                        residue_number = int(resonance[1:j])
                        residue_name = resonance[0]
                        final = j
                nucleus_type = resonance[final:].split("i")[0]
                if i == len(chem_shifts):
                    chem_shift = chem_shifts[0]
                else:
                    chem_shift = chem_shifts[i]
                resonance = Resonance(residue_number, residue_name,
                    nucleus_type, chem_shift)
                resonances.append(resonance)
                #print(peak_name, residue_number, residue_name,nucleus_type, resonance_chem_shift)
                #self.resonances_list.append([residue_number, residue_name, nucleus_type, resonance_chem_shift])
                #for character in resonance[1:]:
                #    number = "".join(character )
        self.resonances = resonances
        return self
    


class Peak:
    """
    The peak in an NMR spectrum, defined by its position in every dimension.
     
    It also holds and uses the object with information on the variable types
    of the parameters.

    Attributes:
        name (str): Name/label of the peak.
        chem_shifts (list[float]): Chemical shifts in ppm, defining
            the position of the peak in every dimension.
        widths (list[float]): The widths of the peak in ppm.
        height (float): The height of the peak in the spectrum.

    Methods:
        __init__(name, chem_shifts, height=None)
        __str__()
        center_peak()
    
    ToDo:
        *?
    """
    def __init__(self,
        name: str,
        chem_shifts: list[float],
        widths: list[float] = None,
        height: float = None,
        ) -> None:
        """
        Class constructor.
        
        Assign the name and chemical shifts of the peak, and other attributes
        but only if provided.

        Required arguments:
            name (str): Name/label of the peak.
            chem_shifts (list[float]): Chemical shifts in ppm, defining
                the position of the peak in every dimension.

        Optional arguments:
            widths (list[float]): The widths of the peak in ppm.
            height (float): The height of the peak in the spectrum.

        Returns:
            None

        ToDo:
            *verify if user's input is an object of the correct class?
        """
        self.name = name
        #chem_shifts_list = []
        chem_shifts_dict = {}
        if isinstance(chem_shifts, list):
            for i in range(len(chem_shifts)):
                #chem_shifts_list.append(chem_shifts[i])
                chem_shifts_dict[i] = float(chem_shifts[i])
        elif isinstance(chem_shifts, dict):
            chem_shifts_dict = chem_shifts #TODO: might need change - to deepcopy the chemical shifts?
        else:
            pass #what when the input is wrong?
        #self.chem_shifts = chem_shifts_list
        self.chem_shifts = chem_shifts_dict
        # Optional attributes (more to come?):
        if widths != None:
            self.widths = widths
        if height != None:
            self.height = height
        #logger.critical("Created Peak object")
        #logger.debug(self.chem_shifts)
        logger.debug(f"Created new Peak named {self.name} with "
            f"chemical shifts: {self.chem_shifts}")
    

    def __str__(self) -> str: # TODO
        peak_info = f"Peak labeled '{self.name}':"
        for attr, value in self.__dict__.items():
            if value != None and attr != "name":
                peak_info += f"\n\t{attr} = {value}"
            #still could be improved, especially regarding printing the shifts etc.
        return peak_info
        

    def center_peak(self):
        pass


    def calc_peak_coords(self,
        spectrum_params: parameters.SpectrumParams
        ) -> list[int]:
        """Calculate coordinates of the peak in the spectrum.
        
        Required arguments:
            spectrum_params (SpectrumParams): Object of the 'SpectrumParams'
                class, containing information about the parameters
                of the spectrum.

        Returns:
            coords (list[int]): Coordinates of the peak.
        """
        coords = []
        #for dim, shift in enumerate(self.chem_shifts):
        for dim, shift in self.chem_shifts.items():
            coord = chemical_shift_from_ppm_to_point(spectrum_params, dim, shift)#ucsf_chemical_shift_from_ppm_to_point(spectrum_params, dim, shift)#chemical_shift_from_ppm_to_point(spectrum_params, dim, shift)
            #coord = new_chemical_shift_from_ppm_to_point(spectrum_params, dim, shift)
            coords.append(coord)
        self.coords = coords
        logger.debug(self.coords)
        return self.coords
    
    def establish_widths(self,
        spectrum: spectra.Spectrum,
        spectrum_real_file: os.PathLike):
        #coords = self.coords
        #point_number = spectrum.calc_point_number_in_spectrum_real_file(coords)
        #spectrum_value = spectrum.read_value_in_spectrum_real_file(point_number, spectrum_real_file)
        widths = []
        for dim in range(0, len(self.coords)):
            widths.append(self.establish_width_in_specified_dimension(dim, spectrum, spectrum_real_file))
            """new_value = spectrum_value
            while new_value > spectrum_value / 2:       
                new_coords = coords
                new_coords[i] = coords - 1"""
        self.widths = widths
                
    def establish_width_in_specified_dimension(self,
        dim: int,
        spectrum: spectra.Spectrum,
        spectrum_real_file: os.PathLike,
        ) -> float:
        coords = self.coords
        point_number = spectrum.calc_point_number_in_spectrum_real_file(coords)
        # TODO: spectrum_value should be a sum of values from a list of spectrum_real_files
        spectrum_value = spectrum.read_value_in_spectrum_real_file(
            point_number, spectrum_real_file)
        compared_value1 = spectrum_value
        new_coords = [coord for coord in coords]
        i = 1
        while abs(compared_value1) > abs(spectrum_value / 2):
            #print(abs(compared_value1), abs(spectrum_value / 2))
            new_coords[dim] = coords[dim] - i
            new_point_number = \
                spectrum.calc_point_number_in_spectrum_real_file(new_coords)
            # TODO: spectrum_value should be a sum of values from a list of spectrum_real_files
            compared_value1 = spectrum.read_value_in_spectrum_real_file(
                new_point_number, spectrum_real_file)
            width1 = (new_coords[dim] - coords[dim]) * 2
            i += 1
        j = 1
        compared_value2 = spectrum_value
        while abs(compared_value2) > abs(spectrum_value / 2):
            new_coords[dim] = coords[dim] + j
            new_point_number = \
                spectrum.calc_point_number_in_spectrum_real_file(new_coords)
            compared_value2 = spectrum.read_value_in_spectrum_real_file(
                new_point_number, spectrum_real_file)
            width2 = (new_coords[dim] - coords[dim]) * 2
            j += 1
        #print(width1, width2)
        width = max(width1, width2)
        #chem_shift_point = fnz/2 - (chem_shift_ppm - center)*reffrq*fnz/sw
        #chem_shift_ppm = (fnz/2 - chem_shift_point)*sw/fnz/reffrq + center
        sw = spectrum.params.get_value("sw", dim)
        fnz = spectrum.params.get_value("fnz", dim)
        reffrq = spectrum.params.get_value("reffrq", dim)
        width_ppm = sw/(fnz*reffrq) * abs(width)
        return width_ppm
    
    def establish_height(self,
        spectrum: spectra.Spectrum,
        spectrum_real_file: os.PathLike
        ) -> float:
        coords = self.coords
        point_number = spectrum.calc_point_number_in_spectrum_real_file(coords)
        spectrum_value = spectrum.read_value_in_spectrum_real_file(
            point_number, spectrum_real_file)
        self.height = spectrum_value
        return self.height



class PeakList:
    """
    The list of peaks in an NMR spectrum.

    Attributes:
        peaks (list[Peak]): List of peaks.

    Methods:
        create_from_peaklist_sparky(peaklist_path)

    ToDo:
        *?
    """
    def __init__(self,
        peaks: list[Peak] = None,
        ) -> None:
        """
        Class constructor.

        Required arguments:
            peaks (list[Peak]): List of peaks.

        Returns:
            None

        ToDo:
            *
        """
        self.peaks = peaks #mismatch with self.peaklist in method create_from_peaklist_sparky
        self.peaklist = peaks


    def create_from_peaklist_sparky(self,
        peaklist_path: os.PathLike,
        ) -> PeakList:
        """
        Create and return a list of peaks based on the list file from Sparky.  

        Required arguments:
            peaklist_path (os.PathLike): Path of the peaklist file.

        Returns:
            self (PeakList): List of created peaks.

        Raises:
            TypeError: description of the error

        ToDo:
            *
        """
        with open(peaklist_path) as file: #TODO: add the path to the properties of the created PeakList object
            file_lines = file.readlines()
        first_line = file_lines[0]
        peaklist_legend = first_line.split()
        dim = 1
        while f"w{dim}" in peaklist_legend:
            dim += 1
        dim = dim - 1
        self.peaklist = [] #mismatch with self.peaks in the constructor
        for line in file_lines[1:]: #[1:]
            if line.strip() == "":
                continue
            line.strip()
            line_contents = line.split()
            #print("line content", line_contents)
            peak_name = line_contents[0]
            #chem_shifts = [float(line_contents[dim])]
            #chem_shifts = []
            chem_shifts = {0: float(line_contents[dim])}
            for index, chem_shift in enumerate(line_contents[1:dim]):
            #for index, chem_shift in enumerate(line_contents[1:dim]):
            #for chem_shift in line_contents[1:dim]:
                #print("shift", peak_name, index, chem_shift)
                #chem_shifts.append(float(chem_shift))
                chem_shifts[index + 1] = float(chem_shift)
            new_peak = Peak(peak_name, chem_shifts)
            #print("Chem shifts", chem_shifts)
            self.peaklist.append(new_peak)
        log_msg = f"Created new peaklist from Sparky file: {peaklist_path}"
        logger.info(log_msg)
        return self



class PeakOverlap:
    #This class will handle everything connected to peak overlap
    # so that it is separated from the Peak class
    def __init__(self, peak: Peak,):
        self.__overlap = []
        self.__overlap_weak = []
        self.__overlap_strong = []
        self.__overlap_weak_acc = []
        self.__overlap_strong_acc = []
        self.peak = peak

    @property
    def overlap(self):
        return self.__overlap

    @property
    def overlap_strong(self):
        return self.__overlap_strong

    @property
    def overlap_weak(self):
        return self.__overlap_weak

    @property
    def overlap_weak_acc(self):
        return self.__overlap_weak_acc

    @property
    def overlap_strong_acc(self):
        return self.__overlap_strong_acc

    def evaluate_overlap_peak_on_peak_1D(self, other, dim: int) -> str:
        """Evaluate whether the peak's position within the specified dimension
            overlaps with the position of the other peak.
        
        Arguments:
    
        Returns:
        str: 'strong', 'weak' or 'no', depending on the overlap detection"""
        if self.peak == other:
            pass
        else:
            self_left_edge = self.peak.chem_shifts[dim] - 0.5*self.peak.widths[dim]
            self_right_edge = self.peak.chem_shifts[dim] + 0.5*self.peak.widths[dim]
            other_left_edge = other.chem_shifts[dim] - 0.5*other.widths[dim]
            other_right_edge = other.chem_shifts[dim] + 0.5*other.widths[dim]
            #check if peak's center is within the linewidth of the other peak:
            if self.peak.chem_shifts[dim] >= other_left_edge and \
                self.peak.chem_shifts[dim] <= other_right_edge:
                return "strong"
            #otherwise check if both peak's linewidths are overlapping:
            elif (self_left_edge <= other_right_edge and 
                self_right_edge >= other_left_edge) or \
                (self_right_edge >= other_left_edge and
                self_left_edge <= other_right_edge):
                return "weak"
            #otherwise no overlap:
            else:
                return "no"
    
    def evaluate_overlap_peak_on_peak(self, other):
        overlap_list = []
        if self.peak == other:
            pass
        else:
            for dim in range(0, len(self.peak.chem_shifts)):
                overlap = self.peak.evaluate_overlap_peak_on_peak_1D(other, dim)
                overlap_list.append(overlap)
            #check if there is at least one dimension preventing from overlap:
            if "no" in overlap_list:
                pass
            #otherwise check if all dimensions indicate strong overlap:
            elif "weak" not in overlap_list:
                self.__overlap_strong.append(other.name)
            #otherwise there is only weak overlap in at least one dimension:
            else:
                self.__overlap_weak.append(other.name)
    
    def evaluate_overlap_cross_section(self, other, cs_dim: int, cs_dim_trim_range):
        overlap_list = []
        if self.peak == other:
            pass
        else:
            for dim in range(0, len(self.peak.chem_shifts)):
                if dim == cs_dim:
                    cs_left_edge = self.peak.chem_shifts[dim] - 0.5*cs_dim_trim_range
                    cs_right_edge = self.peak.chem_shifts[dim] + 0.5*cs_dim_trim_range
                    other_left_edge = other.chem_shifts[dim] - 0.5*other.widths[dim]
                    other_right_edge = other.chem_shifts[dim] + 0.5*other.widths[dim]
                    if (other.chem_shifts[dim] >= cs_left_edge and
                        other.chem_shifts[dim] <= cs_right_edge):
                        #overlap_list.append("strong")
                        overlap_list.append("strong_acc") 
                    #this condition has to be wrong, includes cases that shouldn't be included - RETHINK, DRAW IT
                    elif (((other_right_edge >= cs_left_edge) 
                        and (other_right_edge <= cs_right_edge))
                        or ((other_left_edge <= cs_right_edge)
                        and (other_left_edge >= cs_left_edge))):
                        #overlap_list.append("weak")
                        overlap_list.append("weak_acc")
                    else:
                        #overlap_list.append("no")
                        overlap_list.append("no_acc")
                else:
                    overlap_list.append(self.evaluate_overlap_peak_on_peak_1D(other, dim))
            if "no" in overlap_list or "no_acc" in overlap_list:
                pass
            #elif "strong" in overlap_list or "strong_acc" in overlap_list:
            #elif "strong_acc" in overlap_list: #CHANGE THIS, IT DETECTS STRONG OVERLAP EVEN THOUGH IT CAN BE WEAK IN SOME OF THE DIMENSIONS
            elif "strong_acc" in overlap_list and "weak" not in overlap_list:# and other.height/self.height*100 > 20:
                self.__overlap_strong_acc.append(other.name)
                self.__overlap_strong_acc.append(f"{other.height/self.peak.height*100:.3f}%") #
            else:
                self.__overlap_weak_acc.append(other.name)
                self.__overlap_weak_acc.append(f"{other.height/self.peak.height*100:.3f}%") #

    def evaluate_overlap(self, other, acc_dim, acc_dim_trim_range):
        if self.peak == other:
            pass
            return None
        else:
            overlap_list = []
            for i in range(0, len(self.peak.chem_shifts)):
                #only after checking each dimension, I can assign the overlap CORRECT THIS! #this seems corrected
                #now have to add checking overlap along the accordion dimension - along the cross section
                #strong overlap check
                if i == acc_dim:
                    if other.chem_shifts[i] >= self.peak.chem_shifts[i] - acc_dim_trim_range/2 and other.chem_shifts[i] <= self.peak.chem_shifts[i] + acc_dim_trim_range/2:
                        #overlap_list.append("strong")
                        overlap_list.append("strong_acc")
                    elif (other.chem_shifts[i] + other.widths[i]/2 >= self.peak.chem_shifts[i] - acc_dim_trim_range/2) \
                        or (other.chem_shifts[i] - other.widths[i]/2 <= self.peak.chem_shifts[i] + acc_dim_trim_range/2):
                        #overlap_list.append("weak")
                        overlap_list.append("weak_acc")
                    else:
                        #overlap_list.append("no")
                        overlap_list.append("no_acc")
                else:
                    if self.peak.chem_shifts[i]>=other.chem_shifts[i]-other.widths[i]/2 and self.peak.chem_shifts[i]<=other.chem_shifts[i]+other.widths[i]/2:
                        #there is strong overlap
                        overlap_list.append("strong")
                        #self.__overlap_strong.append(other)
                        #self.__overlap.append(other)
                    #weak overlap check
                    elif (self.peak.chem_shifts[i]-self.peak.widths[i]/2<=other.chem_shifts[i]+other.widths[i]/2 and self.peak.chem_shifts[i]+self.peak.widths[i]/2>=other.chem_shifts[i]-other.widths[i]/2) \
                        or (self.peak.chem_shifts[i]+self.peak.widths[i]/2>=other.chem_shifts[i]-other.widths[i]/2 and self.peak.chem_shifts[i]-self.peak.widths[i]/2<=other.chem_shifts[i]+other.widths[i]/2):
                        #there is weak overlap
                        overlap_list.append("weak")
                        #self.__overlap_weak.append(other)
                        #self.__overlap.append(other)
                    else:
                        #there is no overlap
                        overlap_list.append("no")
                        pass
            #peak on peak overlap
            if "no" in overlap_list:
                pass
            elif "strong" in overlap_list:
                self.__overlap_strong.append(other.name)
            else:
                self.__overlap_weak.append(other.name)
            
            #cross section overlap:
            if "no" in overlap_list or "no_acc" in overlap_list:
                pass
            #elif "strong" in overlap_list or "strong_acc" in overlap_list:
            elif "strong_acc" in overlap_list:
                self.__overlap_strong_acc.append(other.name)
            else:
                self.__overlap_weak_acc.append(other.name)
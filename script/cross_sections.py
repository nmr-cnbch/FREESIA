""" #TODO: fix the documentation of the module
This module can be used for extracting and processing cross sections 
of an NMR spectrum.

Classes:
    CrossSection: for cross sections of the spectrum.
    Profiles: for amplitude profiles obtained from cross sections.
"""
# TODO: clean up the code of this module, raw values of cs should be preserved
# any operation should modify the "modifiable" attribute of the cross section, except for addition
from __future__ import annotations
from copy import deepcopy
import logging
import os
import math
from statistics import median
import time

import matplotlib.pyplot as plt
import numpy as np
#import numpy.fft as fft #this is how it's always been
import scipy.fft as fft
import scipy.signal as ssig

import parameters
import peaks
import spectra

logger = logging.getLogger(__name__)



class Data2D:
    def __init__(self):
        pass


class CrossSection(Data2D):
    """
    Cross section of an NMR spectrum.

    Attributes:
        orig_values (list[float]): Intensities of the subsequent 
            spectral points of the cross section at the moment of its 
            creation, before any processing.
        orig_length (int): Original length of the cross section at 
            the moment of its creation, before any processing.
        values (list[float]): Intensities of the subsequent spectral 
            points of the cross section. Values are after processing, 
            if any was applied.
        length (int): Length of the cross section - after processing, 
            if any was applied.
        slicing_point_coords (list[float]): Slicing point coordinates.
        dim (int): Dimension along which the spectrum was sliced.
        name (str): Name of the cross section or of the peak 
            through which the spectrum was sliced.
        spectrum_path (os.PathLike): Path of the spectrum file 
            from which the cross section was obtained.
        centered_vals (remove?)
        peak_offset (int): Offset of the peak from the center of
            the cross section. Negative value means the peak is 
            on the left from the center, positive value - on the right.
        orig_peak_pos
        cur_peak_pos
        y_values (?)
        x_values (?)
        overlap_factor
        main_overlap_factor
        interfering_peaks
        corr_peak_index
        ft: Values after performing Fourier transformation

    Methods:
        __init__
        create_from_spectrum_real
        add_cross_sections
        center_on_peak
        trim
        weighting
        new_weighting
        zerofilling
        inverse_Fourier_transformation
        create_text_file
        check_overlap
        check_weighted_spectrum_overlap
        correct_peak_pos
        check_weighted_spectrum_overlap_v3
        check_weighted_spectrum_overlap_v2
        find_neighbouring_cross_sections
        
    TODO: add a method of saving parameters into a file
    TODO: should other parameters set at creation also be private attributes?
        like slicing_point_coords or dim etc.
    """
    
    
    def __init__(self, 
        input_values: list[float] | CrossSection,
        spectrum_file_path: os.PathLike = None, 
        dim: int = None, 
        slicing_point_coords: list[float] = None, 
        name: str = None, 
        ) -> CrossSection:
        """
        Create cross section and assign main attributes.

        The constructor should be used on its own only for copying
        already existing cross sections. To construct a new one,
        it is best to use one of the class methods instead.

        Required arguments:
            input_values (list[float] | CrossSection): List of 
                intensities of the spectral points in the cross section.
                Or the CrossSection object which will be copied.
        
        Optional arguments:
            spectrum_file_path (os.PathLike): Path of the spectrum file 
                from which the cross section is obtained.
            dim (int): Number of the dimension along which the spectrum
                is sliced.
            slicing_point_coords (list[int]): Slicing point coordinates.
            name (str): Name of the cross section or of the peak 
                through which the spectrum was sliced.


        TODO: option to take the Peak object instead of the coords of the point; although taking just coords may be more optimal (then coords are calculated only once; although peak coords could be calculated once anyway, and then the method would check if they exist within the Peak object before attempting to calculate them again)
        """
        if isinstance(input_values, CrossSection):
            for key, value in input_values.__dict__.items():
                self.__dict__[key] = value
        else:
            self.__orig_values = input_values
            self.__orig_length = len(self.__orig_values)
            self.slicing_point_coords = slicing_point_coords
            self.dim = dim
            self.name = name
            self.spectrum_path = spectrum_file_path


    @classmethod
    def create_from_spectrum_real(cls,
        spectrum_real_file_path: os.PathLike,
        spectrum: spectra.Spectrum,
        cs_dim: int,
        slicing_point_coords: list[int],
        name: str = None,
        ) -> CrossSection:
        """
        Create cross section of NMR spectrum from spectrum_real file.  

        Required arguments:
            spectrum_real_file_path (os.PathLike): 
                Path of the spectrum_real file.
            spectrum (Spectrum): Object of the 'Spectrum' class 
                containing information about the measured spectrum, 
                especially its params.
            cs_dim (int): Number of the dimension along which 
                the spectrum is sliced to obtain the cross section.
            slicing_point_coords (list[int]): Slicing point coordinates.
            name (str): Name of the cross section or of the peak 
                through which the spectrum was sliced.
        
        Returns:
            CrossSection_obj (CrossSection): New instance of the
                CrossSection class.
        """
        values = []
        points_dir = []
        point_numbers = []
        flip_flag = False

        spectrum_params = spectrum.params
        try:
            dim = spectrum.params["dim"]
        except KeyError:
            dim = spectrum.params.establish_dimensionality()

        for i in range(0, dim):
            if f"points_dir{i}" in spectrum.params:
                points_dir.append(spectrum.params.get_value("points_dir", i))
            else:
                points_dir.append(spectrum.params.get_value("fnz", i))

        for i in range(0, points_dir[cs_dim]):
            indexes = [slicing_point_coords[j] for j in range(0, dim)]
            indexes[cs_dim] = i

            point_number = \
                spectrum.calc_point_number_in_spectrum_real_file(
                indexes
            )
            point_numbers.append(point_number)
        
        with open(spectrum_real_file_path, "rb") as file:
            values = spectrum.read_many_values_in_spectrum_real_file(
                point_numbers, spectrum_real_file_path, file
            )
                
        if values[slicing_point_coords[cs_dim]] < 0:
            flip_flag = True
            for k in range(0, len(values)):
                values[k] = values[k] * (-1)

        CrossSection_obj = cls(values, 
            spectrum_file_path = spectrum_real_file_path,
            dim = cs_dim,
            slicing_point_coords = slicing_point_coords, 
            name = name, 
        )
        return CrossSection_obj


    @classmethod
    def create_from_ucsf(cls,
        ucsf_file_path: os.PathLike,
        spectrum: spectra.UCSF_Spectrum,
        cs_dim: int,
        slicing_point_coords: list[int],
        name: str = None,
        ) -> CrossSection:
        """
        Create cross section of NMR spectrum from ucsf file.  

        Required arguments:
            ucsf_file_path (os.PathLike): 
                Path of the ucsf file.
            spectrum (Spectrum): Object of the 'UCSF_Spectrum' class 
                containing information about the measured spectrum, 
                especially its params.
            cs_dim (int): Number of the dimension along which 
                the spectrum is sliced to obtain the cross section.
            slicing_point_coords (list[int]): Slicing point coordinates.
            name (str): Name of the cross section or of the peak 
                through which the spectrum was sliced.
        
        Returns:
            CrossSection_obj (CrossSection): New instance of the
                CrossSection class.
        """
        values = []
        points_dir = []
        point_numbers = []
        flip_flag = False

        params = spectrum.params
        try:
            dim = params["dim"]
        except KeyError:
            spectrum.read_ucsf_header(ucsf_file_path)
            dim = params["dim"]

        '''for i in range(0, dim):
            if f"points_dir{i}" in spectrum.params:
                points_dir.append(spectrum.params.get_value("points_dir", i))
            else:
                points_dir.append(spectrum.params.get_value("fnz", i))'''

        values = []
        #for i in range(0, points_dir[cs_dim]):
        for i in range(0, params.get_value("fnz", cs_dim)):
            indexes = [slicing_point_coords[j] for j in range(0, dim)]
            indexes[cs_dim] = i

            point_number = \
                spectrum.calc_point_in_ucsf(
                indexes
            )
            point_numbers.append(point_number)
            #value = spectrum.read_value_from_ucsf(point_number, ucsf_file_path)
            
            #values.append(value)

        with open(ucsf_file_path, "rb") as file:
            values = spectrum.read_many_values_from_ucsf(
                point_numbers, ucsf_file_path, file
            )
        
            
        '''with open(spectrum_real_file_path, "rb") as file:
            values = spectrum.read_many_values_in_spectrum_real_file(
                point_numbers, spectrum_real_file_path, file
            )'''
                
        if values[slicing_point_coords[cs_dim]] < 0:
            flip_flag = True
            for k in range(0, len(values)):
                values[k] = values[k] * (-1)
        
        CrossSection_obj = cls(values, 
            spectrum_file_path = ucsf_file_path,
            dim = cs_dim,
            slicing_point_coords = slicing_point_coords, 
            name = name, 
        )
        return CrossSection_obj
    

    @classmethod
    def create_from_any_spectrum(cls,
        spectrum_file_path: os.PathLike,
        spectrum: spectra.UCSF_Spectrum | spectra.Spectrum,
        cs_dim: int,
        slicing_point_coords: list[int],
        name: str = None,
        spectrum_format: str = "ucsf",
        ) -> CrossSection:
        """
        Create cross section of NMR spectrum from a spectrum file.  

        Required arguments:
            spectrum_file_path (os.PathLike): 
                Path of the spectrum file.
            spectrum (Spectrum): Object of the 'UCSF_Spectrum' 
                or 'Spectrum' class containing information about 
                the measured spectrum, especially its params.
            cs_dim (int): Number of the dimension along which 
                the spectrum is sliced to obtain the cross section.
            slicing_point_coords (list[int]): Slicing point coordinates.

        Optional arguments:
            name (str): Name of the cross section or of the peak 
                through which the spectrum was sliced. Defaults to None.
            spectrum_format (str): Format of the spectrum file, either
                "ucsf" or "spectrum_real". Defaults to "ucsf".
        
        Returns:
            CrossSection_obj (CrossSection): New instance of the
                CrossSection class.
        """
        #if str(spectrum_file_path).endswith("ucsf"):
        if spectrum_format == "ucsf":
            create_cs_func = cls.create_from_ucsf
        elif spectrum_format == "spectrum_real":
            create_cs_func = cls.create_from_spectrum_real
        else:
            raise ValueError(f"{spectrum_format} is an unknown spectrum format")
        CrossSection_obj = create_cs_func(spectrum_file_path, spectrum,
            cs_dim, slicing_point_coords, name, 
        )
        return CrossSection_obj    
    
    
    @property
    def orig_values(self):
        return self.__orig_values    


    @property
    def orig_length(self):
        return self.__orig_length
    
    
    @property
    def values(self):
        try:
            values = self._values
        except AttributeError:
            values = self.__orig_values
        return values
    

    @values.setter
    def values(self, vals):
        self._values = vals
        self._length = len(vals)


    @property
    def length(self):
        try:
            length = self._length
        except AttributeError:
            length = self.__orig_length
        return length


    def add_cross_sections(self,
        cross_sections: list[CrossSection],
        ) -> CrossSection:
        """
        Add together several cross sections.

        Useful for situations when the spectrum was measured multiple times and
        the spectra should have been combined together into a single spectrum.

        Required arguments:
            cross_sections (list[CrossSection]): List of objects of the
                'CrossSection' class.

        Returns:
            self (CrossSection): The sum of the original cross sections and
                all the cross sections from the list.

        Raises:
            TypeError: description of the error

        TODO:
            *update docstring, maybe method too?
        """
        if len(cross_sections) > 0:
            for i in range(0, self.length):
                for j in range(0, len(cross_sections)):
                    self.values[i] += cross_sections[j].values[i]
        return self
    

    def center_on_peak(self,
        peak: peaks.Peak,
        center: int = None,
        spectrum_params: spectra.SpectrumParams = None,
        ) -> list[float]:
        """
        Shift the cross section to center it on the specific peak or point.  

        Required arguments:
            peak (Peak): Object of the 'Peak' class.
        
        Optional arguments:
            center (int): Point number of the cross section on which 
                it should be centered. Defaults to None. When no value 
                is provided, then the center corresponds to 
                the coordinate of the peak.
            spectrum_params (SpectrumParams): Object of the 
                'SpectrumParams' class, containing information about 
                the parameters of the spectrum. Defaults to None. 
                It is unused if the 'center' argument is specified. 
                Otherwise it may be necessary in order to calculate 
                the coordinates of the peak but only when they were not 
                calculated for that peak before.

        Returns:
            self.values (list[float]): Values of the centered cross
                section.

        Raises:
            TypeError: description of the error

        TODO:
            *add proper input management, check if center is within the len
                of cross section
        """
        if center == None:
            try:
                offset = self.peak_offset
            except AttributeError:

                if not hasattr(peak, 'coords'):
                    if spectrum_params == None:
                        exit("needs spectrum_params object")
                    else:
                        peak_coords = peak.calc_peak_coords(spectrum_params)
                point_coords = peak.coords
                central_point = point_coords[self.dim]
                offset = int(central_point - self.length/2)
                self.peak_offset = offset

        else:
            central_point = center
            offset = int(central_point - self.length/2)
            self.peak_offset = offset
        central_point = peak.coords[self.dim] #TODO: remove this line later
        centered_cs = [self.values[i] for i in range(0, self.length)]
        if offset < 0:
            for j in range(0, abs(offset)):
                centered_cs.insert(0, centered_cs.pop(-1))
        elif offset > 0:
            for j in range(0, offset):
                centered_cs.insert(-1, centered_cs.pop(0))
        self.values = centered_cs
        self.centered_vals = centered_cs
        self.orig_peak_pos = central_point #TODO: also define elsewhere?
        self.cur_peak_pos = central_point - offset #TODO: also define elsewhere?
        return self.values


    def trim(self,
        remaining_points_count: int,
        central_point: int = None,
        ) -> list[float]:
        """
        Trim the cross section to leave a specific amount of points remaining.

        Trimming occurs around the central point, by default the middle of the
        cross section, for example a cross section with 100 points and trimmed
        to 50 points will have the first 25 and the last 25 points removed,
        the 50 points in the middle are unchanged.
        Removed points are replaced by zeroes, the length of the cross section
        remains the same.

        Required arguments:
            remaining_points_count (int): Number of points to leave unchanged.

        Optional arguments:
            central_point (int): Point number of the cross section around which
                it will be trimmed. Defaults to None. When no value is provided,
                then the middle point of the cross section is used.

        Returns:
            trimmed_cs (list[float]): List of values of the trimmed cross
                section.

        Raises:
            TypeError: description of the error

        ToDo:
            *check and control for where the central point actually is,
                if it's on the left from the middle, on the right,
                or actually in the middle (when cross section has odd number
                of points)
            *might still be necessary to check for a situation with
                the central point being to close to the edge of the cross section
        """
        if not isinstance(remaining_points_count, int):
            print("not integer!")
            exit() #something else so that the program doesn't continue
        if central_point == None:
            central_point = self.length/2
        else:
            if not isinstance(central_point, int):
                print("not integer!")
                exit()
        
        left_edge = int(central_point - remaining_points_count/2)
        right_edge = int(central_point + remaining_points_count/2)
        #if right_edge - left_edge != remaining_points_count:
            #print(f"Wrong: {right_edge - left_edge}")
        trimmed_cs = []
        for i in range(0, self.length):
            if i >= left_edge and i < right_edge:
                trimmed_cs.append(self.values[i])
            else:
                trimmed_cs.append(0)
        self.values = trimmed_cs
        trimmed_msg = (f"Trimmed around the point {central_point} to contain "
            f"only {remaining_points_count} points")
        #self.history.append(trimmed_msg)
        return self.values
    
    #TODO: change name to weighting
    def new_weighting(self,
        weighting_values: list[float],
        ) -> list[float]:
        """
        Multiply the values of the cross section by weights.

        Required arguments:
            weighting_values (list[float]): List of weights by which 
                the cross section should be multiplied. The list has to 
                have the same length as the cross section.
        
        Returns:
            self.values (list[float]): Weighted values of cross section.

        Raises:
            ValueError: If the list of weights has different length
                than the cross section.
        """
        if len(weighting_values) != self.length:
            raise ValueError("The list of weighting values has different "
                "length than the cross section.")
        else:
            new_values = [value * weight for value,weight in \
                zip(self.values, weighting_values)]
            self.values = new_values
        return self.values


    def zerofilling(self, 
        minimal_zerofilling_factor: int
        ) -> list[float]:
        """
        Add zeroes to the beginning and end of the cross section.

        The number of zeroes added is adjusted so that the final length 
        of the cross section is equal to the nearest power of 2.

        Example:
        If the length of the original cross section is 300 and 
        minimal_zerofilling_factor is set to 2, then the minimum
        amount of points in the final cross section should be 600.
        The nearest power of 2 for that number is 1024 (2**10), 
        so 724 zeroes get added to the original 300 points 
        of the cross section (362 zeroes at the beginning and 
        362 zeroes at the end of the original values).

        Required arguments:
            minimal_zerofilling_factor (int): How many times the length 
                of the cross section should be increased by zerofilling.
        
        Returns:
            zerofilling_factor (float): How many times the length 
                of the cross section got increased by zerofilling.
        """
        #TODO: log the length of the cross section before and after zerofilling, 
            #how many zeroes got added etc.
        #TODO: take into account odd number of points in cs
        length = self.length
        two_power = 2**1
        i = 1
        while two_power < minimal_zerofilling_factor*length:
            i += 1
            two_power = 2**i
            
        new_length = two_power
        zeroes_amount = int((new_length - length) / 2)
        zeroes_list = [0 for j in range(0, zeroes_amount)]
        new_values = zeroes_list + self.values + zeroes_list
        zerofilling_factor = len(new_values) / length
        self.values = new_values
        self.zerofilling_factor = zerofilling_factor
        return zerofilling_factor


    def inverse_Fourier_transformation(self,
        mode: str = "r"
        ) -> list[float]:
        """
        Calculate the inverse Fourier transformation of the cross section.

        Optional arguments:
            mode (str): Determines the type of the result. Defaults to "r"
                which corresponds to keeping only the real part
                of the Fourier transform.

        Returns:
            ift_results (list[float]): List of values of the inverse Fourier
                transform of the cross section.
        """
        cs = self.values
        shifted_data = fft.fftshift(cs)
        ift_results = fft.ifft(shifted_data)
        if mode == "r":
            ift_results = np.real(ift_results)
        self.ft = ift_results

        """textfile_path = "/nhome/bartek/Documents/Doktorat/Python_scripts/classes/COdecProfiles_v4_testing/cs.txt"
        with open(textfile_path, "w") as file:
            for i, value in enumerate(cs):
                file.write(f"{i}, {value}\n")
        
        textfile_path2 = "/nhome/bartek/Documents/Doktorat/Python_scripts/classes/COdecProfiles_v4_testing/shifted_data.txt"
        with open(textfile_path2, "w") as file:
            for i, value in enumerate(shifted_data):
                file.write(f"{i}, {value}\n")
        
        textfile_path3 = "/nhome/bartek/Documents/Doktorat/Python_scripts/classes/COdecProfiles_v4_testing/ift_results.txt"
        with open(textfile_path3, "w") as file:
            for i, value in enumerate(ift_results):
                file.write(f"{i}, {value}\n")"""
        #exit()

        return ift_results
       
    
    def create_text_file(self, 
        file_path: os.PathLike, 
        x_values: list[float] = None, 
        reverse: bool = False,
        ) -> None:
        """
        Save the values of the cross section into a text file.

        The output file consists of two columns, the first has values 
        for the x axis, the second column has values for the y axis.

        Required arguments:
            file_path (os.PathLike): Path of the text file to create.
        
        Optional arguments:
            x_values (list[float]): Values to use for the x axis.
                If not provided, the values will correspond to 
                the index of the cross section value.
            reverse (bool): Determines if the values in the output file 
                are reversed (for True) or not (for False). 
        """
        #TODO: add checking if x_values and y_values have equal lengths
        try:
            y_values = self.values
        except:
            #TODO: maybe skip this except by putting values as property
                #for the profile class, so that it can take values
                # from the y_values attribute
            y_values = self.y_values
        
        if x_values == None:
            try:
                x_values = self.x_values
            except:
                x_values = [i for i in range(0, len(y_values))]
        
        if reverse:
            x_values = list(reversed(x_values))
            y_values = list(reversed(y_values))
            
        with open(file_path, "w") as file:
            for k in range(0, len(y_values)):
                file.write(f"{x_values[k]}   {y_values[k]}\n")
    

    def check_overlap(self):
        #self.values
        #self.slicing_point_coords
        #self.spectrum_path
        #self.cs_dim
        main_peak = self.slicing_point_coords[self.dim]
        main_peak_height = self.values[main_peak]
        interfering_peaks, peaks_properties = ssig.find_peaks(self.values,
            height = 0.05 * main_peak_height)
        interfering_peak_heights = []
        for interfering_peak in interfering_peaks:
            interfering_peak_height = self.values[interfering_peak]
            interfering_peak_heights.append(interfering_peak_height)
        return interfering_peaks, interfering_peak_heights
    
    def check_weighted_spectrum_overlap(self):
        main_peak = self.slicing_point_coords[self.dim]
        all_peaks_indexes, _ = ssig.find_peaks(self.values)
        all_peaks_indexes = all_peaks_indexes.tolist()
        if main_peak in all_peaks_indexes:
            all_peaks_indexes.remove(main_peak)
        else:
            all_peaks_heights = [self.values[peak_index] for peak_index in all_peaks_indexes]
            highest_value = max(all_peaks_heights)
            #print(f"Main peak pos: {main_peak}, detected max at: {self.values.index(highest_value)}")
        all_peaks_heights = [self.values[peak_index] for peak_index in all_peaks_indexes]
        #for i, peak_height in enumerate(all_peaks_heights):
        highest_value = max(all_peaks_heights)
        highest_value_index = self.values.index(highest_value)
        ratio = highest_value / self.values[main_peak] * 100
        
        #print("\n", self.name, highest_value_index, highest_value, ratio)
        
        return highest_value_index, highest_value
    
    def correct_peak_pos(self):
        peak_indexes, _ = ssig.find_peaks(self.values)
        try:
            orig_peak_index = self.cur_peak_pos
        except AttributeError:
            orig_peak_index = self.slicing_point_coords[self.dim]
        #print(self.slicing_point_coords)
        orig_peak_height = self.values[orig_peak_index]
        corr_peak_index, corr_peak_height = orig_peak_index, orig_peak_height
        for peak_index in peak_indexes:
            peak_height = self.values[peak_index]
            if (abs(orig_peak_index - peak_index) == 1 and
                abs(peak_height) - abs(orig_peak_height)) > 0:
                corr_peak_index = peak_index
                corr_peak_height = peak_height
                break
        return corr_peak_index, corr_peak_height
    
    def check_weighted_spectrum_overlap_v3(self):
        #peaks_dict = {}
        #orig_peak_index = self.slicing_point_coords[self.cs_dim]
        #orig_peak_height = self.values[orig_peak_index]
        peak_indexes, _ = ssig.find_peaks(self.values)
        corr_peak_index, corr_peak_height = self.correct_peak_pos()
        
        peak_indexes = peak_indexes.tolist()
        if corr_peak_index in peak_indexes:
            peak_indexes.remove(corr_peak_index)
            # it should detect the peak around the middle now, not around the original peak position
        else: 
            #print("PEAK NOT IN?")
            pass #TODO: verify this part?
        #max_height = max(peak_indexes)
        #if orig_peak_index in peak_indexes:
            
        #calculate overall_overlap_factor and main_overlap_factor
        overlap_factor = 0
        interfering_peaks = []
        #print(len(peak_indexes))
        if len(peak_indexes) > 0:
            all_peaks_heights = [self.values[i] for i in peak_indexes]
            #all_peaks_indexes = [peaks[i][0] for i in range(0, len(peaks))]
            all_peaks_indexes = peak_indexes
            highest_value = max(all_peaks_heights)
            highest_value_index = self.values.index(highest_value)
            #noise, _ = calc_noise(self.values)
            noise, _, _, _ = calc_noise(self.orig_values)
            #return highest_value_index, highest_value
            
            
            #for peak in peaks:
            for index in all_peaks_indexes:
                #peak_index, peak_height = peak[0], peak[1]
                peak_height = self.values[index]
                if peak_height >= 2*noise:
                    interfering_peaks.append(index)
                    overlap_factor += peak_height / corr_peak_height
        #else:
        if len(interfering_peaks) > 0:
            highest_value = max([self.values[interfering_peaks[i]] for i in range(0, len(interfering_peaks))])
        else:
            highest_value = 0
        main_overlap_factor = highest_value / corr_peak_height
        
        self.overlap_factor = overlap_factor
        self.main_overlap_factor = main_overlap_factor
        self.interfering_peaks = interfering_peaks
        self.corr_peak_index = corr_peak_index
        #return all_peaks_indexes, all_peaks_heights, noise
        #print(f"{self.name}, overlap factor: {overlap_factor}, main overlap factor: {main_overlap_factor}, interfering peaks:{interfering_peaks}")
        return interfering_peaks, corr_peak_index, overlap_factor, main_overlap_factor


    def check_weighted_spectrum_overlap_v2(self):
        #peaks_dict = {}
        orig_peak_index = self.slicing_point_coords[self.dim]
        orig_peak_height = self.values[orig_peak_index]
        peak_indexes, _ = ssig.find_peaks(self.values)
        #peak_indexes = peak_indexes.tolist()
        #max_height = max(peak_indexes)
        #if orig_peak_index in peak_indexes:
            
        peaks = []
        main_peak_index = orig_peak_index
        main_peak_height = orig_peak_height
        for peak_index in peak_indexes:
            peak_height = self.values[peak_index]
            # check if the main peak needs to be moved
            #print(abs(peak_height) - abs(orig_peak_height))
            if (abs(orig_peak_index - peak_index) <= 1 and
                abs(peak_height) - abs(orig_peak_height)) > 0:
                main_peak_index = peak_index
                main_peak_height = peak_height
            else:
                peak_tuple = (peak_index, self.values[peak_index])
                peaks.append(peak_tuple)
        #print(self.name, peaks)
        peaks.remove((main_peak_index, main_peak_height))
        overlap_factor = 0
        interfering_peaks = []
        #print(len(peaks))
        if len(peaks) > 0:
            all_peaks_heights = [peaks[i][1] for i in range(0, len(peaks))]
            all_peaks_indexes = [peaks[i][0] for i in range(0, len(peaks))]
            highest_value = max(all_peaks_heights)
            highest_value_index = self.values.index(highest_value)
            #noise, _ = calc_noise(self.values)
            #return highest_value_index, highest_value
            
            
            for peak in peaks:
                peak_index, peak_height = peak[0], peak[1]
                if peak_height >= 2*noise:
                    interfering_peaks.append(peak)
                    overlap_factor += peak_height / main_peak_height
        #else:

        #return all_peaks_indexes, all_peaks_heights, noise
        #print(f"{self.name}, overlap factor: {overlap_factor}, interfering peaks:{interfering_peaks}")
        return interfering_peaks, main_peak_index, overlap_factor
        # find all peaks
        # check if one of them is close to the original peak index
        # change the main index of the peak
        # create a list of peaks to verify
        # calculate noise
        # evaluate if an interfering peak is higher than noise*2
        # calculate the ratio of overlapping peak heights
        # assign to group? - maybe only for the finally selected cross section
    
    def find_neighbouring_cross_sections(self, 
        coords, 
        spectrum_real_paths, 
        spectrum,
        ) -> list[CrossSection]:
        cs_dim = self.dim
        dim = len(coords)
        all_coords_to_check = [coords]
        for i in range(0, dim):
            if i == cs_dim:
                continue
            else:
                coords_to_check = [coord for coord in coords]
                coords_to_check[i] = coords_to_check[i] - 1
                all_coords_to_check.append(coords_to_check)

                coords_to_check = [coord for coord in coords]
                coords_to_check[i] = coords_to_check[i] + 1
                all_coords_to_check.append(coords_to_check)

        #print(f"{self.name}, verifying coordinates: {all_coords_to_check}")
        cross_sections_to_compare = []
        for coords_to_check in all_coords_to_check:
            #cs_to_check = CrossSection()
            #cs_to_check = CrossSection.create_from_spectrum_real(
            cs_to_check = CrossSection.create_from_any_spectrum(
                spectrum_real_paths[0],
                spectrum,
                cs_dim,
                coords_to_check,
                self.name
            )
            #cs_to_check2 = CrossSection()
            #cs_to_check2 = CrossSection.create_from_spectrum_real(
            '''cs_to_check2 = CrossSection.create_from_any_spectrum(
                spectrum_real_paths[1],
                spectrum,
                cs_dim,
                coords_to_check,
                self.name
            )
            cs_to_check.add_cross_sections([cs_to_check2])'''
            cross_sections_to_compare.append(cs_to_check)
        return cross_sections_to_compare
        

class Profile(CrossSection):
    def __init__(self,
        cross_section: CrossSection,
        x_list: list[float],
        length: int = None,
        ):
        self.x_values = x_list #TODO: add a check to ensure the lengths of both lists (x and y points) are equal
        if length != None:
            #self.y_values = cross_section.ift[:length]
            self.y_values = cross_section.ft[:length]
        else:
            #self.y_values = cross_section.ift
            self.y_values = cross_section.ft
        self.name = cross_section.name
    
    def rescale_percentage(self): #more suitable for Data2D

        """textfile_path4 = "/nhome/bartek/Documents/Doktorat/Python_scripts/classes/COdecProfiles_v4_testing/profile.txt"
        with open(textfile_path4, "w") as file:
            file.write(self.name)
            for i, value in enumerate(self.y_values):
                file.write(f"{i}, {value}\n")"""
        #exit()

        np_array = np.array(self.y_values)
        maximum = np.max(np_array)
        rescaled_array = np_array/maximum * 100
        self.y_values = list(rescaled_array)

    def establish_minima(self, number_of_minima = 1): #more suitable for more general Data2D
        y_values = np.array(self.y_values)
        min_peaks_coord, _ = ssig.find_peaks(-y_values)
        global_minimum = np.max(y_values)

        for i in range(0, len(min_peaks_coord)):
            if y_values[min_peaks_coord[i]] < global_minimum:
                global_minimum = y_values[min_peaks_coord[i]]
                global_minimum_index = min_peaks_coord[i]

        if number_of_minima > 1:
            pass #TODO: implement this, utilising already found min_peaks_coord
        self.global_minimum_index = global_minimum_index
        self.global_minimum = global_minimum
        global_minimum_shift = self.x_values[global_minimum_index]
        self.chemical_shift = global_minimum_shift
        return global_minimum_shift, global_minimum #TODO: if above is implemented, then it should return lists of indexes and minima
    
    def verify_quality(self):
        #add a ton of logging
        values = np.array(self.y_values)
        min_peaks_coords, _ = ssig.find_peaks(-values)
        peak_heights = [self.y_values[coord] for coord in min_peaks_coords]
        #print("QUALITY CHECK")
        def sec_tuple(x):
            return x[1]
        peak_heights = []
        for coord in min_peaks_coords:
            #print(coord, self.y_values[coord])
            profilepeak = (coord, self.y_values[coord])
            peak_heights.append(profilepeak)
        peak_heights.sort(key=sec_tuple)
        #print(peak_heights)
        main_minimum = peak_heights[0]
        second_minimum = peak_heights[1]
        height_difference = abs(main_minimum[1] - second_minimum[1])
        #print(self.name, height_difference)
        if height_difference >= (100 - main_minimum[1])/2:
            self.quality = "g"
            
        else:
            self.quality = "b"
        return self.quality
    
    def calc_noise(self):
        # Estimating noise and its standard deviation
        #log the results and calculation etc.
        ft_results_zeroed_peaks = []
        ft_results_noise = []
        noise_squares_sum = 0.0
        values = self.y_values
        global_minimum_index = self.global_minimum_index
        global_minimum_start_index = global_minimum_index - 40
        global_minimum_end_index = global_minimum_index + 40
        #second_minimum_start_index = second_minimum_index - 40
        #second_minimum_end_index = second_minimum_index + 40
        
        for k in range(0, len(values)):
            if k >= global_minimum_start_index \
                and k <= global_minimum_end_index:
                #or k >= second_minimum_start_index \
                #and k <= second_minimum_end_index:
                ft_results_zeroed_peaks.append(0)
            else:
                ft_results_zeroed_peaks.append(values[k])
                ft_results_noise.append(ft_results_zeroed_peaks[k])
        if len(ft_results_noise) != 0:
            noise_average = sum(ft_results_noise) / float(len(ft_results_noise))
        
            for m in range(0, len(ft_results_noise)):
                noise_squares_sum += (ft_results_noise[m] - noise_average)**2
            noise_standard_deviation = math.sqrt(
                noise_squares_sum/len(ft_results_noise))
        
            minimum_test_value = noise_average - 6*noise_standard_deviation
        else:
            minimum_test_value = 0
        return noise_standard_deviation, minimum_test_value

"""
#spectrum_real_file_path, - Spectrum class?
output_directory, - only for printing the results
#dim, - number of dimensions (from Spectrum or Params?)
#acc_dim, dim_number - cross section will be across this dimension
#points_dir_list,measured_numbers_of_points - list of numbers of points in each dimension of the spectrum (from Params? - points_dir or fnz)
centered_peak_coords,peak_coords - coordinates of the peak (from Peak?)
    requires function PeakPpmToNumber
peak_name,
peak_number=j

requires implementation of:
    PointNumberInSpectrumFile
"""
class Plot2D: #maybe Plot2D(Data2D)?
    #modify so that there's a decorator(?) that handles adding stuff like labels etc.
    def __init__(self, decorators: list = []):
        #self.id = 0#fig_id
        self.x_list = [0,1,2,3,4,5]
        self.y_list = [x**2 for x in self.x_list]
        #self.plot = plt.figure()
        #plt.plot(x_list, y_list)
        #plt.show()
        #plt.close()
        #self.x_list = 0#x_list
        #self.y_list = 0#y_list
        self.decorators = decorators
        #plt.plot here?
    
    def create_plot_file(self):
        plt.plot(self.x_list, self.y_list)
        
        pass
    
    def change_property(self):
        #xlabel, ylabel, title, color, label, format, bbox_inches?
        pass
    
    def show_plot(self): #lub zamiast tego __call__?
        pass
#should return matplotlib.pyplot.figure(args)


#CreatePlotImageFile - na później
"""fig = plt.figure(layout = "tight")
fig_plot = fig.add_subplot()
fig_plot.set_title("la")
fig_plot.set_xlabel("Chemical shift [ppm]")
fig_plot.set_ylabel("Intensity")
print(fig.get_children())
"""
#fig.plot(x_values, y_values, color = "b")
#fig.savefig(img_path)




















#HEAVILY WORK IN PROGRESS
import figures

def check_cs_overlap(spectrum, accordion_dim, coords, peak, num, output_dir_cs):
    weighted_spectra = [
        "/nhome/bartek/Documents/Doktorat/Research/MBP/MBP_deut_2025_05/151_withNweighting_files/spectrum_real_part1",
        "/nhome/bartek/Documents/Doktorat/Research/MBP/MBP_deut_2025_05/151_withNweighting_files/spectrum_real_part2",
    ]
    #cross_section_for_overlap = cross_sections.CrossSection()
    cross_section_for_overlap = CrossSection()
    cross_section_for_overlap.create_from_spectrum_real(weighted_spectra[0], spectrum, accordion_dim, coords, peak.name)
    more_cs = []
    for spectrum_file in weighted_spectra[1:]:
        #additional_cs = cross_sections.CrossSection()
        additional_cs = CrossSection()
        additional_cs.create_from_spectrum_real(spectrum_file,
            spectrum, accordion_dim, coords, peak.name)
        more_cs.append(additional_cs)
    cross_section_for_overlap.add_cross_sections(more_cs)
    #value_index, value, noise = cross_section_for_overlap.check_weighted_spectrum_overlap_v2()
    interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = cross_section_for_overlap.check_weighted_spectrum_overlap_v3()
    if len(interfering_peaks) > 0:
        #peak_indexes = [peak[0] for peak in interfering_peaks]
        #peak_heights = [peak[1] for peak in interfering_peaks]
        peak_indexes = interfering_peaks
        peak_heights = [cross_section_for_overlap.values[index] for index in interfering_peaks]
    else:
        peak_indexes = None
        peak_heights = None
    
    cs_file_name = (f"CS_p{num}_{peak.name}")
    
    # text output
    #if args.text == 'cs' or args.text == 'both':
    cs_txt_name = cs_file_name + ".txt"


    # figures
    #if args.figures == 'individual' or args.figures == 'both':
    x_values = [i for i in range(0, cross_section_for_overlap.length)]
    y_values = cross_section_for_overlap.values
    xlim = (min(x_values), max(x_values))
    xticks = None #[len(x_values)/2]
    xticks_labels = None #["0"]
    yticks = [0]

    figures.create_2Dplots(
        list_of_x_values = [x_values],
        list_of_y_values = [y_values],
        fig_file_dir = output_dir_cs,
        fig_file_basename = cs_file_name,
        transparent = True,
        xticks = xticks,
        xticks_labels = xticks_labels,
        yticks = yticks,
        xlim = xlim,
        xlabel = r"$\omega_1$",
        ylabel = "intensity",
        list_of_scatter_x_values = [peak_indexes],#[0],#[value_index],
        list_of_scatter_y_values = [peak_heights],#[noise],#[value, noise],
        plot_settings = figures.plot_settings_cross_section_for_paper
    )


def check_cs_overlap_v3(spectrum, accordion_dim, coords, peak, num, output_dir_cs, weighting_values = None):
    weighted_spectra = [
        "/nhome/bartek/Documents/Doktorat/Research/MBP/MBP_deut_2025_05/151_withNweighting_files/spectrum_real_part1",
        "/nhome/bartek/Documents/Doktorat/Research/MBP/MBP_deut_2025_05/151_withNweighting_files/spectrum_real_part2",
    ]
    #cross_section_for_overlap = cross_sections.CrossSection()
    cross_section_for_overlap = CrossSection()
    cross_section_for_overlap.create_from_spectrum_real(weighted_spectra[0], spectrum, accordion_dim, coords, peak.name)
    more_cs = []
    for spectrum_file in weighted_spectra[1:]:
        #additional_cs = cross_sections.CrossSection()
        additional_cs = CrossSection()
        additional_cs.create_from_spectrum_real(spectrum_file,
            spectrum, accordion_dim, coords, peak.name)
        more_cs.append(additional_cs)
    cross_section_for_overlap.add_cross_sections(more_cs)
    
    if weighting_values != None:
        cross_section_for_overlap.new_weighting(weighting_values)
    #value_index, value, noise = cross_section_for_overlap.check_weighted_spectrum_overlap_v2()
    cross_sections_to_compare = cross_section_for_overlap.find_neighbouring_cross_sections(coords, weighted_spectra, spectrum)
    main_factors = []
    for cs in cross_sections_to_compare:
        interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = cs.check_weighted_spectrum_overlap_v3()
        main_factors.append(main_overlap_factor)
    lowest_overlap = min(main_factors)
    ind = main_factors.index(lowest_overlap)
    cross_section_for_overlap = cross_sections_to_compare[ind]
    interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = cs.check_weighted_spectrum_overlap_v3()
    #print(f"Final peak: {main_overlap_factor}, {ind}, {lowest_overlap}, {main_factors}")

    final_cross_section = cross_sections_to_compare[ind]
    if len(interfering_peaks) == 0:
        final_cross_section = cross_section_for_overlap
    
    if len(interfering_peaks) > 0:
        #peak_indexes = [peak[0] for peak in interfering_peaks]
        #peak_heights = [peak[1] for peak in interfering_peaks]
        peak_indexes = interfering_peaks
        peak_heights = [cross_section_for_overlap.values[index] for index in interfering_peaks]
    else:
        peak_indexes = None
        peak_heights = None
    
    cs_file_name = (f"CS_p{num}_{peak.name}")
    
    # text output
    #if args.text == 'cs' or args.text == 'both':
    cs_txt_name = cs_file_name + ".txt"


    # figures
    #if args.figures == 'individual' or args.figures == 'both':
    x_values = [i for i in range(0, cross_section_for_overlap.length)]
    y_values = cross_section_for_overlap.values
    xlim = (min(x_values), max(x_values))
    xticks = None #[len(x_values)/2]
    xticks_labels = None #["0"]
    yticks = [0]

    figures.create_2Dplots(
        list_of_x_values = [x_values],
        list_of_y_values = [y_values],
        fig_file_dir = output_dir_cs,
        fig_file_basename = cs_file_name,
        transparent = True,
        xticks = xticks,
        xticks_labels = xticks_labels,
        yticks = yticks,
        xlim = xlim,
        xlabel = r"$\omega_1$",
        ylabel = "intensity",
        list_of_scatter_x_values = [peak_indexes],#[0],#[value_index],
        list_of_scatter_y_values = [peak_heights],#[noise],#[value, noise],
        plot_settings = figures.plot_settings_cross_section_for_paper
    )
    return final_cross_section


def calc_noise(values):
    # Estimating noise and its standard deviation
    #log the results and calculation etc.
    #ft_results_zeroed_peaks = []
    ft_results_noise = [value for value in values]
    noise_squares_sum = 0.0
    #values = self.y_values
    #global_minimum_index = self.global_minimum_index
    #global_minimum_start_index = global_minimum_index - 40
    #global_minimum_end_index = global_minimum_index + 40

    """global_minimum_start_index = -1
    global_minimum_end_index = -1



    #second_minimum_start_index = second_minimum_index - 40
    #second_minimum_end_index = second_minimum_index + 40
    
    for k in range(0, len(values)):
        if k >= global_minimum_start_index \
            and k <= global_minimum_end_index:
            #or k >= second_minimum_start_index \
            #and k <= second_minimum_end_index:
            ft_results_zeroed_peaks.append(0)
        else:
            ft_results_zeroed_peaks.append(values[k])
            ft_results_noise.append(ft_results_zeroed_peaks[k])
    """
    if len(ft_results_noise) != 0:
        noise_average = sum(ft_results_noise) / float(len(ft_results_noise))
    
        for m in range(0, len(ft_results_noise)):
            noise_squares_sum += (ft_results_noise[m] - noise_average)**2
        noise_standard_deviation = math.sqrt(
            noise_squares_sum/len(ft_results_noise))
    
        minimum_test_value = noise_average - 6*noise_standard_deviation
    else:
        minimum_test_value = 0
    """print(noise_average, noise_squares_sum, noise_standard_deviation, len(ft_results_noise))
    for x in ft_results_noise:
        print(x)
    #print(ft_results_noise)
    exit()"""
    noise_median = median(ft_results_noise)
    return noise_standard_deviation, minimum_test_value, noise_median, noise_average


def select_best_cross_section(cross_sections, peak = None, overlap_data=None, noise_values=None):
    # select the cross section based on the least overlap and best signal to noise ratio
    main_overlap_factors = []
    noise_values = []
    sino_values = []
    for cs in cross_sections:
        try:
            main_overlap_factor = cs.main_overlap_factor
            interfering_peaks = cs.interfering_peaks
            main_peak_index = cs.main_peak_index
            overlap_factor = cs.overlap_factor
        except AttributeError:
            interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = cs.check_weighted_spectrum_overlap_v3()
        main_overlap_factors.append(abs(main_overlap_factor))

        try:
            noise = cs.noise
        except AttributeError:
            #noise, _ = calc_noise(cs.values)
            new_cs = deepcopy(cs)
            new_cs.center_on_peak(peak)
            #noise, _ = calc_noise(cs.orig_values)
            vals_for_noise = new_cs.centered_vals[0:int(cs.length/3)] + \
                new_cs.centered_vals[int(cs.length*2/3) : ]
            noise, _, _, _ = calc_noise(vals_for_noise)
        noise_values.append(noise)
        #sino = cs.values[cs.corr_peak_index] / noise

        #previous sino calc
        #sino = new_cs.values[int(new_cs.length/2)] / noise

        #attempt at new sino calc:
        signal_noise, _, _, _ = calc_noise(new_cs.values[int(cs.length/3) : int(cs.length*2/3)])
        sino = signal_noise / noise
        sino_values.append(abs(sino))
        cs.sino = sino
    lowest_overlap = min(main_overlap_factors)
    ind = main_overlap_factors.index(lowest_overlap)
    
    for i, overlap in enumerate(main_overlap_factors):
        if overlap == lowest_overlap:
            # comparing only noise
            """if noise_values[i] < noise_values[ind]:
                ind = i"""
            # comparing signal to noise ratio
            if sino_values[i] > sino_values[ind]:
                ind = i
    
    final_cs = cross_sections[ind]
    #print(cs.corr_peak_index)
    #print("NOISE:", noise_values, noise_values[ind], min(noise_values))
    #print("SINO:", sino_values, sino_values[ind], max(sino_values))
    #time.sleep(1)
    #try:
    interfering_peaks = final_cs.interfering_peaks
    #main_peak_index = final_cs.main_peak_index
    main_peak_index = final_cs.corr_peak_index
    overlap_factor = final_cs.overlap_factor
    main_overlap_factor = final_cs.main_overlap_factor
    #except AttributeError:
    #    interfering_peaks, main_peak_index, overlap_factor, main_overlap_factor = final_cs.check_weighted_spectrum_overlap_v3()
    method_str = "WIP"
    #print(f"Final peak, {method_str}: {main_overlap_factor}, {ind}, {lowest_overlap}, {main_overlap_factors}")

            #final_cross_section = cross_sections_to_compare[ind]
            #if len(interfering_peaks) == 0:
            #    final_cross_section = cross_section_for_overlap
    final_coords = final_cs.slicing_point_coords

    r"""
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
    )"""
    return final_cs

def calc_sino(cross_section):
    cs = cross_section
    noise_values = []
    vals_for_noise = cs.centered_vals[0:int(cs.length/3)] + \
                cs.centered_vals[int(cs.length*2/3) : ]
    noise, _, _, _ = calc_noise(vals_for_noise)
        #noise_values.append(noise)
        #sino = cs.values[cs.corr_peak_index] / noise

        #previous sino calc
    sino = cs.values[int(cs.length/2)] / noise

        #attempt at new sino calc:
    signal_noise, _, _, _ = calc_noise(cs.values[int(cs.length/3) : int(cs.length*2/3)])
    #sino = abs(signal_noise / noise)
    #sino_values.append(abs(sino))
    cs.sino = sino
    #print("SINONOISE", cs.values[int(cs.length/2)], noise, cs.sino)
    #time.sleep(2)

def calc_sino_profile(profile: Profile):
    y_vals = [-val for val in profile.y_values]
    min_pos = profile.global_minimum_index
    offset = int(min_pos - len(y_vals)/2)
    centered_vals = [y_val for y_val in y_vals]
    if offset < 0:
        for j in range(0, abs(offset)):
            centered_vals.insert(0, centered_vals.pop(-1))
    elif offset > 0:
        for j in range(0, offset):
            centered_vals.insert(-1, centered_vals.pop(0))
    length = len(centered_vals)

    vals_for_noise = centered_vals[0:int(length/5*2)] + \
        centered_vals[int(length*3/5) : ]
    
    noise, test_val, noise_median, noise_average = calc_noise(vals_for_noise)

    peak_indexes, _ = ssig.find_peaks(y_vals)
    peak_indexes = list(peak_indexes)
    #print(peak_indexes, min_pos)
    peak_indexes.remove(min_pos)
    heights = [y_vals[ind] for ind in peak_indexes]
    if len(heights) > 0:
        second_min = max(heights)
        second_min_pos = y_vals.index(second_min)
        sino2 = abs((second_min - noise_average) / noise)
    else:
        sino2 = 0

    sino = abs((y_vals[min_pos] - noise_average) / noise)
    profile.sino = sino
    return sino, sino2


def assign_to_group(cross_section:CrossSection, profile:Profile, peak:peaks.Peak, additional_data):
    # assign based on the overlap and signal to noise in the cross section but also based on the noise of the profile
    sino = cross_section.sino
    overlap_factor = cross_section.overlap_factor
    main_overlap_factor = cross_section.main_overlap_factor
    profile.global_minimum_index
    profile.values = profile.y_values
    #profile._length = len(profile.values)
    profile.dim = 2
    profile.center_on_peak(peak, profile.global_minimum_index)
    #profile_noise = calc_noise
    profile.values = [-val for val in profile.values]
    #profile_sino = calc_sino(profile)
    profile_sino, profile_sino2 = calc_sino_profile(profile)

    if sino < 15:
    #if False:
        group = "D"
    elif overlap_factor >= 0.75:
        group = "C"
    elif overlap_factor == 0:
        group = "A"
    elif overlap_factor < 0.75:
        group = "B"

    '''#if sino < 20:
    if False:
        group = "D"
    elif main_overlap_factor >= 0.75:
        group = "C"
    elif main_overlap_factor == 0:
        group = "A"
    elif main_overlap_factor < 0.75:
        group = "B"'''

    if profile_sino <3:
        group = "D"
    #if profile_sino2 >3:
    #    group = "E"
    #if profile.sino < 2:
    #    group = "E"
#sprawdzanie czy jest podobne minimum w profilu
#sprawdzanie szumu w profilu
#usunąć przeskalowywanie profili
#może warunek że jeśli sino jest przynajmniej 15 i profile_sino jest ok, to wtedy profil jest dobry
# przykład: G24CAi_1-K25Ni-HNi
#jeśli profile_sino2 też jest wysokie to odrzucić (albo porównać samą wysokość pików?):
# sprawdzić czy to pomoże dla K26CAi_1-F27Ni-HNi
#warunek także na main_overlap_factor:
# to pomoże dla T31CA-G32N-G32HN
    return group



def determine_result_validity(cross_section:CrossSection, profile:Profile, peak:peaks.Peak, additional_data):
    sino = cross_section.sino
    overlap_factor = cross_section.overlap_factor
    main_overlap_factor = cross_section.main_overlap_factor
    profile.global_minimum_index
    profile.values = profile.y_values
    profile.dim = cross_section.dim#2
    #profile.center_on_peak(peak, profile.global_minimum_index)
    profile.values = [-val for val in profile.values]
    profile_sino, profile_sino2 = calc_sino_profile(profile)

    if sino < 15:
        validity = False
    elif overlap_factor >= 0.75:
        validity = False
    elif profile_sino < 3:
        validity = False
    else:
        validity = True
    return validity

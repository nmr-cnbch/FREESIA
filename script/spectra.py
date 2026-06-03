"""
This module can be used for parameters of an NMR spectrum and reading 
spectrum_real files.

Classes:
    SpectrumParams: for parameters of an NMR spectrum
    Spectrum: for reading and handling files with spectral data
"""

from __future__ import annotations
import logging
import math
import os
import struct
from typing import Any

from parameters import ParamsTypes, Params

logger = logging.getLogger(__name__)



class SpectrumParams(Params):
    """
    The parameters of a measured NMR spectrum.
     
    Inherits all attributes and methods of the 'Params' class 
    without changes (except for small change to constructor).
    Additionally, it contains new methods specific to parameters
    of an NMR spectrum.
    
    Class attributes:
        spectrum_params_types (dict[str:Any]): Dictionary that assigns
            specific data types to names of spectral parameters.
            Used in the constructor if no ParamsTypes object is passed.
        bioref_values (dict[str:float]): Defines the ratio of 
            the gyromagnetic ratios of a specific nucleus 
            (like 13C or 15N) to the proton (1H)

    Methods:
        __init__(params=None, params_types=None)
        establish_dimensionality()
        check_reffrq()
        calc_center()
        calc_measured_points()
    """
    
    spectrum_params_types = {
        "bioref": int,
        "fnz": int,
        "fp_mult": float,
        "ni": int,
        "np": int,
        "point_dir_begin": int,
        "reffrq": float,
        "sf": float,
        "sw": float,
        "tmax": float,
    }

    bioref_values = {
        "13C": 0.25144953,
        "15N": 0.10132912,
        }
    
    
    def __init__(self,
        params: dict[str: Any] = None,
        params_types: ParamsTypes = None,
        ) -> SpectrumParams:
        """
        Construct the attributes: 'params' and 'params_types'.

        Optional arguments:
            params (dict[str: Any]): Dictionary with parameters
                and their values. If not provided, an empty dictionary
                will be created instead.
            params_types (ParamsTypes): Object of the 'ParamsTypes' 
                class, containing information about the desired
                data types for the parameters. If not provided,
                a new object of 'ParamsTypes' is created using 
                the class attribute 'spectrum_params_types'.

        Raises:
            TypeError: If 'params' argument was provided but it is not
                a dictionary. Or if 'params_types' argument was provided
                but it is not an instance of the ParamsTypes class.
        """
        # using constructor of the parent class
        super().__init__(params, params_types)
        if len(self.params_types._types) == 0:
            self.params_types._types = SpectrumParams.spectrum_params_types


    def establish_dimensionality(self
        ) -> int:
        """
        Establish the dimensionality of the spectrum.

        The dimensionality is based on the presence of 'name' parameters
        inside of the 'self.params' attribute. 
        It defaults to 1 if only one 'name' parameter exists or 
        if no 'name' parameters exist at all.
        If 'dim' was already a defined parameter it is compared to 
        the detected dimensionality. If they are different, the user 
        is asked to input the dimensionality of the experiment.

        For example with 'name', 'name1' and 'name2' present 
        in the parameters, the dimensionality is 3. 

        Returns:
            dim (int): Dimensionality of the spectrum.
        """
        established_dim = 1
        i = 1
        # establish dim based on parameters
        while f"name{i}" in self.params:
            established_dim += 1
            i += 1
        dim = established_dim
        log_msg = ("Based on the parameters, it was established that the "
            f"dimensionality is {dim}.")
        
        # compare established dim with the previously existing dim value
        if "dim" in self.params:
            old_dim = self.params["dim"]
            if established_dim != old_dim:
                input_str = (f"Dimensionality was already defined "
                    f"as {old_dim}, but based on the parameters it got "
                    f"established as {established_dim}. Enter the "
                    "dimensionality value you want to use (integer number): ")
                dim = input(input_str)
                while not isinstance(dim, int):
                    try:
                        dim = int(dim)
                    except ValueError:
                        dim = input(f"Your input: '{dim}' cannot be converted "
                            "to an integer, enter an integer number "
                            "(or type 'quit' to end the program here): ")
                        if dim == "quit" or dim == "'quit'":
                            logger.critical("User's input stopped the program "
                                "while establishing dimensionality.")
                            exit("You typed 'quit', program was stopped.")
                log_msg = ("Based on the parameters, it was established that "
                    f"the dimensionality is {established_dim} but 'dim' "
                    f"parameter was already defined as {old_dim}. "
                    f"User got asked for input and decided dimensionality is "
                    f"{dim}.")
            
        self.params["dim"] = dim
                            
        logger.info(log_msg)
        return dim
    

    def check_reffrq(self
        ) -> list[float]:
        """
        Calculate and replace reffrq values if bioref param equals 1.

        The bioref parameter decides how reffrq values are obtained.
        If bioref equals 1, the reffrq values are calculated based on
        the reffrq for 1H dimension. Otherwise they are obtained 
        directly from the parameters without any change.

        Returns:
            reffrq_values (list[float]): Reference frequencies 
                for all dimensions.

        Raises:
            KeyError: If a nucleus present in the parameters 
                is not defined in the 'bioref_values' dictionary.
        """
        try:
            dim = self["dim"]
        except KeyError:
            dim = self.establish_dimensionality()
        
        bioref_values = SpectrumParams.bioref_values
        reffrq_values = [self.get_value('reffrq', i) for i in range(0, dim)]
        
        try:
            bioref_flag = self["bioref"]
        except KeyError:
            bioref_flag = 0

        if bioref_flag != 1:
            log_msg1 = ("Parameter bioref is not equal 1 "
                "and therefore reference frequency values (reffrq) "
                "were obtained directly from the parameters. "
                f"They are equal: {reffrq_values}")
            logger.info(log_msg1)
            
        else:
            # find reffrq for 1H
            for i in range(0, dim):
                if self.get_value("name", i) == "1H":
                    reffrq_1H = self.get_value("reffrq", i)
                    break

            # calculate reffrq for dimensions other than 1H
            for i in range(0, dim):
                dim_nucleus = self.get_value("name", i)
                
                if dim_nucleus == "1H":
                    continue
                
                elif dim_nucleus in bioref_values:
                    old_reffrq = reffrq_values[i]
                    reffrq_value = bioref_values[dim_nucleus] * reffrq_1H
                    self.set_value("reffrq", i, reffrq_value)

                    log_msg2 = (f"Reference frequency (reffrq) for "
                        f"the dimension {i}, {dim_nucleus} was calculated "
                        f"based on reffrq for 1H, and it is now equal "
                        f"{reffrq_value}. Original value was {old_reffrq}."
                        )
                    logger.info(log_msg2)
                    debug_msg = (f"Dim {i}, {dim_nucleus}: "
                        f"reffrq={reffrq_value} = "
                        f"{bioref_values[dim_nucleus]}*{reffrq_1H}, "
                        f"previous value: {old_reffrq}"
                        )
                    logger.debug(debug_msg)

                else:
                    raise KeyError(f"Unable to calculate the reference "
                        f"frequency for {dim_nucleus} nucleus (dim {i}). "
                        "Its gyromagnetic ratio is not defined in the program."
                        )
        return reffrq_values


    def calc_center(self
        ) -> list[float]:
        """
        Calculate and add the center values for every dimension.
        
        Calculated based on sf and reffrq values.

        Returns:
            center_values (list[float]): Center values for all dimensions.
        """
        center_values = []
        try:
            dim = self["dim"]
        except KeyError:
            dim = self.establish_dimensionality()
            
        for i in range(0, dim):
            dim_nucleus = self.get_value("name", i)
            sf = self.get_value("sf", i)
            reffrq = self.get_value("reffrq", i)
            center = (sf - reffrq)*(10**6)/reffrq
            self.set_value("center", i, center)
            center_values.append(center)

            log_msg = (f"Center point for the dimension "
                f"{i}, {dim_nucleus} was calculated based on "
                f"sf and reffrq values, and it is equal {center}. "
                )
            logger.info(log_msg)
            debug_msg = (f"Dim {i}, {dim_nucleus}: center = {center} = "
                f"({sf}-{reffrq})*(10^6)/{reffrq}."
                )
            logger.debug(debug_msg)
        
        return center_values
    

    def calc_measured_points(self
        ) -> list[int]:
        """
        Calculate and add the mp values for every dimension.

        Parameter 'mp' stands for the number of points measured.
        It is calculated based on sw and tmax.

        Returns:
            mp_values (list[int]): Number of measured points 
                for all dimensions.
        """
        mp_values = []
        try:
            dim = self["dim"]
        except KeyError:
            dim = self.establish_dimensionality()
            
        for i in range(0, dim):
            dim_nucleus = self.get_value("name", i)
            try:
                np_check = self.get_value("np", i)
            except KeyError: 
                np_check = None
            
            # mp equals half of np if np exists
            if np_check != None:
                mp = int(np_check / 2)

            # otherwise mp is calculated from sw and tmax
            else:
                sw = self.get_value("sw", i)
                tmax = self.get_value("tmax", i)
                mp = round(sw * tmax)
                
                log_msg = (f"Number of measured points for the dimension "
                        f"{i}, {dim_nucleus} was calculated based on "
                        f"sw and tmax values, and it is equal {mp}.")
                logger.info(log_msg)
                debug_msg = (f"Dim {i}, {dim_nucleus}: mp = {mp} = "
                            f"round({sw}*{tmax}).")
                logger.debug(debug_msg)
            
            self.set_value("mp", i, mp)
            mp_values.append(mp)
        
        return mp_values



class Spectrum:
    """
    For reading and handling files with spectral data.

    Methods:
        __init__(spectrum_params)
        calc_point_number_in_spectrum_real_file(spectrum_point_coords)
        read_value_in_spectrum_real_file(
            point_number, spectrum_real_file)
        read_many_values_in_spectrum_real_file(point_numbers,
            spectrum_real_file, file)
    """
    def __init__(self,
        spectrum_params: SpectrumParams
        ) -> Spectrum:
        """
        Construct the attribute: 'params'.

        Required arguments:
            spectrum_params (SpectrumParams): Object of the 
                'SpectrumParams' class with the parameters 
                of the spectrum.

        Raises:
            TypeError: If 'spectrum_params' argument was provided
                but it is not an instance of the SpectrumParams class.
        """
        if isinstance(spectrum_params, SpectrumParams):
            self.params = spectrum_params
        else:
            raise TypeError(f"'spectrum_params' argument was provided but "
                f"its type is {type(spectrum_params)}, it should be "
                "an instance of SpectrumParams class instead.")


    def calc_point_number_in_spectrum_real_file(self,
        spectrum_point_coords: list[int],
        ) -> int:
        """
        Calculate the point number in the spectrum_real file
        which corresponds to the provided coordinates in the spectrum.
        
        Required arguments:
            spectrum_point_coords (list[int]): Coordinates of the point.

        Returns:
            point_number (int): Point number in the spectrum file.
        """
        spectrum_params = self.params
        dim = spectrum_params["dim"]
        points_dir = []
        point_number = 0

        for i in range(0, dim):
            if f"points_dir{i}" in spectrum_params.params:
                points_dir.append(spectrum_params.get_value("points_dir", i))
            else:
                points_dir.append(spectrum_params.get_value("fnz", i))
        
        for i in range(0, dim):
            multiplier = 1
            j = i + 1
            while j < dim:
                multiplier = multiplier * points_dir[j]
                j += 1
            #print(i, multiplier)
            point_number += spectrum_point_coords[i] * multiplier
        return point_number
    
    
    def read_value_in_spectrum_real_file(self,
        point_number: int,
        spectrum_real_file: os.PathLike,
        ) -> float:
        """
        Read the value (intensity) from the spectrum_real file.
        
        Required arguments:
            point_number (int): Point number in the spectrum file.
            spectrum_real_file (os.PathLike): Path of the spectrum file.

        Returns:
            spectrum_value (float): Intensity of the spectrum
                at the specified point of the spectrum file.
        
        Raises:
            EOFError: If the point number is too high and leads beyond 
                the file size.
        """
        filesize = os.path.getsize(spectrum_real_file)
        byte_number = point_number * 4
        EOFError_msg = (f"Tried to read byte number {byte_number} from file "
            f"{spectrum_real_file} which contains only {filesize} bytes - "
            "reached the end of file. Check if the parameters file has "
            "correct parameters, especially if 'point_dir_begin' parameters "
            "are correct.")
        if byte_number > filesize:
            raise EOFError(EOFError_msg)
        with open(spectrum_real_file, "rb") as file:
            file.seek(byte_number, 0)
            content = file.read(4)
        spectrum_value = float(struct.unpack(">f", content)[0])
        #try:
            #spectrum_value = float(struct.unpack(">f", content)[0])
        #except struct.error as file_error:
            #raise EOFError(EOFError_msg)
        return spectrum_value


    def read_many_values_in_spectrum_real_file(self,
        point_numbers: list[int],
        spectrum_real_file: os.PathLike,
        file: os.BufferedReader,
        ) -> list[float]:
        """
        Read the values (intensities) from the spectrum_real file.
        
        Required arguments:
            point_numbers (list[int]): Point numbers to read from 
                the spectrum file.
            spectrum_real_file (os.PathLike): Path of the spectrum file.
            file (os.BufferedReader): Open spectrum_real file 
                (in binary mode).

        Returns:
            spectrum_values (list[float]): Intensities of the spectrum
                at the specified points of the spectrum file.
        
        Raises:
            EOFError: If the point number is too high and leads beyond 
                the file size.
        """
        spectrum_values = []
        filesize = os.path.getsize(spectrum_real_file)

        for point_number in point_numbers:
            byte_number = point_number * 4
            EOFError_msg = (f"Tried to read byte number {byte_number} from "
                f"file {spectrum_real_file} which contains only "
                f"{filesize} bytes - reached the end of file. "
                "Check if the parameters file has correct parameters, "
                "especially if 'point_dir_begin' parameters are correct.")
            if byte_number > filesize:
                raise EOFError(EOFError_msg)
            
            file.seek(byte_number, 0)
            content = file.read(4)
            spectrum_value = float(struct.unpack(">f", content)[0])
            #try:
                #spectrum_value = float(struct.unpack(">f", content)[0])
            #except struct.error as file_error:
                #raise EOFError(EOFError_msg)
            spectrum_values.append(spectrum_value)
        return spectrum_values



class UCSF_Spectrum:
    """
    For reading and handling spectra in ucsf sparky format.

    Methods:
        __init__(spectrum_params)
        read_ucsf_header(ucsf_file_path)
        read_ucsf_dimension_header(ucsf_file_path, dim)
        calc_tile_and_point_number(peak_in_points)
        calc_point_in_ucsf(tile_nums, point_in_tile)
        read_value_from_ucsf(point_in_ucsf)
        read_many_values_from_ucsf(points_in_ucsf)
    """
    def __init__(self,
        spectrum_params: SpectrumParams
        ) -> Spectrum:
        """
        Construct the attribute: 'params'.

        Required arguments:
            spectrum_params (SpectrumParams): Object of the 
                'SpectrumParams' class with the parameters 
                of the spectrum.

        Raises:
            TypeError: If 'spectrum_params' argument was provided
                but it is not an instance of the SpectrumParams class.
        """
        if isinstance(spectrum_params, SpectrumParams):
            self.params = spectrum_params
        else:
            raise TypeError(f"'spectrum_params' argument was provided but "
                f"its type is {type(spectrum_params)}, it should be "
                "an instance of SpectrumParams class instead.")


    def read_ucsf_header(self,
        ucsf_file_path: os.PathLike
        ) -> None:
        """
        Getting parameters of the spectrum from ucsf file.
        """
        with open(ucsf_file_path, 'rb') as file:
            # Checking file type:
            file.seek(0)
            file_type_raw = file.read(10)
            file_type = file_type_raw.decode('utf-8')
            file_type = file_type.strip()
            file_type = file_type.replace("\0", "")
            if file_type != "UCSF NMR":
                raise ValueError("The ucsf file seems wrong, its beginning is: "
                    f"{file_type}, while it should be 'UCSF NMR'.")
            
            # Dimensionality: 
            file.seek(10)
            dimensionality_raw = file.read(1)
            dim = int.from_bytes(dimensionality_raw, byteorder='big')
            #self.params["dim"] = dim
            self.params.set_value("dim", 0, dim)

            # Number of data components:
            file.seek(11)
            data_components_raw = file.read(1)
            data_components = int.from_bytes(data_components_raw, 
                byteorder='big')
            if data_components != 1:
                raise ValueError("Unsupported number of components detected: "
                    f"{data_components}. Expected '1' for real data only.")

            # Format version:
            file.seek(13)
            format_version_raw = file.read(1)
            format_version = int.from_bytes(format_version_raw, byteorder='big')
            if format_version != 2:
                raise ValueError("Unsupported format version detected: "
                    f"{format_version}. Expected: '2'.")
            
            #file.seek(14)
            #remaining_header = file.read(166)
            #remaining_header = remaining_header.decode()
            #print(remaining_header)

    
    def read_ucsf_dimension_header(self,
        ucsf_file_path: os.PathLike,
        dim: int,
    ) -> tuple[str, int, int, float, float, float]:
        """
        Getting parameters of the spectrum from ucsf file.
        """
        sp = 180 + dim*128 #starting point in the file
        params = self.params
        dimensionality = self.params["dim"]
        if dim == dimensionality - 1:
            params_dim = 0
        else:
            params_dim = dim + 1
        with open(ucsf_file_path, 'rb') as file:
            # Nucleus/name:
            file.seek(sp)
            nucleus = file.read(6)
            nucleus = nucleus.decode()
            nucleus = nucleus.replace("\0", "")
            params.set_value("name", params_dim, nucleus)

            # Number of points:
            file.seek(sp + 8)
            fnz = file.read(4)
            fnz = int.from_bytes(fnz, byteorder='big')
            params.set_value("fnz", params_dim, fnz)

            # Tile size:
            file.seek(sp + 16)
            tile_size = file.read(4)
            tile_size = int.from_bytes(tile_size, byteorder='big')
            params.set_value("tile_size", params_dim, tile_size)

            # Spectrometer frequency [MHz]:
            file.seek(sp + 20)
            reffrq = file.read(4)
            reffrq = struct.unpack('>f', reffrq)[0]
            params.set_value("reffrq", params_dim, reffrq)

            # Spectral width [Hz]:
            file.seek(sp + 24)
            sw = file.read(4)
            sw = struct.unpack('>f', sw)[0]
            params.set_value("sw", params_dim, sw)

            # Center of data [ppm]:
            file.seek(sp + 28)
            center = file.read(4)
            center = struct.unpack('>f', center)[0]
            params.set_value("center", params_dim, center)
        
        
        #print(nucleus, fnz, tile_size, reffrq, sw, center)
        return nucleus, fnz, tile_size, reffrq, sw, center
    

    def read_ucsf_all_headers(self,
        ucsf_file_path: os.PathLike,
        ) -> None:
        self.read_ucsf_header(ucsf_file_path)
        dim = self.params["dim"]
        for i in range(0, dim):
            dim_output = self.read_ucsf_dimension_header(ucsf_file_path, i)
    
    @classmethod
    def read_ucsf_file(cls,
        ucsf_file_path: os.PathLike,
        ) -> UCSF_Spectrum:
        UCSF_Spectrum_obj = cls(SpectrumParams())
        UCSF_Spectrum_obj.read_ucsf_all_headers(ucsf_file_path)
        params = UCSF_Spectrum_obj.params
        dim = params["dim"]
        tile_sizes = []
        fnz_list = []
        numbers_of_tiles = []
        for i in range(1, dim):
            tile_size = params.get_value("tile_size", i)
            tile_sizes.append(tile_size)
            fnz = params.get_value("fnz", i)
            fnz_list.append(fnz)
            if fnz % tile_size == 0:
                number_of_tiles = int(fnz/tile_size)
            else:
                number_of_tiles = int(fnz/tile_size + 1)
            numbers_of_tiles.append(number_of_tiles)
        tile_size0 = params.get_value("tile_size", 0)
        tile_sizes.append(tile_size0)
        fnz0 = params.get_value("fnz", 0)
        fnz_list.append(fnz0)
        if fnz0 % tile_size0 == 0:
            number_of_tiles0 = int(fnz0/tile_size0)
        else:
            number_of_tiles0 = int(fnz0/tile_size0 + 1)
        numbers_of_tiles.append(number_of_tiles0)
        full_tile_size = math.prod(tile_sizes)
        UCSF_Spectrum_obj.tile_sizes = tile_sizes
        UCSF_Spectrum_obj.fnz_list = fnz_list
        UCSF_Spectrum_obj.numbers_of_tiles = numbers_of_tiles
        UCSF_Spectrum_obj.full_tile_size = full_tile_size
        return UCSF_Spectrum_obj


    def calc_tile_and_point_number(peak_in_points):
        pass


    def calc_point_in_ucsf(self,
        spectrum_point_coords: list[int],
        ) -> int:
        """
        Calculate the point number in the spectrum_real file
        which corresponds to the provided coordinates in the spectrum.
        
        Required arguments:
            spectrum_point_coords (list[int]): Coordinates of the point.

        Returns:
            point_number (int): Point number in the spectrum file.
        """
        params = self.params
        tile_sizes = self.tile_sizes
        fnz_list = self.fnz_list
        numbers_of_tiles = self.numbers_of_tiles
        full_tile_size = self.full_tile_size

        spectrum_point_coords.append(spectrum_point_coords.pop(0))
        #print(spectrum_point_coords)
        #exit()
        tile_number = 0
        remaining_tile_sizes = 0
        for i, point_coord in enumerate(spectrum_point_coords):
            tile_mult = 1
            tile_numbers = numbers_of_tiles[i+1:]
            for num in tile_numbers:
                tile_mult = tile_mult * num
            #tile_mult = [tile_mult * tile_number for tile_number in tile_numbers]
            tile_num = point_coord // tile_sizes[i]
            tile_number += tile_num * tile_mult

            rem = 1
            rems = tile_sizes[i+1:]
            for x in rems:
                rem = rem * x
            point_within_tile = point_coord % tile_sizes[i]
            remaining_tile_sizes += point_within_tile * rem
            #print(tile_num, point_within_tile)
        
        tile_pos = full_tile_size * tile_number

        ppos = tile_pos + remaining_tile_sizes
        #print(ppos)
        return ppos
    
    
    def read_value_from_ucsf(self,
        point_in_ucsf: int,
        ucsf_file_path: os.PathLike,
        ) -> float:
        """
        Read the value (intensity) from the spectrum_real file.
        
        Required arguments:
            point_number (int): Point number in the spectrum file.
            spectrum_real_file (os.PathLike): Path of the spectrum file.

        Returns:
            spectrum_value (float): Intensity of the spectrum
                at the specified point of the spectrum file.
        
        Raises:
            EOFError: If the point number is too high and leads beyond 
                the file size.
        """
        dim = self.params["dim"]
        filesize = os.path.getsize(ucsf_file_path)
        #byte_number = point_number * 4
        byte_number = 180 + 128*dim + 4 * point_in_ucsf
        EOFError_msg = (f"Tried to read byte number {byte_number} from file "
            f"{ucsf_file_path} which contains only {filesize} bytes - "
            "reached the end of file.") 
        '''Check if the parameters file has "
            "correct parameters, especially if 'point_dir_begin' parameters "
            "are correct.")'''
        if byte_number > filesize:
            raise EOFError(EOFError_msg)
        with open(ucsf_file_path, 'rb') as file:
            file.seek(byte_number)
            value_binary = file.read(4)
        value = struct.unpack('>f', value_binary)[0]
        return value


    def read_many_values_from_ucsf(self,
        points_in_ucsf: list[int],
        ucsf_file_path: os.PathLike,
        file: os.BufferedReader,
        ) -> list[float]:
        """
        Read the values (intensities) from the spectrum_real file.
        
        Required arguments:
            point_numbers (list[int]): Point numbers to read from 
                the spectrum file.
            spectrum_real_file (os.PathLike): Path of the spectrum file.
            file (os.BufferedReader): Open spectrum_real file 
                (in binary mode).

        Returns:
            spectrum_values (list[float]): Intensities of the spectrum
                at the specified points of the spectrum file.
        
        Raises:
            EOFError: If the point number is too high and leads beyond 
                the file size.
        """
        spectrum_values = []
        dim = self.params["dim"]
        filesize = os.path.getsize(ucsf_file_path)

        for point_number in points_in_ucsf:
            byte_number = 180 + 128*dim + 4 * point_number
            EOFError_msg = (f"Tried to read byte number {byte_number} from "
                f"file {ucsf_file_path} which contains only "
                f"{filesize} bytes - reached the end of file. "
                "Check if the peak list contains chemical shifts outside of "
                "the appropriate range.")
            if byte_number > filesize:
                raise EOFError(EOFError_msg)
            
            file.seek(byte_number, 0)
            value_binary = file.read(4)
            spectrum_value = float(struct.unpack(">f", value_binary)[0])
            #try:
                #spectrum_value = float(struct.unpack(">f", content)[0])
            #except struct.error as file_error:
                #raise EOFError(EOFError_msg)
            spectrum_values.append(spectrum_value)
        return spectrum_values


if __name__ == "__main__":
    ucsf_file_path = "/nhome/bartek/Documents/Doktorat/Python_scripts/classes/COdecProfiles_v5/ucsf_read_test_files/ucsf_sum_new.ucsf"
    params = SpectrumParams()
    ucsf = UCSF_Spectrum.read_ucsf_file(ucsf_file_path)#UCSF_Spectrum(params)
    params = ucsf.params
    
    ucsf.read_ucsf_header(ucsf_file_path)
    for d in range(0, ucsf.params["dim"]):
        ucsf.read_ucsf_dimension_header(ucsf_file_path, dim=d)
    #print(params)
    for param in params.params:
        print(param, params.params[param])
    
    obj = UCSF_Spectrum.read_ucsf_file(ucsf_file_path)
    print(obj, obj.tile_sizes, obj.fnz_list)
    shifts = [134,160,256]
    ppos = obj.calc_point_in_ucsf(shifts)
    val = obj.read_value_from_ucsf(ppos, ucsf_file_path)
    print(ppos, val)
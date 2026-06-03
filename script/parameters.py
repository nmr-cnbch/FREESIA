"""
This module can be used for parameters of any object or process.

Classes:
    ParamsTypes: for types of parameters
    Params: for parameters of any object or process etc.
"""

from __future__ import annotations
import logging
import os
from typing import Any

logger = logging.getLogger(__name__)



class ParamsTypes:
    """
    The data types of expected parameters like int, str etc.

    This class can be used to store the types of parameters, for example
    so that the values of parameters can be easily converted to 
    appropriate types. Might be especially useful if reading parameters 
    from files - to convert strings into other data types.

    Attributes:
        _types (dict[str:Any]): Dictionary that assigns a specific 
            data type to the names of parameters.

    Methods:
        __init__(params_types)
        __contains__(param)
        __getitem__(param)
        __setitem__(param, param_type)
    """

    
    def __init__(self,
        params_types: dict[str:Any],
        ) -> ParamsTypes:
        """
        Assign the provided dictionary to '_types' attribute.

        Required arguments:
            params_types (dict[str:Any]): Dictionary that assigns 
                a specific data type to the name of the parameter.
                For example: param_name: float

        Returns:
            None

        Raises:
            TypeError: If provided input 'params_types' is not a dict.
        """
        if isinstance(params_types, dict):
            self._types = params_types
        else:
            raise TypeError("The types of parameters should be a dictionary!")
        
    
    def __contains__(self, 
        param: str,
        ) -> bool:
        """
        Behaviour for 'in' operator.
        
        Checks if the parameter has its type assigned.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            bool
        """
        return param in self._types


    def __getitem__(self,
        param: str,
        ) -> Any:
        """
        Return the data type assigned to the name of the parameter.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            param_type (Any): Data type associated with the parameter.
                If no type was associated, it defaults to str.
        """
        try:
            param_type = self._types[param]
        except KeyError:
            param_type = str
            msg = (f"The desired data type of parameter '{param}' is unknown "
                "so str type was used instead.")
            logger.warning(msg)
        
        return param_type


    def __setitem__(self,
        param: str,
        param_type: Any,
        ) -> None:
        """
        Assign the provided data type to the provided parameter name.

        Required arguments:
            param (str): Name of the parameter.
            param_type (Any): Data type to be associated with
                the parameter.

        Returns:
            None
        """
        #if param in self._types:
        if param in self:
            previous_type = self._types[param]
            param_type_existed_msg = (f"The parameter {param} was already "
                f"assigned the data type: {previous_type.__name__}. ")
        else:
            param_type_existed_msg = ""
        self._types[param] = param_type
        msg = (f"{param_type_existed_msg}Data type of the parameter {param} "
            f"is now: {param_type.__name__}.")
        logger.debug(msg)



class Params:
    """
    For storing and handling parameters and their values.
    
    Created specifically for multi-dimensional situations, 
    when parameters should share the same name but apply to different 
    dimensions, for example 'velocity', 'velocity1' and 'velocity2'
    can be names for velocities for dimensions that correspond to numbers
    0, 1 and 2, respectively.

    Attributes:
        params (dict[str:Any]): Dictionary that assigns values
            to names of the parameters.
        params_types (ParamsTypes): Object of the 'ParamsTypes' class,
            containing information about the desired data types
            for the parameters.

    Methods:
        __init__(params=None, params_types=None)
        __contains__(param)
        __getitem__(param)
        __setitem__(param, param_value)
        get_param_corename_number(param)
        get_value(param_core_name, dim)
        set_value(param_core_name, dim, param_value)
        correct_param_type(param)
        read_parameters_file(file_path)
        
    TODO: add a method of saving parameters into a file
    """


    def __init__(self,
        params: dict[str: Any] = None,
        params_types: ParamsTypes = None,
        ) -> Params:
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
                an empty dictionary.

        Raises:
            TypeError: If 'params' argument was provided but it is not
                a dictionary. Or if 'params_types' argument was provided
                but it is not an instance of the ParamsTypes class.
        """
        if isinstance(params, dict):
            self.params = params
        elif params is None:
            self.params = {}
        else:
            raise TypeError(f"'params' argument was provided but "
                f"its type is {type(params_types)}, it should be "
                "a dictionary instead.")
        
        if isinstance(params_types, ParamsTypes):
            self.params_types = params_types
        elif params_types is None:
            self.params_types = ParamsTypes({})
        else:
            raise TypeError(f"'params_types' argument was provided but "
                f"its type is {type(params_types)}, it should be "
                "an instance of ParamsTypes class instead.")
        

    def __contains__(self, 
        param: str,
        ) -> bool:
        """
        Behaviour for 'in' operator.
        
        Checks if the parameter exists in the 'params' attribute.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            bool
        """
        return param in self.params


    def __getitem__(self,
        param: str,
        ) -> Any:
        """
        Return the value of the specified parameter.

        The value is obtained from the dictionary in 'params' attribute.
        If the data type for this parameter is specified within
        the 'params_types' attribute then the value is converted
        to that type, otherwise it is returned as is, without conversion.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            param_value (Any): Value of the specified parameter.

        Raises:
            KeyError: If the specified parameter does not exist.
        """
        if param in self: 
            param_raw_value = self.params[param]
            if param in self.params_types: 
                param_type = self.params_types[param]
                param_value = param_type(param_raw_value)
            else:
                param_value = param_raw_value
            
        else:
            msg = (f"Parameter '{param}' does not exist "
                "and therefore its value is unknown.")
            logger.warning(msg)
            raise KeyError(msg)
        return param_value


    def __setitem__(self,
        param: str,
        param_value: Any,
        ) -> None: 
        """
        Assign the provided value to the parameter.

        Required arguments:
            param (str): Name of the parameter.
            param_value (Any): The value of the parameter.

        Returns:
            None
        """
        self.params[param] = param_value


    def get_param_corename_number(self,
        param: str,
        ) -> tuple[str, int]:
        """
        Return the core name and number of the parameter.

        Parameters can be built of their core names and a number.
        The number at the end can for example determine the order 
        or dimension with which the specific parameter is associated.
        For example, parameter 'phase2' has the core name 'phase' 
        and the number 2, which may indicate it is the phase for 
        the dimension associated with that number.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            param_core_name (str): Core name of the parameter 
                after separating it from the number.
            param_number (int): Number of the parameter. Equals zero 
                if the parameter had no number at the end.
        """
        param_core_name = param
        param_number = 0
        for i in range(len(param) - 1, -1, -1):
            if param[i:].isdigit():
                param_core_name = param[0:i]
                param_number = int(param[i:])
                continue
                
        """for index, char in enumerate(param):
            if param[index:].isdigit():
                param_core_name = param[0:index]
                param_number = int(param[index:])
                break
            """
        return param_core_name, param_number
    
    
    def get_value(self,
        param_core_name: str,
        dim: int,
        ) -> Any:
        """
        Return the value of the parameter for the specified dimension.

        The value is obtained from the dictionary in 'params' attribute.
        If the data type for this parameter is specified within
        the 'params_types' attribute then the value is converted
        to that type, otherwise its type remains unchanged.        

        Required arguments:
            param_core_name (str): Core name of the parameter 
                after separating it from the number.
            dim (int): Number/dimension associated with the parameter.

        Returns:
            param_value (Any): Value of the specified parameter.
        """
        if dim == 0:
            param = param_core_name
        else:
            param = param_core_name + str(dim)

        param_raw_value = self[param]
        # check if parameter's type is defined:
        if param_core_name in self.params_types: 
            param_type = self.params_types[param_core_name]
            param_value = param_type(param_raw_value)
        else:
            param_value = param_raw_value
        return param_value


    def set_value(self,
        param_core_name: str,
        dim: int,
        param_value: Any,
        ) -> None: 
        """
        Assign the provided value to the parameter.

        Required arguments:
            param_core_name (str): Core name of the parameter.
            dim (int): Number/dimension of the parameter.
            param_value (Any): Value of the parameter.

        Returns:
            None
        """
        if dim == 0:
            param = param_core_name
        else:
            param = param_core_name + str(dim)
        self.params[param] = param_value


    def correct_param_type(self,
        param: str,
        ) -> Any:
        """
        Correct the parameter's value by changing its type.

        The type of the parameter's value is corrected according to 
        the type associated with the parameter's core name, 
        as specified within the 'params_types' attribute.
        If no desired type can be found for the specified parameter, 
        its value remains unchanged.

        Required arguments:
            param (str): Name of the parameter.

        Returns:
            param_value_corrected (Any): Parameter value after 
                adjusting its type (if possible).
        """
        param_value = self[param]
        param_core_name, param_number = self.get_param_corename_number(param)
        if param_core_name in self.params_types:
            param_type = self.params_types[param_core_name]
            param_value_corrected = param_type(param_value)
            self[param] = param_value_corrected
        else:
            param_value_corrected = param_value
            msg = f"The desired type of parameter '{param}' is unknown"
            logger.info(msg)
        return param_value_corrected

        """
        if param in self.params:
            param_value = self.params[param]
            param_without_number, param_number = self.get_param_corename_number(param)
            if param_without_number in self.params_types._types:
                param_type = self.params_types[param_without_number]
                param_value_corrected = param_type(param_value)
                self.params[param] = param_value_corrected
            else:
                msg = f"The desired type of parameter '{param}' is unknown"
                #print(msg)
                logger.info(msg)
        else:
            msg = f"Parameter '{param}' does not exist"
            #print(msg)
            logger.warning(msg)
            #exit()"""            


    def read_parameters_file(self,
        file_path: os.PathLike,
        ) -> dict[str:Any]:
        """
        Read and save the parameters from the file.

        Parameters are saved into the 'params' attribute. If possible, 
        their type is adjusted according to the contents 
        of the 'params_types' attribute, using correct_param_type method.
        
        The following structure of the parameters file is required:
        1. Every parameter and its value are located within the same 
        file line.
        2. The line begins with the parameter name, followed by 
        empty space and the value of the parameter.
        3. There can be many values for one parameter name, each value
        separated from each other by spaces. Then the values are put 
        into a list and such a list becomes the value associated 
        with the parameter name.
        4. The file can have empty lines - they are ignored.

        Example of the input file contents:
        Name tulip
        Colour red
        Amount 12
        Name1 freesia
        Colour1 purple orange
        Amount1 25 32
        
        If the 'Amount' parameter is assigned type int within the
        'params_types' attribute, then the dictionary created from 
        the above input file will look like this:
        {
        'Name': 'tulip'
        'Colour': 'red'
        'Amount': 12
        'Name1': 'freesia'
        'Colour1': ['purple', 'orange']
        'Amount1': [25, 32]
        }

        Required arguments:
            file_path (os.PathLike): Path of the file with 
                the parameters and their values.

        Returns:
            self.params (dict[str:Any]): Dictionary with parameters 
                and their values.

        Raises:
            FileNotFoundError: If provided file path does not lead
                to an existing file.
        """ 
        file_path = os.path.abspath(file_path)
        if not os.path.isfile(file_path):
            error_msg = (f"This path was searched for an input file but "
                f"no file was found: {file_path}")
            logger.critical(error_msg)
            raise FileNotFoundError(error_msg)
        else:
            with open(file_path) as file:
                lines = file.readlines()
            for line in lines:
                line = line.strip()
                if line == "":
                    continue
                words = line.split()
                param = words[0]
                if len(words) > 2:
                    param_value = words[1:]
                else:
                    param_value = words[1]
                self.params[param] = param_value
                self.correct_param_type(param)
        msg = f"Finished reading the parameters file: {file_path}"
        logger.info(msg)
        return self.params
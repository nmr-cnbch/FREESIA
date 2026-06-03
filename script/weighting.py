"""
This module is used for weighting functions and calculating weights.

Functions:
    find_float_in_string_ending
    select_weighting_function_from_str
    weight_func_vals

Classes:
    Cos_To_Nth_Power
    One_Minus_Sine_To_Nth_Power
    Gauss
    Exp
"""

from __future__ import annotations
from abc import ABC, abstractmethod

import logging
import math

logger = logging.getLogger(__name__)



def find_float_in_string_ending(
    string: str,
    ) -> tuple[float|False, str]:
    """
    Find the number located at the end of the string.

    Required arguments:
        string (str): String which ends with a number.
    
    Returns:
        num (float|False): Number found at the end of the input string. 
            If no number is found, False is returned.
        shortened_string (str): Original string without the number 
            found at the end.
    """
    shortened_string = string
    for i in range(0, len(string)):
        try:
            #print(f"try {i}")
            num = float(string[i:])
        except:
            #print(f"except {i}")
            num = False
            continue
        else:
            #print(f"else {i}")
            shortened_string = string[0:i]
            break
    return num, shortened_string



def find_float_inside_brackets(
    string: str,
    ) -> tuple[float|False, str]:
    """
    Find the number located inside the brackets [] within the string.

    Required arguments:
        string (str): String which contains a number inside the brackets.
    
    Returns:
        num (float|False): Number found within the brackets [] within 
            the input string. If no number is found, False is returned.
        shortened_string (str): Modified string with the brackets
            and the number removed.
    """
    shortened_string = string
    start = string.find("[")
    end = string.find("]")
    #print(string, start, end)
    if start != -1 and end != -1:
        num = float(string[start+1:end])
        shortened_string = string[0:start] + string[end+1:]
    else:
        num = False
        shortened_string = string
    return num, shortened_string



class WeightingFunction(ABC):
    @abstractmethod
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def __call__(self, x):
        pass



class Cos_To_Nth_Power(WeightingFunction):
    def __init__(self, n, endpoint):
        self.power = n
        self.endpoint = endpoint

    def __call__(self, x):
        result = (math.cos(math.pi/self.endpoint * x)) ** self.power
        return result



class One_Minus_Sine_To_Nth_Power(WeightingFunction):
    def __init__(self, n, endpoint):
        self.power = n
        self.endpoint = endpoint

    def __call__(self, x):
        result = 1 - (math.sin(math.pi/self.endpoint * x)) ** self.power
        return result



class Gauss(WeightingFunction):
    def __init__(self, sigma):
        self.sigma = sigma

    def __call__(self, x):
        result = math.exp(-0.5*(x)**2 / (self.sigma**2))
        return result



class Exp(WeightingFunction):
    def __init__(self, n):
        self.n = n

    def __call__(self, x):
        result = math.exp(-abs(x)/self.n)
        return result



def select_weighting_function_from_str(
    weighting_function_str: str,
    endpoint: int = None,
    ) -> tuple[WeightingFunction, float]:
    """
    Return the weighting function based on the input string.

    Required arguments:
        weighting_function_str (str): String which starts with the name 
            of one of the programmed functions, and ends with the value 
            of the function's main parameter. For example 'cos2' for 
            the function cosine to the power of 2.
    
    Optional arguments:
        endpoint (int): Used by some of the functions to determine 
            the point at which the function should be equal to or 
            approach zero. Defaults to None.
    
    Returns:
        weighting_function (WeightingFunction): Instance of the desired 
            weighting function class, that can be used to calculate 
            the weights.
        param (float): Parameter of the function, extracted from 
            the input string.
    
    Raises:
        ValueError: If the function is not recognized.
    """
    #param, function = find_float_in_string_ending(weighting_function_str)
    param, function = find_float_inside_brackets(weighting_function_str)
    #print(param, function)
    if function == "c":
        weighting_function = Cos_To_Nth_Power(param, endpoint)
    elif function == "s":
        weighting_function = One_Minus_Sine_To_Nth_Power(param, endpoint)
    elif function == "g":
        if not param: param = 5
        sigma = (endpoint/2)/param
        weighting_function = Gauss(sigma)
    elif function == "e":
        if not param: param = 50
        weighting_function = Exp(param)
    else:
        raise ValueError(f"{function} is an unknown weighting function")
    return weighting_function, param



def weight_func_vals(
    weighting_function: WeightingFunction,
    number_of_points: int,
    ) -> list[float]:
    """
    Returns the list of calculated weights.

    Required arguments:
        weighting_function (WeightingFunction): Instance of the desired 
            weighting function class.
        number_of_points (int): Number of weights to calculate.
    
    Returns:
        weight_function_vals (list[float]): List of calculated weights.
    """
    weight_function_vals = []
    for i in range(0, number_of_points):
        point_number = i - number_of_points/2
        value = weighting_function(point_number)
        weight_function_vals.append(value)
    return weight_function_vals

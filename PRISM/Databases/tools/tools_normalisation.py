import math
import numpy as np

def normalisation_sin(value):
    """
    Function performes a sin normalisation for the provided value.

    :param value: value, for which a sin normalisation is undertaken for
    :return normalised_value: sin normalised value
    """

    normalised_value = math.sin(value)
    
    return normalised_value

def normalisation_cos(value):
    """
    Function performes a cos normalisation for the provided value.

    :param value: value, for which a cos normalisation is undertaken for
    :return normalised_value: cos normalised value
    """

    normalised_value = math.cos(value)
    
    return normalised_value

def normalisation_natural_log(value):
    """
    Function performes a natural log normalisation for the provided value.

    :param value: value, for which a natural log normalisation is undertaken for
    :return normalised_value: natural log normalised value
    """

    normalised_value = np.log(value)
    
    return normalised_value

def normalisation_radial_basis_function(value, mean, variance):
    """
    Function performes a radial basis function normalisation for the provided value.

    :param value: value, for which a radial basis function normalisation is undertaken for
    :param mean: mean of all the group of values, which are normalised with this function
    :param variance: variance of all the group of values, which are normalised with this function
    :return normalised_value: radial basis function normalised value
    """
    
    normalised_value = math.exp(-((mean-value)**2/variance))
    
    return normalised_value

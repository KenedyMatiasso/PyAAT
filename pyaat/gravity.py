"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2020 Kenedy Matiasso Portella
Distributed under MIT License

This file implements the gravity models.
"""

import numpy as np
from numpy import radians, sin, cos, array, sqrt
from abc import ABC, abstractmethod

from pyaat.constants import GRAVITY, MU_EARTH, R_EARTH, J2, J3, J4


class metaGravity(ABC):
    """
    Metaclass for the gravity models. The interface for all models is exactly the same.
    """
    def __init__(self, altitude = 0.0, latitude = radians(45), longitude = 0.0):
        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        self._latitude = latitude
        self._altitude = altitude
        self._longitude = longitude

    def set_altitude(self, entry):
        self._altitude = entry

    def get_altitude(self):
        return self._altitude
    
    def set_longitude(self, entry):
        self._longitude = entry

    def get_longitude(self):
        return self._longitude

    def set_latitude(self, entry):
        self._latitude = entry
    
    def get_latitude(self):
        return self._latitude

    def get_gravity(self):
        return self._gravity

    @property
    @abstractmethod
    def _gravity(self):
        pass

class Earth_VerticalConstant(metaGravity):
    """
    Vertical constant gravity
    
    """
    def __init__(self, altitude = 0.0, latitude = radians(45), longitude = 0.0):
        super().__init__(altitude, latitude, longitude)

    @property
    def _gravity(self):
        return array([0, 0, GRAVITY])      
        
class Earth_NewtonGravity(metaGravity):
    """
    Newton gravity
    g = -G*M*\vec(r)/r³
    
    Reference
    ---------
    A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.
    
    """
    def __init__(self, altitude = 0.0, latitude = radians(45), longitude = 0.0):
        super().__init__(altitude, latitude, longitude)
    
    @property
    def _gravity(self):
        gr = MU_EARTH/(R_EARTH+self._altitude)**2
        return array([0.0, 0.0, gr])

class Earth_highOrder(metaGravity):
    """
    Description
    ----------
    Implements a high order gravitational model based on spheric harmonics.
    
    Reference
    ---------
    A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.

    """
    def __init__(self, altitude = 0.0, latitude = radians(45), longitude = 0.0):
        super().__init__(altitude, latitude, longitude)
        
    @property
    def _gravity(self):

        radius = R_EARTH + self._altitude
        co_latitude = radians(90) - self._latitude

        gphi = 3*MU_EARTH*R_EARTH**2/(radius)**4*sin(co_latitude)*cos(co_latitude)*(J2+0.5*J3*(R_EARTH/radius)/cos(co_latitude)*(5*cos(co_latitude)**2-1)+5/6*J4*(R_EARTH/radius)**2*(7*cos(co_latitude)**2-1))
        P2 = 1/2*(3*cos(co_latitude)**2-1)
        P3 = 1/2*(5*cos(co_latitude)**3-3*cos(co_latitude))
        P4 = 1/8*(35*cos(co_latitude)**4-30*cos(co_latitude)**2+3)
        gr = MU_EARTH/radius**2*(1-3*J2*(R_EARTH/radius)**2*P2-4*J3*(R_EARTH/radius)**3*P3-5*J4*(R_EARTH/radius)**4*P4)
        return array([0, gphi, gr])
"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2020 Kenedy Matiasso Portella
Distributed under MIT License

Version control:
21/10/2020                          KMP: File creation

This file implements the gravity models that can be used by PyAAT.
"""

import numpy as np
from pyaat.constants import GRAVITY, MU_EARTH, R_EARTH, J2, J3, J4
from numpy import radians, sin, cos, array, sqrt

# this class deppends on altitude, longiude and latidude
# It also deppends from de attitude
     
class VerticalConstant(object):
    """
    Vertical constant gravity
    
    """
    def __init__(self):
        self._magnitude = GRAVITY
        self._versor = np.array([0, 0, 1], dtype=float)
        self._altitude = 0.0
        self._latitude = 0.0
        self._longitude = 0.0
        
    @property
    def _gravity(self):
        return self._magnitude*self._versor
    
    
class NewtonGravity(object):
    """
    Newton gravity
    g = -G*M*\vec(r)/r³
    
    Reference
    ---------
    A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.
    
    """
    def __init__(self):
        self._versor = np.array([0, 0, 1], dtype=float)
        self._altitude = 0.0
        self._altitude = 0.0
        self._latitude = 0.0
        self._longitude = 0.0
        
    @property
    def _magnitude(self):
        return MU_EARTH/(R_EARTH+self._altitude)**2

    @property
    def _gravity(self):
        return self._magnitude*self._versor

class HighOrder(object):
    """
    Description
    ----------
    Implements a high order gravitational model based on spheric harmonics.
    
    
    
    Reference
    ---------
    A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.

    """
    def __init__(self):
        self._altitude = 0.0
        self._latitude = 0.0
        self._longitude = 0.0
    
    @property
    def phi(self):
        """
        Returns
        -------
        phi : float
            co-latitude angle in rad.
            
        """
        phi = radians(90)-self._latitude
        return phi
    
    @property
    def r(self):
        return R_EARTH+self._altitude
    
    @property
    def _gravity(self):
        gphi = 3*MU_EARTH*R_EARTH**2/(self.r)**4*sin(self.phi)*cos(self.phi)*(J2+0.5*J3*(R_EARTH/self.r)/cos(self.phi)*(5*cos(self.phi)**2-1)+5/6*J4*(R_EARTH/self.r)**2*(7*cos(self.phi)**2-1))
        P2 = 1/2*(3*cos(self.phi)**2-1)
        P3 = 1/2*(5*cos(self.phi)**3-3*cos(self.phi))
        P4 = 1/8*(35*cos(self.phi)**4-30*cos(self.phi)**2+3)
        gr = MU_EARTH/self.r**2*(1-3*J2*(R_EARTH/self.r)**2*P2-4*J3*(R_EARTH/self.r)**3*P3-5*J4*(R_EARTH/self.r)**4*P4)
        return array([0,gphi,gr])
    
    @property
    def _magnitude(self):
        return sqrt(self._gravity[0]**2+self._gravity[1]**2+self._gravity[2]**2)
    
    @property
    def _versor(self):
        return self._gravity/self._magnitude
"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Stores information about the aircraft
"""

from numpy import array
from constants import RHO_SEA
from numpy.linalg import inv

class airplane(object):
    def __init__(self):
        self._mass = 45e3
        self.Ixx = 0.554e6
        self.Iyy = 2.53e6
        self.Izz = 3.01e6
        self.Ixz = 0.106e6
        self.Izy = 0.0
        self.Iyx = 0.0
        self.inertia = array([[self.Ixx, -self.Iyx, -self.Ixz],
                              [-self.Iyx, self.Iyy, -self.Izy],
                              [-self.Ixz, -self.Izy, self.Izz]])
        
        self.invInertia = inv(self.inertia)
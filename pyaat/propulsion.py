"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

propulsion model for turbo-propeller engine
"""

from numpy import array, cross
from pyaat.constants import RHO_SEA

class SimpleModel(object):
    def __init__(self):
        # control
        self.delta_p = 0.0
        
        # States
        self.TAS = 200.0
        self._rho = RHO_SEA
    
        # parameters
        self.Fmaxi = 70e3
        self.nrho = 0.775
        self.rhoi = 0.41271
        self.nv = 0
        self.Vi = 200.0
        self.position = array([0,0,1.42])
        self.attitude = array([0.,0.,0.])
        self.Pmax = self.Fmaxi*self.Vi
        
    @property
    def Forces(self):
        Fx = self.delta_p*self.Fmaxi*(self._rho/self.rhoi)**(self.nrho)
        return array([Fx, 0.0, 0.0])
    
    @property
    def Moments(self):
        return cross(self.Forces, self.position)
    
    @property
    def Power(self):
        return self.Fmaxi/self.Vi**self.nv*(self._rho/self.rhoi)**(self.nrho)
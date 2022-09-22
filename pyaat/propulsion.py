"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

propulsion model for turbo-propeller engine
"""

from abc import ABC, abstractmethod
from numpy import array, cross, radians, zeros
from constants import RHO_SEA
from tools import body2earth

class metaPropulsion(ABC):
    def __init__(self, name, position = array([0.0, 0.0, 0.0]), attitude = array([0.0, 0.0, 0.0])):
        self._delta_p = 0.0
        self._attitude = attitude
        self._position = position
        self._name = name
        self._rho = RHO_SEA
        self._TAS = 0.0

    def __len__(self):
        return 1000
    
    def set_position(self, entry):
        self._position = entry

    def set_attitude(self, entry):
        self._attitude = entry

    def set_delta_p(self, entry):
        self._delta_p = entry

    def set_airDensity(self, entry):
        self._rho = entry

    def set_airTemperature(self,entry):
        self._airTemperature = entry

    def set_TAS(self, entry):
        self._TAS = entry   

    def get_forces(self):
        return self._forces

    def get_moments(self):
        return self._moments
    
    def get_power(self):
        self._power

    @property
    @abstractmethod
    def _forces(self):
        pass

    @property
    @abstractmethod
    def _moments(self):
        pass

    @property
    @abstractmethod
    def _power(self):
        pass

class JetModel(metaPropulsion):
    def __init__(self, name, Fmaxi = 0.0 , nrho =0.0, rhoi = 0.0, nv = 0.0, Vi =0.0,
                 position = array([0,0,0]), attitude = array([0.,0.,0.])):
        super().__init__(name, position, attitude)

        # parameters
        self._Fmaxi = Fmaxi
        self._nrho = nrho
        self._rhoi = rhoi
        self._nv = nv
        self._Vi = Vi
        self._Pmaxi = self._Fmaxi*self._Vi
        
    def __len__(self):
        return 1000
        
    @property
    def _forces(self):
        Fx = self._delta_p*self._Fmaxi*(self._rho/self._rhoi)**(self._nrho)
        F_prop =  array([Fx, 0.0, 0.0])
        phiprop = self._attitude[0]
        thetaprop = self._attitude[1]
        psiprop = self._attitude[2]
        
        return body2earth(F_prop, psiprop, thetaprop, phiprop)
    
    @property
    def _moments(self):
        return cross(self._position, self._forces)
    
    @property
    def _power(self):
        return self._Fmaxi*self.Vi*(self._rho/self._rhoi)**(self._nrho)
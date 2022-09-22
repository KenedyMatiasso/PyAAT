"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

propulsion model for turbo-propeller engine
"""

from numpy import array, cross, radians, zeros
from constants import RHO_SEA
from tools import body2earth

class JetModel(object):
    def __init__(self, name, Fmaxi = 0.0 , nrho =0.0, rhoi = 0.0, nv = 0.0, Vi =0.0,
                 position = array([0,0,0]), attitude = array([0.,0.,0.])):
        # control
        self.delta_p = 0.0
        self.name = name
        # States
        self.TAS = 0.0
        self._rho = RHO_SEA
    
        # parameters
        self.Fmaxi = Fmaxi
        self.nrho = nrho
        self.rhoi = rhoi
        self.nv = nv
        self.Vi = Vi
        self.position = position
        self.attitude = attitude
        self.Pmaxi = self.Fmaxi*self.Vi
        
    def __len__(self):
        return 1000
        
    @property
    def Forces(self):
        Fx = self.delta_p*self.Fmaxi*(self._rho/self.rhoi)**(self.nrho)
        F_prop =  array([Fx, 0.0, 0.0])
        phiprop = self.attitude[0]
        thetaprop = self.attitude[1]
        psiprop = self.attitude[2]
        
        return body2earth(F_prop, psiprop, thetaprop, phiprop)
    
    @property
    def Moments(self):
        return cross(self.position, self.Forces)
    
    @property
    def Power(self):
        return self.Fmaxi*self.Vi*(self._rho/self.rhoi)**(self.nrho)
        
class Electric(object):
    def __init__(self, name, MaxSpeed =0.0, ct =0.0,  cd =0.0, sense =  'clockwise',
                 position = array([0,0,0]), attitude = array([0.,0.,0.])):
        # control
        self.delta_p = 0.0
        self.name = name
        # States
        self.TAS = 0.0
        self._rho = RHO_SEA
        self.sense = sense
    
        # parameters
        self.MaxSpeed = MaxSpeed
        self.ct = ct
        self.cd = cd
        self.position = position
        self.attitude = attitude
        
    def __len__(self):
        return 1000
        
    @property
    def Forces(self):
        Fx = (self.delta_p*self.MaxSpeed)**2*self.ct
        F_prop =  array([Fx, 0.0, 0.0])
        phiprop = self.attitude[0]
        thetaprop = self.attitude[1]
        psiprop = self.attitude[2]
        
        return body2earth(F_prop, psiprop, thetaprop, phiprop)
    
    @property
    def Moments(self):
        MD = (self.delta_p*self.MaxSpeed)**2*self.cd
        
        if self.sense == 'counter-clockwise':
            MD = -MD
        elif self.sense == 'clockwise':
            MD = MD
        else:
            MD = 0.0
            print('Incorrect sense')
        MD_vec = array([MD, 0.0, 0.0])
        phiprop = radians(self.attitude[0])
        thetaprop = radians(self.attitude[1])
        psiprop = radians(self.attitude[2])
        return cross(self.position, self.Forces) + body2earth(MD_vec, psiprop, thetaprop, phiprop)
    
    
class NoEngine(object):
    def __init__(self):
        # control
        self.delta_p = 0.0
        
        # States
        self.TAS = 0.0
        self._rho = RHO_SEA
    
        # parameters
        self.Fmaxi = 0.0
        self.nrho = 0.0
        self.rhoi = 0.0
        self.nv = 0.0
        self.Vi = 0.0
        self.position = array([0,0,0.0])
        self.attitude = array([0.,0.,0.])
        self.Pmax = self.Fmaxi*self.Vi
        
    @property
    def Forces(self):      
        return array([0.0, 0.0, 0.0])
    
    @property
    def Moments(self):
        return cross(self.position, self.Forces)
    
    @property
    def Power(self):
        return self.Fmaxi/self.Vi**self.nv*(self._rho/self.rhoi)**(self.nrho)
    
class Cluster(object):
    def __init__(self, propulsion_list):
        # control
        self.props = propulsion_list
        self.number = len(propulsion_list)
        self.dp = zeros(len(propulsion_list))
            
    @property
    def Forces(self):
        total_trhust = array([0.0, 0.0, 0.0])
        for i in range(len(self.props)):
            propeller = self.props[i]
            propeller.delta_p = self.dp[i]
            indForce = propeller.Forces
            total_trhust  = total_trhust + indForce
        return total_trhust
    
    @property
    def Moments(self):
        total_Moment = array([0.0, 0.0, 0.0])
        for propeller in self.props:
            indMoment = propeller.Moments
            total_Moment  = total_Moment + indMoment
        return total_Moment   
    
    @property
    def Power(self):
        return self.Fmaxi/self.Vi**self.nv*(self._rho/self.rhoi)**(self.nrho)
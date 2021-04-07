"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Aerodynamic model for a convencional aicraft
"""
from numpy import array, radians
from pyaat.constants import RHO_SEA
from numpy.linalg import inv

class Aircraft(object):
    def __init__(self):
        # Relevant states
        self.alpha= 0.0
        self.beta= 0.0
        self.TAS= 0.0
        self.alpha_d= 0.0
        self.beta_d= 0.0
        self.p= 0.0
        self.q= 0.0
        self.r= 0.0
        
        #Control action
        self.delta_e= 0.0
        self.delta_a= 0.0
        self.delta_r= 0.0
        
        # Aerodynamic coeficients
        self.CL0= 0.382
        self.CD0= 0.0252
        self.CY0= 0.0
        self.Cl0= 0.0
        self.Cm0= 0.0622
        self.Cn0 = 0.0
        self.CLa= 6.29
        self.CDa= 0.2010
        self.Cma= -3.63
        
        self.CYb= 0.785
        self.Clb= -0.121
        self.Cnb= 0.174
        
        self.CYp= -0.0794
        self.Clp= -0.522
        self.Cnp= -0.0587
        
        self.CLq= 14.6
        self.CDq= 0.281
        self.Cmq= -45.5
        
        self.CYr= 0.572
        self.Clr= 0.254
        self.Cnr= -0.277
        
        self.CLad= 4.04
        self.CDad= 0.0
        self.Cmad= -16.5
        
        self.CYbd= 0.0214
        self.Clbd= 0.0035
        self.Cnbd= 0.0098
    
        # Control coefficient
        self.CLde= 0.3891
        self.CDde= 0.0126
        self.Cmde= -1.5980
        
        self.CYda= 0.0094
        self.Clda= -0.1784
        self.Cnda= 0.0080
        
        self.CYdr= -0.3030
        self.Cldr= -0.0464
        self.Cndr= -0.1594
        
        # Geometric parameters
        self.Cbarw= 3.666
        self.bw= 28.42
        self.Sw= 95.0
        
        # Environment
        self._rho= RHO_SEA
        self.Vinf = 224.6
        
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
        
        self.CLmax = 2.4
        self.qmax = 23052.05 # considering 700km/h at sea level
        self.k = self.CDa/(2*self.CLa*self.CL0) #linear esimation TODO:imporve it
        
    @property
    def qdyn(self):
        return 0.5*self.TAS*self.TAS*self._rho
               
    @property
    def CL(self):
        return self.CL0 + self.CLa*self.alpha + self.Cbarw/(self.Vinf)*(self.CLq*self.q + self.alpha_d*self.CLad) + self.CLde*self.delta_e
    
    @property 
    def CD(self):
        return self.CD0 + self.CDa*self.alpha + self.Cbarw/(self.Vinf)*(self.CDq*self.q + self.CDad*self.alpha_d) + self.CDde*self.delta_e  
    
    @property 
    def CY(self):
        return self.CY0 + self.CYb*self.beta + self.bw/(self.Vinf)*(self.CYp*self.p + self.CYr*self.r + self.CYbd*self.beta_d) + self.CYda*self.delta_a + self.CYdr*self.delta_r
    
    @property
    def Cm(self):
        return self.Cm0 + self.Cma*self.alpha + self.Cbarw/(self.Vinf)*(self.Cmq*self.q + self.Cmad*self.alpha_d) + self.Cmde*self.delta_e
    
    @property
    def Cl(self):
        return self.Cl0 + self.Clb*self.beta + self.bw/(self.Vinf)*(self.Clp*self.p + self.Clr*self.r + self.Clbd*self.beta_d) + self.Clda*self.delta_a + self.Cldr*self.delta_r
    
    @property
    def Cn(self):
        return self.Cn0 + self.Cnb*self.beta + self.bw/(self.Vinf)*(self.Cnp*self.p + self.Cnr*self.r + self.Cnbd*self.beta_d) + self.Cnda*self.delta_a + self.Cndr*self.delta_r
    
    @property
    def Forces(self):
        return self.qdyn*self.Sw*array([-self.CD, -self.CY, -self.CL])
    
    @property
    def Moments(self):
        return self.qdyn*self.Sw*array([self.bw, self.Cbarw, self.bw])*array([self.Cl, self.Cm, self.Cn])

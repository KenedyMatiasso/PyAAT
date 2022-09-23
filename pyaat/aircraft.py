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
    def __init__(self, mass=1, Ixx=1, Iyy=1, Izz=1, Iyz=0, Izx=0,
                 Ixy=0, Cbarw =1, bw=1, Sw=1):

        # Relevant states
        self.number_aero_controls = 0
        self.controlsNames = []
        self.coefControls = []
        self.controlsLimits = []
        self.alpha= 0.0
        self.beta= 0.0
        self.TAS= 0.0
        self.alpha_d= 0.0
        self.beta_d= 0.0
        self.p= 0.0
        self.q= 0.0
        self.r= 0.0
        
        self.u = 0
        self.v = 0
        self.w = 0
        
        self.coefU = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefV = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefW = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        self.coef0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefALPHA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefBETA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefP = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefQ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.coefR = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        self.U = None
        self.ALPHAmax = radians(180)
                
        # Geometric parameters
        self.Cbarw= Cbarw
        self.bw= bw
        self.Sw= Sw
        
        # Environment
        self._rho= RHO_SEA
        self.Vinf = 224.6
        
        self._mass = mass
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.Ixz = Izx
        self.Izy = Iyz
        self.Iyx = Ixy
        

        self.CLmax = 2.4
        self.qmax = 23052.05 # considering 700km/h at sea level

    @property
    def k(self):
        return self.coefALPHA[0]/(2*self.coefALPHA[2]*self.coef0[2])
    
    @property
    def inertia(self):
        return array([[self.Ixx + 0.0001, -self.Iyx, -self.Ixz],
                      [-self.Iyx, self.Iyy - 0.0004, -self.Izy],
                      [-self.Ixz, -self.Izy, self.Izz + 0.0003]])
                    
    @property
    def invInertia(self):
        return inv(self.inertia)
        
    @property
    def qdyn(self):
        return 0.5*self.TAS*self.TAS*self._rho
    
    def set_control_surface(self, Name = 'Name of the surface', coeficients = None, limits = (radians(-20), radians(20))):
        self.number_aero_controls  += 1
        self.controlsNames.append(Name)
        self.coefControls.append(coeficients)
        self.controlsLimits.append(limits)
               
    @property
    def CL(self):
        CL = self.coef0[2] + self.coefALPHA[2]*self.alpha + self.coefBETA[2]*self.beta +\
        self.Cbarw/(self.Vinf)*(self.coefP[2]*self.p + self.coefQ[2]*self.q + self.coefR[2]*self.r) +\
        self.coefU[2]*self.u + self.coefV[2]*self.v + self.coefW[2]*self.w
        if CL>self.CLmax:
            CL = self.CLmax
        
        if self.number_aero_controls!=0:
            for i in range(self.number_aero_controls):
                CL += self.coefControls[i][2]*self.U[i]
        if self.alpha>= self.ALPHAmax:
            CL = 0
        return CL
    
    @property 
    def CD(self):      
        CD = self.coef0[0] + self.coefALPHA[0]*self.alpha + self.coefBETA[0]*self.beta +\
        self.Cbarw/(self.Vinf)*(self.coefP[0]*self.p + self.coefQ[0]*self.q + self.coefR[0]*self.r)  +\
        self.coefU[0]*self.u + self.coefV[0]*self.v + self.coefW[0]*self.w
        if self.number_aero_controls!=0:
            for i in range(self.number_aero_controls):
                CD += self.coefControls[i][0]*self.U[i]
        return  CD
    
    @property 
    def CY(self):       
        CY = self.coef0[1] + self.coefALPHA[1]*self.alpha + self.coefBETA[1]*self.beta +\
        self.Cbarw/(self.Vinf)*(self.coefP[1]*self.p + self.coefQ[1]*self.q + self.coefR[1]*self.r) +\
        self.coefU[1]*self.u + self.coefV[1]*self.v + self.coefW[1]*self.w
        if self.number_aero_controls!=0:
            for i in range(self.number_aero_controls):
                CY += self.coefControls[i][1]*self.U[i]
        return CY
    
    @property
    def Cm(self):
        Cm = self.coef0[4] + self.coefALPHA[4]*self.alpha + self.coefBETA[4]*self.beta +\
        self.Cbarw/(self.Vinf)*(self.coefP[4]*self.p + self.coefQ[4]*self.q + self.coefR[4]*self.r)+\
        self.coefU[4]*self.u + self.coefV[4]*self.v + self.coefW[4]*self.w
        if self.number_aero_controls!=0:
            for i in range(self.number_aero_controls):
                Cm += self.coefControls[i][4]*self.U[i]
        return Cm
    
    @property
    def Cl(self):
        Cl = self.coef0[3] + self.coefALPHA[3]*self.alpha + self.coefBETA[3]*self.beta +\
        self.bw/(self.Vinf)*(self.coefP[3]*self.p + self.coefQ[3]*self.q + self.coefR[3]*self.r) +\
        self.coefU[3]*self.u + self.coefV[3]*self.v + self.coefW[3]*self.w
        if self.number_aero_controls!=0:     
            for i in range(self.number_aero_controls):
                Cl += self.coefControls[i][3]*self.U[i]
        return Cl
    
    @property
    def Cn(self):
        Cn = self.coef0[5] + self.coefALPHA[5]*self.alpha + self.coefBETA[5]*self.beta +\
        self.bw/(self.Vinf)*(self.coefP[5]*self.p + self.coefQ[5]*self.q + self.coefR[5]*self.r) +\
        self.coefU[5]*self.u + self.coefV[5]*self.v + self.coefW[5]*self.w
        if self.number_aero_controls!=0:
            for i in range(self.number_aero_controls):
                Cn += self.coefControls[i][5]*self.U[i]
        return Cn
    
    @property
    def Forces(self):
        return self.qdyn*self.Sw*array([-self.CD, self.CY, -self.CL])
    
    @property
    def Moments(self):
        return self.qdyn*self.Sw*array([self.bw, self.Cbarw, self.bw])*array([self.Cl, self.Cm, self.Cn])

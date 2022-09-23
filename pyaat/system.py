"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This is the system file of PyAAT.
"""
#import sys 
#sys.setrecursionlimit(10**4)
import numpy as np

from pyaat.tools import computeTAS, earth2body, aero2body, body2earth, body2euler
from pyaat.tools import trimmerClimb, linearization, trimmerGeneral
from pyaat.tools import trimmerPullUp, trimmerCurve, trimmerCruize


from numpy import array, cross, arange, radians, tan, sqrt, sin, copy, transpose
from scipy.integrate import odeint
from pyaat.pyaatcontrol import equilibrium

class system(object):
    def __init__(self, aircraft = None, atmosphere = None, propulsion = None,
                 gravity = None, wind = None, control = None, extraDynamics = None):
        
        self.aircraft = aircraft
        self.atmosphere = atmosphere
        if type(propulsion)==list:
            self.propulsion = propulsion
        else:
            self.propulsion = [propulsion]
        self.gravity = gravity
        eq = equilibrium()
        self.control = [eq]
        self.Xe = []
        self.Ue = []
        self.U = []
        self.wind = wind
        self.windval = []
        self.time = None
        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        self.controlAction = []
        self.states = None
        self.feedback = None
        self.actuators = []
        self.sensors = []
        self.statesNames = ['xo', 'yo', 'zo', 'u', 'v', 'w', 'phi', 'theta', 'psi','p', 'q', 'r']
        self.XpControls = []
        
        if extraDynamics != None:
            for anyDynamic in extraDynamics:
                if anyDynamic.type == 'actuator':
                    self.actuators.append(anyDynamic)
                    self.statesNames.append(anyDynamic.name)
                    
            for anyDynamic in extraDynamics:
                if anyDynamic.type == 'sensor':
                    self.sensors.append(anyDynamic)
                    self.statesNames.append(anyDynamic.name)
                    
        if control != None:
            self.control.extend(control)
            for cont in control:
                if cont.number_of_states > 0:
                    for val in cont.names:
                        self.statesNames.append(val)
                        
    @property
    def controlsNames(self):
        controlsNames = []
        for i in self.aircraft.controlsNames:
            controlsNames.append(i)
        if len(self.propulsion)!=1000:
            for anyPropeller in self.propulsion: 
                controlsNames.append(anyPropeller._name + ' thrust')
                controlsNames.append(anyPropeller._name + ' psip')
                controlsNames.append(anyPropeller._name + ' thetap')
                
        else:
            controlsNames.append(self.propulsion._name + ' thrust')
            controlsNames.append(self.propulsion._name + ' psip')
            controlsNames.append(self.propulsion._name + ' thetap')
            
        return controlsNames
    
    @property
    def controlsLimits(self):
        controlsLimits = []
        for i in self.aircraft.controlsLimits:
            controlsLimits.append(i)
            
        if len(self.propulsion)!=1000:
            for anyPropeller in self.propulsion: 
                controlsLimits.append((0, 1))
                controlsLimits.append((-radians(180), radians(180)))
                controlsLimits.append((-radians(90), radians(90)))
                
        else:
            controlsLimits.append((0, 1))
            controlsLimits.append((-radians(180), radians(180)))
            controlsLimits.append((-radians(90), radians(90)))
            
        return controlsLimits
    
    def dynamics(self, t, X, U):
        #print('X', X)
        #print(U)
        
        # State space
        z = X[2]
        
        u = X[3]
        v = X[4]
        w = X[5]
        uvw = array([u, v, w])
        
        phi = X[6]
        theta = X[7]
        psi = X[8]
        
        p = X[9]
        q = X[10]
        r = X[11]
        pqr = array([p, q, r])
        
        # Aircraft
        m = self.aircraft._mass
        Inertia = self.aircraft.inertia
        InvInertia = self.aircraft.invInertia
        
        #gravity model
        self.gravity.set_altitude(-z)
        g = self.gravity.get_gravity()
        F_gravity_earth = m*g
        F_gravity_body = earth2body(F_gravity_earth, psi, theta, phi)
        
        #Atmospheric model
        self.atmosphere.set_altitude(-z)
        rho = self.atmosphere.get_airDensity()

        #Aerodynamic model
        self.aircraft._rho = rho
        
        #print('uvw', uvw)
        if t>0.1:
            if self.wind != None:
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw                
                uvw = uvw +  earth2body(winduvw, psi, theta, phi)
            #print('uvw', uvw)
        alpha, beta, TAS = computeTAS(uvw)
       
        self.aircraft.U = U
                
        self.aircraft.alpha = alpha
        self.aircraft.beta = beta
        self.aircraft.TAS = TAS
        
        self.aircraft.u = u
        self.aircraft.v = w
        self.aircraft.w = w
        
        self.aircraft.p = p
        self.aircraft.q = q
        self.aircraft.r = r
        
        F_aero_aero = self.aircraft.Forces
        M_aero_aero = self.aircraft.Moments
        
        F_aero_body = aero2body(F_aero_aero, alpha, beta)
        M_aero_body = aero2body(M_aero_aero, alpha, beta)
        
        #Propulsion model
        F_prop_body = np.zeros(3)
        M_prop_body = np.zeros(3)
        
        if len(self.propulsion)!=1000:
            for i in range(len(self.propulsion)):
                if U[self.aircraft.number_aero_controls + 3*i]>=1:
                    self.propulsion[i].set_delta_p(1.0)
                elif U[self.aircraft.number_aero_controls + 3*i]<=0:
                    self.propulsion[i].set_delta_p(0.0)
                else:
                    self.propulsion[i].set_delta_p(U[self.aircraft.number_aero_controls + 3*i])
                    
                self.propulsion[i].set_airDensity(rho)
                self.propulsion[i].set_attitude(array([U[self.aircraft.number_aero_controls + 3*i + 1], U[self.aircraft.number_aero_controls + 3*i +2], 0.]))
                F_prop_body += self.propulsion[i].get_forces()
                M_prop_body += self.propulsion[i].get_moments()
        else:
            self.propulsion.delta_p = U[self.aircraft.number_aero_controls + 1]
            self.propulsion._rho = rho
            self.propulsion.attitude = array([U[self.aircraft.number_aero_controls + 1], U[self.aircraft.number_aero_controls + 2], 0.])
            F_prop_body = np.array(self.propulsion.get_forces())
            M_prop_body = np.array(self.propulsion.get_moments())
        
        """
        Kinematics
        """
        
        # Linear Kinematics
        xyz_d = body2earth(uvw, psi, theta, phi)
        
        x_d = xyz_d[0]
        y_d = xyz_d[1]
        z_d = xyz_d[2]
        
        # Angular Kinematics
        phithetapsi_d = body2euler(pqr, theta, phi)
        phi_d = phithetapsi_d[0]
        theta_d = phithetapsi_d[1]
        psi_d = phithetapsi_d[2]
        
        """
        Dynamics
        """
        # Linear Dynamics
        uvw_d = 1/m*(F_prop_body + F_aero_body + F_gravity_body) - cross(pqr, uvw)
        u_d = uvw_d[0]
        v_d = uvw_d[1]
        w_d = uvw_d[2]
            
        # Angular Dynamics
        Moments = M_aero_body + M_prop_body
        
        #print('Moments', Moments)
        pqr_d = InvInertia.dot(Moments - cross(pqr, Inertia.dot(pqr)))
        p_d = pqr_d[0]
        q_d = pqr_d[1]
        r_d = pqr_d[2]
        
        #print('pqr_d', pqr_d)
        return array([
            x_d,
            y_d,
            z_d,
            u_d,
            v_d,
            w_d,
            phi_d,
            theta_d,
            psi_d,
            p_d,
            q_d,
            r_d
            ])
    
    def AugmentedDynamics(self, t, X, U):
        #Estados da aeronave
        XAircraft = X[0:12]
        Uapp = copy(U)
        
        #Computa ação de controle real
        if len(self.actuators) !=0:
            for i in range(len(self.actuators)):
                Uapp[self.controlsNames.index(self.actuators[i].command)] = X[12+i]
        
        #Derivada dos estados da aeronave
        Xp =  list(self.dynamics(t, XAircraft, Uapp))
        
        #Derivada dos estados dos atuadores
        if len(self.actuators) !=0:
            for i in range(len(self.actuators)):
                self.actuators[i].ddelta = U[self.controlsNames.index(self.actuators[i].command)]
                self.actuators[i].integral = X[12+i]
                Xp.append(self.actuators[i].derivative)
                
        if len(self.sensors) !=0:
            for i in range(len(self.sensors)):
                self.sensors[i].ddelta = X[self.statesNames.index(self.sensors[i].state)]
                self.sensors[i].integral = X[12 + len(self.actuators) + i]
                self.t = t
                Xp.append(self.sensors[i].derivative)
                
        if self.control != None:
            Xp.extend(self.XpControls)

        return Xp
    
    def simularaviao(self, X, t):
        Ue = copy(self.Ue)
        phi = X[6]
        theta = X[7]
        psi = X[8]   
        
        Xcont = copy(X)
        
        if t>0.1:
            if self.wind != None:
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw  
                    winduvw = earth2body(winduvw, psi, theta, phi)
                    
                Xcont[3] = Xcont[3] + winduvw[0]
                Xcont[5] = Xcont[4] + winduvw[1]
                Xcont[5] = Xcont[5] + winduvw[2]
        
        for cont in self.control:
            if cont.type == 'open-loop':
                cont.statesNames = self.statesNames
                cont.controlsNames = self.controlsNames
                cont.t = t
                cont.Xe = self.Xe
                cont.Ue = Ue
                cont.X = Xcont
                self.U = cont.U
                cont.controlsNames =  self.controlsNames
                Ue = copy(self.U)
        
        return self.dynamics(t, X, self.U)
    
    def simularDinamicaAumentada(self, X, t):
        self.XpControls = []
        Ue = copy(self.Ue)
        phi = X[6]
        theta = X[7]
        psi = X[8]
        Xcont = copy(X)
        if t>0.1:
            if self.wind != None:
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw
                    winduvw = earth2body(winduvw, psi, theta, phi)
                   
                Xcont[3] = Xcont[3] + winduvw[0]
                Xcont[5] = Xcont[4] + winduvw[1]
                Xcont[5] = Xcont[5] + winduvw[2]
                
        for cont in self.control:
            cont.statesNames = self.statesNames
            cont.controlsNames = self.controlsNames
            cont.t = t
            cont.Xe = self.Xe
            cont.Ue = Ue
            cont.X = Xcont
            self.U = cont.U
            Ue = copy(self.U)
            if cont.number_of_states > 0:
                self.XpControls.extend(cont.Xp_local)
                
        if len(self.propulsion)!=1000:
            for i in range(len(self.propulsion)):
                if self.U[self.aircraft.number_aero_controls + 3*i]>=1:
                    self.U[self.aircraft.number_aero_controls + 3*i] = 1.0
                elif self.U[self.aircraft.number_aero_controls + 3*i]<=0:
                    self.U[self.aircraft.number_aero_controls + 3*i] = 0.0
                else:
                    pass
        
        #print(self.XpControls)        
        return self.AugmentedDynamics(t, X, self.U)
    
    def ExternalSimulation(self, t, X, U_external):
        if t == 0.0:
            pass
            #self.simulate(X, U_external, T0 = 0.0, TF = 0.1)
            
        self.XpControls = []
        Ue = copy(U_external)
        phi = X[6]
        theta = X[7]
        psi = X[8]
        Xcont = copy(X)
        if t>0.1:
            if self.wind != None:
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw
                    winduvw = earth2body(winduvw, psi, theta, phi)
                   
                Xcont[3] = Xcont[3] + winduvw[0]
                Xcont[5] = Xcont[4] + winduvw[1]
                Xcont[5] = Xcont[5] + winduvw[2]
        
        for cont in self.control:
            cont.statesNames = self.statesNames
            cont.controlsNames = self.controlsNames
            cont.t = t
            cont.Xe = self.Xe
            cont.Ue = Ue
            cont.X = Xcont
            self.U = cont.U
            Ue = copy(self.U)
            if cont.number_of_states > 0:
                self.XpControls.extend(cont.Xp_local)
        
        #print(self.XpControls)        
        return array(self.AugmentedDynamics(t, X, self.U))
    
    def trimmer(self, condition = 'cruize', **kwargs):
        
        if condition == 'general':
            free = kwargs.get('free')
            desired = kwargs.get('desired')
            dependent = kwargs.get('dependent')
            fixed = kwargs.get('fixed')
            return trimmerGeneral(self.dynamics, free, desired, fixed, dependent , self.controlsNames, self.controlsLimits)
        
        zo = kwargs.get('zo')
        u = kwargs.get('u')
        zod = kwargs.get('zod')
        thetad = kwargs.get('thetad')
        psid = kwargs.get('psid')
        beta = kwargs.get('beta')

        # Check if HE is valid
        if zo == None:
            self._flagError = 'INVALID'
            self._error = 'You shall define an equilibrium z position zo'
            
        elif (type(zo) != int) and type(zo) != float:
            self._flagError = 'INVALID'
            self._error = 'The equilibrium altitude shall be an integer or a float'
                
        elif zo > 0:
            self._flagError = 'INVALID'
            self._error = 'The equilibrium altitude shall be positive'
        
        # Check if VE is valid
        if u == None:
            self._flagError = 'INVALID'
            self._error = 'You shall define an equilibrium velocity u'
            
        elif (type(u) != int) and type(u) != float:
            self._flagError = 'INVALID'
            self._error = 'The equilibrium velocity shall be an integer or a float'
                
        elif u < 0:
            self._flagError = 'INVALID'
            self._error = 'The equilibrium velocity shall be positive'
            
        # If cruize condition
        if condition == 'cruize':
            if self._flagError == 'VALID':
                return trimmerCruize(self, zo, u)
   
        # If climb condition
        elif condition == 'climb':
            # Check dH
            if zod == None:
                self._flagError = 'INVALID'
                self._error = 'You shall define a climb rate -zod'
                    
            elif (type(zod) != int) and type(zod) != float:
                self._flagError = 'INVALID'
                self._error = 'The equilibrium altitude shall be an integer or a float'
                
            if self._flagError == 'VALID':
                return trimmerClimb(self, zo, u, zod)
            else:
                raise Exception(self._error)
            
        elif condition == 'pullUp':
            if thetad == None:
                self._flagError = 'INVALID'
                self._error = 'You shall define a value for theta rate thetad'
                                    
            elif thetad < 0:
                self._flagError = 'INVALID'
                self._error = 'The theta rate dTH shall be positive'       
            
            if self._flagError == 'VALID':

                return trimmerPullUp(self, zo, u, thetad)
            else:
                raise Exception(self._error)
                
        elif condition == 'turn':
            if psid == None:
                self._flagError = 'INVALID'
                self._error = 'You shall define a value for turn rate psid'
                
            if self._flagError == 'VALID':
                return trimmerCurve(self, zo, u, psid, beta)
            else:
                raise Exception(self._error)
                
    def propagate(self, Xe, Ue, T0 = 0.0, TF = 10.0, dt = 0.01,
                  state = None, control = None):
        self.Xe = copy(Xe)
        tcritL = []
        self.Ue = copy(Ue)
          
        if state != None:
            for value in state:
                if value == 'xo':
                    Xe[0] = Xe[0] + state['xo']
                    
                if value == 'yo':
                    Xe[1] = Xe[1] + state['yo']
                    
                if value == 'zo':
                    Xe[2] = Xe[2] + state['zo']
                    
                if value == 'altitude':
                    Xe[2] = Xe[2] - state['altitude']
                
                if value == 'u':
                    Xe[3] = Xe[3] + state['u']
                    
                if value == 'v':
                    Xe[4] = Xe[4] + state['v']
                    
                if value == 'w':
                    Xe[5] = Xe[5] + state['w']
                    
                if value == 'phi':
                    Xe[6] = Xe[6] + state['phi']
                    
                if value == 'theta':
                    Xe[7] = Xe[7] + state['theta']
                    
                if value == 'psi':
                    Xe[8] = Xe[8] + state['psi']
                    
                if value == 'p':
                    Xe[9] = Xe[9] + state['p']
                    
                if value == 'q':
                    Xe[10] = Xe[10] + state['q']
                    
                if value == 'r':
                    Xe[11] = Xe[11] + state['r']               
                
                if value == 'alpha':
                    Xe[5] = Xe[5] + Xe[3]*tan(state['alpha'])
                    
                if value == 'beta':                 
                    V = sqrt(Xe[3]**2 + Xe[4]**2+ Xe[5]**2)
                    Xe[4] = Xe[4] + V*sin(state['beta'])

            
        if control != None:
            self.control.extend(control)
            for i in control:
                if type(i.t_init) == float:
                    tcritL.append(i.t_init)
                elif type(i.t_init) == int:
                    tcritL.append(float(i.t_init))
            
        self.time = arange(T0, TF, dt)
        
        solution = odeint(self.simularaviao, Xe, self.time, tcrit=tcritL)
        
        # Create control vector (plot purpuses only)      
        for i in range (0, len(self.time)):
            
            t = self.time[i]
            X = solution[i]          
            Ue = copy(self.Ue)
            
            for cont in self.control:
                cont.t = t
                cont.Xe = self.Xe
                cont.Ue = Ue
                cont.X = X
                U = cont.U
                Ue = copy(U)
            
            self.controlAction.append(U)
        self.controlAction = transpose(array(self.controlAction))
        
        self.states = transpose(solution)
        
        if self.wind != None:
            for i in range (0, len(self.time)):
                t = self.time[i]
                X = solution[i]        
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw    
                self.windval.append(list(winduvw))
            self.windval = array(self.windval)
        return solution, self.controlAction

    def simulate(self, Xe, Ue, T0 = 0.0, TF = 10.0, dt = 0.01,
                  perturbation = True, state = {'beta':0., 'alpha':0.},
                  control = None, extraDynamics  = None):
        tcritL = []
        self.Xe = copy(Xe)
        ## X =[x, y, z, u, v, w, phi, theta, psi, p, q, r]
        self.Ue = copy(Ue)
        
        if perturbation ==True:
            for value in state:
                if value == 'xo':
                    Xe[0] = self.Xe[0] + state['xo']
                    
                if value == 'yo':
                    Xe[1] = self.Xe[1] + state['yo']
                    
                if value == 'zo':
                    Xe[2] = self.Xe[2] + state['zo']
                    
                if value == 'altitude':
                    Xe[2] = self.Xe[2] - state['altitude']
                
                if value == 'u':
                    Xe[3] = self.Xe[3] + state['u']
                    
                if value == 'v':
                    Xe[4] = self.Xe[4] + state['v']
                    
                if value == 'w':
                    Xe[5] = self.Xe[5] + state['w']
                    
                if value == 'phi':
                    Xe[6] = self.Xe[6] + state['phi']
                    
                if value == 'theta':
                    Xe[7] = self.Xe[7] + state['theta']
                    
                if value == 'psi':
                    Xe[8] = self.Xe[8] + state['psi']
                    
                if value == 'p':
                    Xe[9] = self.Xe[9] + state['p']
                    
                if value == 'q':
                    Xe[10] = self.Xe[10] + state['q']
                    
                if value == 'r':
                    Xe[11] = self.Xe[11] + state['r']             
                
                if value == 'alpha':
                    Xe[5] = self.Xe[5] + self.Xe[3]*tan(state['alpha'])
                    
                if value == 'beta':
                    V = sqrt(self.Xe[3]**2 + self.Xe[4]**2+ self.Xe[5]**2)
                    Xe[4] = self.Xe[4] + V*sin(state['beta'])
            
            # Cria a flag para reduzir o passo de tempo
            if control != None:
                self.control.extend(control)
                for i in control:
                    if i.type == 'open-loop':
                        if type(i.t_init) == float:
                            tcritL.append(i.t_init)
                        elif type(i.t_init) == int:
                            tcritL.append(float(i.t_init))
     
            Xe = list(Xe)
            if extraDynamics != None:
                for anyDynamic in extraDynamics:
                    if anyDynamic.type == 'actuator':
                        self.actuators.append(anyDynamic)
                    elif anyDynamic.type == 'sensor':
                        self.sensors.append(anyDynamic)
                                    
            for actuator in self.actuators:
                Xe.append(self.Ue[self.controlsNames.index(actuator.command)])
     
            for sensor in self.sensors:
                Xe.append(self.Xe[self.statesNames.index(sensor.state)])

            condInicial = []
            for cont in self.control:
                if cont.number_of_states > 0:
                    for val in cont.names:
                        condInicial.append(0.0)
                        Xe.append(0.0)
        self.time = arange(T0, TF, dt)
        #print(np.degrees(Xe[7]))
        solution = odeint(self.simularDinamicaAumentada, Xe, self.time, tcrit=tcritL)
                
        for i in range (0, len(self.time)):
            
            t = self.time[i]
            X = solution[i]          
            Ue = copy(self.Ue)
            
            for cont in self.control:
                cont.statesNames = self.statesNames
                cont.controlsNames = self.controlsNames
                cont.t = t
                cont.Xe = self.Xe
                cont.Ue = Ue
                cont.X = X
                U = cont.U
                Ue = copy(U)
                
                if len(self.propulsion)!=1000:
                    for i in range(len(self.propulsion)):
                        if U[self.aircraft.number_aero_controls + 3*i]>=1:
                            U[self.aircraft.number_aero_controls + 3*i] = 1.0
                        elif U[self.aircraft.number_aero_controls + 3*i]<=0:
                            U[self.aircraft.number_aero_controls + 3*i] = 0.0
                        else:
                            pass
                
            self.controlAction.append(U)
        self.controlAction = transpose(array(self.controlAction))
        self.states = transpose(solution)
                
        if self.wind != None:
            for i in range (0, len(self.time)):
                t = self.time[i]
                X = solution[i]        
                winduvw = array([0, 0, 0])
                for eachwind in self.wind:
                    eachwind.X = X
                    eachwind.t = t
                    winduvw = winduvw + eachwind.uvw    
                self.windval.append(list(winduvw))
            self.windval = array(self.windval)
            
        return solution, self.controlAction
    
    def Linearize(self, Xe, Ue):
        return linearization(self.dynamics, Xe, Ue)
    
    def FullLinearize(self, Xe, Ue):
        from pyaat.tools import fullLinearization
        A, B = fullLinearization(self.ExternalSimulation, Xe, Ue)
        return A, B
    

"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This is the system file of PyAAT.
"""
from tools import computeTAS, earth2body, aero2body, body2earth, body2euler
from tools import trimmer, trimmerClimb, linearization
from tools import trimmerPullUp, trimmerCurve, lateroMatrix, modesMatrix
from tools import longMatrix, body2aero

from numpy import array, cross, arange, radians, tan, sqrt, sin, cos
from scipy.integrate import odeint


class system(object):
    def __init__(self, aircraft = None, atmosphere = None, propulsion = None, gravity = None):
        self.aircraft = aircraft
        self.atmosphere = atmosphere
        self.propulsion = propulsion
        self.gravity = gravity
        self.Xe = None
        self.Ue = None
        self.time = None
        
    def dynamics(self, t, X, U):
        # State space
        z = X[2]
        
        u = X[3]
        v = X[4]
        w = X[5]
        uvw =array([u, v, w])
        
        phi = X[6]
        theta = X[7]
        psi = X[8]
        
        p = X[9]
        q = X[10]
        r = X[11]
        pqr =array([p, q, r])
        
        #Controls
        delta_p = U[0]
        delta_e = U[1]
        delta_a = U[2]
        delta_r = U[3]
        
        # Aircraft
        m = self.aircraft._mass
        Inertia = self.aircraft.inertia
        InvInertia = self.aircraft.invInertia
        
        
        #gravity model
        self.gravity._altitude = -z
        g = self.gravity._gravity
        F_gravity_earth = m*g
        F_gravity_body = earth2body(F_gravity_earth, psi, theta, phi)
        #Atmospheric model
        self.atmosphere._altitude = -z
        rho = self.atmosphere._rho

        #Aerodynamic model
        self.aircraft._rho = rho
        alpha, beta, TAS = computeTAS(uvw)
       
        self.aircraft.delta_r = delta_r
        self.aircraft.delta_e = delta_e
        self.aircraft.delta_a = delta_a
        self.aircraft.alpha = alpha
        self.aircraft.beta = beta
        self.aircraft.TAS = TAS
        self.aircraft.p = p
        self.aircraft.q = q
        self.aircraft.r = r
        
        F_aero_aero = self.aircraft.Forces
        M_aero_aero = self.aircraft.Moments
        F_aero_body = aero2body(F_aero_aero, alpha, beta)
        M_aero_body = aero2body(M_aero_aero, alpha, beta)
        
        #Propulsion model
        self.propulsion.delta_p = delta_p
        self.propulsion._rho = rho
        F_prop_body = self.propulsion.Forces
        M_prop_body = self.propulsion.Moments
        
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
        pqr_d = InvInertia.dot(Moments + cross(pqr, Inertia.dot(pqr)))
        p_d = pqr_d[0]
        q_d = pqr_d[1]
        r_d = pqr_d[2]
        
        # Obtain alpha_d and beta_d to use on next iteration of the aeodynamic model
        # TODO: It slows down hugly the simulation. Solve it.
        """
        if t!=0:
            Valphabeta_d = array([1., 1./(TAS),1./(TAS*cos(beta))])*body2aero(uvw_d, alpha, beta)
            beta_d = Valphabeta_d[1]
            alpha_d = Valphabeta_d[2]
            # Send information to aerodynamic model
            self.aircraft.alpha_d = alpha_d
            self.aircraft.beta_d = beta_d
            
        """
        
        return array([x_d, y_d, z_d, u_d, v_d, w_d, phi_d, theta_d, psi_d, p_d, q_d, r_d])
    
    def simularaviao(self, X, t):
        Xp = self.dynamics(t, X, self.Ue)
        return Xp
    
    def trimmer(self, condition = 'cruize', HE = 5000.0, VE =150.0, dH = 0.0, dTH = 5., dPS = 2., BTA = 0.0):
        BTA = radians(BTA)
        dTH = radians(dTH)
        dPS = radians(dPS)
        if condition == 'cruize':
            return trimmer(self.dynamics, HE, VE)
        
        elif condition == 'climb':
            return trimmerClimb(self.dynamics, HE, VE, dH)
        
        elif condition == 'pullUp':
            return trimmerPullUp(self.dynamics, HE, VE, dTH)
        
        elif condition == 'curve':
            return trimmerCurve(self.dynamics, HE, VE, dPS, BTA)
    
    def propagate(self, Xe, Ue, T0 = 0.0, TF=10.0, dt = 0.01, perturbation = False, state = {'beta':0., 'alpha':0.}):
        if perturbation == False:
            self.Ue = Ue

        elif perturbation ==True:
            for value in state:
                if value == 'x':
                    Xe[0] = Xe[0] + state['x']
                    
                if value == 'y':
                    Xe[1] = Xe[1] + state['y']
                    
                if value == 'z':
                    Xe[2] = Xe[2] + state['z']
                    
                if value == 'altitude':
                    Xe[2] = Xe[2] - state['altitude']
                
                if value == 'u':
                    Xe[3] = Xe[3] + state['u']
                    
                if value == 'v':
                    Xe[4] = Xe[4] + state['v']
                    
                if value == 'w':
                    Xe[5] = Xe[5] + state['w']
                    
                if value == 'phi':
                    Xe[6] = Xe[6] + radians(state['phi'])
                    
                if value == 'theta':
                    Xe[7] = Xe[7] + radians(state['theta'])
                    
                if value == 'psi':
                    Xe[8] = Xe[8] + radians(state['psi'])
                    
                if value == 'p':
                    Xe[9] = Xe[9] + radians(state['p'])
                    
                if value == 'q':
                    Xe[10] = Xe[10] + radians(state['q'])
                    
                if value == 'r':
                    Xe[11] = Xe[11] + radians(state['r'])                  
                
                if value == 'alpha':
                    Xe[5] = Xe[5] + Xe[3]*tan(radians(state['alpha']))
                    
                if value == 'beta':
                    V = sqrt(Xe[3]**2 + Xe[4]**2+ Xe[5]**2)
                    Xe[4] = Xe[4] + V*sin(radians(state['beta']))
            self.Ue = Ue
            
        self.time = arange(T0, TF, dt)
        solution = odeint(self.simularaviao, Xe, self.time)
        
        # Create control vector (plot purpuses only)
        dp = []
        de = []
        da = []
        dr = []
        
        for t in self.time:
            dp.append(Ue[0])
            de.append(Ue[1])
            da.append(Ue[2])
            dr.append(Ue[3])
        
        control = array([dp, de, da, dr])    

        return solution, control
    
    def Linearize(self, Xe, Ue):
        return linearization(self.dynamics, Xe, Ue)
    
    def LinearModes(self, Xe, Ue):
        A, B = linearization(self.dynamics, Xe, Ue)
        return modesMatrix(A, B)
    
    def LinearLatero(self, Xe, Ue):
        A, B = linearization(self.dynamics, Xe, Ue)
        return lateroMatrix(A, B)
    
    def LinearLong(self, Xe, Ue):
        A, B = linearization(self.dynamics, Xe, Ue)
        return longMatrix(A, B)

    
    

"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This is the main file of PyAAT.
"""
from atmosphere import atmosISA
from aerodynamic import Aircraft
from propulsion import SimpleModel
from gravity import NewtonGravity
from airplane import airplane

from tools import computeTAS, earth2body, aero2body, body2earth, body2euler
from tools import trimmer, printInfo, trimmerClimb
from post_processing import plotter

from numpy import array, cross, arange, radians
from scipy.integrate import odeint

def dynamics(X,U):
    # State space
    x = X[0]
    y = X[1]
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
    
    #gravity model
    gravity._altitude = -z
    g = gravity._gravity
    F_gravity_earth = m*g
    F_gravity_body = earth2body(F_gravity_earth, psi, theta, phi)
    
    #Atmospheric model
    atmosphere._altitude = -z
    rho = atmosphere._rho
    
    #Aerodynamic model
    aerodynamic._rho = rho
    alpha, beta, TAS = computeTAS(uvw)
   
    aerodynamic.delta_r = delta_r
    aerodynamic.delta_e = delta_e
    aerodynamic.delta_a = delta_a
    aerodynamic.alpha = alpha
    aerodynamic.beta = beta
    aerodynamic.TAS = TAS
    aerodynamic.p = p
    aerodynamic.q = q
    aerodynamic.r = r
    
    F_aero_aero = aerodynamic.Forces
    M_aero_aero = aerodynamic.Moments
    F_aero_body = aero2body(F_aero_aero, alpha, beta)
    M_aero_body = aero2body(M_aero_aero, alpha, beta)
    
    #Propulsion model
    propulsion.delta_p = delta_p
    propulsion._rho = rho
    F_prop_body = propulsion.Forces
    M_prop_body = propulsion.Moments
    
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
    # TODO: It slows down strongly the simulation.
    
    '''
    Valphabeta_d = array([1., 1./(TAS),1./(TAS*cos(beta))])*body2aero(uvw_d, alpha, beta)
    beta_d = Valphabeta_d[1]
    alpha_d = Valphabeta_d[2]
    # Send information to aerodynamic model
    aerodynamic.alpha_d = alpha_d
    aerodynamic.beta_d = beta_d
    '''
    
    return array([x_d, y_d, z_d, u_d, v_d, w_d, phi_d, theta_d, psi_d, p_d, q_d, r_d])

def simularaviao(X,t):
    Xp = dynamics(X,Ue)
    return Xp

atmosphere = atmosISA()
propulsion = SimpleModel()
aerodynamic = Aircraft()
gravity = NewtonGravity()
airplane = airplane()

m = airplane._mass
Inertia = airplane.inertia
InvInertia = airplane.invInertia

HE = 10000
UE = 200
dH = 5

Xe, Ue = trimmer(dynamics,HE, UE)
#Xe, Ue = trimmerClimb(dynamics,HE, UE, dH)


sol = list(dynamics(Xe,Ue))
printInfo(Xe,Ue, frame ='aero')
printInfo(Xe,Ue, frame='controls')

T0=0
TF=180
dt=0.01
time = arange(T0, TF, dt)

Xe[5] = Xe[5] + 0.0
Xe[4] = Xe[4] + 0.0

Ue[1] = Ue[1] + radians(0)
Ue[2] = Ue[2] + radians(0)

solution = odeint(simularaviao, Xe, time)

# Create control vector (plot purpuses only)
dp = []
de = []
da = []
dr = []

for t in time:
    dp.append(Ue[0])
    de.append(Ue[1])
    da.append(Ue[2])
    dr.append(Ue[3])

control = array([dp, de, da, dr])    

pltr = plotter()
pltr.states = solution
pltr.time = time
pltr.control = control

pltr.LinVel(frame = 'body')
pltr.LinVel(frame = 'aero')
pltr.LinPos()
pltr.Attitude()
pltr.AngVel()
pltr.Controls()
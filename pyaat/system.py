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

from scipy.optimize import least_squares
from numpy import array, degrees, radians, cross, cos


def dynamic(X,U):
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
    print(F_prop_body)
    
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
    
    return array([x_d, y_d, z_d, u_d, v_d, w_d, phi_d, theta_d, psi_d, p_d, q_d, r_d])


def obj(Z):
    xe = 0.0
    ye = 0.0
    ze = -HE
    ue = UE
    ve = 0.0
    we = Z[0]   
    phie = 0.0
    psie = 0.0
    thetae = Z[1]
    pe = 0.0
    qe = 0.0
    re = 0.0
    
    delta_ae = 0.0
    delta_re = 0.0
    delta_pe = Z[2]
    delta_ee = Z[3]

    Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
    Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
    
    sol = dynamic(Xe, Ue)
    Zp = array([sol[2], sol[3], sol[5], degrees(sol[10])])
    return Zp

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

wg = 5
thetag = radians(3)
delta_pg = 0.3
delta_eg = radians(-5)
Zg = array([wg , thetag, delta_pg, delta_eg])

wlim = [-20,20]
thetalim = [radians(-20), radians(20)]
delta_p_lim = [0,1]
delta_e_lim = [radians(-30), radians(30)]

boundary1 = (wlim[0], thetalim[0], delta_p_lim[0],delta_e_lim[0])
boundary2 = (wlim[1], thetalim[1], delta_p_lim[1],delta_e_lim[1])


root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
we = root.x[0]
thetae = root.x[1]
delta_ee = root.x[3]
delta_pe = root.x[2]
ue = UE

'''
root = array([5.8911, radians(1.6879), radians(-1.1753), 0.3459])
alphae = radians(1.6879)
Ve = 200
we = root[0]
thetae = root[1]
delta_ee = root[2]
delta_pe = root[3]
ue = Ve*cos(alphae)
'''

xe = 0.0
ye = 0.0
ze = -HE
phie = 0.0
psie = 0.0
ve = 0.0
pe = 0.0
qe = 0.0
re = 0.0
delta_ae = 0.0
delta_re = 0.0


Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
Ue = array([delta_pe, delta_ee, delta_ae, delta_re])

sol = list(dynamic(Xe,Ue))

alpha, beta, TAS = computeTAS([ue,ve,we])

print('V')
print(TAS)
print('-------------')
print('alpha')
print(degrees(alpha))
print('-------------')
print('beta')
print(degrees(beta))
print('-------------')
print('dp')
print(Ue[0]*100)
print('-------------')
print('de')
print(degrees(Ue[1]))
print('-------------')
print('da')
print(degrees(Ue[2]))
print('-------------')
print('dr')
print(degrees(Ue[3]))

print('-------------------------------------------------------')
print('-------------------------------------------------------')
print(sol)

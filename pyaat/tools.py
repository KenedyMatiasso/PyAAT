"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of tools that can be used by PyAAT.
"""

from numpy import cos, sin, array, transpose, tan, around
from numpy import arctan2, arcsin, sqrt, degrees, radians, zeros, copy
from numpy import absolute, append, arange

import numpy as np
import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['axes.formatter.useoffset'] = False

from scipy.optimize import fsolve

def C1(theta):
    """
    Parameters
    ----------
    theta : float
        Angle 'theta' to be rotated around the X axis in rad .

    Returns
    -------
    array
        Cossine matix (3X3) of a 'theta' rotation around the X axis.

    """
    return array([[1, 0, 0],
                [0, cos(theta), sin(theta)],
                [0, -sin(theta), cos(theta)]])    
    
def C2(theta):
    """
    Parameters
    ----------
    theta : float
        Angle 'theta' to be rotated around the X axis in rad .

    Returns
    -------
    array
        Cossine matix (3X3) of a 'theta' rotation around the Y axis.
    """
    return array([[cos(theta), 0, -sin(theta)],
                [0, 1, 0],
                [sin(theta), 0, cos(theta)]])
    
def C3(theta):
    """
    Parameters
    ----------
    theta : float
        Angle 'theta' to be rotated around the X axis in rad.

    Returns
    -------
    array
        Cossine matix (3X3) of a 'theta' rotation around the Z axis.
    """
    return array([[cos(theta), sin(theta), 0],
                [-sin(theta), cos(theta), 0],
                [0,0,1]])

def aero2body(vector,alpha,beta):
    """
    Parameters
    ----------
    vector : array or list
        Vector to be converted form wind frame to body frame.
    alpha : float
        Attack angle.
    beta : float
        Sideslip angle.

    Returns
    -------
    array
        Vector rotated from aerodynamic frame to body fixed frame.

    """
    return transpose(array([[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta)],
                  [-cos(alpha)*sin(beta), cos(beta), -sin(alpha)*sin(beta)],
                  [-sin(alpha), 0, cos(alpha)]])).dot(vector)

def body2aero(vector,alpha,beta):
    """
    Parameters
    ----------
    vector : array or list
        Vector to be converted from body fixed frame to wind frame.
    alpha : float
        Attack angle.
    beta : float
        Sideslip angle.

    Returns
    -------
    array
        Vector rotated from body fixed to frame aerodynamic frame.

    """
    return array([[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta)],
                  [-cos(alpha)*sin(beta), cos(beta), -sin(alpha)*sin(beta)],
                  [-sin(alpha), 0, cos(alpha)]]).dot(vector)


def earth2body(vector, psi, theta, phi):
    """
    Parameters
    ----------
    vector : array or list
        Vector to be rotated from NED frame to body frame.
    psi : float
        Yaw angle in rad.
    theta : float
        Pitch angle in rad.
    phi : float
        Roll angle in rad.

    Returns
    -------
    TYPE
        Vector rotated from NED frame to body fixed frame.
        
    """
    return array([[cos(psi)*cos(theta), sin(psi)*cos(theta),-sin(theta)],
                  [sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi), sin(phi)*sin(psi)*sin(theta)+cos(theta)*cos(psi), sin(phi)*cos(theta)],
                  [cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi), cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi), cos(phi)*cos(theta)]]).dot(vector)

def body2earth(vector,psi,theta,phi):
    """
    Parameters
    ----------
    vector : array or list
        Vector to be rotated from body frame to NED frame.
    psi : float
        Yaw angle in rad.
    theta : float
        Pitch angle in rad.
    phi : float
        Roll angle in rad.

    Returns
    -------
    array
        Vector rotated from  body fixed frameto NED frame.
        
    """
    return transpose(array([[cos(psi)*cos(theta), sin(psi)*cos(theta),-sin(theta)],
                  [sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi), sin(phi)*sin(psi)*sin(theta)+cos(theta)*cos(psi), sin(phi)*cos(theta)],
                  [cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi), cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi), cos(phi)*cos(theta)]])).dot(vector)

def body2euler(vector,theta, phi):
    """
    Parameters
    ----------
    vector : array
        array containing p,q and r anular speed.
    theta : float
        pitch angle in rad.
    phi : float
        roll angle in rad.

    Returns
    -------
    array
        returns the derivativative of the Euler's angles phi_dot, theta_dot
        and psi_dot.
        
    Reference
    ---------
    B.L.Stevens, F.L.Lewis and E.N.Johnson, Aircraft control and Simulation:
        Dynamics, control design and autonomous systems, 3th edition, Willey,2016. 

    """
    return array([[1,sin(phi)*tan(theta), cos(phi)*tan(theta)],
                  [0,cos(phi), -sin(phi)],
                  [0, sin(phi)/cos(theta), cos(phi)/cos(theta)]]).dot(vector)


def computeTAS(uvw, wind=[0,0,0]):
    """
    Parameters
    ----------
    uvw : array or list
        array containg the body fixed velocity of the aircraft in m/s.
    wind : array or list
        array containing the inertial components of the wind.

    Returns
    -------
    alpha : float
        Attack angle in rad.
    beta : float
        Slideslip angle in rad.
    TAS : float
        True air speed in m/s.

    """
    u= uvw[0]
    v= uvw[1]
    w= uvw[2]
    
    TAS= sqrt(u**2+v**2+w**2)
    alpha = arctan2(w,u)
    if TAS > 0:
        beta = arcsin(v/TAS)
    else:
        beta = 0
    return alpha, beta, TAS


def trimmerGeneral(dynamic, free, desired, fixed, dependent, controlsNames, controlsLimits):
    from scipy.optimize import least_squares
    states = ['xo', 'yo', 'zo', 'u', 'v', 'w', 'phi', 'theta', 'psi', 'p', 'q', 'r']
    dstates = ['xod', 'yod', 'zod', 'ud', 'vd', 'wd', 'phid', 'thetad', 'psid','pd', 'qd', 'rd']
    
    zplen = len(desired)
    
    statesLimits = [(-10000, 10000),
                    (-10000, 10000),
                    (0, 100000),
                    (-1000,1000),
                    (-1000, 1000),
                    (-1000, 1000),
                    (-radians(180), radians(180)),
                    (-radians(180), radians(180)),
                    (-radians(180), radians(180)),
                    (-radians(180),radians(180)),
                    (-radians(180),radians(180)),
                    (-radians(180),radians(180))]    
    
    def obj(Z):
        Zp = []
        Xguess = zeros(12)
        Uguess = zeros(len(controlsNames))

        for variable in fixed:
            try:
                Xguess[states.index(variable)] = fixed[states[states.index(variable)]]
            except:
                pass
            
        for variable in fixed:
            try:
                Uguess[controlsNames.index(variable)] = fixed[controlsNames[controlsNames.index(variable)]]
            except:
                pass
        
        for variable in free:
            try:
                Xguess[states.index(variable)] = Z[Zlabels.index(variable)]
            except:
                pass
        
        for variable in free:
            try:
                if type(free[variable]) != str:
                    Uguess[controlsNames.index(variable)] = Z[Zlabels.index(variable)]
                else:
                    Zp.append(Z[Zlabels.index(variable)] - Z[Zlabels.index(free[variable])])
            except:
                pass
            
        if dependent !=None:
            for variable in dependent:
                Uguess[controlsNames.index(variable)] = Z[Zlabels.index(dependent[variable])]
            
        sol = dynamic(0, Xguess, Uguess)
        
        for variable in desired:
            try:
                if dstates.index(variable) == 'pd' or dstates.index(variable) == 'qd' or dstates.index(variable) == 'rd':
                    Zp.append(1000*(sol[dstates.index(variable)] - desired[variable]))
                elif dstates.index(variable) == 'phid' or dstates.index(variable) == 'thetad' or dstates.index(variable) == 'psid':
                    Zp.append(10*(sol[dstates.index(variable)] - desired[variable]))
                elif dstates.index(variable) == 'ud' or dstates.index(variable) == 'vd' or dstates.index(variable) == 'wd':
                    Zp.append(50*(sol[dstates.index(variable)] - desired[variable]))
                else:
                    Zp.append(1*(sol[dstates.index(variable)] - desired[variable]))
            except:
                pass
            
        if len(Zp) < len(Z):
            for i in range(len(desired) < len(free)):
                Zp.append(0.0)
                
        return Zp
    
    Zg = []
    boundary =[]
    Zlabels = []
    
    for variable in free:
        try:
            Zg.append(free[states[states.index(variable)]])
            boundary.append(statesLimits[states.index(variable)])
            Zlabels.append(states[states.index(variable)])
        except:
            pass
    for variable in free:
        try:
            if type(free[variable]) != str:
                Zg.append(free[controlsNames[controlsNames.index(variable)]])
                boundary.append(controlsLimits[controlsNames.index(variable)])
                Zlabels.append(controlsNames[controlsNames.index(variable)])
            else:
                #Zg.append(free[controlsNames[controlsNames.index(free(variable))]])
                Zg.append(free[controlsNames[controlsNames.index(free[variable])]])
                boundary.append(controlsLimits[controlsNames.index(variable)])
                Zlabels.append(controlsNames[controlsNames.index(variable)])
                zplen +=1
        except:
            pass
        
    if len(desired) > len(free):
        for i in range(len(desired) > len(free)):
            Zg.append(0.0)
            boundary.append((0, 1))
            Zlabels.append('token')
        
    boundary = array(boundary)
    boundary1 = boundary[:,0]
    boundary2 = boundary[:,1]
    
    root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
    
    Xe = zeros(12)
    Ue = zeros(len(controlsNames))
    
    for variable in fixed:
        try:
            Xe[states.index(variable)] = fixed[states[states.index(variable)]]
        except:
            pass
        
    for variable in fixed:
        try:
            Ue[controlsNames.index(variable)] = fixed[controlsNames[controlsNames.index(variable)]]
        except:
            pass
        
    for variable in free:
        try:
            Xe[states.index(variable)] = root.x[Zlabels.index(variable)]
        except:
            pass
    for variable in free:
        try:
            Ue[controlsNames.index(variable)] = root.x[Zlabels.index(variable)]
        except:
            pass
    if dependent !=None:
        for variable in dependent:
            Ue[controlsNames.index(variable)] = root.x[Zlabels.index(dependent[variable])]
    Ue = array(Ue)
    Xe = array(Xe)
    return Xe, Ue


def trimmerCruize(system, zo, u):   
    #print(system.aircraft.number_aero_controls)
    varFixas = {'xo': 0.0,
                'yo':0.0,
                'zo': zo,
                'psi': 0.0,
                'phi':0.0,
                'u':u}
    
    varLivres = {'w':0.0,
             'theta':0.0,
             'p': 0.0,
             'q':0.0,
             'r':0.0}
    
    varDep = {}
    
    for i in range(system.aircraft.number_aero_controls):
        varLivres[system.controlsNames[i]] = 0.0
        
    if len(system.controlsNames)>system.aircraft.number_aero_controls:
        varFixas[system.controlsNames[1 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[0]
        varFixas[system.controlsNames[2 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[1]
        varLivres[system.controlsNames[system.aircraft.number_aero_controls]] = 0.5
        
        for l in range(len(system.propulsion)-1):
            i = l+1
            varDep[system.controlsNames[3*i + system.aircraft.number_aero_controls]] = system.controlsNames[system.aircraft.number_aero_controls]
            varFixas[system.controlsNames[3*i + 1 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[0]
            varFixas[system.controlsNames[3*i + 2 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[1]
            
    varDesejadas = {'zod':0.0,
                'ud':0.0,
                'vd':0.0,
                'wd':0.0,
                'phid':0.0,
                'psid': 0.0,
                'thetad':0.0,
                'pd':0.0,
                'qd':0.0,
                'rd':0.0}
    
    Xe, Ue = trimmerGeneral(system.dynamics, varLivres, varDesejadas, varFixas, varDep, system.controlsNames, system.controlsLimits)

    return Xe, Ue


def trimmerClimb(system, zo, u, zod):
    varFixas = {'xo': 0.0,
                'yo':0.0,
                'zo': zo,
                'psi': 0.0,
                'phi':0.0,
                'u':u}
    
    varLivres = {'w':0.0,
             'theta':0.0,
             'p': 0.0,
             'q': 0.0,
             'r': 0.0,
             'v': 0.0}
    
    varDep = {}
    
    for i in range(system.aircraft.number_aero_controls):
        varLivres[system.controlsNames[i]] = 0.0
        
    if len(system.controlsNames)>system.aircraft.number_aero_controls:
        varFixas[system.controlsNames[1 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[0]
        varFixas[system.controlsNames[2 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[1]
        varLivres[system.controlsNames[system.aircraft.number_aero_controls]] = 0.5
        
        for l in range(len(system.propulsion)-1):
            i = l+1
            varDep[system.controlsNames[3*i + system.aircraft.number_aero_controls]] = system.controlsNames[system.aircraft.number_aero_controls]
            varFixas[system.controlsNames[3*i + 1 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[0]
            varFixas[system.controlsNames[3*i + 2 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[1]
            
    varDesejadas = {'zod':0.0 + zod,
                'ud':0.0,
                'vd':0.0,
                'wd':0.0,
                'phid':0.0,
                'psid': 0.0,
                'thetad':0.0,
                'pd':0.0,
                'qd':0.0,
                'rd':0.0}
    
    Xe, Ue = trimmerGeneral(system.dynamics, varLivres, varDesejadas, varFixas, varDep, system.controlsNames, system.controlsLimits)

    return Xe, Ue

def trimmerPullUp(system, zo, u, thetad):
    varFixas = {'xo': 0.0,
                'yo':0.0,
                'zo':zo,
                'psi': 0.0,
                'phi':0.0,
                'u':u}
    
    varLivres = {'w':0.0,
             'theta':0.0,
             'p': 0.0,
             'q': 0.0,
             'r': 0.0,
             'v': 0.0}
    
    varDep = {}
    
    for i in range(system.aircraft.number_aero_controls):
        varLivres[system.controlsNames[i]] = 0.0
        
    if len(system.controlsNames)>system.aircraft.number_aero_controls:
        varFixas[system.controlsNames[1 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[0]
        varFixas[system.controlsNames[2 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[1]
        varLivres[system.controlsNames[system.aircraft.number_aero_controls]] = 0.5
        
        for l in range(len(system.propulsion)-1):
            i = l+1
            varDep[system.controlsNames[3*i + system.aircraft.number_aero_controls]] = system.controlsNames[system.aircraft.number_aero_controls]
            varFixas[system.controlsNames[3*i + 1 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[0]
            varFixas[system.controlsNames[3*i + 2 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[1]
            
    varDesejadas = {'zod':0.0,
                'ud':0.0,
                'vd':0.0,
                'wd':0.0,
                'phid':0.0,
                'psid': 0.0,
                'thetad':0.0 +thetad,
                'pd':0.0,
                'qd':0.0,
                'rd':0.0}
    
    Xe, Ue = trimmerGeneral(system.dynamics, varLivres, varDesejadas, varFixas, varDep, system.controlsNames, system.controlsLimits)

    return Xe, Ue

def trimmerCurve(system, zo, u, psid, beta):
    varFixas = {'xo': 0.0,
                'yo':0.0,
                'zo':zo,
                'psi': 0.0,
                'u': u,
                'v': u*sin(beta)}
    
    varLivres = {'w':0.0,
             'theta':0.0,
             'p': 0.0,
             'q': 0.0,
             'r': 0.0,
             'phi':0.0}
    
    varDep = {}
    
    for i in range(system.aircraft.number_aero_controls):
        varLivres[system.controlsNames[i]] = 0.0
        
    if len(system.controlsNames)>system.aircraft.number_aero_controls:
        varFixas[system.controlsNames[1 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[0]
        varFixas[system.controlsNames[2 + system.aircraft.number_aero_controls]] = system.propulsion[0].attitude[1]
        varLivres[system.controlsNames[system.aircraft.number_aero_controls]] = 0.5
        
        for l in range(len(system.propulsion)-1):
            i = l+1
            varDep[system.controlsNames[3*i + system.aircraft.number_aero_controls]] = system.controlsNames[system.aircraft.number_aero_controls]
            varFixas[system.controlsNames[3*i + 1 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[0]
            varFixas[system.controlsNames[3*i + 2 + system.aircraft.number_aero_controls]] = system.propulsion[i].attitude[1]
            
    varDesejadas = {'zod':0.0,
                'ud':0.0,
                'vd':0.0,
                'wd':0.0,
                'phid':0.0,
                'psid': 0.0 + psid,
                'thetad':0.0,
                'pd':0.0,
                'qd':0.0,
                'rd':0.0}
    
    Xe, Ue = trimmerGeneral(system.dynamics, varLivres, varDesejadas, varFixas, varDep, system.controlsNames, system.controlsLimits)

    return Xe, Ue

def printInfo(X, U, frame = 'body'):
    """
    Print states and control

    Parameters
    ----------
    X : array
        States.
    U : array
        controls.
    frame : TYPE, optional
        DESCRIPTION. The default is 'body'.

    Returns
    -------
    None.

    """
    x = X[0]
    y = X[1]
    z = X[2]
    
    u = X[3]
    v = X[4]
    w = X[5]
    
    phi = X[6]
    theta = X[7]
    psi = X[8]
    
    p = X[9]
    q = X[10]
    r = X[11]
    
    #Controls
    delta_p = U[3]
    delta_e = U[0]
    delta_a = U[1]
    delta_r = U[2]
    
    alpha, beta, TAS = computeTAS([u,v,w])
    if frame== 'aero':
        print('--------------------------------')
        print('------------ STATES ------------')
        print('------------- AERO -------------')
        print('--------------------------------')

        print('V')
        print(TAS)
        print('-------------')
        print('alpha')
        print(degrees(alpha))
        print('-------------')
        print('beta')
        print(degrees(beta))
        print('-------------')
        print('phi')
        print(degrees(phi))
        print('-------------')
        print('theta')
        print(degrees(theta))
        print('-------------')
        print('psi')
        print(degrees(psi))
        print('-------------')
        print('p')
        print(degrees(p))
        print('-------------')
        print('q')
        print(degrees(q))
        print('-------------')
        print('r')
        print(degrees(r))
        print('-------------')
        print('x0')
        print(x)
        print('-------------')
        print('y0')
        print(y)
        print('-------------')
        print('H')
        print(-z)
        
    elif frame =='controls':
        print('--------------------------------')
        print('----------- CONTROLS -----------')
        print('--------------------------------')
        print('delta_p')
        print(delta_p*100)
        print('-------------')
        print('delta_e')
        print(degrees(delta_e))
        print('-------------')
        print('delta_a')
        print(degrees(delta_a))
        print('-------------')
        print('delta_r')
        print(degrees(delta_r))
    else:
        print('--------------------------------')
        print('------------ STATES ------------')
        print('------------- BODY -------------')
        print('--------------------------------')
        print('x')
        print(x)
        print('-------------')
        print('y')
        print(y)
        print('-------------')
        print('z')
        print(z)
        print('-------------')
        print('u')
        print(u)
        print('-------------')
        print('v')
        print(v)
        print('-------------')
        print('w')
        print(w)
        print('-------------')
        print('phi')
        print(degrees(phi))
        print('-------------')
        print('theta')
        print(degrees(theta))
        print('-------------')
        print('psi')
        print(degrees(psi))
        print('-------------')
        print('p')
        print(degrees(p))
        print('-------------')
        print('q')
        print(degrees(q))
        print('-------------')
        print('r')
        print(degrees(r))

def linearization(dynamics, Xe, Ue):
    eps = 0.02
    A = zeros((len(Xe), len(Xe)))
    B = zeros((len(Xe), len(Ue)))
    
    for i in range(0,len(Xe)):
        Xm = copy(Xe)
        Xm[i] += eps
        sol1 = dynamics(0, Xm, Ue)
        Xm = copy(Xe)
        Xm[i] -= eps
        sol2 = dynamics(0, Xm,Ue)
        A[:,i] = (sol1 -sol2)/(2*eps)
        
    for i in range(0,len(Ue)):
        Um = copy(Ue)
        Um[i] += eps
        sol1= dynamics(0, Xe, Um)
        
        Um = copy(Ue)
        Um[i] -= eps
        sol2 = dynamics(0, Xe, Um)
        
        B[:,i] = (sol1 -sol2)/(2*eps)
    return A, B

def fullLinearization(dynamics, Xe, Ue):
    eps = 0.02
    A = zeros((len(Xe), len(Xe)))
    B = zeros((len(Xe), len(Ue)))
    
    for i in range(0,len(Xe)):
        Xm = copy(Xe)
        Xm[i] += eps
        sol1 = dynamics(0, Xm, Ue)
        Xm = copy(Xe)
        Xm[i] -= eps
        sol2 = dynamics(0, Xm,Ue)
        A[:,i] = (sol1 -sol2)/(2*eps)
        
    for i in range(0,len(Ue)):
        Um = copy(Ue)
        Um[i] += eps
        sol1= dynamics(0, Xe, Um)
        
        Um = copy(Ue)
        Um[i] -= eps
        sol2 = dynamics(0, Xe, Um)
        
        B[:,i] = (sol1 -sol2)/(2*eps)
    return A, B

def quaternion2DCM(quaternion):
    qs = quaternion[0]
    qx = quaternion[1]
    qy = quaternion[2]
    qz = quaternion[3]
    
    return array([[qs**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qz*qs), 2*(qx*qz + qy*qs)],
               [2*(qx*qy + qz*qs), qs**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qx*qs)],
               [2*(qx*qz - qy*qs), 2*(qy*qz + qx*qs), qs**2 - qx**2 - qy**2 + qz**2]])

def DCM2quaternion(DCM):
    qs = absolute(1/2*sqrt(DCM[0,0] + DCM[1,1] + DCM[2,2] + 1))
    qx = (DCM[1,2] - DCM[2,1])/(4*qs)
    qy = (DCM[2,0] - DCM[0,2])/(4*qs)
    qz = (DCM[0,1] - DCM[1,0])/(4*qs)
    return array([qs, qx, qy,qz])

def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)



def rpm2radps(vel):
    return vel*0.10472

class RootFinder:
    def __init__(self, start, stop, step=0.01, root_dtype="float64", xtol=1e-9):

        self.start = start
        self.stop = stop
        self.step = step
        self.xtol = xtol
        self.roots = array([], dtype=root_dtype)

    def add_to_roots(self, x):

        if (x < self.start) or (x > self.stop):
            return  # outside range
        if any(abs(self.roots - x) < self.xtol):
            return  # root already found.

        self.roots = append(self.roots, x)

    def find(self, f, *args):
        current = self.start

        for x0 in arange(self.start, self.stop + self.step, self.step):
            if x0 < current:
                continue
            x = self.find_root(f, x0, *args)
            if x is None:  # no root found.
                continue
            current = x
            self.add_to_roots(x)

        return self.roots

    def find_root(self, f, x0, *args):

        x, _, ier, _ = fsolve(f, x0=x0, args=args, full_output=True, xtol=self.xtol)
        if ier == 1:
            return x[0]
        return None

class plotter(object):
    def __init__(self):
        # State space
        self.states = array([0,0,0,0,0,0,0,0,0,0,0,0])
        self.time = None
        self.control = array([0,0,0,0])

    @property
    def x(self):
        return self.states[:,0]
    
    @property
    def y(self):
        return self.states[:,1]
    
    @property
    def z(self):
        return self.states[:,2]
    
    @property
    def H(self):
        return -self.states[:,2]   
    
    @property
    def u(self):
        return self.states[:,3]
    
    @property
    def v(self):
        return self.states[:,4]
    
    @property
    def w(self):
        return self.states[:,5]
    
    @property
    def phi(self):
        return self.states[:,6]
    
    @property
    def theta(self):
        return self.states[:,7]
    
    @property
    def psi(self):
        return self.states[:,8] 
    
    @property
    def p(self):
        return self.states[:,9]
    
    @property
    def q(self):
        return self.states[:,10]
    
    @property
    def r(self):
        return self.states[:,11]
    
    @property
    def uvw(self):
        return array([self.u, self.v, self.w])
    
    @property
    def TAS(self):
        alpha, beta, TAS = computeTAS(self.uvw)
        return TAS
    
    @property
    def alpha(self):
        alpha, beta, TAS = computeTAS(self.uvw)
        return alpha
    
    @property
    def beta(self):
        alpha, beta, TAS = computeTAS(self.uvw)
        return beta
    
    @property
    def delta_p(self):
        return self.control[0,:]
    
    @property
    def delta_e(self):
        return self.control[1,:]

    @property
    def delta_a(self):
        return self.control[2,:]
    
    @property
    def delta_r(self):
        return self.control[3,:]

    def LinVel(self, frame = 'body'):
        if frame == 'aero':
            plt.figure()
            plt.subplot(311)
            plt.title("Linear Velocity")
            plt.plot(self.time, around(self.TAS, decimals=4), color = 'red')
            plt.ylabel('V [m/s]')
            plt.grid()
            plt.subplot(312)
            plt.plot(self.time, around(degrees(self.alpha), decimals=4), color = 'red')
            plt.ylabel(r'$\alpha$ [deg]')
            plt.grid()
            plt.subplot(313)
            plt.plot(self.time, around(degrees(self.beta), decimals=4), color = 'red')
            plt.ylabel(r'$\beta$ [deg]')
            plt.grid()
            plt.show()
        else:
            plt.figure()
            plt.subplot(311)
            plt.title("Velocity")
            plt.plot(self.time, around(self.u, decimals=4), color = 'red')
            plt.ylabel('u [m/s]')
            plt.grid()
            plt.subplot(312)
            plt.plot(self.time, around(self.v, decimals=4), color = 'red')
            plt.ylabel('v [m/s]')
            plt.grid()
            plt.subplot(313)
            plt.plot(self.time, around(self.w, decimals=4), color = 'red')
            plt.ylabel('w [m/s]')
            plt.grid()
            plt.show()
            
    def LinPos(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Linear position")
        plt.plot(self.time, around(self.x, decimals=4), color = 'red')
        plt.ylabel('x [m]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, around(self.y, decimals=4), color = 'red')
        plt.ylabel('y [m]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time,around(-self.z, decimals=4), color = 'red')
        plt.ylabel('H [m]')
        plt.grid()
        plt.show()
        
    def LinPos3D(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(around(self.x, decimals=4), around(self.y, decimals=4), around(-self.z, decimals=4), 'red')
        
    def Attitude(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Attitude")
        plt.plot(self.time, around(degrees(self.phi), decimals=4), color = 'red')
        plt.ylabel('$\Phi$ [deg]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, around(degrees(self.theta), decimals=4), color = 'red')
        plt.ylabel('$\Theta$ [deg]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time, around(degrees(self.psi), decimals=4), color = 'red')
        plt.ylabel('$\Psi$ [deg]')
        plt.grid()
        plt.show()
        
    def linPos2D(self):
        plt.figure()
        plt.subplot(211)
        plt.title("Attitude")
        plt.plot(around(self.x, decimals=4), around(self.y, decimals=4), color = 'red')
        plt.ylabel('Y position]')
        plt.grid()
        plt.subplot(212)
        plt.plot(around(self.x, decimals=4), around(-self.z, decimals=4), color = 'red')
        plt.ylabel('X position')
        plt.grid()
        plt.show()
        
    def AngVel(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Angular Speed")
        plt.plot(self.time, around(degrees(self.p), decimals=4), color = 'red')
        plt.ylabel('$p$ [deg/s]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, around(degrees(self.q), decimals=4), color = 'red')
        plt.ylabel('$q$ [deg/s]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time, around(degrees(self.r), decimals=4), color = 'red')
        plt.ylabel('$r$ [deg/s]')
        plt.grid()
        plt.show()
        
    def Controls(self):
        plt.figure()
        plt.subplot(411)
        plt.title("Controls")
        plt.plot(self.time, around(self.delta_p, decimals=4), color = 'red')
        plt.ylabel('$\delta_p$ [%]')
        plt.grid()
        plt.subplot(412)
        plt.plot(self.time, around(degrees(self.delta_e), decimals=4), color = 'red')
        plt.ylabel('$\delta_e$ [deg]')
        plt.grid()
        plt.subplot(413)
        plt.plot(self.time, around(degrees(self.delta_a), decimals=4), color = 'red')
        plt.ylabel('$\delta_a$ [deg]')
        plt.grid()
        plt.subplot(414)
        plt.plot(self.time, around(degrees(self.delta_r), decimals=4), color = 'red')
        plt.ylabel('$\delta_r$ [deg]')
        plt.grid()
        plt.show()
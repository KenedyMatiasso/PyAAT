"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of tools that can be used by PyAAT.
"""

from numpy import cos, sin, array, transpose, tan
from numpy import arctan2, arcsin, sqrt, degrees, radians

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


def earth2body(vector,psi,theta,phi):
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
    beta = arcsin(v/TAS)
    
    return alpha, beta, TAS
    
def trimmer(dynamic, HE, UE):
    """
    Trimmer for horizontal steady flight (cruize)

    Parameters
    ----------
    dynamic : function
        Function f(X,U) containing the dynamic and kinematics of the body.
    HE : float
        Altitude.
    UE : float
        aerodynamic speed.

    Returns
    -------
    tuple
        tuple containing the states and controls at trimmed condition.

    """
    from scipy.optimize import least_squares
    
    def obj(Z):
        xe = 0.0
        ye = 0.0
        ze = -HE
        ve = 0.0
        we = Z[0]   
        phie = 0.0
        psie = 0.0
        thetae = Z[1]
        pe = 0.0
        qe = 0.0
        re = 0.0
        
        ue = sqrt(UE**2 - we**2 -ve**2 )
    
        delta_ae = 0.0
        delta_re = 0.0
        delta_pe = Z[2]
        delta_ee = Z[3]
    
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(Xe, Ue)
        Zp = array([sol[2], sol[3], sol[5], degrees(sol[10])])
        return Zp
    
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
    
    ue = sqrt(UE**2 - we**2 -ve**2 )
    
    Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
    Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
    return (Xe, Ue)

def trimmerClimb(dynamic, HE, UE, dH):
    """
    Trimmer for constant climb

    Parameters
    ----------
    dynamic : function
        Function f(X,U) containing the dynamic and kinematics of the body.
    HE : float
        Altitude.
    UE : float
        aerodynamic speed.

    Returns
    -------
    tuple
        tuple containing the states and controls at trimmed condition.

    """
    from scipy.optimize import least_squares
    
    def obj(Z):
        xe = 0.0
        ye = 0.0
        ze = -HE
        ve = 0.0
        we = Z[0]   
        phie = 0.0
        psie = 0.0
        thetae = Z[1]
        pe = 0.0
        qe = 0.0
        re = 0.0
        
        ue = sqrt(UE**2 - we**2 -ve**2 )
    
        delta_ae = 0.0
        delta_re = 0.0
        delta_pe = Z[2]
        delta_ee = Z[3]
    
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(Xe, Ue)
        Zp = array([sol[2]+dH, sol[3], sol[5], degrees(sol[10])])
        return Zp
    
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
    
    ue = sqrt(UE**2 - we**2 -ve**2 )
    
    Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
    Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
    return (Xe, Ue)

    
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
    delta_p = U[0]
    delta_e = U[1]
    delta_a = U[2]
    delta_r = U[3]
    
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

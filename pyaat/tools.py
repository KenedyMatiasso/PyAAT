"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of tools that can be used by PyAAT.
"""

from numpy import cos, sin, array, transpose, tan
from numpy import arctan2, arcsin,sqrt

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
    u =uvw[0]
    v=uvw[1]
    w=uvw[2]
    
    TAS= sqrt(u**2+v**2+w**2)
    alpha = arctan2(w,u)
    beta = arcsin(u/TAS)
    
    return alpha,beta,TAS
    
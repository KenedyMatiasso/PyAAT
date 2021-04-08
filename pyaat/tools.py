"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of tools that can be used by PyAAT.
"""

from numpy import cos, sin, array, transpose, tan, around
from numpy import arctan2, arcsin, sqrt, degrees, radians, zeros, copy
from numpy import absolute, append, arange

import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['axes.formatter.useoffset'] = False

from scipy.optimize import fsolve


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
        we = Z[0]  
        ve = Z[1]
        ue = sqrt(UE**2 - we**2 -ve**2 )
        phie = 0.0
        psie = 0.0
        thetae = Z[2]
        pe = Z[3]
        qe = Z[4]
        re = Z[5]
    
        delta_pe = Z[6]
        delta_ee = Z[7]
        delta_ae = Z[8]
        delta_re = Z[9]
        
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(0, Xe, Ue)
        Zp = array([sol[2], sol[3], sol[4], sol[5], degrees(sol[6]),degrees(sol[7]), degrees(sol[8]), degrees(sol[9]), degrees(sol[10]), degrees(sol[11])])
        return Zp
    
    wg = 0.0
    vg = 0.0
    thetag = radians(0)
    pg = 0.0
    qg = 0.0
    rg = 0.0
    delta_pg = 0.5
    delta_eg = radians(0.)
    delta_ag = 0.0
    delta_rg = 0.0
    Zg = array([wg, vg, thetag, pg, qg, rg, delta_pg, delta_eg, delta_ag, delta_rg])
    
    wlim = [-20, 20]
    vlim = [-50, 50]
    thetalim = [radians(-20), radians(20)]
    plim = [radians(-20), radians(20)]
    qlim = [radians(-20), radians(20)]
    rlim = [radians(-20), radians(20)]
    delta_p_lim = [0,1]
    delta_e_lim = [radians(-40), radians(40)]
    delta_a_lim = [radians(-40), radians(40)]
    delta_r_lim = [radians(-40), radians(40)]
    
    boundary1 = (wlim[0], vlim[0], thetalim[0], plim[0], qlim[0], rlim[0], delta_p_lim[0],delta_e_lim[0], delta_a_lim[0], delta_r_lim[0])
    boundary2 = (wlim[1], vlim[1], thetalim[1], plim[1], qlim[1], rlim[1], delta_p_lim[1],delta_e_lim[1], delta_a_lim[1], delta_r_lim[1])

    root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
    we = root.x[0]
    ve = root.x[1]
    thetae = root.x[2]
    pe = root.x[3]
    qe = root.x[4]
    re = root.x[5]
    delta_pe = root.x[6]
    delta_ee = root.x[7]
    delta_ae = root.x[8]
    delta_re = root.x[9]
    xe = 0.0
    ye = 0.0
    ze = -HE
    phie = 0.0
    psie = 0.0
    
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
        we = Z[0]  
        ve = Z[1]
        ue = sqrt(UE**2 - we**2 -ve**2 )
        phie = 0.0
        psie = 0.0
        thetae = Z[2]
        pe = Z[3]
        qe = Z[4]
        re = Z[5]
    
        delta_pe = Z[6]
        delta_ee = Z[7]
        delta_ae = Z[8]
        delta_re = Z[9]
        
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(0, Xe, Ue)
        Zp = array([sol[2]+dH, sol[3], sol[4], sol[5], degrees(sol[6]),degrees(sol[7]), degrees(sol[8]), degrees(sol[9]), degrees(sol[10]), degrees(sol[11])])
        return Zp
    
    wg = 0.0
    vg = 0.0
    thetag = radians(0)
    pg = 0.0
    qg = 0.0
    rg = 0.0
    delta_pg = 0.5
    delta_eg = radians(0.)
    delta_ag = 0.0
    delta_rg = 0.0
    Zg = array([wg, vg, thetag, pg, qg, rg, delta_pg, delta_eg, delta_ag, delta_rg])
    
    wlim = [-20, 20]
    vlim = [-50, 50]
    thetalim = [radians(-20), radians(20)]
    plim = [radians(-20), radians(20)]
    qlim = [radians(-20), radians(20)]
    rlim = [radians(-20), radians(20)]
    delta_p_lim = [0,1]
    delta_e_lim = [radians(-40), radians(40)]
    delta_a_lim = [radians(-40), radians(40)]
    delta_r_lim = [radians(-40), radians(40)]
    
    boundary1 = (wlim[0], vlim[0], thetalim[0], plim[0], qlim[0], rlim[0], delta_p_lim[0],delta_e_lim[0], delta_a_lim[0], delta_r_lim[0])
    boundary2 = (wlim[1], vlim[1], thetalim[1], plim[1], qlim[1], rlim[1], delta_p_lim[1],delta_e_lim[1], delta_a_lim[1], delta_r_lim[1])

    root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
    we = root.x[0]
    ve = root.x[1]
    thetae = root.x[2]
    pe = root.x[3]
    qe = root.x[4]
    re = root.x[5]
    delta_pe = root.x[6]
    delta_ee = root.x[7]
    delta_ae = root.x[8]
    delta_re = root.x[9]
    xe = 0.0
    ye = 0.0
    ze = -HE
    phie = 0.0
    psie = 0.0
    
    ue = sqrt(UE**2 - we**2 -ve**2 )
    
    Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
    Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
    return (Xe, Ue)

def trimmerPullUp(dynamic, HE, UE, dTH):
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
        we = Z[0]  
        ve = Z[1]
        ue = sqrt(UE**2 - we**2 -ve**2 )
        phie = 0.0
        psie = 0.0
        thetae = Z[2]
        pe = Z[3]
        qe = Z[4]
        re = Z[5]
    
        delta_pe = Z[6]
        delta_ee = Z[7]
        delta_ae = Z[8]
        delta_re = Z[9]
        
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(0, Xe, Ue)
        Zp = array([sol[2], sol[3], sol[4], sol[5], degrees(sol[6]),10*degrees(sol[7]-dTH), degrees(sol[8]), degrees(sol[9]), degrees(sol[10]), degrees(sol[11])])
        return Zp
    
    wg = 0.0
    vg = 0.0
    thetag = radians(0)
    pg = 0.0
    qg = dTH
    rg = 0.0
    delta_pg = 0.5
    delta_eg = radians(0.)
    delta_ag = 0.0
    delta_rg = 0.0
    Zg = array([wg, vg, thetag, pg, qg, rg, delta_pg, delta_eg, delta_ag, delta_rg])
    
    wlim = [-20, 20]
    vlim = [-50, 50]
    thetalim = [radians(-40), radians(40)]
    plim = [radians(-40), radians(40)]
    qlim = [radians(-40), radians(40)]
    rlim = [radians(-40), radians(40)]
    delta_p_lim = [0,1]
    delta_e_lim = [radians(-40), radians(40)]
    delta_a_lim = [radians(-40), radians(40)]
    delta_r_lim = [radians(-40), radians(40)]
    
    boundary1 = (wlim[0], vlim[0], thetalim[0], plim[0], qlim[0], rlim[0], delta_p_lim[0],delta_e_lim[0], delta_a_lim[0], delta_r_lim[0])
    boundary2 = (wlim[1], vlim[1], thetalim[1], plim[1], qlim[1], rlim[1], delta_p_lim[1],delta_e_lim[1], delta_a_lim[1], delta_r_lim[1])

    root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
    we = root.x[0]
    ve = root.x[1]
    thetae = root.x[2]
    pe = root.x[3]
    qe = root.x[4]
    re = root.x[5]
    delta_pe = root.x[6]
    delta_ee = root.x[7]
    delta_ae = root.x[8]
    delta_re = root.x[9]
    xe = 0.0
    ye = 0.0
    ze = -HE
    phie = 0.0
    psie = 0.0
    
    ue = sqrt(UE**2 - we**2 -ve**2 )
    
    Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
    Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
    return (Xe, Ue)

def trimmerCurve(dynamic, HE, UE, dPS, BTA):
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
        we = Z[0]  
        ve = UE*sin(BTA)
        ue = sqrt(UE**2 - we**2 -ve**2 )
        phie = Z[1]
        psie = 0.0
        thetae = Z[2]
        pe = Z[3]
        qe = Z[4]
        re = Z[5]
    
        delta_pe = Z[6]
        delta_ee = Z[7]
        delta_ae = Z[8]
        delta_re = Z[9]
        
        Xe = array([xe, ye, ze, ue, ve, we, phie, thetae, psie, pe, qe, re])
        Ue = array([delta_pe, delta_ee, delta_ae, delta_re])
        
        sol = dynamic(0, Xe, Ue)
        Zp = array([sol[2], sol[3], sol[4], sol[5], degrees(sol[6]),degrees(sol[7]), degrees(sol[8]-dPS), degrees(sol[9]), degrees(sol[10]), degrees(sol[11])])
        return Zp
    
    wg = 0.0
    vg = 0.0
    thetag = radians(0)
    pg = 0.0
    qg = 0.0
    rg = 0.0
    delta_pg = 0.5
    delta_eg = radians(0.)
    delta_ag = 0.0
    delta_rg = 0.0
    Zg = array([wg, vg, thetag, pg, qg, rg, delta_pg, delta_eg, delta_ag, delta_rg])
    
    wlim = [-20, 20]
    vlim = [-50, 50]
    thetalim = [radians(-20), radians(20)]
    plim = [radians(-20), radians(20)]
    qlim = [radians(-20), radians(20)]
    rlim = [radians(-20), radians(20)]
    delta_p_lim = [0,1]
    delta_e_lim = [radians(-40), radians(40)]
    delta_a_lim = [radians(-40), radians(40)]
    delta_r_lim = [radians(-40), radians(40)]
    
    boundary1 = (wlim[0], vlim[0], thetalim[0], plim[0], qlim[0], rlim[0], delta_p_lim[0],delta_e_lim[0], delta_a_lim[0], delta_r_lim[0])
    boundary2 = (wlim[1], vlim[1], thetalim[1], plim[1], qlim[1], rlim[1], delta_p_lim[1],delta_e_lim[1], delta_a_lim[1], delta_r_lim[1])

    root = least_squares(obj, Zg, bounds = (boundary1, boundary2))
    we = root.x[0]
    ve = UE*sin(BTA)
    thetae = root.x[2]
    pe = root.x[3]
    qe = root.x[4]
    re = root.x[5]
    delta_pe = root.x[6]
    delta_ee = root.x[7]
    delta_ae = root.x[8]
    delta_re = root.x[9]
    xe = 0.0
    ye = 0.0
    ze = -HE
    phie = root.x[1]
    psie = 0.0
    
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

def linearization(dynamics, Xe, Ue):
    dX = array([0.01, 0.01, 0.01, 0.01, radians(0.01), radians(0.01), radians(0.01), radians(0.01), radians(0.01), radians(0.01)])
    dU = array([0.001, radians(0.01), radians(0.01), radians(0.01)])
    A = zeros((10,10))
    B = zeros((10,4))
    for i in range(0,len(dX)):
        Xm = copy(Xe)
        Xm[i+2] = Xm[i+2] + dX[i]
        
        sol = dynamics(0, Xm,Ue)
        A[:,i] = 1/dX[i]*sol[2:12]
        
    for i in range(0,len(dU)):
        Um = copy(Ue)
        Um[i] = Um[i] + dU[i]
        sol= dynamics(0, Xe, Um)
        B[:,i] = 1/dU[i]*sol[2:12]
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
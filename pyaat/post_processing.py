"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Plots
"""
from numpy import array, degrees
from tools import computeTAS
import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['axes.formatter.useoffset'] = False

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
            plt.plot(self.time, self.TAS, color = 'red')
            plt.ylabel('V [m/s]')
            plt.grid()
            plt.subplot(312)
            plt.plot(self.time, degrees(self.alpha), color = 'red')
            plt.ylabel(r'$\alpha$ [deg]')
            plt.grid()
            plt.subplot(313)
            plt.plot(self.time,degrees(self.beta), color = 'red')
            plt.ylabel(r'$\beta$ [deg]')
            plt.grid()
            plt.show()
        else:
            plt.figure()
            plt.subplot(311)
            plt.title("Velocity")
            plt.plot(self.time, self.u, color = 'red')
            plt.ylabel('u [m/s]')
            plt.grid()
            plt.subplot(312)
            plt.plot(self.time, self.v, color = 'red')
            plt.ylabel('v [m/s]')
            plt.grid()
            plt.subplot(313)
            plt.plot(self.time,self.w, color = 'red')
            plt.ylabel('w [m/s]')
            plt.grid()
            plt.show()
            
    def LinPos(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Linear position")
        plt.plot(self.time, self.x, color = 'red')
        plt.ylabel('x [m]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, self.y, color = 'red')
        plt.ylabel('y [m]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time,-self.z, color = 'red')
        plt.ylabel('H [m]')
        plt.grid()
        plt.show()
        
    def LinPos3D(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(self.x, self.y, -self.z, 'red')
        
    def Attitude(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Attitude")
        plt.plot(self.time, degrees(self.phi), color = 'red')
        plt.ylabel('$\Phi$ [deg]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, degrees(self.theta), color = 'red')
        plt.ylabel('$\Theta$ [deg]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time, degrees(self.psi), color = 'red')
        plt.ylabel('$\Psi$ [deg]')
        plt.grid()
        plt.show()
        
    def AngVel(self):
        plt.figure()
        plt.subplot(311)
        plt.title("Angular Speed")
        plt.plot(self.time, degrees(self.p), color = 'red')
        plt.ylabel('$p$ [deg/s]')
        plt.grid()
        plt.subplot(312)
        plt.plot(self.time, degrees(self.q), color = 'red')
        plt.ylabel('$p$ [deg/s]')
        plt.grid()
        plt.subplot(313)
        plt.plot(self.time, degrees(self.r), color = 'red')
        plt.ylabel('$r$ [deg/s]')
        plt.grid()
        plt.show()
        
    def Controls(self):
        plt.figure()
        plt.subplot(411)
        plt.title("Controls")
        plt.plot(self.time, self.delta_p, color = 'red')
        plt.ylabel('$\delta_p$ [%]')
        plt.grid()
        plt.subplot(412)
        plt.plot(self.time, degrees(self.delta_e), color = 'red')
        plt.ylabel('$\delta_e$ [deg]')
        plt.grid()
        plt.subplot(413)
        plt.plot(self.time, degrees(self.delta_a), color = 'red')
        plt.ylabel('$\delta_a$ [deg]')
        plt.grid()
        plt.subplot(414)
        plt.plot(self.time, degrees(self.delta_r), color = 'red')
        plt.ylabel('$\delta_r$ [deg]')
        plt.grid()
        plt.show()
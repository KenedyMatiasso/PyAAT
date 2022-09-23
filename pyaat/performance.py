"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains methods and functions for performance analysis
"""

from numpy import linspace, sqrt
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from pyaat.tools import RootFinder
from pyaat.constants import GRAVITY

class Envelope():
    def __init__(self, system = None, plot = True, limits = [0,3000]):
        # import classes
        self.atmosphere = system.atmosphere
        self.propulsion  = system.propulsion      
        
        # Import parameters
        self.qmax = system.aircraft.qmax
        self.mass = system.aircraft._mass
        self.S = system.aircraft.Sw
        self.CD0 = system.aircraft.coef0[0] + 0.05
        self.k = system.aircraft.k + 0.08
        self.CLmax = system.aircraft.CLmax
        
        # Create outputs 
        self.lim_prop_min = []
        self.lim_prop_max = []
        self.lim_aero = []
        self.lim_est = []
        self.h_prop_min = []
        self.h_prop_max = []
        self.limits = limits
        
        # Define variables
        self.g = GRAVITY
                
        self.celling = fsolve(self.findCeling , 10000)
        self.hlist = linspace(0, self.celling , 50)
        self.compute_limits()
        
        if plot == True:
            self.plot_results()
    @property
    def rho_min(self):
        Pmaxtotal = 0
        rhoin = 0
        nrhon = 0
        if type(self.propulsion) ==list:
            for prop in self.propulsion:
                Pmaxtotal += prop._Pmaxi
                rhoin += prop._Pmaxi*prop._rhoi
                nrhon += prop._Pmaxi*prop._nrho
                
            nrho = nrhon/Pmaxtotal
            rhoi = rhoin/Pmaxtotal
            
            term = rhoi**nrho/Pmaxtotal
                      
        return (4/3*term*sqrt(2*(self.mass*self.g)**3/self.S*
            sqrt(3*self.k**3*self.CD0)))**( 2/(2*nrho + 1))
        
    def findCeling(self, h):
        self.atmosphere.set_altitude(h)
        rho = self.atmosphere.get_airDensity()
        return self.rho_min - rho
    
    def aeroLim(self, rho):
        return sqrt(2*self.mass*self.g/(rho*self.S*self.CLmax))

    def estLim(self,rho):
        return sqrt(2*self.qmax)/rho
    
    def propLim(self, rho):
        def f(V):
            if type(self.propulsion) ==list:
                Pmax = 0
                for prop in self.propulsion:
                    Pmax += prop._Pmaxi*(rho/prop._rhoi)**prop._nrho
            else:
                Pmax = self.propulsion._Pmaxi*(self.Vi/V)*(rho/self.propulsion._rhoi)**self.propulsion._nrho
                
            return (0.5*rho*V**3*self.S*self.CD0 + 
                (2*self.k*(self.mass*self.g)**2)/(rho*V*self.S) - Pmax)

        r = RootFinder(0, 500, 1)
        args = ()
        roots = r.find(f, *args) 
        return roots
    
    def compute_limits(self):
        for h in self.hlist:
            self.atmosphere.set_altitude(h)
            rho = self.atmosphere.get_airDensity()   
            self.lim_aero.append(self.aeroLim(rho))
            
            try:
                self.lim_prop_min.append(min(self.propLim(rho)))
                self.h_prop_min.append(h)
            except:
                pass
            
            try:
                self.lim_prop_max.append(max(self.propLim(rho)))
                self.h_prop_max.append(h)
            except:
                pass
            
            self.lim_est.append(self.estLim(rho))
            
    def plot_results(self):
        plt.figure(figsize=(8,4))
        plt.plot(self.lim_aero, self.hlist, label ='Aerodynamic limit', color ='blue')
        plt.plot(self.lim_prop_min, self.h_prop_min, label ='Propulsive limit', color ='red')
        plt.plot(self.lim_prop_max, self.h_prop_max, color ='red')
        plt.plot(self.lim_est, self.hlist, label ='Structural limit', color ='green')
        plt.scatter(self.lim_prop_max[-1], self.celling, color = 'k', label = f"Service celling: {round(self.celling[0],2)} m")
        plt.xlim(self.limits)
        plt.grid()
        plt.xlabel(r'Velocity $[m/s]$')
        plt.ylabel(r'Altitude $[m]$')
        plt.legend(loc = 'upper right')
        #plt.show()
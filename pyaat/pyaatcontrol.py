#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:51:00 2021

@author: ydor9e
"""
import numpy as np
import copy

class doublet():
    def __init__(self, command, amplitude=1, t_init=1, T=1):
        self.controlsNames = None
        self.name = 'doublet'
        self.type = 'open-loop'
        self.number_of_states = 0
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.command = command
        self.amplitude = amplitude
        self.T = T
        self.t_init = t_init
          
    @property
    def U(self):
        U = self.Ue
        if self.controlsNames!=None:
            if self.t >= self.t_init and self.t <= (self.T/2 + self.t_init):
                try:
                    U[self.controlsNames.index(self.command)] += self.amplitude
                except:
                     pass
            elif self.t >= (self.T/2 + self.t_init) and self.t <= (self.T + self.t_init):
                try:
                    U[self.controlsNames.index(self.command)] -= self.amplitude
                except:
                     pass    
        return U

class equilibrium(object):
    def __init__(self):
        self.controlsNames = None
        self.name = 'equilibrium'
        self.number_of_states = 0
        self.type ='open-loop'
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        
    @property
    def U(self):            
        return self.Ue
    
class actuator():
    def __init__(self, name, command, tau = 0.00001):
        self.type ='actuator'
        self.name = name
        self.number_of_states = 1
        self.tau = tau
        self.command = command
        self.integral = 0
        self.control_action = 0
        self.ddelta = 0
        self.limits = None
    
    @property
    def derivative(self):
        if type(self.limits) ==list:
            if self.ddelta<=self.limits[0]:    
                return -1/self.tau*self.integral + 1/self.tau*self.limits[0]
            elif self.ddelta>=self.limits[1]:
                return -1/self.tau*self.integral + 1/self.tau*self.limits[1]
            else:
                return -1/self.tau*self.integral + 1/self.tau*self.ddelta
        else:
            return -1/self.tau*self.integral + 1/self.tau*self.ddelta
        
class sensor():
    def __init__(self, name, state, tau):
        self.t = 0
        self.type ='sensor'
        self.name = name
        self.number_of_states = 1
        self.tau = tau
        self.noise = 'off'
        self.state = state
        self.integral = 0
        self.control_action = 0
        self.ddelta = 0
        self.noiseNorm = 1
        self.noiseStep = 0.001
    
    @property
    def derivative(self):
        if self.noise =='on':
            #randomValue =  random.uniform(-self.noiseNorm, self.noiseNorm)
            randomValue = self.noiseNorm*np.sin(200*self.t)
            return -1/self.tau*self.integral + 1/self.tau*(self.ddelta + randomValue)
        else:
            return -1/self.tau*self.integral + 1/self.tau*self.ddelta

class step():
    def __init__(self, command, amplitude =1, t_init=1):
        self.type ='open-loop'
        self.name = 'step'
        self.number_of_states = 0
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.controlsNames = None
        self.command = command
        self.amplitude = amplitude
        self.t_init = t_init
    
    @property
    def U(self):
        U = self.Ue
        if self.controlsNames!=None:
            if self.t >= self.t_init:
                try:
                    U[self.controlsNames.index(self.command)] += self.amplitude
                except:
                     pass
        return U
    
class loop():
    def __init__(self, name, Xreference, Xmeasure, command):
        self.controlsNames = None
        self.statesNames = None
        self.type ='closed-loop'
        self.name = name
        self.loopControls = []
        self.Xreference = Xreference
        self.Xmeasure = Xmeasure
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.command = command
        self.Xp_local = []
            
    @property
    def names(self):
        names =[]
        for cont in self.loopControls:
            if cont.number_of_states > 0:
                names.append(cont.name)
        return names
    
    @property
    def number_of_states(self):
        number_of_states = 0
        for cont in self.loopControls:
            number_of_states += cont.number_of_states
        return number_of_states
    
    @property
    def U(self):
        Uapp = np.zeros(np.shape(self.Ue))
        self.Xp_local = []
        self.loopControls[0].error = self.X[self.statesNames.index(self.Xmeasure)] - self.Xe[self.statesNames.index(self.Xreference)]
        self.loopControls[0].statesNames = self.statesNames
        self.loopControls[0].X = self.X
        self.loopControls[0].Xe = self.Xe
        self.loopControls[0].command = self.controlsNames.index(self.command)
        if self.loopControls[0].number_of_states !=0:
            self.Xp_local.append(self.loopControls[0].derivative)
        Uapp[self.controlsNames.index(self.command)] = self.loopControls[0].U
        for l in range(len(self.loopControls) -1):
            i = l+1
            self.loopControls[i].statesNames = self.statesNames
            self.loopControls[i].X = self.X
            self.loopControls[i].Xe = self.Xe
            self.loopControls[i].command = self.controlsNames.index(self.command)
            self.loopControls[i].error = self.loopControls[i-1].U
            if self.loopControls[i].number_of_states !=0:
                self.Xp_local.append(self.loopControls[i].derivative)
            Uapp[self.controlsNames.index(self.command)] = self.loopControls[i].U
        
        Uapp += self.Ue
        return Uapp

class mimoloop():
    def __init__(self):
        self.controlsNames = None
        self.statesNames = None
        self.type ='closed-loop'
        self.name = 'loop Name'
        self.loopControls = []
        self.Xreference = None
        self.Xmeasure = None
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.command = 0
        self.Xp_local = []
    
    @property
    def names(self):
        names =[]
        for cont in self.loopControls:
            if cont.number_of_states > 0:
                if type(cont.name) == list:
                    names.extend(cont.name)
                else:
                    names.append(cont.name)
        return names
    
    @property
    def number_of_states(self):
        number_of_states = 0
        for cont in self.loopControls:
            number_of_states += cont.number_of_states
        return number_of_states
    
    @property
    def U(self):
        Uapp = np.zeros(len(self.controlsNames))
        self.Xp_local = []
        error = []
        for x in range(len(self.Xmeasure)):          
            state_measurement = self.X[self.statesNames.index(self.Xmeasure[x])]
            state_reference = self.Xe[self.statesNames.index(self.Xreference[x])]
            error.append(state_measurement - state_reference)
        
        self.loopControls[0].error = np.array(error)
        self.loopControls[0].statesNames = self.statesNames
        self.loopControls[0].X = self.X
        self.loopControls[0].Xe = self.Xe
        if self.loopControls[0].number_of_states !=0:
            self.Xp_local.extend(self.loopControls[0].derivative)
            
        U = self.loopControls[0].U
        for val in range(len(U)):
            Uapp[self.controlsNames.index(self.command[val])] += self.loopControls[0].U[val]
            
        Uapp += self.Ue
        return Uapp

class P():
    def __init__(self, name, gain):
        self.controlsNames = None
        self.number_of_states = 0
        self.type = 'compensator'
        self.name = name
        self.gain = gain
        self.error = None
        self.command = 0
        self.Xmeasure = None
        self.Xreference = None
        self.X = None
        
    @property
    def U(self):
        if self.Xmeasure == None:
            return self.gain*self.error
        else:
            if self.Xreference != None:
                error2 = self.X[self.statesNames.index(self.Xmeasure)] - self.Xe[self.statesNames.index(self.Xreference)]
                return self.gain*(self.error - error2)
            else:
                error2 = self.X[self.statesNames.index(self.Xmeasure)]
                return self.gain*(self.error - error2) 
    @property
    def derivative(self):
        return 0.0
    
class LowPassFilter():
    def __init__(self, name, tau):
        self.type ='filter'
        self.name = name
        self.number_of_states = 1
        self.tau = tau
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.gain = None
        self.error = None
        self.coomand = 0
        
    @property
    def U(self):
        return self.X[self.statesNames.index(self.name)]
    
    @property
    def derivative(self):
        
        integral = self.X[self.statesNames.index(self.name)]
        ddelta = self.error
        return -1/self.tau*integral + 1/self.tau*ddelta
    
class HighPassFilter():
    def __init__(self, name, tau):
        self.type ='filter'
        self.name = name
        self.number_of_states = 1
        self.tau = tau
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.gain = None
        self.error = None
    
    @property
    def derivative(self):      
        return -1/self.tau*self.X[self.statesNames.index(self.name)] + 1/self.tau*self.error
    
    @property
    def U(self):
        return -self.X[self.statesNames.index(self.name)] + self.error
    
class leadlag():
    def __init__(self):
        self.type ='compensator'
        self.name = 'compensator Name'
        self.number_of_states = 1
        self.pole = None
        self.zero = None
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.gain = None
        self.error = None
    
    @property
    def derivative(self):
        return -self.pole*self.X[self.statesNames.index(self.name)] + (self.pole -self.zero)*self.error
    
    @property
    def U(self):
        return -self.X[self.statesNames.index(self.name)] + self.error
        
class PI():
    def __init__(self, name, zero, Kp):
        self.type ='compensator'
        self.name = name
        self.number_of_states = 1
        self.zero = zero
        self.Kp = Kp
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.gain = None
        self.error = None
    
    @property
    def derivative(self):
        return  self.error
    
    @property
    def U(self):
        return -self.zero *self.X[self.statesNames.index(self.name)] + self.Kp*self.error
    
class PImimo():
    def __init__(self):
        self.type ='compensator'
        self.name = 'compensator Name'
        self.number_of_states = 1
        self.zero = None
        self.Kp = None
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.gain = None
        self.error = None
    
    @property
    def derivative(self):
        return self.error
    
    @property
    def U(self):
        return -self.zero *self.X[self.statesNames.index(self.name)] + self.Kp*self.error

class lqr():
    def __init__(self):
        self.type ='compensator'
        self.name = 'compensator Name'
        self.number_of_states = 0
        self.K = None
        self.controlsNames = None
        self.statesNames = None
        self.X = None
        self.error = None
    
    @property
    def derivative(self):
        return  np.zeros(np.shape(self.K)[1])
    
    @property
    def U(self):
        return self.K.dot(np.array(self.error))
    
    
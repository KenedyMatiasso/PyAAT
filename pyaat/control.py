#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:51:00 2021

@author: ydor9e
"""
from numpy import radians, array

class doublet(object):
    def __init__(self):
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.command = 'elevator'
        self.amplitude = 1
        self.T = 1
        self.t_init = 1
    
    @property
    def U(self):
        if self.command == 'thrust':
            if self.t >= self.t_init and self.t <= (self.T/2 + self.t_init):
                delta_p = self.Ue[0] + self.amplitude
            elif self.t >= (self.T/2 + self.t_init) and self.t <= (self.T + self.t_init):
                delta_p = self.Ue[0] - self.amplitude
            else:
                delta_p = self.Ue[0]
                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
       
        elif self.command == 'elevator':

            if self.t >= self.t_init and self.t <= (self.T/2 + self.t_init):
                delta_e = self.Ue[1] + radians(self.amplitude)
                
            elif self.t >= (self.T/2 + self.t_init) and self.t <= (self.T + self.t_init):
                delta_e = self.Ue[1] - radians(self.amplitude)
            else:
                delta_e = self.Ue[1]
                
            delta_p = self.Ue[0]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
            
        elif self.command == 'aileron':
            if self.t >= self.t_init and self.t <= (self.T/2 + self.t_init):
                delta_a = self.Ue[2] + radians(self.amplitude)
                
            elif self.t >= (self.T/2 + self.t_init) and self.t <= (self.T + self.t_init):
                delta_a = self.Ue[2] - radians(self.amplitude)
                
            else:
                delta_a = self.Ue[2]

            delta_p = self.Ue[0]    
            delta_e = self.Ue[1]
            delta_r = self.Ue[3]
            
        elif self.command == 'rudder':
            if self.t >= self.t_init and self.t <= (self.T/2 + self.t_init):
                delta_r = self.Ue[3] + radians(self.amplitude)
            elif self.t >= (self.T/2 + self.t_init) and self.t <= (self.T + self.t_init):
                delta_r = self.Ue[3] - radians(self.amplitude)
            else:
                delta_r = self.Ue[3]
                
            delta_p = self.Ue[0]                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            
        else:
            delta_p = self.Ue[0]                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
        return array([delta_p, delta_e, delta_a, delta_r])


class equilibrium(object):
    def __init__(self):
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        
    @property
    def U(self):            
        return array([self.Ue[0], self.Ue[1], self.Ue[2], self.Ue[3]])
    
class step(object):
    def __init__(self):
        self.t = 0.0
        self.Xe = None
        self.X = None
        self.Ue = None
        self.command = 'elevator'
        self.amplitude = 1
        self.t_init = 1
    
    @property
    def U(self):
        if self.command == 'thrust':
            if self.t >= self.t_init:
                delta_p = self.Ue[0] + self.amplitude
            else:
                delta_p = self.Ue[0]
                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
       
        elif self.command == 'elevator':

            if self.t >= self.t_init:
                delta_e = self.Ue[1] + radians(self.amplitude)
            else:
                delta_e = self.Ue[1]
                
            delta_p = self.Ue[0]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
            
        elif self.command == 'aileron':
            if self.t >= self.t_init:
                delta_a = self.Ue[2] + radians(self.amplitude)                
            else:
                delta_a = self.Ue[2]

            delta_p = self.Ue[0]    
            delta_e = self.Ue[1]
            delta_r = self.Ue[3]
            
        elif self.command == 'rudder':
            if self.t >= self.t_init:
                delta_r = self.Ue[3] + radians(self.amplitude)
            else:
                delta_r = self.Ue[3]
                
            delta_p = self.Ue[0]                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            
        else:
            delta_p = self.Ue[0]                
            delta_e = self.Ue[1]
            delta_a = self.Ue[2]
            delta_r = self.Ue[3]
            
        return array([delta_p, delta_e, delta_a, delta_r])
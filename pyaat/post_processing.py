"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Plots
"""
from numpy import array

class states(object):
    def __init__(self, X = array([0,0,0,0,0,0,0,0,0,0,0,0])):
        self.V= X[:0]
        self.alpha= X[:1]
        
    def plotLatero(self):
        

class Controls(object):
    def __init__(self, U = array([0,0,0,0])):
        self.de = U[:1]
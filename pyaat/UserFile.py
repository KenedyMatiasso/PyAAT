"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This is the user file.
"""

from atmosphere import atmosISA
from aircraft import Aircraft
from propulsion import SimpleModel
from gravity import NewtonGravity
from system import system

from tools import printInfo
from tools import plotter

from numpy.linalg import eig

#from numpy import radians, tan


atm = atmosISA()
prop = SimpleModel()
airc = Aircraft()
grav = NewtonGravity()

Mysys = system(atmosphere =atm, propulsion = prop, aircraft = airc, gravity=grav)

Xe, Ue = Mysys.trimmer(condition='cruize', HE =5000., VE= 150, dPS= 2, BTA =2, dTH = 10.,dH =5.)


printInfo(Xe,Ue, frame ='aero')
printInfo(Xe,Ue, frame='controls')

solution, control = Mysys.propagate(Xe, Ue, TF =100, perturbation=False, state={'alpha': 5, 'beta':2})

pltr = plotter()
pltr.states = solution
pltr.time = Mysys.time
pltr.control = control

pltr.LinVel(frame = 'body')
pltr.LinVel(frame = 'aero')
pltr.LinPos()
pltr.Attitude()
pltr.AngVel()
pltr.Controls()
pltr.LinPos3D()
pltr.linPos2D()

A, B = Mysys.LinearModes(Xe, Ue)
Ald, Bld = Mysys.LinearLatero(Xe,Ue)
Al, Bl =Mysys.LinearLong(Xe,Ue)
wld, vld = eig(Ald)
wl, vl = eig(Al)

print('--------------------------------')
print('Eigenvalues')
print(wl)
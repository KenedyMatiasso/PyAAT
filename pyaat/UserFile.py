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


atm = atmosISA()
prop = SimpleModel()
airc = Aircraft()
grav = NewtonGravity()

Mysys = system(atmosphere =atm, propulsion = prop, aircraft = airc, gravity=grav)

Xe, Ue = Mysys.trimmer(condition='curve', HE =5000., VE = 150, dPS = 2, BTA =2)

printInfo(Xe,Ue, frame ='aero')
printInfo(Xe,Ue, frame='controls')

solution, control = Mysys.propagate(Xe, Ue, TF =100)

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

'''
#sol = list(dynamics(Xe,Ue))




Xe[5] = Xe[5] + 0.0
Xe[4] = Xe[4] + 0.0
Xe[3] = Xe[3] + 0.0

Ue[1] = Ue[1] + radians(0)
Ue[2] = Ue[2] + radians(0)

solution = odeint(simularaviao, Xe, time)









A, B = linearization(dynamics, Xe, Ue)
Af, Bf = modesMatrix(A,B)
An, Bn = lateroMatrix(A,B)
Al, Bl = longMatrix(A, B)
#w, v = eig(A)

'''
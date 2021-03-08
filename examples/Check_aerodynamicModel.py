"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Example using aerodynamic model
"""

from Aerodynamic import Aircraft
from numpy import radians
from atmosphere import atmosISA

atmosphere = atmosISA()
atmosphere._altitude = 10000
rho= atmosphere._rho


aircraft = Aircraft()
aircraft.TAS = 200.0
aircraft.alpha = radians(3.)
aircraft.beta = radians(0.0)
aircraft.p = radians(0.)
aircraft.r = radians(0.)
aircraft.q = radians(0.)
aircraft.delta_a = radians(0.)
aircraft.delta_r = radians(0.)
aircraft.delta_e = radians(-4.55)


aircraft._rho = rho
print(aircraft._rho)
print(aircraft.AeroForces)
print(aircraft.AeroMoments)

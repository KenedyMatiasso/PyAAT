"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Before you run this example, move this file to forder PyAAt/pyaat

Example using propulsion
"""

from propulsion import SimpleModel
from atmosphere import atmosISA
from numpy import arange
import matplotlib.pyplot as plt


prop = SimpleModel()
atmos = atmosISA()

# Thust as function of altitude
hlist = arange(0,12000,10)
prop.delta_p = 0.3
thrust = []

for h in hlist:
    atmos._altitude = h
    prop.rho = atmos._rho
    thrust.append(prop.Forces[0])
    
plt.figure()
plt.plot(thrust,hlist, color='r', linestyle = '-')
plt.xlabel('Thrust [N]')
plt.ylabel('Altitude')
plt.title('Thrust with altitude for a turbofan engine')
plt.legend()
plt.grid()
plt.show()

# Thust as function of altitude and troatle

troatlelist = arange(0,1,0.1)

thrust = []
thrust2 = []
thrust3 = []

h1 = 4000
h2 = 8000
h3 = 12000

atmos._altitude = h1
prop.rho = atmos._rho

for tr in troatlelist:
    prop.delta_p = tr
    thrust.append(prop.Forces[0])
    
atmos._altitude = h2
prop.rho = atmos._rho

for tr in troatlelist:
    prop.delta_p = tr
    thrust2.append(prop.Forces[0])
    
atmos._altitude = h3
prop.rho = atmos._rho

for tr in troatlelist:
    prop.delta_p = tr
    thrust3.append(prop.Forces[0])


plt.figure()
plt.plot(thrust,troatlelist*100, color='r', label = '4000 meters', linestyle = '-')
plt.plot(thrust2,troatlelist*100, color='k', label = '8000 meters', linestyle = '-')
plt.plot(thrust3,troatlelist*100, color='b', label = '12000 meters', linestyle = '-')

plt.xlabel('Thrust [N]')
plt.ylabel('troatle [%]')
plt.title('Thrust with troatle for a turbofan engine')
plt.legend()
plt.grid()
plt.show()
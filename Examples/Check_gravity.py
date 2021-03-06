"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

Before you run this example, move this file to forder PyAAt/pyaat

Code for checking gravity models
"""

from gravity import VerticalConstant, NewtonGravity, HighOrder
from numpy import arange, radians, degrees
import matplotlib.pyplot as plt

# Compair three gravity at Seal level and 45deg latitude
g1 = VerticalConstant()
g2 = NewtonGravity()
g3 = HighOrder()

g2._altitude=0
g3._latitude = radians(45)
print('Vertical constant')
print(g1._gravity)
print('-------------------------')
print('Newton gravity (sphere)')
print(g2._gravity)
print('-------------------------')
print('High order model')
print(g3._gravity)

# Example 3.1 from A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.

"""
Example 3.1. Construct a model of the earth’s gravity using the ﬁrst four
Jeﬀery’s constants in the series expansion of gravitational potential. Compare
the acceleration due to gravity with that of the spherical earth model (R =
Re = 6378.14 km) for a trajectory in which the latitude (in degrees) varies
with altitude, h = r − R e (in kilometers), as follows:
lat = h − 100, (0 ≤ h ≤ 200 km)
"""

hlist = arange(0,200,1)
g = []
gb = []
latlist = []
for h in hlist:
    lat = radians(h-100)
    g3._altitude = h*1000
    g2._altitude = h*1000
    g3._latitude = lat
    g.append(g3._gravity[2])
    gb.append(g2._gravity[2])
    latlist.append(degrees(lat))

plt.figure()
plt.plot(latlist,gb, label ='Newton gravity', color = 'r', linestyle='-')
plt.plot(latlist,g, label='High Order model', color='k', linestyle = '--')
plt.ylabel('$-g_r$ [$m/s^2$]')
plt.xlabel('Latitude $\delta$ [deg]')
plt.title('Gravity as function of latitude and altitude')
plt.legend()
plt.grid()
plt.show()
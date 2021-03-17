"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of constants that can be used by PyAAT.
"""

from numpy import pi
# TODO: implement orbital parameters, moon proprieties, distances between celestial objects, etc.

########### Air constants at sea level ###########
# Source: M. V. Cook., Flight Dynamics Principles, 2nd edition, Oxford: Elsivier Ltda, 2007.
RHO_SEA = 1.225             # Sea level air density [kg/m**3]
SOUND_VEL_SEA = 340.29      # Speed of sound at sea level [m/s]

# Source: U.S Standard Atmosphere 1976, NASA, 1976.
# >https://ntrs.nasa.gov/search.jsp?R=19770009539<
GAMMA_AIR = 1.4             # specific heat at 20 degree Celsius
TEMP_SEA = 288.15           # temperature at sea level [K]
DYNAMIC_VISC = 1.7894e-5    # dynamic Viscosity [N*s/(m**2)]
KINEMATIC_VISC = 1.4607e-5  # kinematic Viscosity [m**2/s]
PRESURE_SEA = 101.325e3     # pressure at seal level [Pascal]
R_AIR = 287.05287           # ideal gas constant [J/(Kg*K)]

########### Earth proprieties ###########
# source: B Wie, Space vehicle Dynamics and Control, 2nd edition, Virginia: AIAA Education Series, 2008.
GRAVITY = 9.80665           # mean constant gravitational field of Earth [m/s**2]
R_EARTH = 6378e3            # mean equatorial radius of the Earth [m]

# source: A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.
G = 6.67259e-11             # universal gravitational constant [m**3/kg/s**2]
MU_EARTH = 3.986004418e14      # G*M_EARTH [m**3/s**2]
M_EARTH = 5.972e24          # Earth's mass [kg]
J2 = 1.08263e-3             # Jeffery’s constant for Earth's 2nd spherical harmonic [null]
J3 = -2.532153e-7           # Jeffery’s constant for Earth's 3th spherical harmonic [null]
J4 = -1.6109876e-7          # Jeffery’s constant for Earth's 4th spherical harmonic [null]

########### Mars proprieties ###########
# source: A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.
R_MARS = 3397e3             # mean equatorial radius of the Mars [m]
MU_MARS = 42828.3719e9      # G*M_MARS on mars [m**3/s**2]
M_MARS = 6.39e23            # Mars's mass [kg]
J2MARS = 0.001955453679     # Jeffery’s constant for Mars's 2nd spherical harmonic [null]
J3MARS = 3.144980942e-5     # Jeffery’s constant for Mars's 3th spherical harmonic [null]
J4MARS = -1.53773961e-5     # Jeffery’s constant for Mars's 4th spherical harmonic [null]
J5MARS = 5.718547184e-6     # Jeffery’s constant for Mars's 5th spherical harmonic [null]

########### Sun proprieties ###########
# source: A. Tewary, Atmospheric ans Space Flight Dynamics, Boston: Birkhauser, 2007.
MU_SUN = 1.32712e20         # G*M_SUN [m**3/(s**2)]
M_SUN = 1.989e30            # Sun's mass [kg]

########### Conversion of units from Imperial Units to Systéme International (SI) ###########
# Source: E.L Houghton et al, Aerodyamics for Engineering Students, 6th edition, Oxford: Elsevier Ltda, 2013.
ft2m = 0.3048               # feet [ft] to meters [m]
in2mm = 25.4                # inch [in] to millimeters [mm]
nauticMile2m = 1853.2       # nautic miles to meters [m]
statuteMile2m = 1609.3      # satute miles to meters [m]
ft_sq2m_sq = 0.0929         # ft**2 to m**2
in_sq2m_sq = 6.4516e-4      # in**2 to m**2
in_sq2mm_sq = 645.16        # in**2 to mm**2
in_cub2m_cub = 1.6387e-5    # in**3 to m**3
in_cub2mm_cub = 16387       # in**3 to mm**3
slug2kg = 14.594            # slug [slug] to kilograms [kg]
slugPft_cub2kgPm_cub = 515.38   # slug/(ft**3) to kg/(m**3)
lbf2N = 4.4482              # pounds force [lbf] to Newton [N]
lbfPft_sq2NPm_sq = 47.880   # lbf/(ft**2) to N/(m**2)
lfbPin_sq2NPm_sq = 68948    # lbf/(in**2) to N/(m**2)
ftlbf2J = 1.3558            # lbf*ft to Joule [J]
hp2W = 745.7                # horse power [hp] to Watt [W]
lbfft2Nm = 1.3558           # lbf*ft to N*m
ftPs2mPs = 0.3048           # feet per second [ft/s] to meters per second [m/s]
milePhr2mPs = 0.44704       # mile per hour to m/s
knot2mPs = 0.21477          # nautical mile per hour to m/s

########### Other conversions ###########
# Source: myself
deg2rad = pi/180            # Degree [^o] to radians [rad]
rpm2radPs = 30/pi           # rotation per minute [RPM] to radians per second [rad/s]

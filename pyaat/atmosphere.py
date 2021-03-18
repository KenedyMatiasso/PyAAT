"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2020 Kenedy Matiasso Portella
Distributed under MIT License

Version control:
Sep 05 2020                          KMP: File creation.
Sep 05 2020                          KMP: Implements atmosISA, atmosCOESA, seaLevel.
Nov 12 2020                          KMP: Implements simplified model;
                                          Computes T, rho and P as @property.

This file implements the atmosphere models that can be used by PyAAT.
"""

from numpy import sqrt, array, interp
from pyaat.constants import R_AIR, GRAVITY, GAMMA_AIR, TEMP_SEA, PRESURE_SEA, RHO_SEA, SOUND_VEL_SEA

class atmosISA(object):
    """
    Implements International Standard Atmosphere (ISA).
    source: AIRBUS, Getting to grips with aircraft performance, Toulouse 2002.

    """
    def __init__(self, _altitude = 0):  
        self._altitude = 0           # altitude [m]     
        self._parameters = array(
            [[0., 1524., 3048., 4572., 6096., 7620., 9144., 10668., 12192.],                    # geometric altitude [m]
             [288.15, 278.25, 268.35, 258.45, 248.55, 238.65, 228.75, 218.85, 216.65],          # base atmospheric temperature [K]
             [101300., 84300., 69700, 57200, 46600, 37600, 30100, 23800, 18800]])      # base atmospheric pressure [Pa]
          # Store information about atmosphere behavior as functions
          
                                        # of geometric altitude

        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        
    @property
    def _temperature(self): # temperature [K]
        return interp(self._altitude, self._parameters[0], self._parameters[1])
    
    @property 
    def _pressure(self): # pressure [Pa]      
        return interp(self._altitude, self._parameters[0], self._parameters[2])
    
    @property
    def _rho(self): # air density [kg/m³]
        return self._pressure / (R_AIR * self._temperature)
        
    @property
    def _soundVel(self):  # velocity of sound [m/s]
        return sqrt(GAMMA_AIR * R_AIR * self._temperature)

class atmosCOESA(object):
    """
    Implements US Standard Atmosphere (COESA) 1976.
    source: U.S Standard Atmosphere 1976, NASA, 1976. >https://ntrs.nasa.gov/search.jsp?R=19770009539<
    Velid up to 86000
    """
    
    def __init__(self, _altitude = 0):
        self._altitude = None           # altitude [m]     
        self._parameters = array([[-1000, -500, 0, 500, 1000, 1500, 2000, 2500,
                                   3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500,
                                   7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500,
                                   11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500,
                                   15000, 15500, 16000, 16500, 17000, 17500, 18000, 18500,
                                   19000, 19500, 20000, 21000, 22000, 23000, 24000, 25000,
                                   26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000,
                                   34000, 35000, 36000, 37000, 38000, 39000, 40000, 41000,
                                   42000, 43000, 44000, 45000, 46000, 47000, 48000, 49000,
                                   50000, 55000, 60000, 65000, 70000, 75000, 80000, 85500],
                                  [294.651, 291.4, 288.15, 284.9, 281.651, 278.402, 275.154,271.906,
                                   268.659, 265.413, 262.166, 258.921, 255.676, 252.431, 249.187, 245.943,
                                   242.7, 239.457, 236.215, 232.974, 229.733, 226.492, 223.252, 220.013,
                                   216.774, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65,
                                   216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65,216.65,
                                   216.65, 216.65, 216.65, 217.581, 218.574, 219.567, 220.56, 221.552,
                                   222.554, 223.536, 224.527, 225.518, 226.509, 227.5, 228.49, 230.973,
                                   233.743, 236.513, 239.282, 242.050, 244.818, 247.031, 250.35, 253.114,
                                   255.878, 258.641, 261.403, 264.164, 266.925, 269.684, 270.65, 270.65,
                                   270.65, 260.771, 247.021, 233.292, 219.585, 208.399, 198.639, 187.920],
                                  [1.1393e3, 1.0747e3, 1.01325e3, 9.5461e2, 8.9876e2, 8.4559e2, 7.9501e2, 7.4691e2,
                                   7.0121e2, 6.5780e2, 6.1660e2, 5.7752e2, 5.4048e2, 5.0539e2, 4.7217e2, 4.4075e2,
                                   4.1105e2, 3.8299e2, 3.5651e2, 3.3154e2, 3.08e2, 2.8584e2, 2.6499e2, 2.4540e2,
                                   2.2699e2, 2.0984e2, 1.9399e2, 1.7932e2, 1.6579e2, 1.5327e2, 1.4170e2, 1.31e2,
                                   1.2111e2, 1.1197e2, 1.0352e2, 9.5717e1, 8.8497e1, 8.1822e1, 7.5652e1, 6.9948e1,
                                   6.4674e1, 5.9799e1, 5.5293e1, 4.7289e1, 4.0475e1, 3.4668e1, 2.9717e1, 2.5492e1,
                                   2.1883e1, 1.8799e1, 1.6161e1, 1.3904e1, 1.1970e1, 1.0313e1, 8.8906, 7.6730,
                                   6.6341, 5.7459, 4.9852, 4.3324, 3.7713, 3.2882, 2.8714, 2.5113,
                                   2.1996, 1.9295, 1.6949, 1.4910, 1.3134, 1.1585, 1.0229, 9.2610e-1,
                                   7.9779e-1, 4.2525e-1, 2.1958e-1, 1.0929e-1, 5.2209e-2, 2.3881e-2, 1.0524e-2, 4.0802e-3]])
          # Store information about atmosphere behavior as functions
          
                                        # of geometric altitude

        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        
    @property
    def _temperature(self): # temperature [K]
        return interp(self._altitude, self._parameters[0], self._parameters[1])
    
    @property 
    def _pressure(self): # pressure [Pa]      
        return 100*interp(self._altitude, self._parameters[0], self._parameters[2])
    
    @property
    def _rho(self): # air density [kg/m³]
        return self._pressure / (R_AIR * self._temperature)
        
    @property
    def _soundVel(self):  # velocity of sound [m/s]
        return sqrt(GAMMA_AIR * R_AIR * self._temperature)

class SimplifiedModel(object):
    """
    Simplified model implemented on M. V. Cook., Flight Dynamics Principles, 2nd edition, Oxford: Elsivier Ltda, 2007.
    valid up to 32000 meters
    """
    
    def __init__(self, _altitude = 0):
        self._pressure = PRESURE_SEA    # pressure [Pa]
        self._altitude = 0.0            # altitude [m]
        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        self._Ir = -0.0065              # Lapse rate [K/m]    
        
    @property
    def _temperature(self): # temperature [K]
        return 288.16 + self._Ir * self._altitude
    
    @property
    def _rho(self): # air density [kg/m³]
        return RHO_SEA *(self._temperature/288.16)**(-(GRAVITY/(self._Ir*R_AIR)+1))
        
    @property
    def _soundVel(self):  # velocity of sound [m/s]
        return sqrt(GAMMA_AIR * R_AIR * self._temperature)

class seaLevel(object):
    """
    Implements atmosphere condition at sea level 20 degrees celsius.
    Source: M. V. Cook., Flight Dynamics Principles, 2nd edition, Oxford: Elsivier Ltda, 2007.
    source: U.S Standard Atmosphere 1976, NASA, 1976. >https://ntrs.nasa.gov/search.jsp?R=19770009539<
    """

    def __init__(self, _altitude = 0):
        self._pressure = PRESURE_SEA    # pressure [Pa]
        self._altitude = 0.0            # altitude [m]
        self._temperature = TEMP_SEA    # temperature [K]
        self._soundVel = SOUND_VEL_SEA  # velocity of sound [m/s]
        self._rho = RHO_SEA      # air density [kg/m³]
        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description

class windTunnel(object):
    # TODO: Implement /atmosphere/windTunnel

    """
    Using this class the user is able to set al atmospheric conditions manually to simulate wind tunnel conditions.

    By default it is considered an open-circuit wind tunnel at sea level.

    """
    def __init__(self):
        self._pressure = None    # pressure [Pa]
        self._altitude = None            # altitude [m]
        self._temperature = None    # temperature [K]
        self._soundVel = None  # velocity of sound [m/s]
        self._rho = None      # air density [kg/m³]
        self._flagError = 'ERROR'       # flag error
        self._error = 'Not implemented yet'                # error description


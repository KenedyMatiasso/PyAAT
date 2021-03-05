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
from constants import R_AIR, GRAVITY, GAMMA_AIR, TEMP_SEA, PRESURE_SEA, RHO_SEA, SOUND_VEL_SEA

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

    """
    
    def __init__(self, _altitude = 0):
        self._altitude = None           # altitude [m]     
        self._parameters = array(
            [[-611., 11019., 20063., 32162., 47350., 51413., 71802., 86000.],       # geometric altitude [m]
             [292.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.87],      # base atmospheric temperature [K]
             [108900., 22632., 5474.9, 868.02, 110.91, 66.94, 3.95, 0.3734]])       # base atmospheric pressure [Pa]
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

class SimplifiedModel(object):
    """
    Simplified model implemented on Cook.
    Implemented only to compair results aginst the book results
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


"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2020 Kenedy Matiasso Portella
Distributed under MIT License

Version control:
Jan 21 2021                 Creation

This file implements the wind models that can be used by PyAAT.
"""
import numpy as np
from abc import ABC, abstractmethod

class metaWind(ABC):
    def __init__(self, MaxSpeed = 1, x_start = 0.0, x_end = 1.0, y_start = 1.0, y_end = 10,
                z_start = 0.0, z_end = 0.0):
        self._flagError = 'VALID'       # flag error
        self._error = ''                # error description
        self._X = None
        self._MaxSpeed = MaxSpeed
        self._x_start = x_start
        self._x_end = x_end
        self._y_start = y_start
        self._y_end = y_end
        self._z_start = z_start
        self._z_end = z_end

    @property
    @abstractmethod
    def _windSpeed(self):
        pass

    def set_states(self, entry):
        self._X = entry

    def get_windSpeed(self):
        return self._windSpeed
    
class windFAR25(metaWind):
    """
    FAR 25.341-1 gust
    U.S. Department of Transportation. Advisory Circular No. 25.241-1. [S.l.], 2014. ANM-
    115, 24 p. Acesso em 14 dez. 2021. Disponível em: <https://www.faa.gov/documentLibrary/
    media/Advisory_Circular/AC_25_341-1.pdf>.
    """
    
    def __init__(self, MaxSpeed = 1, x_start = 0.0, x_end = 1.0, y_start = 1.0, y_end = 10,
                z_start = 0.0, z_end = 0.0):
        super().__init__(MaxSpeed, x_start, x_end, y_start, y_end, z_start, z_end)

    @property
    def _windSpeed(self):
        Hw = (self._x_end - self._x_start)/2
        Vw = 0
        if self._x_start < self._X[0]:
            if self._x_end > self._X[0]:
                if self._y_start < self._X[1]:
                    if self._y_end > self._X[1]:                      
                        Vw = self._MaxSpeed/2*(1 - np.cos(np.pi*(self._X[0]- self._x_start)/Hw))
        return np.array([0, 0, Vw])
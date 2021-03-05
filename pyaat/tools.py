from numpy import cos, sin, array, radians, degrees

def C1(theta):
    return array([[1, 0, 0],
                [0, cos(theta), sin(theta)],
                [0, -sin(theta), cos(theta)]])    
    
def C2(theta):
    return array([[cos(theta), 0, -sin(theta)],
                [0, 1, 0],
                [sin(theta), 0, cos(theta)]])
    
def C3(theta):
    return array([[cos(theta), sin(theta), 0],
                [-sin(theta), cos(theta), 0],
                [0,0,1]])    

def earth2body(vector, psi, tht, phi):        
    return C1(phi).dot(C2(tht)).dot(C3(psi)).dot(vector)

def body2earth(vector, psi,tht, phi):
    return C3(phi).dot(C2(tht)).dot(C1(psi)).dot(vector)

def aero2body(vector,alpha,beta):
    pass

def body2aero(vector, psi, tht, phi):
    pass



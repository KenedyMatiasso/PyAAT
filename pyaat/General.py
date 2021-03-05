"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains most of constants that can be used by PyAAT.
"""
from numpy import sin, cos, tan, interp, array, pi
from scipy.optimize import fsolve


R_AIR = 287.05287
H=10000
Vinf =224.6
m=45e3
Ixx = 554e3
Iyy = 2.53e6
Izz=3.01e6
Ixz=106e3
Cbarw=3.666
bw=28.42
Sw=95
zt=1.42
alphaf=0
Fmaxi =70e3
Vi=200
nv=0
rhoi=0.41271
nrho=0.775
g=9.80665


def ISA(_altitude):
    _parameters = array(
            [[0., 1524., 3048., 4572., 6096., 7620., 9144., 10668., 12192.],                    # geometric altitude [m]
             [288.15, 278.25, 268.35, 258.45, 248.55, 238.65, 228.75, 218.85, 216.65],          # base atmospheric temperature [K]
             [101300., 84300., 69700, 57200, 46600, 37600, 30100, 23800, 18800]])      # base atmospheric pressure [Pa]
    
    _pressure = interp(_altitude,_parameters[0], _parameters[2])
    _temperature= interp(_altitude, _parameters[0], _parameters[1])
    
    return _pressure / (R_AIR * _temperature)

def propulsion(delta_p,V,rho):
    F = delta_p*Fmaxi*(rho/rhoi)**(nrho)*(V/Vi)**nv
    Mt =F*zt
    return F,Mt


def modeloAerodinamico(V,rho,alpha,beta,p,q,r,delta_e,delta_a,delta_r):
    Vinf=224.6
    CL0=0.382
    CD0=0.0252
    Cm0=0.0622
    CDa=0.2013
    CLa=6.29
    Cma=-3.63
    CLad =4.04
    Cmad =-16.5
    CDq =0.281
    CLq =14.6
    Cmq =-45.5
    Cyb =-0.785
    Clb=-0.121
    Cnb =0.174
    CYbd =0.0214
    Clbd = 0.0035
    Cnbd =0.0098
    Cyp =-0.0794
    Clp=-0.0522
    Cnp=-0.0587
    Cyr =0.572
    Clr =0.254
    Cnr =-0.277
    
    CDde =0.0126
    CLde =0.3891
    Cmde =-1.5980
    Cyda =0.0094
    Clda =-0.1784
    Cnda =0.0080
    Cydr =-0.3030
    Cldr =-0.0464
    Cndr =-0.1594
    
    qdyn =0.5*rho*V**2
    
    L = qdyn*Sw*(CL0+CLa*alpha+CLq*q*Cbarw/(2*Vinf)+CLde*delta_e)
    D = qdyn*Sw*(CD0+CDa*alpha+CDq*q*Cbarw/(2*Vinf)+CDde*delta_e)
    Y = qdyn*Sw*(Cyb*beta+(Cyp*p+Cyr*r)*bw/(2*Vinf))+Cyda*delta_a+Cydr*delta_r
    
    La = qdyn*Sw*bw*(Clb*beta+(Clp*p+Clr*r)*bw/(2*Vinf)+Clda*delta_a+Cldr*delta_r)
    Ma = qdyn*Sw*Cbarw*(Cm0+Cma*alpha+(Cmq*q)*Cbarw/(2*Vinf)+Cmde*delta_e)
    Na = qdyn*Sw*bw*(Cnb*beta+(Cnp*p+Cnr*r)*bw/(2*Vinf)+Cnda*delta_a+Cndr*delta_r)
    
    return L,D,Y,La,Ma,Na,CLad,Cmad,CYbd,Clbd,Cnbd


def dynamics(t,X,U):
    # States
    V = X[0]
    alpha = X[1]
    beta = X[2]
    phi = X[3]
    theta = X[4]
    psi = X[5]
    p = X[6]
    q = X[7]
    r = X[8]
    x0 = X[9]
    y0 = X[10]
    H = X[10]
    
    #Controls
    delta_p = U[0]
    delta_e = U[1]
    delta_a = U[2]
    delta_r = U[3]
    
    # Atmospheric model
    rho = ISA(H)
    
    # Propulsive model
    F,Mf = propulsion(delta_p,V,rho)
    
    # Aerodynamic Model   
    L,D,Y,La,Ma,Na,CLad,Cmad,CYbd,Clbd,Cnbd = modeloAerodinamico(V,rho,alpha,beta,p,q,r,delta_e,delta_a,delta_r)

    #Kinematics
    u = V*cos(alpha)*cos(beta)
    v = V*sin(beta)
    w = V*cos(beta)*sin(alpha)
    
    #Linear kinematics
    x0_d = cos(theta)*cos(psi)*u + (cos(psi)*sin(theta)*sin(phi)-cos(phi)*sin(psi))*v + (cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))*w
    y0_d = cos(theta)*sin(psi)*u + (cos(phi)*cos(psi)+sin(theta)*sin(phi)*sin(psi))*v + (-cos(psi)*sin(phi)+cos(phi)*sin(theta)*sin(psi))*w
    H_d = u*sin(theta)- v*sin(phi)*cos(theta)-w*cos(phi)*cos(theta)
    
    # Angular kinematics
    phi_d = p + (cos(phi)*r+q*sin(phi)*tan(theta))
    theta_d = cos(phi)*q-r*sin(phi)
    psi_d = 1/cos(theta)*(cos(phi)*r+q*sin(phi))
    
    # Linear  Dynamics
    V_d = 1/m*(-D+F*cos(alpha+alphaf)*cos(beta)+m*g*(cos(beta)*cos(theta)*cos(phi)*sin(alpha)-cos(alpha)*cos(beta)*sin(theta)+cos(theta)*sin(beta)*sin(phi)))
    alpha_d = (1/(m*V*cos(beta)))*(-L-F*sin(alpha+alphaf)+m*g*(cos(alpha)*cos(theta)*cos(phi)+sin(alpha)*sin(theta)))+q-(p*cos(alpha)+r*sin(alpha))*tan(beta)
    beta_d = (1/(m*V))*(-Y-F*cos(alpha+alphaf)*sin(beta)+m*g*(-cos(theta)*cos(phi)*sin(alpha)*sin(beta)+cos(alpha)*sin(beta)*sin(theta)+cos(beta)*cos(theta)*sin(phi)))-r*cos(alpha)+p*sin(alpha)
    
    
    qinf =Vinf**2*0.5*rhoi
    
    #Correção aerodynâmica não-estacionária
    alpha_d = (1+(qinf*Sw*(Cbarw/(2*Vinf))*CLad)/(m*V*cos(beta)))**(-1)*alpha_d
    beta_d = (1+(qinf*Sw*(bw/(2*Vinf))*CYbd)/(m*V))**(-1)*beta_d
    
    La = La + qinf*Sw*bw*(bw/(2*Vinf))*Clbd*beta_d
    Ma = Ma + qinf*Sw*Cbarw*(Cbarw/(2*Vinf))*Cmad*alpha_d
    Na = Na + qinf*Sw*bw*(bw/(2*Vinf))*Cnbd*beta_d
    
    Ma = Ma+Mf
    
    p_d = (Ixx*Izz-Ixz**2)**(-1)*(Izz*La+Ixz*Na+Ixx*Ixz*p*q-Ixz*Iyy*p*q+Ixz*Izz*p*q+Ixz**2*q*r+Iyy*Izz*q*r-Izz**2*r*q)
    q_d = Iyy**(-1)*(Ma-Ixz*p**2-Ixx*p*r-Izz*p*r+Ixz*r**2)
    r_d = (Ixx*Izz-Ixz**2)**(-1)*(Ixz*La+Ixx*Na+Ixx**2*p*q-Ixz**2*p*q-Ixx*Iyy*p*q-Ixx*Ixz*q*r+Ixz*Iyy*q*r-Ixz*Izz*q*r)
    
    #Output
    Xdot = [V_d, alpha_d, beta_d, phi_d, theta_d, psi_d, p_d, q_d, r_d, x0_d, y0_d, H_d]
    return Xdot

def obj(z):
    Ve =VE
    alphae =z[0]
    beatae = z[1]
    phie = 0.
    thetae =z[2]
    psie = 0.
    pe =z[3]
    qe =z[4]
    re =z[5]
    x0e = 0.
    y0e = 0.
    He = HE
    
    deltape =z[6]
    deltaee =z[7]
    deltaae =z[8]
    deltare =z[9]
    
    Xe = [Ve, alphae,beatae, phie,thetae, psie,pe,qe,re,x0e,y0e,He]
    Ue = [deltape, deltaee, deltaae, deltare]
    
    sol = dynamics(0,Xe,Ue)
    
    equilvar = [sol[0],sol[1],sol[2],sol[3],sol[4],sol[5],sol[6],sol[7],sol[8],sol[11]]
    return equilvar


#Xg = [200,2*pi/180,0,0,2*pi/180,0,0,0,0,0,0,5000]
#Ug = [0.3,2*pi/180,0,0]


VE = 200.
HE = 10000.
Zg = [pi/180., pi/180, pi/180, pi/180,pi/180,pi/180, 0.3, pi/180.,pi/180.,pi/180.]


root = fsolve(obj, Zg)

Xg = [VE,root[0],root[1],0.,root[2],0.,root[3],root[4],root[5],0.,0,HE]
Ug = [root[6],root[7],root[8],root[9]]
sol =dynamics(0,Xg,Ug)
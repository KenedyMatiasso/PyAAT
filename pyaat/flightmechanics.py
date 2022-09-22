"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains methods and functions for flight mechanics analysis
"""
from numpy import zeros, array, dot, arange, sqrt
from numpy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import StrMethodFormatter

def modesMatrix(A,B):
    Ai = zeros((10,10))
    Af = zeros((10,10))
    Bf = zeros((10, np.shape(B)[1]))
    
    Ai[:,0] = A[:,1]
    Ai[:,1] = A[:,3]
    Ai[:,2] = A[:,5]
    Ai[:,3] = A[:,8]
    Ai[:,4] = A[:,0]
    Ai[:,5] = A[:,2]
    Ai[:,6] = A[:,4]
    Ai[:,7] = A[:,6]
    Ai[:,8] = A[:,7]
    Ai[:,9] = A[:,9]
    
    Af[0] = Ai[1]
    Af[1] = Ai[3]
    Af[2] = Ai[5]
    Af[3] = Ai[8]
    Af[4] = Ai[0]
    Af[5] = Ai[2]
    Af[6] = Ai[4]
    Af[7] = Ai[6]
    Af[8] = Ai[7]
    Af[9] = Ai[9]

    Bf[0] = B[1]
    Bf[1] = B[3]
    Bf[2] = B[5]
    Bf[3] = B[8]
    Bf[4] = B[0]
    Bf[5] = B[2]
    Bf[6] = B[4]
    Bf[7] = B[6]
    Bf[8] = B[7]
    Bf[9] = B[9]
    
    return Af, Bf

def cut_A_andB(A,B):
    Anew = zeros((10,10))
    Bnew = zeros((10, np.shape(B)[1]))
    
    Anew = A[2:12,2:12]
    Bnew = B[2:12,:]
    return Anew, Bnew
    
def lateroMatrix(A,B):
    A, B = cut_A_andB(A, B)

    Al = zeros((5,5))
    Bl = zeros((5,np.shape(B)[1]))
    Am, Bm =  modesMatrix(A,B)
    Al = Am[5:,5:]
    Bl = Bm[5:,:]
    return Al, Bl

def longMatrix(A,B):
    A, B = cut_A_andB(A,B)
    
    Al = zeros((5,5))
    Bl = zeros((5, np.shape(B)[1]))
    Am, Bm =  modesMatrix(A,B)
    Al = Am[:5,:5]
    Bl = Bm[:5,:]
    return Al, Bl

def shortPeriodMatrix(A,B):
    A, B = cut_A_andB(A,B)

    Asp = array([[A[3,3], A[3,8]],
                 [A[8,3], A[8,8]]])
    
    Bsp = array([B[3],
                 B[8]])
    return Asp, Bsp

def phugoidMatrix(A,B):
    A, B = cut_A_andB(A,B)

    Asp = array([[A[3,3], A[3,8]],
                 [A[8,3], A[8,8]]])
    
    Bsp = array([B[0],B[8]])
       
    Aspph = array([[A[3,1], A[3,5], A[3,0]],
                  [A[8,1], A[8,5], A[8,0]]])
    
    Aphph = array([[A[1,1], A[1,5], A[1,0]],
                [A[5,1], A[5,5], A[5,0]],
                [A[0,1], A[0,5], A[0,0]]])
    
    Bphph = array([B[1], B[5], B[0]])
    
    Aphsp =array([[A[1,3], A[1,8]],
                  [A[5,3], A[5,8]],
                  [A[0,3], A[0,8]]])
        
    Aph = Aphph - Aphsp.dot(inv(Asp).dot(Aspph))
    Bph = Bphph - Aphsp.dot(inv(Asp).dot(Bsp))
    return Aph, Bph

def rollMatrix(A,B):
    A, B = cut_A_andB(A,B)

    Ar = array([[A[7,7]]])
    Br = array([[B[7,2]]])
    return Ar, Br

def dutchRollMatrix(A,B):
    A, B = cut_A_andB(A,B)

    Adr = array([[A[2,2], A[2,9]],
                 [A[9,2], A[9,9]]])
    
    Bdr = array([B[2], B[7]])
    
    return Adr, Bdr

def spiralMatrix(A,B):
    A, B = cut_A_andB(A,B)

    Abprbpr = array([[A[2,2], A[2,7], A[2,9]],
                  [A[7,2], A[7,7], A[7,9]],
                  [A[9,2], A[9,7], A[9,9]]])
    
    Abprf = array([[A[2,4]],
                   [A[7,4]],
                   [A[9,4]]])
    
    BbprDr = array([B[2], B[7], B[9]])
    
    Afbrp = array([A[4,2], A[4,7], A[4,9]])
    
    As = array([A[5,5] - Afbrp.dot(inv(Abprbpr).dot(Abprf))])
    Bs = array([B[5] - Afbrp.dot(inv(Abprbpr).dot(BbprDr))])
    
    return As, Bs
    

def shortPeriodQuality(omega, zeta, aircraft = None, phase = None):
    FL =array([[0.32805384412649785, 5.24932301009952],
    [0.29722785662292733, 4.9395072413178855],
    [0.2733177598440804, 4.6578416572352905],
    [0.25381554044619625, 4.4513188033157665],
    [0.23805698905324912, 4.188394818710652],
    [0.226609247013785, 3.9536210188045753],
    [0.21786027204249464, 3.6906340583947648],
    [0.21363029478006532, 3.446388697462435],
    [0.21579076870505184, 3.20206776556447],
    [0.22234070444176265, 2.8731325424771805],
    [0.23134790765476965, 2.6569240097954365],
    [0.2455319932135751, 2.4218731163486984],
    [0.2684347334461723, 2.1773506618757077],
    [0.31140897879115953, 1.9608650356533026],
    [0.36664031703455435, 1.8289055344936473],
    [0.4648937652833816, 1.7249450761017613],
    [0.5953045992926286, 1.696127347872947],
    [0.7363337328241915, 1.7237737261344188],
    [0.928905381889186, 1.845329624358322],
    [1.0720439465671996, 1.9671122354791288],
    [1.2130237787581462, 2.0701532471224553],
    [1.293493466178002, 2.145157430515196]])
    
    SL = array([[0.5969664582403355, 4.74981152241309],
    [0.5434699631115381, 4.656090929864746],
    [0.4972234266749034, 4.552961752094846],
    [0.44601452187083274, 4.4122989947262266],
    [0.398111427657013, 4.25285685239731],
    [0.36968754790271247, 4.112105928902116],
    [0.3416178259854679, 3.8961996800829124],
    [0.32357737549246396, 3.6990224355802503],
    [0.30953371769328025, 3.5112159908163285],
    [0.30203864987340345, 3.2106072846812257],
    [0.31119408952254307, 2.9380480019576485],
    [0.327016849905606, 2.731210269014645],
    [0.3453499055464702, 2.514963950850083],
    [0.3889333514174703, 2.2703659254114568],
    [0.4401503928770601, 2.147903175599934],
    [0.515655629423155, 2.0253522596618376],
    [0.6378453596858475, 1.9872267074989791],
    [0.7811780277075797, 2.033690256203566],
    [0.9661769051763085, 2.1552965350712254],
    [1.0985883055842474, 2.3428888620991817],
    [1.249088220958795, 2.596253119551468],
    [1.3450455469253035, 2.82156795359223],
    [1.4199323665036854, 3.1221010887616982],
    [1.4989686894555954, 3.4414262040524033],
    [1.5590928119512668, 3.7701850948865445],
    [1.5900690969170355, 3.9204705552126877]])
    
    TL = array([
    [0.6093526385727405, 3.678618274858799],
    [0.5519951276721128, 3.62249423771387],
    [0.5176548882024024, 3.54749005432113],
    [0.48786067735233457, 3.4630772857068317],
    [0.4552686704456923, 3.3411057471719365],
    [0.43334683650072225, 3.2002918478720472],
    [0.42071767977690383, 3.0688235579890213],
    [0.41869243316474186, 2.8997083320588226],
    [0.42498938293243493, 2.768126685727344],
    [0.43782231808765465, 2.655299234034284],
    [0.4532906690243831, 2.50487522693781],
    [0.481040695544303, 2.3919722042791163],
    [0.5284480799922592, 2.3353569558575593],
    [0.5833936791727925, 2.3069170824569207],
    [0.6252119099687357, 2.3067407502037724],
    [0.690186498624594, 2.3440728072274633],
    [0.7619033990712971, 2.4001968443723927],
    [0.8245833580609042, 2.484559232342934],
    [0.896819880238679, 2.606492985395011],
    [0.9850722397656022, 2.7377975381858297],
    [1.0503503375831607, 2.9067616221847565],
    [1.087214044940544, 3.0758012771493197],
    [1.097903089585066, 3.22611192779734],
    [1.0924004708311115, 3.338876403685704],
    [1.055192805617888, 3.3953405101759904],
    [0.9894101514375416, 3.5176521180562426],
    [0.9231707013678775, 3.611788350915578],
    [0.857145480128069, 3.6589572286327594],
    [0.7802659975739912, 3.6873845068724598],
    [0.706787811874387, 3.697032400151862],
    [0.653022964175021, 3.6972339227268884]])
    
    plt.figure(figsize=(8,4))
    plt.plot(FL[:,0],FL[:,1], linestyle ='-.', color = 'b')
    plt.text(0.2,1.2, 'Unacceptable', fontsize=9)
    plt.plot(SL[:,0],SL[:,1], linestyle ='-.', color = 'b')
    plt.text(0.25,4, 'Poor', fontsize=9)
    plt.text(0.4,4, 'Acceptable', fontsize=9)
    plt.plot(TL[:,0],TL[:,1], linestyle ='-.', color = 'b')
    plt.text(0.45,3, 'Satisfactory', fontsize=9)
    
    plt.scatter(zeta, omega, color = 'r')
    plt.title('Short Period Quality, Aircraft Class I, II, III or IV, Phase A, B or C')
    plt.xscale('log')
    plt.xlim([0.1, 4])
    plt.ylim([0.5, 6])
    plt.xlabel(r'Damping ration $\zeta_{sp}$')
    plt.ylabel(r'Undamped natural frequency $\omega_{sp}$ [rad/s]')
    plt.grid(axis='both', which='minor')
    plt.grid()
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    
def dutchRollQuality(omega, zeta, aircraft = None, phase = None):
    if phase =='A':
        if aircraft == 'I' or aircraft == 'IV':
            L2_x = arange(0.19, 0.35/1.0, 0.02)
            L1_x = [0.19, 0.19]
            L1_y = [20, 0.35/0.19]
            L2_y = 0.35/L2_x
            L3_x =[0.35/1.0, 50]
            L3_y = [1.0, 1.0]
            title = 'Dutch Roll Quality, Aircraft class I or IV, Flight phase A'
        else:
            L2_x = arange(0.19, 0.35/0.5, 0.02)
            L1_x = [0.19, 0.19]
            L1_y = [20, 0.35/0.19]
            L2_y = 0.35/L2_x
            L3_x =[0.35/0.5, 50]
            L3_y = [0.5, 0.5]

            title = 'Dutch Roll Quality, Aircraft class II or III, Flight phase A'
    elif phase == 'B':
        L2_x = arange(0.08, 0.15/0.5, 0.02)
        L1_x = [0.08, 0.08]
        L1_y = [20, 0.15/0.08]
        L2_y = 0.15/L2_x
        L3_x =[0.15/0.5, 50]
        L3_y = [0.5, 0.5]
        title = 'Dutch Roll Quality, Aircraft class I, II, III or IV, Flight phase B'
    else:
        if aircraft == 'I' or aircraft == 'IV':
            L2_x = arange(0.08, 0.15/1.0, 0.02)
            L1_x = [0.08, 0.08]
            L1_y = [20, 0.15/0.08]
            L2_y = 0.15/L2_x
            L3_x =[0.15/1.0, 50]
            L3_y = [1.0, 1.0]
            title = 'Dutch Roll Quality, Aircraft class I or IV, Flight phase C'
        else:
            L2_x = arange(0.08, 0.1/0.5, 0.02)
            L1_x = [0.08, 0.08]
            L1_y = [20, 0.1/0.08]
            L2_y = 0.1/L2_x
            L3_x =[0.1/0.5, 50]
            L3_y = [0.5, 0.5]

            title = 'Dutch Roll Quality, Aircraft class II or III, Flight phase C'
    
    L2_x2 = arange(0.02, 0.05/0.5, 0.01)
    L1_x2 = [0.02, 0.02]
    L1_y2 = [20, 0.05/0.02]
    L2_y2 = 0.05/L2_x2
    L3_x2 =[0.1, 50]
    L3_y2 = [0.5, 0.5]
    
    L1_x3 = [0.0, 0.0]
    L1_y3 = [20, 0.4]
    L3_x3 =[0, 50]
    L3_y3 = [0.4, 0.4]
        
    plt.figure(figsize=(8,4))
    plt.plot(L2_x, L2_y, linestyle ='-.', color = 'b')
    plt.plot(L1_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L3_x, L3_y, linestyle ='-.', color = 'b')
    plt.text(0.5,2, 'level 1', fontsize=9)
    
    plt.plot(L2_x2, L2_y2, linestyle ='-.', color = 'b')
    plt.plot(L1_x2, L1_y2, linestyle ='-.', color = 'b')
    plt.plot(L3_x2, L3_y2, linestyle ='-.', color = 'b')
    plt.text(0.05,2.5, 'level 2', fontsize=9)
    
    
    plt.plot(L1_x3, L1_y3, linestyle ='-.', color = 'b')
    plt.plot(L3_x3, L3_y3, linestyle ='-.', color = 'b')
    plt.text(0.02, 1, 'level 3', fontsize=9)
    
    plt.scatter(zeta, omega, color = 'r')
    
    plt.xscale('log')
    plt.xlim([0.005, 1])
    plt.ylim([0, 4])
    plt.title(title)
    plt.xlabel(r'Damping ration $\zeta_{dr}$')
    plt.ylabel(r'Natural frequency $\omega_{dr}$ [rad/s]')
    
    plt.grid(axis='both', which='minor')
    plt.grid()
    
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    
def phugoidQuality(omega, zeta, aircraft = None, phase = None):
    L1_x = [0.04, 0.04]
    L1_y = [-0.5, 0.5]
    
    L2_x = [0.0, 0.0]
    L2_y = [-0.5, 0.5]
    
    L3_y2 = arange(0.1, 0.5, 0.01)

    L3_x2 = -1/(55*L3_y2)
    
    plt.figure(figsize=(8,4))
    plt.plot(L1_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L2_x, L2_y, linestyle ='-.', color = 'b')
    plt.plot(L3_x2, L3_y2, linestyle ='-.', color = 'b')
    plt.plot(L3_x2, -L3_y2, linestyle ='-.', color = 'b')

    
    plt.text(0.1, 0, 'level 1', fontsize=9)
    plt.text(0.01, 0, 'level 2', fontsize=9)
    plt.text(-0.04, 0, 'level 3', fontsize=9)
    plt.text(-0.12, -0.4, 'Unacceptable', fontsize=9)
    plt.text(-0.12, 0.4, 'Unacceptable', fontsize=9)
    
    plt.scatter(zeta, omega, color = 'r')
    
    #plt.xscale('log')
    plt.xlim([-0.15, 0.15])
    plt.ylim([-0.5, 0.5])
    plt.title('Phugoid Quality, Aircraft Class I, II, III or IV, Phase A, B or C ')
    plt.xlabel(r'Damping ration $\zeta_{ph}$')
    plt.ylabel(r'Natural frequency $\omega_{ph}$ [rad/s]')
    plt.grid()
    
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    
def spiralQuality(T, aircraft = None, phase = None):
    if phase == 'A' or  phase == 'C':
        L1_x = [17.3,  17.3]
        title = 'Spiral Quality, Aircraft Class I, II, III or IV, Phase A or C'
    else:
        L1_x = [28.9,  28.9]
        title = 'Spiral Quality, Aircraft Class I, II, III or IV, Phase B'
        
    L1_y = [-0.5, 0.5]
    L2_x = [11.5, 11.5]
    L3_x = [7.2, 7.2]
    L4_x = [0.0, 0.0]

    plt.figure(figsize=(8,3))
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.plot(L1_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L2_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L3_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L4_x, L1_y, linestyle ='-.', color = 'b')
    
    plt.text(40, 0.4, 'level 1', fontsize=9)
    plt.text(20, 0.2, 'level 2', fontsize=9)
    plt.text(7, 0, 'level 3', fontsize=9)
    plt.text(-3, -0.4, 'Unacceptable', fontsize=9)
    plt.text(-9, 0.4, 'level 1', fontsize=9)
   
    plt.scatter(T, 0.0, color = 'r')
    
    #plt.xscale('log')
    if T<=-10:
        plt.xlim([T-1, -10])
    elif T>=50:
        plt.xlim([-10, T+1])
    else:
        plt.xlim([-10, 30])
    
    plt.ylim([-0.5, 0.5])
    plt.title(title)
    plt.xlabel(r'Time constant $T_s$')
    #plt.ylabel(r'Natural frequency $\omega_{ph}$ [rad/s]')
    plt.grid()
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    
def rollQuality(T, aircraft = None, phase = None):
    if phase == 'B':
        title = 'Roll Quality, Aircraft Class I, II, III or IV, Phase B'
        L1_x = [1.4,  1.4]
        L2_x = [3.0, 3.0]
    else:
        if aircraft == 'I' or aircraft == 'IV':
            title = 'Roll Quality, Aircraft Class I or IV, Phase A or C'
            L1_x = [1.0,  1.0]
            L2_x = [1.4, 1.4]
        else:
            title = 'Roll Quality, Aircraft Class II or III, Phase A or C'
            L1_x = [1.4,  1.4]
            L2_x = [3.0, 3.0]
    
    L1_y = [-0.5, 0.5]
    plt.figure(figsize=(8,3))
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.plot(L1_x, L1_y, linestyle ='-.', color = 'b')
    plt.plot(L2_x, L1_y, linestyle ='-.', color = 'b')
    
    plt.text(0, 0, 'level 1', fontsize=9)
    plt.text(2, 0, 'level 2', fontsize=9)
    plt.text(3.5, 0, 'level 3', fontsize=9)
    plt.scatter(T, 0.0, color = 'r')
    #plt.xscale('log')
    plt.xlim([-5, 5])
    plt.ylim([-0.5, 0.5])
    plt.title(title)
    plt.xlabel(r'Time constant $T_r$')
    #plt.ylabel(r'Natural frequency $\omega_d$ [rad/s]')
    plt.grid()
    
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    
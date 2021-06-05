"""
Python Aerospace Analysis Toolbox - PyAAT
Copyright (c) 2021 Kenedy Matiasso Portella
Distributed under MIT License

This file contains methods and functions for flight mechanics analysis
"""
from numpy import zeros, array, dot
from numpy.linalg import inv

def modesMatrix(A,B):
    Ai = zeros((10,10))
    Af = zeros((10,10))
    Bf = zeros((10,4))
    
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
def lateroMatrix(A,B):
    Al = zeros((5,5))
    Bl = zeros((5,2))
    Am, Bm =  modesMatrix(A,B)
    Al = Am[5:,5:]
    Bl = Bm[5:,2:]
    return Al, Bl

def longMatrix(A,B):
    Al = zeros((5,5))
    Bl = zeros((5,2))
    Am, Bm =  modesMatrix(A,B)
    Al = Am[:5,:5]
    Bl = Bm[:5,:2]
    return Al, Bl

def shortPeriodMatrix(A,B):
    Asp = array([[A[3,3], A[3,8]],
                 [A[8,3], A[8,8]]])
    
    Bsp = array([[B[3,0], B[3,1]],
                 [B[8,0], B[8,1]]])
    return Asp, Bsp

def phugoidMatrix(A,B):
    Asp = array([[A[3,3], A[3,8]],
                 [A[8,3], A[8,8]]])
    
    Bsp = array([[B[3,0], B[3,1]],
                 [B[8,0], B[8,1]]])
    
    Aspph = array([[A[3,1], A[3,5], A[3,0]],
                  [A[8,1], A[8,5], A[8,0]]])
    
    Aphph = array([[A[1,1], A[1,5], A[1,0]],
                [A[5,1], A[5,5], A[5,0]],
                [A[0,1], A[0,5], A[0,0]]])
    
    Bphph = array([[B[1,0], B[1,1]],
                  [B[5,0], B[5,1]],
                  [B[0,0], B[0,1]]])
    
    Aphsp =array([[A[1,3], A[1,8]],
                  [A[5,3], A[5,8]],
                  [A[0,3], A[0,8]]])
    
    Aph = Aphph - Aphsp.dot(inv(Asp).dot(Aspph))
    Bph = Bphph - Aphsp.dot(inv(Asp)).dot(Bsp)
    return Aph, Bph
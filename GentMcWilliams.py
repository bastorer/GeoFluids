## Linear Stability of a Barotropic QG Vortex

# This is an attempted re-write of the code
#    qg_BTvortex_stab_SpecFD_loop.
# provided by Francis.

import timeit
import scipy
import time
import sys

import scipy.sparse as sp
import scipy.linalg as spalg
import numpy as np
import numpy.linalg as nlg
import matplotlib.pyplot as plt

from scipy.sparse.linalg import eigs
from scipy.interpolate import interp1d
from scipy.misc import factorial
from cheb import cheb
from FiniteDiff import FiniteDiff

class Parameters:
    ## Class to hold parameter values
    H        = 2.4e3
    L        = 200e3
    f0       = 8e-5
    g        = 9.81
    N        = np.sqrt(5)*1e-3
    Lr       = 6.25
    Nr       = 201
    N2       = 100
    Nt       = 40
    kts      = np.arange(1,2,1)
    kzs      = np.arange(0.4,1,0.1)
    nmodes   = 1
    printout = False

    def display(self):
        print 'H = {0}'.format(self.H)
        print 'L = {0}'.format(self.L)
        print 'f0 = {0}'.format(self.f0)
        print 'g = {0}'.format(self.g)
        print 'N = {0}'.format(self.N)
        print 'Lr = {0}'.format(self.Lr)
        print 'Nr = {0}'.format(self.Nr)
        print 'N2 = {0}'.format(self.N2)
        print 'Nt = {0}'.format(self.Nt)
        print 'kts = {0}'.format(self.kts)
        print 'kzs = {0}'.format(self.kzs)
        print 'nmodes = {0}'.format(self.nmodes)

class Geometry:
    ## Class to hold geometric values
    def __init__(self, method, params):
        self.method = method
        if method == 'cheb':
            Dr, r = cheb(params.Nr)
            self.r = r*params.Lr
            self.Dr = Dr/params.Lr
            self.Dr2 = np.dot(self.Dr,self.Dr)
        elif method == 'FD':
            self.r = np.arange(params.Lr, -params.Lr-2*params.Lr/(params.Nr), -2*params.Lr/(params.Nr))
            self.Dr = FiniteDiff(self.r, 8, True, True)
            self.Dr2 = np.dot(self.Dr, self.Dr)


def Build_Laplacian(params, geom):
    
    D1d = geom.Dr2[1:params.N2+1, 1:params.N2+1]
    D2d = geom.Dr2[np.arange(1,params.N2+1,1),:][:,np.arange(params.Nr-1,params.N2,-1)]

    E1d = geom.Dr[1:params.N2+1, 1:params.N2+1]
    E2d = geom.Dr[np.arange(1,params.N2+1,1),:][:,np.arange(params.Nr-1,params.N2,-1)]

    if sp.issparse(geom.Dr):
        R = sp.spdiags(np.transpose(1.0/geom.r[1:params.N2+1]), np.array([0]), params.N2, params.N2)
    else:
        R = np.diag(1.0/np.ravel(geom.r[1:params.N2+1]))
        
    Lap = D1d + D2d + np.dot(R, E1d + E2d)

    return Lap

def Print_npArray(fp, arr):
    for ii in xrange(0,arr.shape[0]):
        for jj in xrange(0,arr.shape[1]):
            if jj == arr.shape[1]-1:
                fp.write('{0:+2.2e}'.format(arr[ii,jj]))
            else:
                fp.write('{0:+2.2e}, '.format(arr[ii,jj]))
        fp.write('\n')

            
def QG_Vortex_Stability():

    ## Initialize parameters
    paramsCheb = Parameters()
    paramsFD   = Parameters()
    paramsFD.Nr = 1001
    paramsFD.N2 = 500

    ## Set-up the geometry

    GeomCheb = Geometry('cheb', paramsCheb)
    GeomFD   = Geometry('FD', paramsFD)

    GeomCheb.Lap = Build_Laplacian(paramsCheb, GeomCheb)
    GeomFD.Lap   = Build_Laplacian(paramsFD, GeomFD)

    ## Set up the profiles
    rin    = GeomCheb.r[1:paramsCheb.N2+1]
    Prsp   = np.ravel(-0.5*np.exp(-rin**2))            # 1/r*Psi_r
    Qrsp   = np.ravel(-2*np.exp(-rin**2)*(rin**2-2))   # 1/r*Q_r
    
    rin    = GeomFD.r[1:paramsFD.N2+1] 
    Prfd   = np.ravel(-0.5*np.exp(-rin**2))            # 1/r*Psi_r
    Qrfd   = np.ravel(-2*np.exp(-rin**2)*(rin**2-2))   # 1/r*Q_r

    kts    = paramsCheb.kts       
    kzs    = paramsCheb.kzs
    nmodes = paramsCheb.nmodes
 
    growthsp = np.zeros([kzs.shape[0], kts.shape[0], nmodes])
    frequysp = np.zeros([kzs.shape[0], kts.shape[0], nmodes])
    growthfd = np.zeros([kzs.shape[0], kts.shape[0], nmodes])
    frequyfd = np.zeros([kzs.shape[0], kts.shape[0], nmodes])

    ## Start solving

    for cntz in xrange(0, kzs.shape[0]):

        kz  = kzs[cntz]
        kz2 = kz**2
  
        for cntt in xrange(0, kts.shape[0]):

            kt  = kts[cntt]
            kt2 = kt**2
    
            # Build A and B for eigen-analysis

            R2invC = np.diag(np.ravel(1/GeomCheb.r[1:paramsCheb.N2+1]**2))
            Bcheb = GeomCheb.Lap - kt2*R2invC - kz2*np.eye(paramsCheb.N2,paramsCheb.N2)
            Acheb = np.dot(np.diag(Prsp),Bcheb) - np.diag(Qrsp)

            R2invF = np.diag(np.ravel(1./GeomFD.r[1:paramsFD.N2+1]**2))
            Bfd = GeomFD.Lap - kt2*R2invF - kz2*np.eye(paramsFD.N2,paramsFD.N2)
            Afd = np.dot(np.diag(Prfd),Bfd) - np.diag(Qrfd)
            
            # Find eigen-space (Direct)
            
            t0 = timeit.timeit()
            eigValCheb, eigVecCheb = spalg.eig(Acheb,Bcheb)
            t1 = timeit.timeit()
            timesp = t1 - t0
            
            ind = (-eigValCheb.imag).argsort()
            eigVecCheb = eigVecCheb[:,ind]
            eigValCheb = eigValCheb[ind]

            omegaCheb = eigValCheb*kt
            growthsp[cntz,cntt,:] = omegaCheb[0:nmodes].imag;
            frequysp[cntz,cntt,:] = omegaCheb[0:nmodes].real;
            
            # Loop over modes 
            for ii in xrange(0,nmodes):
            
                grow = omegaCheb[ii].imag
                freq = omegaCheb[ii].real
      
                # Find Eigenvalues (Indirect)
                sig0 = eigValCheb[ii]

                X = np.hstack([np.array([paramsCheb.Lr]),\
                              np.ravel(GeomCheb.r[1:paramsCheb.N2+1]),\
                              np.array([0])])[::-1]
                Y = np.hstack([np.array([0]), eigVecCheb[:,ii], np.array([0])])[::-1]
                
                Xnew = np.ravel(GeomFD.r[1:paramsFD.N2+1])[::-1]
                
                interp_fcn = interp1d(X, Y, kind='cubic')
                chebvec = interp_fcn(Xnew)


                tmp = chebvec
                tmp[tmp==0] = 1
                tmp = tmp.conj()
                T = np.diag(np.ravel(tmp))
                Tinv = nlg.inv(T)
                
                t0 = timeit.timeit()
                try:
                    sig1, vec1 = eigs(np.dot(Afd, Tinv), 1, np.dot(Bfd,Tinv),\
                                      sigma=sig0,v0=np.dot(T,chebvec), maxiter=500)
                    vec1 = np.dot(Tinv, vec1)
                    plt.plot(Xnew[::-1], vec1.real, '-b', Xnew[::-1], vec1.imag, '-r')
                    plt.show()
                except:
                    sig1 = [np.nan+1j*np.nan];
                
                t1 = timeit.timeit()
                timefd = t1 - t0
                omegafd = kt*sig1[0]
                
                growfd = omegafd.imag
                freqfd = omegafd.real

                growthfd[cntz,cntt,ii] = growfd;
                frequyfd[cntz,cntt,ii] = freqfd;
        
                # Display the results                
                if paramsCheb.printout:
                    print '----------'
                    print 'kz = {0:4f}, kt = {1:2d}'.format(kz, kt)
                    print 'eig : growth rate = {0:4e}, frequency = {1:4e}, cputime = {2:4e}'\
                      .format(grow, freq, timesp)
                    print 'eigs: growth rate = {0:4e}, frequency = {1:4e}, cputime = {2:4e}'\
                      .format(growfd, freqfd, timefd)

    for ii in xrange(0, kts.shape[0]):
        plt.plot(kzs, 4*np.ravel(growthfd[:,ii,0]), '-o', kzs, 4*np.ravel(growthsp[:,ii,0]), '-*')
        plt.show()

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()

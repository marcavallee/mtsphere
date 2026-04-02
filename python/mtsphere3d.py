
import numpy as np
from scipy.constants import pi, mu_0, epsilon_0
from layeredearthspherecorrection import MTLayeredEarthCorrection,MTSphereReflectionCorrection
from layeredearthfunctions import PlaneWaveSphericalHarmonics
import sphericalfunctions as sphfun
import halfspacesphere as hssphere

def mtsphere3d(nw, nf, nlyr, nterms, nx, ny, freq, pxy, thk, res, depth,
               radius, sres, Zhat, E, htarg):

    dpthl = np.zeros(nlyr+1)
    for jz in range(1,nlyr):
        dpthl[jz+1] = dpthl[jz] + thk[jz]
    zs = depth
    sxlyr = 0
    for jz in range(nlyr,0,-1):
        if zs > dpthl[jz]:
            sxlyr = jz
            break
    impedance = np.zeros((nf,nx,ny,3,2),dtype=complex)
    Es = np.zeros((nf,nx, ny, 3, 2), dtype=complex)
    Hs = np.zeros((nf,nx, ny, 3, 2), dtype=complex)
    Et = np.zeros((nf,nx, ny, 3, 2), dtype=complex)
    Ht = np.zeros((nf,nx, ny, 3, 2), dtype=complex)
    for jf in range(nf):
        nw.write('\nFrequency: {:15.7g}\n\n'.format(freq[jf]))
        w = 2. * pi * freq[jf]
        yhatl = np.zeros(nlyr+1,dtype=complex)
        zhatl = np.zeros(nlyr+1,dtype=complex)
        for jz in range(nlyr+1):
            if jz == 0:
                yhatl[jz] = 1j * w * epsilon_0
            else:
                yhatl[jz] = 1 / res[jz] + 1j * w * epsilon_0
            zhatl[jz] = 1j * w * mu_0
        yhats = 1 / sres + 1j * w * epsilon_0
        zhats = 1j * w * mu_0
        rtm, rte = sphfun.ReflectionCoefficients(nterms, radius, yhatl[sxlyr], zhatl[sxlyr], yhats, zhats)
        AI = MTLayeredEarthCorrection(nw, nterms, nlyr, thk, zs, radius, sxlyr, yhatl, \
                                    zhatl, yhats, zhats, rtm, rte, htarg)
        nw.write('\n Reflection coefficients:\n')
        for n in range(1,nterms+1):
            nw.write('{:5}{:15.7g}{:15.7g}{:15.7g}{:15.7g}\n'.\
                format(n,rtm[n].real,rtm[n].imag,rte[n].real,rte[n].imag))
        psiai, psifi = PlaneWaveSphericalHarmonics(nterms, yhatl[sxlyr], zhatl[sxlyr], zs - dpthl[sxlyr], E[jf, sxlyr, :] )
        for id in range(2):
            nw.write('\n Plane wave spherical harmonics:\n')
            if id == 0:
                nw.write(' X directed propagation\n')
            else:
                nw.write(' Y directed propagation\n')
            for n in range(1,nterms+1):
                for m in [-1,1]:
                    nw.write('{:5}{:5}{:15.7g}{:15.7g}{:15.7g}{:15.7g}\n'.\
                             format(n,m,psiai[id,n,m].real,psiai[id,n,m].imag,psifi[id,n,m].real,psifi[id,n,m].imag))
            psiac, psifc = MTSphereReflectionCorrection(nterms, AI, psiai[id,:,:], psifi[id,:,:])
            nw.write('\nCorrected potentials\n')
            for n in range(1,nterms+1):
                for m in [-1,1]:
                    nw.write('{:5}{:5}{:15.7g}{:15.7g}{:15.7g}{:15.7g}\n'.format(n,m,\
                        psiac[n,m].real,psiac[n,m].imag,psifc[n,m].real,psifc[n,m].imag))
            psiar = np.zeros((nterms+1,3),dtype=complex)
            psifr = np.zeros((nterms+1,3),dtype=complex)
            for n in range(nterms+1):
                psiar[n,:] = rtm[n] * psiac[n,:]
                psifr[n,:] = rte[n] * psifc[n,:]
            Es[jf,:,:,:,id], Hs[jf,:,:,:,id] = hssphere.MTSphere_SingleSource(nterms, psiar, psifr,\
                                    depth, nx, ny, pxy, nlyr, thk, yhatl, zhatl, htarg)
            nw.write('\nInduced fields\n')
            for jx in range(nx):
                for jy in range(ny):
                    nw.write('{:15.7g}{:5}{:5}'.format(freq[jf],jx+1,jy+1))
                    for ic in range(3):
                        nw.write('{:15.7g}{:15.7g}'.format(Es[jf,jx,jy,ic,id].real,Es[jf,jx,jy,ic,id].imag))
                    for ic in range(3):
                        nw.write('{:15.7g}{:15.7g}'.format(Hs[jf,jx,jy,ic,id].real,Hs[jf,jx,jy,ic,id].imag))
                    nw.write('\n')
        impedance[jf,:,:,:,:], Et[jf,:,:,:,:], Ht[jf,:,:,:,:] = getimpedance(nx, ny, Zhat[jf,1], E[jf, 1, :],
                                             Es[jf,:,:,:,:], Hs[jf,:,:,:,:])

    return impedance, Es, Hs, Et, Ht

def getimpedance(nx, ny, Zhat, E, Es, Hs):

    impedance = np.zeros((nx,ny,3,2), dtype=complex)
    Et = np.zeros((nx,ny,3,2),dtype=complex)
    Ht = np.zeros((nx,ny,3,2),dtype=complex)
    Ep = E[0] + E[1]
    Hp = Ep / Zhat
    Et = Es
    Ht = Hs
    Et[:,:,0,0] = Et[:,:,0,0] + Ep
    Et[:,:,1,1] = Et[:,:,1,1] + Ep
    Ht[:,:,0,1] = Ht[:,:,0,1] - Hp
    Ht[:,:,1,0] = Ht[:,:,1,0] + Hp
    for ix in range(nx):
        for iy in range(ny):
            A = np.zeros(2, dtype=complex)
            B = np.zeros(2, dtype=complex)
            for id in range(2): # x and y directions
                A[id] = Ht[ix,iy,0,id]
                B[id] = Ht[ix,iy,1,id]
            den = A[0] * B[1] - B[0] * A[1]
            for ic in range(3):
                C = np.zeros(2, dtype=complex)
                if ic < 2:
                    for id in range(2):
                        C[id] = Et[ix,iy,ic,id]
                else:
                    for id in range(2):
                        C[id] = Ht[ix,iy,2,id]
                if abs(den) > 0:
                    impedance[ix,iy,ic,0] = ( C[0] * B[1] - B[0] * C[1] ) / den
                    impedance[ix,iy,ic,1] = ( A[0] * C[1] - C[0] * A[1] ) / den

    return impedance, Et, Ht



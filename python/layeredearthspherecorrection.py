
import numpy as np
import scipy
from sphericalfunctions import (SphericalGeometryInitialisation,
                                RadialFields,
                                MTSphericalHarmonicAnalysis)
from hssphere_mt import MTSphere_MultipleSources

def MTLayeredEarthCorrection(write0, nterms, nlyr, thk, depth, radius, sxlyr, yhatl, zhatl, yhats, zhats, rtm, rte, htarg):

    nlat = 4 * nterms + 1
    nlon = 4 * nterms
    center = np.zeros(3)
    center[2] = depth
    ps = SphericalGeometryInitialisation(center, radius, nlat, nlon)
    E, H = MTSphere_MultipleSources(nterms,nlat,nlon,depth,depth,nlyr,thk,radius,yhatl,zhatl,rtm,rte,htarg)
    AI = np.zeros((3,2 * nterms, 2 * nterms), dtype=complex)
    write0(' Spherical Harmonic Analysis\n',1)
    for m in [-1,1]:
        if m < 0:
            write0('\nSphere retroaction, negative order\n',1)
        else:
            write0('\nSphere retroaction, positive order\n',1)
        icr = -1
        C = np.zeros((2*nterms, 2*nterms), dtype=complex)
        for ic in range(2):
            for n in range(1, nterms + 1):
                icr += 1
                Er, Hr = RadialFields(nlat, nlon, ps, center, E[ic,n,m,:,:,:], H[ic,n,m,:,:,:])
                psia, psif = MTSphericalHarmonicAnalysis(nterms, nlat, nlon, radius, \
                                yhatl[sxlyr], zhatl[sxlyr], Er, Hr, False)
                ihr = -1
                for ih in range(2):
                    for nh in range(1, nterms + 1):
                        ihr += 1
                        if ih == 0:
                            C[ihr, icr] = psia[nh, m]
                        else:
                            C[ihr, icr] = psif[nh, m]
                        write0('{:5d}{:5d}{:5d}{:5d}{:15.7g}{:15.7}\n'. \
                                 format(ic+1,n,ih+1,nh,C[ihr, icr].real,C[ihr, icr].imag),1)
        A = np.zeros((2*nterms,2*nterms),dtype=complex)
        for ir in range(2*nterms):
            for ih in range(2*nterms):
                if ir == ih:
                    A[ir, ih] = 1. - C[ir, ih]
                else:
                    A[ir, ih] = - C[ir, ih]
        lu, piv = scipy.linalg.lu_factor(A)
        Ar = scipy.linalg.lu_solve((lu, piv), np.eye(A.shape[0]))
        if m < 0:
            write0('\nSphere correction, negative order\n',1)
        else:
            write0('\nSphere correction, positive order\n',1)
        icr = -1
        for ic in range(2):
            for n in range(1,nterms+1):
                icr += 1
                ihr = -1
                for ih in range(2):
                    for nh in range(1,nterms+1):
                        ihr += 1
                        write0('{:5d}{:5d}{:5d}{:5d}{:15.7g}{:15.7}\n'. \
                                 format(ic+1,n,ih+1,nh,Ar[ihr, icr].real,Ar[ihr, icr].imag),1)
        AI[m,:,:] = Ar

    return AI

def MTSphereReflectionCorrection(nterms, AI, psiai, psifi):

    b = np.zeros(2*nterms, dtype=complex)
    x = np.zeros(2*nterms, dtype=complex)
    psiac = np.zeros((nterms + 1, 3), dtype=complex)
    psifc = np.zeros((nterms + 1, 3), dtype=complex)
    for m in [-1,1]:
        i = -1
        for ic in range(2):
            for n in range(1,nterms + 1):
                i += 1
                if ic == 0:
                    b[i] = psiai[n, m]
                else:
                    b[i] = psifi[n, m]
        x = np.matmul(AI[m,:,:], b)
        i = -1
        for ic in range(2):
            for n in range(1,nterms+1):
                i += 1
                if ic == 0:
                    psiac[n,m] = x[i]
                else:
                    psifc[n,m] = x[i]

    return psiac, psifc


import numpy as np
import scipy
import sphericalfunctions as sphfun
import hssphere_mt as hssphere_mt

def MTLayeredEarthCorrection(nw, nterms, nlyr, thk, depth, radius, sxlyr, yhatl, zhatl, yhats, zhats, rtm, rte, htarg):

    nlat = 4 * nterms + 1
    nlon = 4 * nterms
    center = np.zeros(3)
    center[2] = depth
    ps = sphfun.SphericalGeometryInitialisation(center, radius, nlat, nlon)
    E, H = hssphere_mt.MTSphere_MultipleSources(nw,nterms,nlat,nlon,depth,depth,nlyr,thk,radius,yhatl,zhatl,rtm,rte,htarg)
    AI = np.zeros((3,2 * nterms, 2 * nterms), dtype=complex)
    nw.write(' Spherical Harmonic Analysis\n')
    for m in [-1,1]:
        if m < 0:
            nw.write('\nSphere retroaction, negative order\n')
        else:
            nw.write('\nSphere retroaction, positive order\n')
        icr = -1
        C = np.zeros((2*nterms, 2*nterms), dtype=complex)
        for ic in range(2):
            for n in range(1, nterms + 1):
                icr += 1
                Er, Hr = sphfun.RadialFields(nlat, nlon, ps, center, E[ic,n,m,:,:,:], H[ic,n,m,:,:,:])
                psia, psif = sphfun.MTSphericalHarmonicAnalysis(nterms, nlat, nlon, radius, \
                                yhatl[sxlyr], zhatl[sxlyr], Er, Hr, False)
                ihr = -1
                for ih in range(2):
                    for nh in range(1, nterms + 1):
                        ihr += 1
                        if ih == 0:
                            C[ihr, icr] = psia[nh, m]
                        else:
                            C[ihr, icr] = psif[nh, m]
                        nw.write('{:5d}{:5d}{:5d}{:5d}{:15.7g}{:15.7}\n'. \
                                 format(ic+1,n,ih+1,nh,C[ihr, icr].real,C[ihr, icr].imag))
        A = np.zeros((2*nterms,2*nterms),dtype=complex)
        for ir in range(2*nterms):
            for ih in range(2*nterms):
                if ir == ih:
                    A[ir, ih] = 1. - C[ir, ih]
                else:
                    A[ir, ih] = - C[ir, ih]
        Ar = scipy.linalg.inv(A)
        if m < 0:
            nw.write('\nSphere correction, negative order\n')
        else:
            nw.write('\nSphere correction, positive order\n')
        icr = -1
        for ic in range(2):
            for n in range(1,nterms+1):
                icr += 1
                ihr = -1
                for ih in range(2):
                    for nh in range(1,nterms+1):
                        ihr += 1
                        nw.write('{:5d}{:5d}{:5d}{:5d}{:15.7g}{:15.7}\n'. \
                                 format(ic+1,n,ih+1,nh,Ar[ihr, icr].real,C[ihr, icr].imag))
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


import math
import numpy as np
import layeredearthfunctions as lefun
import sphericalfunctions as sphfun
import empymod

def MTSphere_SingleSource(nterms, psia, psif, psiat, psift, nx, ny, nd, ps, nlyr, thk, yhat, zhat, yhats, zhats, dpthl, zs, radius, htarg):

    # Function estimating electric and magnetic field associated to spherical functions in a layered earth.
    # Computation of the coefficients of the electric and magnetic fields from the spherical potential coefficients.
    # Hankel domain ransformation
    # Computation of the cylindrical coefficients for singular functions
    # Vertical potential computation
    # Application of progagation coefficients
    # Computation of induced fields in frequency domain

    # Input
    # nterms: degree number
    # psia: transverse magnetic spherical potential
    # psif: transverse magnetic spherical potential
    # alt: altitude of the source

    E = np.zeros((nx,ny,nd,3),dtype=complex)
    H = np.zeros((nx,ny,nd,3),dtype=complex)
    ht, htarg = empymod.model.check_hankel('dlf',htarg,0)

    sxlyr = 0                # Identify layer containing loop or GW source
    if nlyr > 0:
        for jz in range(nlyr,0,-1):
            if zs > dpthl[jz]:
                sxlyr = jz
                break
    for jd in range(nd):
        zr = ps[0,0,jd,2]
        rxlyr = 0
        for jz in range(nlyr, 0, -1):
            if zr > dpthl[jz]:
               rxlyr = jz
               break
        for jx in range(nx):
            for jy in range(ny):
                if abs(zr-zs) < radius:
                    dd = ps[jx,jy,jd,2] - zs
                    r = np.sqrt(ps[jx,jy,jd,0]**2+ps[jx,jy,jd,1]**2+dd**2)
                    if r < 0.001:
                        r = 0.001
                        phi = 0.0
                    else:
                        phi = math.atan2( ps[jx,jy,jd,1], ps[jx,jy,jd,0])
                    theta = math.acos(dd/r)
                    if r <= radius:
                        singular = False
                        Enm, Hnm = sphfun.MTFieldSphericalCoefficients(nterms,psiat,psift,yhats,zhats)
                        k = np.sqrt(-yhats*zhats)
                        E[jx,jy,jd,:],H[jx,jy,jd,:] = sphfun.MTFieldsFromSphericalCoefficients(Enm,Hnm,nterms,k,r,theta,phi,singular)
                    else:
                        singular = True
                        Enm, Hnm = sphfun.MTFieldSphericalCoefficients(nterms,psia,psif,yhat[rxlyr],zhat[rxlyr])
                        k = np.sqrt(-yhat[rxlyr]*zhat[rxlyr])
                        E[jx,jy,jd,:],H[jx,jy,jd,:] = sphfun.MTFieldsFromSphericalCoefficients(Enm,Hnm,nterms,k,r,theta,phi,singular)
                else:
                    rho = np.sqrt(ps[jx, jy, jd, 0] ** 2 + ps[jx, jy, jd, 1] ** 2)
                    if rho < 0.001:
                        rho = 0.001
                        phi = 0.
                    else:
                        phi = math.atan2(ps[jx, jy, jd, 1], ps[jx, jy, jd, 0])
                    Eh, Hh = MTSphere_Hnk_SingleSource(nterms, sxlyr, rxlyr, psia, psif, nlyr, thk, dpthl, \
                                            yhat, zhat, zs, zr, rho, htarg)
                    for m in range(-2,3):
                        for ic in range(3):
                                emphi = np.exp ( 1j * m * phi)
                                E[jx,jy,jd,ic] = E[jx,jy,jd,ic] + Eh[m,ic] * emphi
                                H[jx,jy,jd,ic] = H[jx,jy,jd,ic] + Hh[m,ic] * emphi

    return E, H

def MTSphere_Hnk_SingleSource(nterms, sxlyr, rxlyr, psia, psif, nlyr, thk, dpthl, \
                 yhat, zhat, zs, zr, rho, htarg):
    off = np.empty((1,))
    off[0] = rho
    ht, htarg = empymod.model.check_hankel('dlf',htarg,0)
    lambd, int_pts = empymod.transform.get_dlf_points(htarg['dlf'], off, htarg['pts_per_dec'])

    # Call the kernel
    nfilt = lambd.shape[1]
    filt = htarg['dlf']
    Eh = np.zeros((5,3),dtype=complex)
    Hh = np.zeros((5,3),dtype=complex)

    for ift in range(0,nfilt):
        lmbda = lambd[0,ift]
        if lmbda > 0.5:
            continue
        An, Fn = lefun.HSSphere_Ker(lmbda, sxlyr, nlyr, thk, dpthl, zs, yhat, zhat)
        Az, Fz = MTSingleSourceCoefficients(nterms,psia,psif,lmbda,
                            yhat[sxlyr],zhat[sxlyr])
        Em, Hm = MTSingleSourceFields(nlyr, dpthl, sxlyr, rxlyr, zs, zr,
                            lmbda, yhat[rxlyr], zhat[rxlyr], Az, Fz, An, Fn)
        Ew = np.zeros((5, 3), dtype=complex)
        Hw = np.zeros((5, 3), dtype=complex)
        Wj = np.zeros(3)
        Wj[0] = filt.j0[ift]
        if rho > 0.0015:
            Wj[1] = filt.j1[ift]
            Wj[2] = 2 * Wj[1] / (rho * lmbda) - Wj[0]
        for m in range(-2,3):
            for ic in range(3):
                if rho <= 0.0015 and m != 0:
                    continue
                Ew[m,ic] = Ew[m,ic] + Em[ic, m] * Wj[abs(m)]
                Hw[m,ic] = Hw[m,ic] + Hm[ic, m] * Wj[abs(m)]
        Ew = Ew * lmbda
        Hw = Hw * lmbda

        Eh = Eh + Ew
        Hh = Hh + Hw

    Eh = Eh / (4 * np.pi * rho)
    Hh = Hh / (4 * np.pi * rho)

    return Eh, Hh

def MTSingleSourceCoefficients(nterms, psia, psif, lmbda, yhat, zhat):

    k = np.sqrt(- yhat * zhat)
    lmbdasq = lmbda * lmbda
    u = np.sqrt(lmbdasq + yhat * zhat)
    psiac = np.zeros((nterms+1,2*nterms+1), dtype=complex)
    psifc = np.zeros((nterms+1,2*nterms+1), dtype=complex)
    for m in [-1,1]:
        psiac[:,m] = psia[:,m]
        psifc[:,m] = psif[:,m]
    
    Eznm, Hznm = sphfun.VerticalFieldSphericalCoefficients(nterms, psiac, psifc, yhat, zhat)
    zeta = SphericalCylindricalConversionCoefficients(nterms + 1, k, lmbda, u)

    Az = np.zeros((2,5),dtype=complex)
    Fz = np.zeros((2,5),dtype=complex)
    for id in range(2):
         for m in [-1,1]:
             for n in range(1,nterms+2):
                Az[id, m] = Az[id, m] + zeta[id, n, m] * Eznm[n, m]
                Fz[id, m] = Fz[id, m] + zeta[id, n, m] * Hznm[n, m]

    Az = Az * yhat / lmbdasq
    Fz = Fz * zhat / lmbdasq
    
    return Az, Fz

def SphericalCylindricalConversionCoefficients(nterms, k, lmbda, u):
    # Developpment of the conversion matrix from spherical to cylindrical coefficients.
    # Earth element is a cylindrical coefficient. The first two index refer to spherical
    # functions and the last index refer to the cylindrical functions. 

    A, B, C = sphfun.SphericalFactorInitialisation(nterms)

    zeta = np.zeros((2, nterms + 1, 2 * nterms + 1), dtype=complex)
    for id in range(2):
        zeta[id, 0, 0] = 1j * np.sqrt(4 * np.pi) / (u * k)
        for n in range(1,nterms+1):
            for m in range(-n,n+1):
                if m == n:
                    zeta[id,n,m] = - lmbda * zeta[id,n-1,m-1] / ( k * B[n,-n] )
                elif m == -n:
                    zeta[id,n,m] = - lmbda * zeta[id,n-1,m+1] / ( k * B[n,-n] )
                else:
                    if n > abs(m) + 1:
                        tz = zeta[id,n-2,m] * A[n-2,m] / A[n-1,m]
                    else:
                        tz = 0
                    if id == 0:
                        zeta[id,n,m] = tz + u * zeta[0,n-1,m] / ( k * A[n-1,m] )
                    else:
                        zeta[id,n,m] = tz - u * zeta[0,n-1,m] / ( k * A[n-1,m] )

    return zeta

def MTSingleSourceFields(nlyr, dpthl, sxlyr, rxlyr, zs, zr, lmbda, yhat, zhat, Az, Fz, An, Fn):

    lmbdasq = lmbda ** 2
    ur = np.sqrt(lmbdasq + yhat * zhat)
    dz = np.array([- ur, ur],dtype=complex)
    Em = np.zeros((3,5),dtype=complex)
    Hm = np.zeros((3,5),dtype=complex)
    CD = lefun.CylindricalDerivatives(1, lmbda)
    for isr in range(2):
        Azisr = np.zeros(3, dtype=complex)
        Fzisr = np.zeros(3, dtype=complex)
        for m in range(-1,2,2):
            Azisr[m] = Az[isr,m]
            Fzisr[m] = Fz[isr,m]
        Azdx = CD.derivate(0, Azisr)
        Azdy = CD.derivate(1, Azisr)
        Fzdx = CD.derivate(0, Fzisr)
        Fzdy = CD.derivate(1, Fzisr)
        for id in range(2):
            if ( id == 0 and rxlyr == 0 ) or ( id == 1 and rxlyr == nlyr ):
                continue
            elif id == 0:
                decay = np.exp ( - ur * (zr - dpthl[rxlyr]) )
            elif rxlyr == 0:
                decay = np.exp(ur * zr)
            else:
                decay = np.exp( ur * (zr - dpthl[rxlyr]))
            rtm = An[id,isr, rxlyr] * decay
            rte = Fn[id,isr, rxlyr] * decay
            Em[0,:] = Em[0,:] + rtm * Azdx * dz[id] / yhat
            Em[1,:] = Em[1,:] + rtm * Azdy * dz[id] / yhat
            Em[2,:] = Em[2,:] + rtm * Az[isr,:] * lmbdasq / yhat
            Hm[0,:] = Hm[0,:] + rtm * Azdy
            Hm[1,:] = Hm[1,:] - rtm * Azdx
            Em[0,:] = Em[0,:] - rte * Fzdy
            Em[1,:] = Em[1,:] + rte * Fzdx
            Hm[0,:] = Hm[0,:] + rte * Fzdx * dz[id] / zhat
            Hm[1,:] = Hm[1,:] + rte * Fzdy * dz[id] / zhat
            Hm[2,:] = Hm[2,:] + rte * Fz[isr,:] * lmbdasq / zhat
        if sxlyr == rxlyr:
            dzrs = zr - zs
            if dzrs == 0:
                dzrs = 0.001
            if zr > zs and isr == 0 :
                rtm = np.exp ( - ur * dzrs )
                rte = np.exp ( - ur * dzrs )
                dzs = - ur
            elif zr < zs and isr == 1:
                rtm = np.exp ( ur * dzrs )
                rte = np.exp ( ur * dzrs )
                dzs = ur
            else:
                continue
            Em[0,:] = Em[0,:] + rtm * Azdx * dzs / yhat
            Em[1,:] = Em[1,:] + rtm * Azdy * dzs / yhat
            Em[2,:] = Em[2,:] + rtm * Az[isr,:] * lmbdasq / yhat
            Hm[0,:] = Hm[0,:] + rtm * Azdy
            Hm[1,:] = Hm[1,:] - rtm * Azdx
            Em[0,:] = Em[0,:] - rte * Fzdy
            Em[1,:] = Em[1,:] + rte * Fzdx
            Hm[0,:] = Hm[0,:] + rte * Fzdx * dzs / zhat
            Hm[1,:] = Hm[1,:] + rte * Fzdy * dzs / zhat
            Hm[2,:] = Hm[2,:] + rte * Fz[isr,:] * lmbdasq / zhat

    return Em, Hm

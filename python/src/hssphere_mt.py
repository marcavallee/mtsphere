
import empymod
import math
import numpy as np
import layeredearthfunctions as lefun
import halfspacesphere as hssphere
import sphericalfunctions as sphfun

def MTSphere_MultipleSources(nw, nterms, nlat, nlon, zs, zr, nlyr, thkd, radius, yhat, zhat, rtm, rte, htarg):
    """
    Computes the vertical electric and magnetic field spherical coefficients
    above and below the receiver layer.The source is represented by its potential
    spherical coefficients.This version is for a plane wave source, with m = -1, 1.

                          InpuT
                             -----
    nterms - maximum number of degrees.
        zs - source elevation positive downward.
        zr = receiver elevation
      nlyr - number of layers
      thkd - vector of layer thicknesses
      yhat - vector of layer admittivity
      zhat - vector of layer impedivity

                             OuTPuT
                             ------
        E - electric field intensity around the sphere.
        H - magnetic field intensity around the sphere.
    """

    E = np.zeros((2, nterms+1, 3, nlat, nlon, 3), dtype=complex)
    H = np.zeros((2, nterms+1, 3, nlat, nlon, 3), dtype=complex)
    
    sxlyr = 0
    rxlyr = 0
    dpthl = np.zeros(nlyr+1)
    if nlyr > 0:
        for jz in range(1,nlyr):
            dpthl[jz+1] = dpthl[jz] + thkd[jz]
        for jz in range(nlyr,0,-1):
            if zs > dpthl[jz]:
                sxlyr = jz
                break
        for jr in range(nlyr, 0, -1):
            if zr > dpthl[jz]:
                rxlyr = jz
                break
    
    for jlat in range(nlat):
        theta = jlat * np.pi / (nlat - 1)
        zlat = radius * math.cos(theta) + zr
        rho = radius * math.sin(theta)
        if rho < 0.01:
            rho = 0.01
        if jlat == 19:
            x = 1
        Et, Ht = MTSphere_Hnk_MultipleSources(nterms, sxlyr, rxlyr, nlyr, thkd,
                            dpthl, yhat, zhat, zs, zlat, rho, rtm, rte, htarg)   #     nw.write('\nMultipleSources: jlat:{:5}\n'.format(jlat+1))
        for ic in range(2):
            for n in range(1, nterms + 1):
                for m in [-1, 1]:
                    for id in range(3):
                        for jlon in range(nlon):
                            phi = jlon * 2. * np.pi / nlon
                            Ei = 0
                            Hi = 0
                            for mp in range(-2, 3):
                                emphi = np.exp(1j * mp * phi)
                                Ei = Ei + Et[ic, n, m, mp, id] * emphi
                                Hi = Hi + Ht[ic, n, m, mp, id] * emphi
                            E[ic, n, m, jlat, jlon, id] = Ei
                            H[ic, n, m, jlat, jlon, id] = Hi
    return E, H

def MTSphere_Hnk_MultipleSources(nterms, sxlyr, rxlyr, nlyr, thk, dpthl, yhat, zhat, zs, zr, rho, rtm, rte, htarg):

    off = np.empty((1,))
    off[0] = rho
    ht, htarg = empymod.model.check_hankel('dlf',htarg,0)
    lambd, int_pts = empymod.transform.get_dlf_points(htarg['dlf'], off, htarg['pts_per_dec'])

    # call the kernel
    nfilt = lambd.shape[1]
    filt = htarg['dlf']
    Et = np.zeros((2,nterms+1,3,5,3),dtype=complex)
    Ht = np.zeros((2,nterms+1,3,5,3),dtype=complex)
    nfilt2 = np.argmax(np.abs(filt.j0))

    def testjump(nterms, Ew, Et, Hw, Ht):
        tol = 1.e-6
        tol2 = 1.e-35
        jump = True
        for ic in range(2):
            for n in range(1,nterms+1):
                for m in [-1,1]:
                    for mp in range(-2,3):
                        for int in range(3):
                            Qr = abs(Et[ic, n, m, mp, int].real)
                            Qi = abs(Et[ic, n, m, mp, int].imag)
                            if Qr > tol2 and abs(Ew[ic, n, m, mp, int].real) > tol * Qr:
                                jump = False
                            if Qi > tol2 and abs(Ew[ic, n, m, mp, int].imag) > tol * Qi:
                                jump = False
                            Qr = abs(Ht[ic, n, m, mp, int].real)
                            Qi = abs(Ht[ic, n, m, mp, int].imag)
                            if Qr > tol2 and abs(Hw[ic, n, m, mp, int].real) > tol * Qr:
                                jump = False
                            if Qi > tol2 and abs(Hw[ic, n, m, mp, int].imag) > tol * Qi:
                                jump = False

        return jump

    def loop(ift):
        lmbda = lambd[0,ift]
        An, Fn = lefun.HSSphere_Ker(lmbda, sxlyr, nlyr, thk, dpthl, zs, yhat, zhat)
        Az, Fz = MTMultipleSourceCoefficients(nterms, lmbda, yhat[sxlyr], zhat[sxlyr], rtm, rte)
        Em, Hm = MTMultipleSourceFields(nterms, nlyr, dpthl, rxlyr, lmbda, zr, yhat[rxlyr], zhat[rxlyr], Az, Fz, An, Fn)
        Wj = np.zeros(3)
        Ew = np.zeros((2, nterms+1, 3, 5, 3), dtype=complex)
        Hw = np.zeros((2, nterms+1, 3, 5, 3), dtype=complex)
        Wj[0] = filt.j0[ift]
        if rho > 0.01:
            Wj[1] = filt.j1[ift]
            Wj[2] = 2 * Wj[1] / ( rho * lmbda ) - Wj[0]
        for m in range(-2,3):
            Ew[:,:,:, m,:] = Ew[:,:,:, m,:] + Em[:,:,:, m,:] * Wj[abs(m)] * lmbda
            Hw[:,:,:, m,:] = Hw[:,:,:, m,:] + Hm[:,:,:, m,:] * Wj[abs(m)] * lmbda
        return Ew, Hw

    for ift in range(nfilt2,nfilt):
        Ew, Hw = loop(ift)
        Et = Et + Ew
        Ht = Ht + Hw
        jump = testjump(nterms,Ew,Et,Hw,Ht)
        if jump and ift - nfilt2 > 50:
            break

    for ift in range(nfilt2, -1, -1):
        Ew, Hw = loop(ift)
        Et = Et + Ew
        Ht = Ht + Hw
        jump = testjump(nterms, Ew, Et, Hw, Ht)
        if jump and ift - nfilt2 < -50:
            break
            
    Et = Et / (rho * 4. * np.pi)
    Ht = Ht / (rho * 4. * np.pi)
    
    return Et, Ht

def MTMultipleSourceCoefficients(nterms, lmbda, yhat, zhat, rtm, rte):

    Az = np.zeros((2, 2, nterms+1, 5), dtype=complex)
    Fz = np.zeros((2, 2, nterms+1, 5), dtype=complex)

    k = np.sqrt(- yhat * zhat)
    lmbdasq = lmbda * lmbda
    u = np.sqrt(lmbdasq + yhat * zhat)
    zeta = hssphere.SphericalCylindricalConversionCoefficients(nterms + 1, k, lmbda, u)
    for ic in range(2):
        for n in range(1,nterms+1):
            for m in [-1,1]:
                psia = np.zeros((nterms+1, 2*nterms+1), dtype=complex)
                psif = np.zeros((nterms+1, 2*nterms+1), dtype=complex)
                if ic == 0:
                    psia[n, m] = rtm[n]
                else:
                    psif[n, m] = rte[n]
                Eznm, Hznm = sphfun.VerticalFieldSphericalCoefficients(nterms, psia, psif, yhat, zhat)
                for id in range(2):
                    Azs = 0.
                    Fzs = 0.
                    for nr in range(1, nterms + 2):
                        Azs = Azs + zeta[id, nr, m] * Eznm[nr, m]
                        Fzs = Fzs + zeta[id, nr, m] * Hznm[nr, m]
                    Az[id, ic, n, m] = Azs
                    Fz[id, ic, n, m] = Fzs
    Az = Az * yhat / lmbdasq
    Fz = Fz * zhat / lmbdasq

    return Az, Fz

def MTMultipleSourceFields(nterms, nlyr, dpthl, rxlyr, lmbda, z, yhat, zhat, Az, Fz, An, Fn):
    
    """" 
    Computation of the spherical coefficients for cylindrical fields
    
               Input
               -----
       nterms - maximum spherical harmonic degree number
       nlyr - number of layers
       depth - layer top depths
       rxlyr - layer number of the receiver
       zs - depth of the source
       zr - depth of the receiver
       lmbda - cylindrical wavenumber
       yhat - layer admittances
       zhat - layer impedances
       Az - Transverse magnetic vertical potential
       Fz - Transverse electric vertical potential
       An - Tm layered earth propagation
       Fn - TE layered earth propagation
    
               Output
               ------
       Em - Vertical electrical field Hankel transform spherical coefficients
       Hm - Vertical magnetic field Hankel transform spherical coefficients
    
    In the various vectors, id is the propagation direction, 1: Downward, 2: upward
    """

    lmbdasq = lmbda ** 2
    u = np.sqrt(lmbdasq + yhat * zhat)
    Em = np.zeros((2,nterms+1,3,5,3), dtype=complex)
    Hm = np.zeros((2,nterms+1,3,5,3), dtype=complex)
    for id in range(2):
        if id == 0 and rxlyr == 0:
            continue
        if id == 1 and rxlyr == nlyr:
            continue
        if id == 0:
            decay = np.exp(- u * (z - dpthl[rxlyr]))
        elif rxlyr == 0:
            decay = np.exp(u * (z - dpthl[1]))
        else:
            decay = np.exp(u * (z - dpthl[rxlyr]))
        rtm = An[id, rxlyr] * decay
        rte = Fn[id, rxlyr] * decay
        if id == 0:
            Dz = - u
        else:
            Dz = u
        for ic in range(2):
            CD = lefun.CylindricalDerivatives(1, lmbda)
            for n in range(1,nterms+1):
                for m in [-1,1]:
                    Azid = np.zeros(3,dtype=complex)
                    Fzid = np.zeros(3,dtype=complex)
                    Azid[m] = Az[id, ic, n, m]
                    Fzid[m] = Fz[id, ic, n, m]
                    AzDx = CD.derivate(0, Azid)
                    AzDy = CD.derivate(1, Azid)
                    FzDx = CD.derivate(0, Fzid)
                    FzDy = CD.derivate(1, Fzid)
                    Em[ic, n, m, :, 0] = Em[ic, n, m, :, 0] + rtm * AzDx * Dz / yhat
                    Em[ic, n, m, :, 1] = Em[ic, n, m, :, 1] + rtm * AzDy * Dz / yhat
                    Em[ic, n, m, m, 2] = Em[ic, n, m, m, 2] + rtm * Az[id, ic, n, m] * lmbdasq / yhat
                    Hm[ic, n, m, :, 0] = Hm[ic, n, m, :, 0] + rtm * AzDy
                    Hm[ic, n, m, :, 1] = Hm[ic, n, m, :, 1] - rtm * AzDx
                    Em[ic, n, m, :, 0] = Em[ic, n, m, :, 0] - rte * FzDy
                    Em[ic, n, m, :, 1] = Em[ic, n, m, :, 1] + rte * FzDx
                    Hm[ic, n, m, :, 0] = Hm[ic, n, m, :, 0] + rte * FzDx * Dz / zhat
                    Hm[ic, n, m, :, 1] = Hm[ic, n, m, :, 1] + rte * FzDy * Dz / zhat
                    Hm[ic, n, m, m, 2] = Hm[ic, n, m, m, 2] + rte * Fz[id, ic, n, m] * lmbdasq / zhat

    return Em,Hm



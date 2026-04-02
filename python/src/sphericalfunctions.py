
import math
import cmath
import numpy as np
import scipy.special as scispe
import pyshtools as pysh

def gauss_legendre(n):
    """
    Calculate nodes and weights for Gauss-Legendre quadrature.

    Parameters:
    -----------
    n : int
        number of nodes

    Returns:
    --------
    x : numpy.ndarray
        node locations
    w : numpy.ndarray
        node weights
    """
    x = np.zeros(n)
    w = np.zeros(n)
    eps = 1e-15
    pi = np.pi

    for i in range((n + 1) // 2):
        # Initial guess
        z = np.cos(pi * (i + 1 - 0.25) / (n + 0.5))

        # newton-Raphson iteration
        for _ in range(100):
            p1 = 1.0
            p2 = 0.0

            # Compute Legendre polynomial
            for j in range(1, n + 1):
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j

            # Compute derivative
            pp = n * (z * p1 - p2) / (z ** 2 - 1.0)
            z1 = z
            z = z1 - p1 / pp

            if abs(z - z1) < eps:
                break

        x[i] = -z
        x[n - i - 1] = z
        w[i] = 2.0 / ((1.0 - z ** 2) * pp ** 2)
        w[n - i - 1] = w[i]

    return x, w

def SphericalGeometryInitialisation(center, radius, nlat, nlon):

    ps = np.zeros((nlat, nlon, 3))
    if nlat == 1:
        ps[0,0,:] = center
    else:
        xgauss, wgauss = gauss_legendre(nlat)
        for ilat in range(nlat):
            theta = math.acos(xgauss[ilat])
            for ilon in range(nlon):
                phi = ilon * 2 * np.pi / nlon
                ps[ilat, ilon, 0] = center[0] + radius * math.sin(theta) * math.cos(phi)
                ps[ilat, ilon, 1] = center[1] + radius * math.sin(theta) * math.sin(phi)
                ps[ilat, ilon, 2] = center[2] + radius * math.cos(theta)
    return ps

def SphericalFactorInitialisation(nterms):
    
    A = np.zeros((nterms + 1, 2 * nterms + 1))
    B = np.zeros((nterms + 1, 2 * nterms + 1))
    C = np.zeros((nterms + 1, 2 * nterms + 1))
    for n in range(nterms + 1):
        for m in range(-n, n + 1):
            A[n, m] = math.sqrt((n + 1 + m) * (n + 1 - m) / ((2 * n + 1) * (2 * n + 3)))
            if m >= 0:
                B[n, m] = math.sqrt((n - m - 1) * (n - m) / ((2 * n - 1) * (2 * n + 1)))
                C[n, m] = math.sqrt((n - m) * (n + m + 1))
            else:
                B[n, m] = - math.sqrt((n - m - 1) * (n - m) / ((2 * n - 1) * (2 * n + 1)))
                C[n, m] = - math.sqrt((n - m) * (n + m + 1))

    return A, B, C

def Ynm(n, m, theta, phi):

    factor = (-1) ** m * math.sqrt( (2 * n + 1)  * math.factorial(n - abs(m) ) / \
                                    (4 * np.pi * math.factorial(n + abs(m) ) ) )
    Pmn_z, Pmn_d_z = scispe.lpmn(abs(m),n,math.cos(theta))
    result = factor * Pmn_z[abs(m),n] * np.exp(1j * m * phi)

    return result

def ReflectionCoefficients(nterms, radius, yhat, zhat, yhats, zhats):

    RTM = np.zeros(nterms+1,dtype=complex)
    RTE = np.zeros(nterms+1,dtype=complex)

    k = np.sqrt(- yhat * zhat)
    ks = np.sqrt(- yhats * zhats)
    ka = k * radius
    ksa = ks * radius
    
    for n in range(1,nterms+1):
        jka   = scispe.spherical_jn(n,ka)
        djka  = scispe.spherical_jn(n,ka,True)
        jksa  = scispe.spherical_jn(n,ksa)
        djksa  = scispe.spherical_jn(n,ksa,True)
        yka   = scispe.spherical_yn(n,ka)
        dyka  = scispe.spherical_yn(n,ka,True)
        hka   = jka   - 1j * yka
        dhka  = djka  - 1j * dyka

        try:
            dajksajksa = 1 + ksa * djksa / jksa
        except:
            dajksajksa = 1 + 1j * ksa
        dajkajka   = 1 + ka  * djka  / jka
        dahkahka   = 1 + ka  * dhka  / hka

        RTMW = ( yhat * dajksajksa - yhats * dajkajka ) / ( yhat * dajksajksa - yhats * dahkahka )
        RTEW = ( zhat * dajksajksa - zhats * dajkajka ) / ( zhat * dajksajksa - zhats * dahkahka )
        RTM[n] = - RTMW * jka / hka
        RTE[n] = - RTEW * jka / hka

    return RTM, RTE

def RadialFields(nlat, nlon, ps, center, E, H):

    Er = np.zeros((nlat,nlon),dtype=complex)
    Hr = np.zeros((nlat,nlon),dtype=complex)
    for jlat in range(nlat):
        for jlon in range(nlon):
            x = ps[jlat, jlon, 0] - center[0]
            y = ps[jlat, jlon, 1] - center[1]
            z = ps[jlat, jlon, 2] - center[2]
            r = math.sqrt(x * x + y * y + z * z)
            theta = math.acos(z / r)
            phi = math.atan2(y, x)
            Er[jlat, jlon] = math.sin(theta) * math.cos(phi) * E[jlat, jlon, 0] + \
                             math.sin(theta) * math.sin(phi) * E[jlat, jlon, 1] + \
                             math.cos(theta)                 * E[jlat, jlon, 2]
            Hr[jlat, jlon] = math.sin(theta) * math.cos(phi) * H[jlat, jlon, 0] + \
                             math.sin(theta) * math.sin(phi) * H[jlat, jlon, 1] + \
                             math.cos(theta)                 * H[jlat, jlon, 2]

    return Er, Hr

def VerticalFieldSphericalCoefficients(nterms, psia, psif, yhat, zhat):
    
    """
    Computation of the spherical field coefficients
    from the potential spherical coefficients
    
              Input
              -----
     nterms: number of degrees
    phia: transverse magnetic spherical potential
    psif: transverse electric spherical potential
    yhat: admittivity
    zhat: impedivity
    
              Output
              -----
    Eznm : electric field spherical coefficients
    Hznm : magnetic field spherical coefficients
    """

    A, B, C = SphericalFactorInitialisation(nterms + 1)
    k = np.sqrt(- yhat * zhat)

    def ComputeField(phi, psi, phat, sign):

        Fnm = np.zeros((nterms + 2, 2 * nterms + 3), dtype=complex)
        for n in range(nterms + 2):
            for m in range(- n, n + 1):
                Tz = 0.
                if n > 0 and n <= nterms:
                    Tz = Tz + 1j * m * phi[n, m] * sign
                if n < nterms:
                    Tz = Tz + k * (n + 2) * A[n, m] * psi[n + 1, m] / phat
                if n > 1 and abs(m) < n:
                    Tz = Tz + k * (n - 1) * A[n - 1, m] * psi[n - 1, m] / phat
                Fnm[n, m] = Tz

        return Fnm

    Eznm = ComputeField(psif, psia, yhat, 1)
    Hznm = ComputeField(psia, psif, zhat, -1)
    return Eznm, Hznm

def MTSphericalHarmonicAnalysis(nterms, nlat, nlon, radius, yhat, zhat, Er, Hr, singular):

    """
    This subroutine computes the spherical coefficients of an electromagnetic field
    from the radial components using spherical harmonic analysis
    
             Input
             -----
    nterms: degree number
    nlat: number of latitudes
    nlon: number of longitudes
    radius: radius of the sphere
    yhat: background admittivity
    zhat: background impedivity
    er: radial electric field matrix
    hr: radial magnetic field matrix
    singular: select singular or regular solution
    
             Output
             -----
    psia: transverse magnetic potential coefficients
    psif: transverse electric potential coefficients
    """
    
    psia = MTSphericalAnalysis(nterms, nlat, nlon, Er)
    psif = MTSphericalAnalysis(nterms, nlat, nlon, Hr)
    
    kr = np.sqrt(- yhat * zhat) * radius
    for n in range(1,nterms+1):
        psia[n,:] = psia[n,:] * radius * yhat / (n * (n + 1))
        psif[n,:] = psif[n,:] * radius * zhat / (n * (n + 1))
        if singular:
            jkr = scispe.spherical_jn(n,kr)
            ykr = scispe.spherical_yn(n,kr)
            sr = jkr - 1j * ykr
        else:
            sr = scispe.spherical_jn(n,kr)
        psia[n,:] = psia[n,:] / sr
        psif[n,:] = psif[n,:] / sr

    return psia, psif
    
def MTSphericalAnalysis(nterms, nlat, nlon, f):
    
    psi = np.zeros((nterms+1, 3), dtype=complex)
    xgauss, wgauss = gauss_legendre(nlat)
    deltalon = 2. * np.pi / nlon
    for n in range(1,nterms+1):
        for m in [-1,1]:
            integral = 0.
            for jlon in range(nlon):
                phi = jlon * 2. * np.pi / nlon
                for jlat in range(nlat):
                    deltalat = wgauss[jlat]
                    theta = math.acos(xgauss[jlat])
                    integral = integral + f[jlat, jlon] * Ynm(n, -m, theta, phi) * deltalon * deltalat
            psi[n, m] = integral

    return psi

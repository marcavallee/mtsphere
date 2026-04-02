
import numpy as np
from scipy.constants import mu_0,epsilon_0

def PlaneWaveImpedance(nw, nlyr, thk, res, nf, freq):

    # Subroutine to compute the coefficients of a plane wave propagating
    # in a layered earth

    # Convention for layer indexes. nlyr is the number of layer below the ground surface.
    # res has nlyr members, representing the resistivity below the surface.
    # thk has nlyr-1 components.

    E = np.zeros((nf,nlyr+1,2),dtype=complex)
    Zhat = np.zeros((nf,nlyr+1),dtype=complex)
    for jf in range(nf):
        w = 2 * np.pi * freq[jf]
        k = np.zeros(nlyr+1, dtype=complex)
        Z = np.zeros(nlyr+1, dtype=complex)
        for jl in range(nlyr+1):
            if jl == 0:
                k[0] = w * np.sqrt( mu_0 * epsilon_0)
            else:
                k[jl] = np.sqrt(-1j * mu_0 * w/res[jl])
            Z[jl] = w * mu_0 / k[jl]
        Zhat[jf, nlyr] = Z[nlyr]
        for jl in range(nlyr-1,0,-1):
            A = Zhat[jf,jl+1] + Z[jl] * np.tanh(1j * k[jl] * thk[jl])
            B = Z[jl] + Zhat[jf,jl+1] * np.tanh(1j * k[jl] * thk[jl])
            Zhat[jf, jl] = Z[jl] * A / B
        for jl in range(nlyr+1):
            if jl == 0:
                E[jf,jl,0] = 1 + 0j
                E[jf,jl,1] = ( Zhat[jf,1] - Z[0] ) / ( Zhat[jf,1] + Z[0] )
            else:
                if jl == 1:
                    Ei = E[jf,0,0] + E[jf,0,1]
                else:
                    Ei = E[jf,jl-1,0] * np.exp(-1j * k[jl-1] * thk[jl-1]) + \
                         E[jf,jl-1,1] * np.exp( 1j * k[jl-1] * thk[jl-1])
                if jl == nlyr:
                    E[jf,jl,0] = Ei
                    E[jf,jl,1] = 0.
                else:
                    den = ( Zhat[jf,jl+1] - Z[jl] ) * np.exp( - 2j * k[jl] * thk[jl] ) + \
                          ( Zhat[jf,jl+1] + Z[jl] )
                    E[jf, jl, 0] = ( Zhat[jf,jl+1] + Z[jl] ) * Ei / den
                    E[jf, jl, 1] = ( Zhat[jf,jl+1] - Z[jl] ) * np.exp( - 2j * k[jl] * thk[jl] ) * Ei / den

    nw.write('\n Layer impedance\n')
    for jf in range(nf):
        for jl in range(1,nlyr+1):
            nw.write(' Frequency: {:15.7g}, Layer: {:5}, Impedance: {:15.7g} {:15.7g}\n'. \
                format(freq[jf],jl,Zhat[jf,jl].real,Zhat[jf,jl].imag))

    nw.write('\n Propagation coefficients\n')
    for jf in range(nf):
        for jl in range(nlyr+1):
            nw.write(' Frequency:, {:15.7g}, Layer: {:5}, Earthward propagation: {:15.7g} {:15.7g}, Sunward propagation: {:15.7g} {:15.7g}\n'. \
            format(freq[jf],jl,E[jf,jl,0].real,E[jf,jl,0].imag, E[jf,jl,1].real, E[jf,jl,1].imag ) )

    return Zhat, E

def apparentresistivity(nf, nx, ny, freq, imp):

    appres = np.zeros((nf,nx,ny,2,2))
    phase = np.zeros((nf,nx,ny,2,2))
    
    for jf in range(nf):
        den = 2 * np.pi * freq[jf] * mu_0
        for jx in range(nx):
            for jy in range(ny):
                for i in range(2):
                    for j in range(2):
                        appres[jf,jx,jy,i,j] = abs ( imp[jf,jx,jy,i,j] ** 2 ) / den
                        phase[jf,jx,jy,i,j] = np.angle( imp[jf,jx,jy,i,j] ) * 180 / np.pi

    return appres, phase

def PlaneWaveSphericalHarmonics(nterms, yhat, zhat, h, E):

    """
    Subroutine to compute the spherical coefficients of the upward and downward
    propagation of a sphere located in a layered earth

             Input
             -----
    nterms: number of spherical coefficients.
    yhat: admittivity
    zhat: impeditivity
    h: elevation of the source
    E: source electric field

              Output
              -----
    psia: transverse magnetic radial coefficients
    psif: transverse electric radial coefficients
    """

    psia = np.zeros((2,nterms+1,3),dtype=complex)
    psif = np.zeros((2,nterms+1,3),dtype=complex)
    k = np.sqrt(- yhat * zhat)
    for ic in range(2): # Downward and upward orientations
        if E[ic] == 0:
            continue
        if ic == 0:
            sic = -1
        else:
            sic = 1
        E0 = E[ic] * np.exp ( sic * 1j * k * h )
        for id in range(2):
            psiat, psift = PlaneWaveFieldssphericalharmonics(E0,nterms,ic,id,k,yhat,zhat)
            psia[id,:,:] = psia[id,:,:] + psiat
            psif[id,:,:] = psif[id,:,:] + psift

    return psia, psif

def PlaneWaveFieldssphericalharmonics(E0, nterms, ic, id, k, yhat, zhat):

    # Kernel of the representation of a plane wave as spherical expansion

    psia = np.zeros((nterms+1,3), dtype=complex)
    psif = np.zeros((nterms+1,3), dtype=complex)
    impedance = np.sqrt(zhat/yhat)
    if ic == 0:
        sic = -1
    else:
        sic = 1
    for n in range(1, nterms + 1):
        for m in [-1,1]:
            An = E0 * (sic * 1j ) ** n * np.sqrt(np.pi*(2*n+1)/(n*(n+1))) / k
            if id == 0:
                Afact = - sic * 1j
                Ffact = m
            else:
                Afact = - sic * m
                Ffact = - 1j
            psia[n, m] = Afact * An * yhat
            psif[n, m] = Ffact * An * zhat / impedance

    return psia, psif

def HSSphere_Ker(lmbda,sxlyr,nlyr,thk,dpthl,zs,yhat,zhat):

    # thk index is n-1 as the thk vector starts with zero.

    if nlyr == 0:
        An = np.zeros((2,1))
        Fn = np.zeros((2,1))
        return An, Fn
            
    s  = np.zeros(nlyr+1,dtype=complex)
    adm = np.zeros(nlyr+1,dtype=complex)
    imp = np.zeros(nlyr+1,dtype=complex)
    admhat = np.zeros(nlyr+1,dtype=complex)
    imphat = np.zeros(nlyr+1,dtype=complex)
    EXP_TOL = 80.0

    lmbsq = lmbda ** 2
    for j in range(nlyr+1):
        s[j] = np.sqrt ( lmbsq + yhat[j] * zhat[j] )
        adm[j] = s[j] / zhat[j]
        imp[j] = s[j] / yhat[j]
    
     # Admittance and impedance computation for layers higher than source layer.
    
    if sxlyr < nlyr:
        admhat[nlyr] = adm[nlyr]
        imphat[nlyr] = imp[nlyr]
        for j in range(nlyr-1,sxlyr,-1):
            Ej = np.tanh ( s[j] * thk[j] )
            admhat[j] = adm[j] * ( admhat[j+1] + adm[j] * Ej ) / \
                                 ( adm[j] + admhat[j+1] * Ej )
            imphat[j] = imp[j] * ( imphat[j+1] + imp[j] * Ej ) / \
                                 ( imp[j] + imphat[j+1] * Ej )

    # Admittance and impedance computation for layers higher than source layer.

    if sxlyr > 0:
        admhat[1] = - adm[0]
        imphat[1] = - imp[0]
        for j in range(1,sxlyr):
            Ei = s[j] * thk[j]
            if abs(Ei) < 100:
                Ej = np.tanh ( Ei )
            else:
                Ej = 0.
            admhat[j+1] = adm[j] * ( admhat[j] - adm[j] * Ej ) / \
                                   ( adm[j] - admhat[j] * Ej )
            imphat[j+1] = imp[j] * ( imphat[j] - imp[j] * Ej ) / \
                                   ( imp[j] - imphat[j] * Ej )

    An = np.zeros((2, nlyr + 1), dtype=complex)
    Fn = np.zeros((2, nlyr + 1), dtype=complex)
    At = np.zeros(nlyr+1,dtype=complex)
    Ft = np.zeros(nlyr+1,dtype=complex)

    if sxlyr == 0:
        ep = s[0] * ( dpthl[0] - zp )
        if ep.real < EXP_TOL:
            Ej = np.exp ( - ep )
        else:
            Ej = 0
        Fn[1,0] = Ej * ( adm[0] - admhat[1] ) / ( adm[0] + admhat[1] )
        An[1,0] = Ej * ( imp[0] - imphat[1] ) / ( imp[0] + imphat[1] )
        Ft[1] = Fn[1,0] + Ej
        At[1] = An[1,0] + Ej
    elif sxlyr == nlyr:
        ep = s[nlyr] * ( zs - dpthl[nlyr] )
        if ep.real < EXP_TOL:
            Ej = np.exp ( - ep )
        else:
            Ej = 0
        if abs ( adm[nlyr] - admhat[nlyr] ) > 0:
            Fn[0,nlyr] = Ej * ( adm[nlyr] + admhat[nlyr] ) / ( adm[nlyr] - admhat[nlyr] )
        if abs ( imp[nlyr] - imphat[nlyr] ) > 0:
            An[0,nlyr] = Ej * ( imp[nlyr] + imphat[nlyr] ) / ( imp[nlyr] - imphat[nlyr] )
        Ft[nlyr] = Fn[0,nlyr] + Ej
        At[nlyr] = An[0,nlyr] + Ej
    else:
        Ep = 2. * s[sxlyr] * thk[sxlyr]
        if Ep.real < EXP_TOL:
            Emh = np.exp(-Ep/2)
            Eph = np.exp(Ep/2)
            E2h = np.exp(-Ep)
            Es1 = np.exp( - s[sxlyr] * ( zs - dpthl[sxlyr] ) )
            Es2 = np.exp( - s[sxlyr] * ( dpthl[sxlyr] - zs ) )
        else:
            Emh = 0
            Eph = 0
            E2h = 0
            Es1 = 0
            Es2 = 0
        Fden = ( admhat[sxlyr] - adm[sxlyr] ) * ( admhat[sxlyr+1] + adm[sxlyr] ) - \
               ( admhat[sxlyr] + adm[sxlyr] ) * ( admhat[sxlyr+1] - adm[sxlyr] ) * E2h
        Aden = ( imphat[sxlyr] - imp[sxlyr] ) * ( imphat[sxlyr+1] + imp[sxlyr] ) - \
               ( imphat[sxlyr] + imp[sxlyr] ) * ( imphat[sxlyr+1] - imp[sxlyr] ) * E2h
        Fn[0,sxlyr] = - ( ( admhat[sxlyr] + adm[sxlyr] ) * \
                          ( ( admhat[sxlyr+1] + adm[sxlyr] ) * Es1 - \
                            ( admhat[sxlyr+1] - adm[sxlyr] ) * Emh * Es2 ) ) / Fden
        An[0,sxlyr] = - ( ( imphat[sxlyr] + imp[sxlyr] ) * \
                          ( ( imphat[sxlyr+1] + imp[sxlyr] ) * Es1 - \
                            ( imphat[sxlyr+1] - imp[sxlyr] ) * Emh * Es2 ) ) / Aden
        Fn[1,sxlyr] = ( ( admhat[sxlyr+1] - adm[sxlyr] ) * \
                          ( ( admhat[sxlyr] + adm[sxlyr] ) * E2h * Es1 - \
                            ( admhat[sxlyr] - adm[sxlyr] ) * Emh * Es2 ) ) / Fden
        An[1,sxlyr] = ( ( imphat[sxlyr+1] - imp[sxlyr] ) * \
                          ( ( imphat[sxlyr] + imp[sxlyr] ) * E2h * Es1 - \
                            ( imphat[sxlyr] - imp[sxlyr] ) * Emh * Es2 ) ) / Aden
        Ft[sxlyr] = Fn[0,sxlyr] + Fn[1,sxlyr] + Es1
        At[sxlyr] = An[0,sxlyr] + An[1,sxlyr] + Es1
        Ft[sxlyr+1] = Fn[0,sxlyr] * Emh + Fn[1,sxlyr] * Eph + Es2
        At[sxlyr+1] = An[0,sxlyr] * Emh + An[1,sxlyr] * Eph + Es2

    for j in range(sxlyr-1, -1, -1):
        if j == 0:
            An[1, 0] = At[1]
            Fn[1, 0] = Ft[1]
        else:
            Ep = 2. * s[j] * thk[j]
            if Ep.real < EXP_TOL:
                Emh = np.exp(-Ep/2)
                E2h = np.exp(-Ep)
            else:
                Emh = 0
                E2h = 0
            At[j] = At[j + 1] * 2 * imp[j] * Emh / (imp[j] * (1. + E2h) - imphat[j] * (1. - E2h))
            Ft[j] = Ft[j + 1] * 2 * adm[j] * Emh / (adm[j] * (1. + E2h) - admhat[j] * (1. - E2h))
            An[0, j] = At[j] * (imp[j] + imphat[j]) / ( 2.0 * imp[j] )
            An[1, j] = At[j] * (imp[j] - imphat[j]) / ( 2.0 * imp[j] )
            Fn[0, j] = Ft[j] * (adm[j] + admhat[j]) / ( 2.0 * adm[j] )
            Fn[1, j] = Ft[j] * (adm[j] - admhat[j]) / ( 2.0 * adm[j] )

    for j in range(sxlyr+1,nlyr+1):
        if j == nlyr:
            An[0, nlyr] = At[nlyr]
            Fn[0, nlyr] = Ft[nlyr]
        else:
            Ep = s[j] * thk[j]
            if Ep.real < EXP_TOL:
                Emh = np.exp(-Ep)
                Eph = np.exp(Ep)
            else:
                Emh = 0
                Eph = 0
            An[0, j] = At[j] * (imp[j] + imphat[j]) / ( 2.0 * imp[j] )
            An[1, j] = At[j] * (imp[j] - imphat[j]) / ( 2.0 * imp[j] )
            Fn[0, j] = Ft[j] * (adm[j] + admhat[j]) / ( 2.0 * adm[j] )
            Fn[1, j] = Ft[j] * (adm[j] - admhat[j]) / ( 2.0 * adm[j] )
            At[j + 1] = An[0, j] * Emh + An[1, j] * Eph
            Ft[j + 1] = Fn[0, j] * Emh + Fn[1, j] * Eph

    return An, Fn

class CylindricalDerivatives:

    """
    Class to compute the cylindrical derivatives.
    It is used the following way.
    CD = lefun.CylindricalDerivatives(nterms,lmbda)
    Azdx = CD.derivate(der, Az)
    If der = 0, x derivative, otherwise der = 1, y derivative.
    The length of Az must be equal or larger than 2 * nterms + 1.
    """

    def __init__(self,n,lmbda):

        self.n = n
        dxpiy = np.zeros((2*(n+1)+1,2*n+1))
        dxmiy = np.zeros((2*(n+1)+1,2*n+1))

        for m in range(-n,n+1):
            if m == 0:
                dxpiy[ 1,0]  = -1
                dxmiy[-1,0]  = -1
            elif m < 0:
                mp = m + 1
                dxpiy[mp,m] = 1
                mp = m - 1
                dxmiy[mp,m] = -1
            elif m > 0:
                mp = m + 1
                dxpiy[mp,m] = -1
                mp = m - 1
                dxmiy[mp,m] =  1

        dxpiy = dxpiy * lmbda
        dxmiy = dxmiy * lmbda
        self.dx =   0.5  * ( dxpiy + dxmiy )
        self.dy = - 0.5j * ( dxpiy - dxmiy )

    def derivate(self,ider,F):

        if len(F) == 2*self.n+1:
            Fc = F
        elif len(F) > 2*self.n+1:
            Fc = np.zeros(2*self.n+1,dtype=complex)
            Fc[0] = F[0]
            for i in range(1,self.n+1):
                Fc[i] = F[i]
                Fc[-i] = F[-i]
        else:
            print('Wrong F size in CylindricalDerivatives.')
            return
        dF = np.zeros(2*(self.n+1)+1,dtype=complex)
        for m in range(-self.n-1,self.n+2):
            for mp in range(-self.n,self.n+1):
                if ider == 0:
                    dF[m] = dF[m] + self.dx[m,mp] * Fc[mp]
                elif ider == 1:
                    dF[m] = dF[m] + self.dy[m,mp] * Fc[mp]

        return dF

 
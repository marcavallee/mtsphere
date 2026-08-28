
import math
import numpy as np
from scipy.constants import pi, mu_0, epsilon_0
from layeredearthspherecorrection import MTLayeredEarthCorrection,MTSphereReflectionCorrection
from layeredearthfunctions import PlaneWaveSphericalHarmonics
import sphericalfunctions as sphfun
import halfspacesphere as hssphere

import numpy as np

def dynamic_cleaning_vectorized(F):
    """Nettoie le bruit numérique sur une grille 3D de tenseurs complexes.

    Version vectorisée NumPy (très rapide).
    Modifie le tableau F 'in-place' (directement en mémoire).
    """
    # 1. Extraction des parties réelles et imaginaires de toutes les composantes
    # F a la forme (NX, NY, ND, 3)
    real_part = np.real(F)
    imag_part = np.imag(F)

    # 2. Calcul du max_local pour chaque cellule (jx, jy, jd)
    # On prend le max absolu parmi les réels et imaginaires des 3 composantes
    # axis=-1 permet de chercher le max sur la dernière dimension (les 3 composantes)
    max_real = np.max(np.abs(real_part), axis=-1)
    max_imag = np.max(np.abs(imag_part), axis=-1)

    # max_local aura la forme (NX, NY, ND)
    max_local = np.maximum(max_real, max_imag)

    # 3. Calcul du seuil (ceiling) pour chaque cellule
    # On ajoute np.newaxis pour pouvoir appliquer le seuil aux 3 composantes
    ceiling_local = (max_local * 1e-12)[..., np.newaxis]

    # 4. Application du filtre (mise à zéro si inférieur au seuil)
    real_cleaned = np.where(np.abs(real_part) < ceiling_local, 0.0, real_part)
    imag_cleaned = np.where(np.abs(imag_part) < ceiling_local, 0.0, imag_part)

    # 5. Reconstruction du tableau complexe d'origine
    F[:] = real_cleaned + 1j * imag_cleaned

def mtsphere3d(nw, nf, nlyr, nterms, nx, ny, nd, freq, pxyd, thk, res, depth,
               radius, sres, E, dpthl, htarg):

    zs = depth
    sxlyr = 0
    for jz in range(nlyr,0,-1):
        if zs > dpthl[jz]:
            sxlyr = jz
            break
    impedance = np.zeros((nf,nx,ny,nd,3,2),dtype=complex)
    Es = np.zeros((nf, nx, ny, nd, 3, 2), dtype=complex)
    Hs = np.zeros((nf, nx, ny, nd, 3, 2), dtype=complex)
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
        rtm, rte, ttm, tte = sphfun.ReflectionCoefficients(nterms, radius, yhatl[sxlyr], zhatl[sxlyr], yhats, zhats)
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
            psiat = np.zeros((nterms + 1, 3), dtype=complex)
            psift = np.zeros((nterms + 1, 3), dtype=complex)
            for n in range(1,nterms+1):
                psiar[n,:] = rtm[n] * psiac[n,:]
                psifr[n,:] = rte[n] * psifc[n,:]
                psiat[n,:] = ttm[n] * psiac[n,:]
                psift[n,:] = tte[n] * psifc[n,:]
            Es[jf,:,:,:,:,id], Hs[jf,:,:,:,:,id] = hssphere.MTSphere_SingleSource(nterms, psiar, psifr, psiat, psift,\
                                    nx, ny, nd, pxyd, nlyr, thk, yhatl, zhatl, yhats, zhats, dpthl, zs, radius, htarg)
            dynamic_cleaning_vectorized(Es[jf,:,:,:,:,id])
            dynamic_cleaning_vectorized(Hs[jf,:,:,:,:,id])
            nw.write('\nInduced fields\n')
            for jx in range(nx):
                for jy in range(ny):
                    for jd in range(nd):
                        nw.write('{:15.7g}{:5}{:5}{:5}'.format(freq[jf],jx+1,jy+1,jd+1))
                        for ic in range(3):
                            nw.write('{:15.7g}{:15.7g}'.format(Es[jf,jx,jy,jd,ic,id].real,Es[jf,jx,jy,jd,ic,id].imag))
                        for ic in range(3):
                            nw.write('{:15.7g}{:15.7g}'.format(Hs[jf,jx,jy,jd,ic,id].real,Hs[jf,jx,jy,jd,ic,id].imag))
                        nw.write('\n')

    return impedance, Es, Hs

def zpinv2(A, tol):
    """Calcule la pseudo-inverse de Moore-Penrose d'une matrice complexe 2x2.

    Équivalent strict de la subroutine Fortran ZPINV2.
    """
    # SVD de la matrice complexe A
    try:
    # U, S, VH correspondent aux sorties de ZGESVD
    # Note : NumPy retourne déjà VH (V conjuguée transposée), et non VT
        U, S, VH = np.linalg.svd(A, full_matrices=True)
    except np.linalg.LinAlgError:
    # Équivalent de IF (info /= 0) THEN
        return np.zeros((2, 2), dtype=complex)

    # Seuil de tolérance (S(i) > tol * smax)
    smax = np.max(S)

    # Inversion des va0leurs singulières sous condition
    S_plus = np.where((S>tol*smax)&(S>0),1.0/S,0.0)

    # Création de la matrice diagonale Sigma+
    W = np.diag(S_plus)

    # Calcul final : Aplus = V * Sigma+ * U^H
    # En NumPy : VH est V^H, donc V est VH.conj().T
    # U^H est U.conj().T
    V = VH.conj().T
    UH = U.conj().T

    Aplus = V @ W @ UH
    Aplus.real = np.where(np.abs(Aplus.real) < tol, 0.0, Aplus.real)
    Aplus.imag = np.where(np.abs(Aplus.imag) < tol, 0.0, Aplus.imag)

    return Aplus

    return Aplus

def getimpedance(nf, nx, ny, nd, nlyr, nterms, pxyd, dpthl, K, Z, depth, radius, E, Es, Hs):

    impedance = np.zeros((nf,nx,ny,nd,3,2), dtype=complex)
    Et = Es.copy()
    Ht = Hs.copy()
    for jf in range(nf):
        for jd in range(nd):
            dr = pxyd[0,0,jd,2]
            rxlyr = 0
            for jl in range(nlyr,0,-1):
                if dr > dpthl[jl]:
                    rxlyr = jl
                    break
            for jx in range(nx):
                for jy in range(ny):
                    r = math.sqrt ( pxyd[jx,jy,jd,0]**2 + pxyd[jx,jy,jd,1]**2 + ( pxyd[jx,jy,jd,2] - depth )**2 )
                    if r > radius:
                        if rxlyr == nlyr:
                            EE = E[jf,nlyr,0] * np.exp ( - 1j * K[jf,rxlyr] * ( dr - dpthl[rxlyr] ) )
                            HH = E[jf,nlyr,0] / Z[jf,rxlyr]
                        elif rxlyr == 0:
                            EE = E[jf,0,0] * np.exp ( - 1j * K[jf,rxlyr] * dr ) \
                               + E[jf,0,1] * np.exp (   1j * K[jf,rxlyr] * dr )
                            HH = ( E[jf,0,0] * np.exp ( - 1j * K[jf,rxlyr] * dr ) \
                                 - E[jf,0,1] * np.exp (   1j * K[jf,rxlyr] * dr ) ) / Z[jf,rxlyr]
                        else:
                            EE = E[jf,rxlyr,0] * np.exp ( - 1j * K[jf,rxlyr] * ( dr - dpthl[rxlyr] ) ) \
                               + E[jf,rxlyr,1] * np.exp (   1j * K[jf,rxlyr] * ( dr - dpthl[rxlyr] ) )
                            HH = ( E[jf,rxlyr,0] * np.exp ( - 1j * K[jf,rxlyr] * ( dr - dpthl[rxlyr] ) ) \
                                 - E[jf,rxlyr,1] * np.exp (   1j * K[jf,rxlyr] * ( dr - dpthl[rxlyr] ) ) ) / Z[jf,rxlyr]
                        Et[jf,jx,jy,jd,0,0] = Et[jf,jx,jy,jd,0,0] + EE
                        Et[jf,jx,jy,jd,1,1] = Et[jf,jx,jy,jd,1,1] + EE
                        Ht[jf,jx,jy,jd,0,1] = Ht[jf,jx,jy,jd,0,1] - HH
                        Ht[jf,jx,jy,jd,1,0] = Ht[jf,jx,jy,jd,1,0] + HH
                    Hmat = np.array([Ht[jf,jx,jy,jd,0,:],Ht[jf,jx,jy,jd,1,:]])
                    Hplus = zpinv2(Hmat,1.0e-12)
                    for ic in range(3):
                        C = np.zeros(2, dtype=complex)
                        if ic < 2:
                            for ip in range(2):
                                C[ip] = Et[jf,jx,jy,jd,ic,ip]
                        else:
                            for ip in range(2):
                                C[ip] = Ht[jf,jx,jy,jd,2,ip]
                        impedance[jf,jx,jy,jd,ic,0] = C[0] * Hplus[0,0] + C[1] * Hplus[1,0]
                        impedance[jf,jx,jy,jd,ic,1] = C[0] * Hplus[0,1] + C[1] * Hplus[1,1]

    return impedance, Et, Ht



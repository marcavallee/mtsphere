    
SUBROUTINE MTSphere_MultipleSources(NW,NTERMS,NLAT,NLON,ZS,ZR,NLYR,THKD,RADIUS,YHAT,ZHAT,RTM,RTE,E,H)
!  Computes the vertical electric and magnetic field spherical coefficients
!  above and below the receiver layer. The source is represented by its potential
!  spherical coefficients. This version is for a plane wave source, with M = -1,1.
!
!                             INPUT
!                             -----
!    NTERMS - maximum number of degrees.
!      PSIA - TM spherical coefficients.
!      PSIF - TE spherical coefficients.
!        ZS - source elevation positive downward.
!        ZR = receiver elevation
!      NLYR - number of layers
!      THKD - vector of layer thicknesses
!      YHAT - vector of layer admittivity
!      ZHAT - vector of layer impedivity
!
!                             OUTPUT
!                             ------
!        E - electric field intensity around the sphere.
!        H - magnetic field intensity around the sphere.
!

use precision

 IMPLICIT NONE

 INTEGER NTERMS,NLAT,NLON,NLYR,JZ,SLYR,RLYR,IC,N,M,MP,INT,JLAT,JLON,NW
 REAL(KIND=QL) THKD(NLYR-1),DPTHL(NLYR),ZS,ZR,RADIUS,ZLAT,LMBDA,THETA, &
                RHO, PHI
 
 COMPLEX(KIND=QL) YHAT(0:NLYR), zhat(0:NLYR), RTE(NTERMS), RTM(NTERMS), EMPHI, &
                   E(2,NTERMS,-1:1,NLAT,NLON,3), H(2,NTERMS,-1:1,NLAT,NLON,3), &
                   ET(2,nterms,-1:1,-2:2,3), HT(2,nterms,-1:1,-2:2,3), EI, HI
  
    SLYR = 0              
    RLYR = 0
    DPTHL = 0._QL
    IF ( NLYR > 0 ) THEN
        DO JZ = 1, NLYR-1
        DPTHL(JZ+1) = DPTHL(JZ) + THKD(JZ)
        END DO
        DO JZ = NLYR,1,-1
        IF (ZS > DPTHL(JZ)) THEN
            SLYR = JZ
            EXIT
        END IF
        END DO
        DO JZ = NLYR,1,-1
        IF (ZR > DPTHL(JZ)) THEN
            RLYR = JZ
            EXIT
        END IF
        END DO
    END IF
 
    do JLAT = 1, NLAT
        THETA = ( JLAT - 1 ) * PI / ( NLAT - 1 )
        ZLAT = RADIUS * COS(THETA) + ZR
        RHO = RADIUS * SIN(THETA)
        if ( RHO < 0.01_QL ) RHO = 0.01_QL
        call MTSPHERE_HNK_MultipleSources(NTERMS, SLYR, RLYR, NLYR, THKD, DPTHL, YHAT, ZHAT, ZS, ZLAT, RHO, RTM, RTE, ET, HT)
        do IC = 1, 2
            do n = 1, NTERMS
                do M = -1, 1, 2
                    do INT = 1, 3
                        do JLON = 1, NLON
                            PHI = ( JLON - 1 ) * 2._QL * PI / NLON
                            EI = ZERO
                            HI = ZERO
                            do MP = -2, 2
                                EMPHI = EXP ( CI * MP * PHI ) 
                                EI = EI + ET(IC,N,M,MP,INT) * EMPHI
                                HI = HI + HT(IC,N,M,MP,INT) * EMPHI
                            end do
                            E(IC,N,M,JLAT,JLON,INT) = EI
                            H(IC,N,M,JLAT,JLON,INT) = HI
                        end do
                    end do
                end do
            end do
        end do
    end do
    
end subroutine MTSphere_MultipleSources

subroutine MTSPHERE_HNK_MultipleSources(NTERMS, SLYR, RLYR, NLYR, THKD, DPTHL, YHAT, ZHAT, ZS, ZR, RHO, RTM, RTE, ET, HT)

use FILT_COEF_Q

implicit NONE

 INTEGER NTERMS,NLYR,I,SLYR,RLYR
 REAL(KIND=QL) DEL_JN,Y,RHO_JN,LMBDA,RHO,THKD(NLYR-1),DPTHL(NLYR),ZS,ZR
 COMPLEX(KIND=QL) ET(2,NTERMS,-1:1,-2:2,3), HT(2,NTERMS,-1:1,-2:2,3), &
    RTM(NTERMS),RTE(NTERMS),YHAT(0:NLYR),ZHAT(0:NLYR), &
    FN(2,0:NLYR), AN(2,0:NLYR)
 LOGICAL JUMP

 DEL_JN = LOG (10.D0) / DBLE (NDEC_JN)
 RHO_JN = -LOG (RHO) - SHFTJN
    
 ET = ZERO
 HT = ZERO
 
 DO I = -50, JNHI             ! Start at I = -50 to pick up low values.
   Y = RHO_JN + DBLE(I) * DEL_JN
   LMBDA = EXP (Y)
   IF ( NLYR > 0 ) THEN
        CALL HSSPHERE_KER (LMBDA,SLYR,NLYR,THKD,DPTHL,ZS,YHAT,ZHAT,AN,FN)
   END IF
   CALL MTSPHERE_JMP_MultipleSource
   IF (JUMP .AND. I > -40) EXIT
 END DO

 JUMP = .FALSE.           ! Finish off the low end 
 DO I = -51, JNLO, -1
   Y = RHO_JN + DBLE(I) * DEL_JN
   LMBDA = EXP (Y)
   IF ( NLYR > 0 ) THEN
        CALL HSSPHERE_KER (LMBDA,SLYR,NLYR,THKD,DPTHL,ZS,YHAT,ZHAT,AN,FN)
   END IF
   CALL MTSPHERE_JMP_MultipleSource
   IF (JUMP .AND. I < -60) EXIT

 END DO

ET = ET / ( RHO * 4._QL * PI )
HT = HT / ( RHO * 4._QL * PI )

contains

   SUBROUTINE MTSPHERE_JMP_MultipleSource
!  ---------------------

!***  Called by: HSMD_HNK

!  Accumulates function calls for the Hankel transformation &
!  checks convergence.
   
     REAL(KIND=QL), PARAMETER :: TOL=1.D-6, TOL2=1.D-35
     integer IC, N, M, INT, MP
     REAL(KIND=QL) QR,QI,WJ(0:2)
     COMPLEX(KIND=QL), ALLOCATABLE :: AZ(:,:,:,:), FZ(:,:,:,:), EM(:,:,:,:,:), HM(:,:,:,:,:), &
                                      EW(:,:,:,:,:), HW(:,:,:,:,:)
     
     ALLOCATE ( AZ (2,2,NTERMS,-2:2), FZ (2,2,NTERMS,-2:2), &
                     EM(2,NTERMS,-1:1,-2:2,3), HM(2,NTERMS,-1:1,-2:2,3), &
                     EW(2,NTERMS,-1:1,-2:2,3), HW(2,NTERMS,-1:1,-2:2,3) )
     CALL MTMultipleSourceCoefficients(NTERMS,LMBDA,YHAT(SLYR),ZHAT(SLYR),RTM,RTE,AZ,FZ)
     CALL MTMultipleSourceFields(NTERMS,NLYR,DPTHL,RLYR,LMBDA,ZR,yhat(RLYR),zhat(RLYR), &
         AZ,FZ,AN,FN,EM,HM)
      
    ! Calcul des champs
  
    EW = ZERO
    HW = ZERO
     
    wj = 0._QL
    wj(0) = WJ0(I)
    if ( RHO .GT. 0.01_QL ) THEN
        WJ(1) = WJ1(I)
        WJ(2) = 2._QL * WJ(1) / ( RHO * LMBDA ) - WJ(0)
    end if
    do M = -2, 2
        EW(:,:,:,M,:) = EW(:,:,:,M,:) + EM(:,:,:,M,:) * WJ(ABS(M)) * LMBDA
        HW(:,:,:,M,:) = HW(:,:,:,M,:) + HM(:,:,:,M,:) * WJ(ABS(M)) * LMBDA
    end DO

    ET = ET + EW
    HT = HT + HW
    
    JUMP = .TRUE.
    do ic = 1, 2
        do n = 1, NTERMS
            do M = -1, 1, 2
                do MP = -2, 2
                    do INT = 1, 3
                        QR = ABS (REAL  (ET(IC,N,M,MP,INT), KIND=QL) )
                        QI = ABS (AIMAG (ET(IC,N,M,MP,INT) ) )
                        IF (QR > TOL2 .AND. ABS (REAL  (EW(IC,N,M,MP,INT))) > TOL * QR) JUMP = .FALSE.
                        IF (QI > TOL2 .AND. ABS (AIMAG (EW(IC,N,M,MP,INT))) > TOL * QI) JUMP = .FALSE.
                        QR = ABS (REAL  (HT(IC,N,M,MP,INT), KIND=QL) )
                        QI = ABS (AIMAG (HT(IC,N,M,MP,INT) ) )
                        IF (QR > TOL2 .AND. ABS (REAL  (HW(IC,N,M,MP,INT))) > TOL * QR) JUMP = .FALSE.
                        IF (QI > TOL2 .AND. ABS (AIMAG (HW(IC,N,M,MP,INT))) > TOL * QI) JUMP = .FALSE.
                    end do
                end do
            end do
        END DO
    END DO
    
   deallocate ( AZ, FZ, EM, HM, EW, HW )
   
   END SUBROUTINE MTSPHERE_JMP_MultipleSource
   
    END SUBROUTINE MTSPHERE_HNK_MultipleSources
 
SUBROUTINE MTMultipleSourceCoefficients(NTERMS,LMBDA,YHAT,ZHAT,RTM,RTE,AZ,FZ)

    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NTERMS, N, M, NP, ID, IC, NM
    REAL(KIND=QL) LMBDA, LMBDASQ
    COMPLEX(KIND=QL) YHAT,ZHAT,K,U, &                      
                     AZ(2,2,NTERMS,-2:2), FZ(2,2,NTERMS,-2:2), &
                     RTM(NTERMS), RTE(NTERMS), AZS, FZS  
    COMPLEX(KIND=QL), ALLOCATABLE :: PSIA(:,:), PSIF(:,:), EZNM(:,:), HZNM(:,:), ZETA(:,:,:)  
    ALLOCATE ( PSIA(NTERMS, -NTERMS:NTERMS), PSIF(NTERMS, -NTERMS:NTERMS), &
                     ZETA(2,0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                     EZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                     HZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1) )
    
    AZ = (0._QL,0._QL)
    FZ = (0._QL,0._QL)
    K = SQRT ( - YHAT * ZHAT )    
    LMBDASQ = LMBDA * LMBDA
    U = SQRT ( LMBDASQ + YHAT * ZHAT )
    CALL SPHERICALCYLINDRICALCONVERSIONCoefficients(NTERMS+1,K,LMBDA,U,ZETA)
    DO IC = 1, 2
        do n = 1, NTERMS
            do M = -1, 1, 2
                PSIA = (0._QL,0._QL)
                psif = (0._QL,0._QL)
                if ( IC == 1 ) then
                    PSIA(n,M) = RTM(n)
                else
                    psif(n,M) = RTE(n)
                end if
                CALL VerticalFieldSphericalCoefficients(NTERMS,PSIA,PSIF,YHAT,ZHAT,EZNM,HZNM)
                DO ID = 1, 2
                    AZS = (0._QL,0._QL)
                    FZS = (0._QL,0._QL)
                    do NP = 1, NTERMS+1
                        AZS = AZS + ZETA(ID,NP,M) * EZNM(NP,M)
                        FZS = FZS + ZETA(ID,NP,M) * HZNM(NP,M)
                    end do
                    AZ(ID,IC,N,M) = AZS
                    FZ(ID,IC,N,M) = FZS
                end do
            end DO
        end do
    END DO

    AZ = AZ * YHAT / LMBDASQ
    FZ = FZ * ZHAT / LMBDASQ 
    
    DEALLOCATE ( PSIA, PSIF, EZNM, HZNM, ZETA )
    
END SUBROUTINE MTMultipleSourceCoefficients


SUBROUTINE MTMultipleSourceFields(NTERMS,NLYR,DPTHL,RLYR,LMBDA,Z,YHAT,ZHAT,AZ,FZ,AN,FN,EM,HM)

! Computation of the spherical coefficients for cylindrical fields

!           Input
!           -----
!   NTERMS - maximum spherical harmonic degree number
!   NLYR - number of layers
!   DEPTH - layer top depths
!   RXLYR - layer number of the receiver
!   ZS - depth of the source
!   ZR - depth of the receiver
!   LMBDA - cylindrical wavenumber
!   YHAT - layer admittances
!   ZHAT - layer impedances
!   AZ - Transverse magnetic vertical potential
!   FZ - Transverse electric vertical potential
!   AN      - TM layered earth propagation
!   FN      - TE layered earth propagation

!           Output
!           ------
!   EM    - Vertical electrical field Hankel transform spherical coefficients
!   HM    - Vertical magnetic field Hankel transform spherical coefficients

!   In the various vectors, ID is the propagation direction, 1 : Downward, 2 : Upward

use precision

    IMPLICIT NONE
    INTEGER NTERMS,NLYR,RLYR,ID,IC,N,M,IL
    REAL(KIND=QL)    LMBDA, lmbdasq, DPTHL(NLYR), Z
    COMPLEX(KIND=QL) EM(2,NTERMS,-1:1,-2:2,3), DZ,  &
                     HM(2,NTERMS,-1:1,-2:2,3), &
                     AZ(2,2,NTERMS,-2:2), &
                     FZ(2,2,NTERMS,-2:2), &
                     YHAT, ZHAT, U, DECAY, RTM, RTE, AN(2,0:NLYR), FN(2,0:NLYR), &
                     AZID(-1:1), FZID(-1:1), AZDX(-2:2), AZDY(-2:2), FZDX(-2:2), FZDY(-2:2)
        
    LMBDASQ = LMBDA ** 2
    U = SQRT ( LMBDASQ + YHAT * ZHAT ) 
    EM = ZERO
    HM = ZERO
    DO ID = 1, 2
        if ( RLYR == 0 .AND. ID == 1 ) CYCLE
        if ( RLYR == NLYR .AND. ID == 2 ) CYCLE
        if ( ID == 1 ) THEN
            DECAY = EXP ( - U * ( Z - DPTHL(RLYR) ) )
        ELSE if ( RLYR == 0 ) then
                DECAY = EXP ( U * ( Z - DPTHL(1) ) )
        else
                DECAY = EXP ( U * ( Z - DPTHL(RLYR) ) )
        end if
        RTM = AN(ID,RLYR) * DECAY
        RTE = FN(ID,RLYR) * DECAY
        if ( ID == 1 ) then
            DZ = - U
        else
            DZ = U
        end IF
        do IC = 1, 2
            do n = 1, NTERMS
                do M = -1, 1, 2
                    AZID = (0._QL,0._QL)
                    FZID = (0._QL,0._QL)
                    AZID(M) = AZ(ID,IC,N,M)
                    FZID(M) = FZ(ID,IC,N,M)
                    call cylindricalderivatives(1, 1, LMBDA, AZID, AZDX)
                    call Cylindricalderivatives(1, 2, LMBDA, AZID, AZDY)
                    call cylindricalderivatives(1, 1, LMBDA, FZID, FZDX)
                    call cylindricalderivatives(1, 2, LMBDA, FZID, FZDY)
                    EM(IC,N,M,:,1) = EM(IC,N,M,:,1) + RTM * AZDX     * DZ  / yhat
                    EM(IC,N,M,:,2) = EM(IC,N,M,:,2) + RTM * AZDY     * DZ  / yhat
                    EM(IC,N,M,M,3) = EM(IC,N,M,M,3) + RTM * AZ(ID,IC,N,M) * LMBDASQ / yhat
                    HM(IC,N,M,:,1) = HM(IC,N,M,:,1) + RTM * AZDY
                    HM(IC,N,M,:,2) = HM(IC,N,M,:,2) - RTM * AZDX
                    EM(IC,N,M,:,1) = EM(IC,N,M,:,1) - RTE * FZDY
                    EM(IC,N,M,:,2) = EM(IC,N,M,:,2) + RTE * FZDX
                    HM(IC,N,M,:,1) = HM(IC,N,M,:,1) + RTE * FZDX     * DZ  / zhat
                    HM(IC,N,M,:,2) = HM(IC,N,M,:,2) + RTE * FZDY     * DZ  / zhat
                    HM(IC,N,M,M,3) = HM(IC,N,M,M,3) + RTE * FZ(ID,IC,N,M) * LMBDASQ / zhat
                end do
            end do
        end do
    end do

END SUBROUTINE MTMultipleSourceFields      
      
SUBROUTINE MTSPHERE_SingleSource(NTERMS,PSIA,PSIF,ZS,NX,NY,PS,NLYR,THKD,YHAT,ZHAT,E,H)
    
!  Compute the electric and magnetic fields associated to a spherical source embedded in a 
!  plane wave field (maximum absolute order number is one).
! Computation of the vertical and electric magnetic fields from spherical potentials.
! Transformation to TE and TM vertical potentials
! Computation of cylindrical coefficients
! Transmission in a layered earth before surface field computation

!                             INPUT
!                             -----
!    NTERMS - maximum number of degrees.
!      PSIA - TM spherical coefficients.
!      PSIF - TE spherical coefficients.
!        ZS - Sphere center elevation (positive downward).
!      NX   - number of X positions.
!      NY   - number of Y positions.
!        PS - measurement point matrix.
!     TXA90 - true for vertical co-planar briadside array
!      ZRXD - vertical receiver offset (QL) from transmitter, (below = +)
!      XRXD - in-line horizontal offset (QL)  of RX J;        (behind = +)
!      YRXD - transverse horizontal offset (QL) of RX J       (left = +)
!      NLYR - number of layers
!      THKD - vector of layer thicknesses
!      YHAT - vector of layer admittivity
!      ZHAT - vector of layer impedivity
!
!                             OUTPUT
!                             ------
!         E - electric field vector at the measurement points.
!         H - magnetic field vector at the measurement points.

USE FILT_COEF_Q

 IMPLICIT NONE

 REAL, PARAMETER :: TWOPI = 2._QL * PI
 INTEGER NTERMS,NLYR,JZ,IX,IY,IXY,NX,NY,NXY,SXLYR,RXLYR,OPTION,I,M,IC
 REAL W
 REAL(KIND=QL) DEL_JN,THKD(NLYR-1),DPTHL(NLYR),ZS,ZR,EPR,EPI,HPR,HPI,RHO,PHI
 REAL(KIND=QL) PS(NX,NY,3)   
 COMPLEX(KIND=QL) E(NX,NY,3), H(NX,NY,3), &
                  YHAT(0:NLYR), zhat(0:NLYR), EMPHI, &
                  PSIA(NTERMS,-1:1), PSIF(NTERMS,-1:1), EH(-2:2,3), HH(-2:2,3)
                  
 E = ZERO
 H = ZERO
 DEL_JN = LOG (10.D0) / DBLE (NDEC_JN)
 
SXLYR = 0                ! Identify layer containing loop or GW source
ZR = 0
DPTHL = 0._QL
RXLYR = 0
IF ( NLYR > 0 ) THEN
    DO JZ = 1, NLYR-1
    DPTHL(JZ+1) = DPTHL(JZ) + THKD(JZ)
    END DO
    DO JZ = NLYR,1,-1
    IF (ZS > DPTHL(JZ)) THEN
        SXLYR = JZ
        EXIT
    END IF
    END DO
END IF
     
    DO IX = 1, NX
        do IY = 1, NY
            RHO = SQRT (PS(IX,IY,1)**2 + PS(IX,IY,2)**2)
            IF ( RHO < 0.01_QL  ) THEN 
                RHO = 0.01_QL
                PHI = 0._QL
            ELSE
                PHI = ATAN2 ( PS(IX,IY,2), PS(IX,IY,1) )
            END IF
            CALL MTSPHERE_HNK_SingleSource(NTERMS,SXLYR,RXLYR,PSIA,PSIF, &
                        NLYR,THKD,DPTHL,YHAT,ZHAT,ZS,ZR,RHO,EH,HH)
            DO M = -2, 2
                DO IC = 1, 3
                    EMPHI = EXP ( CI * M * PHI )
                    E(IX,IY,IC) = E(IX,IY,IC) + EH(M,IC) * EMPHI
                    H(IX,IY,IC) = H(IX,IY,IC) + HH(M,IC) * EMPHI
                END DO
            end DO
        END DO
    END DO
    
END SUBROUTINE MTSPHERE_SingleSource
    
SUBROUTINE MTSPHERE_HNK_SingleSource (NTERMS,SXLYR,RXLYR,PSIA,PSIF, &
                        NLYR,THKD,DPTHL,YHAT,ZHAT,ZS,ZR,RHO,Eh,Hh)
!-------------------------------------------------------------------------------------------

!***  Called by: HSPHERE_FD

!  Computes transform integrals HLYR(3) which are used to compute vertical
!  and horizontal frequency-domain magnetic field components at the RX from
!  VMD and HMD sources.  It evaluates the Hankel transform integral using a
!  15 points per decade filter coefficient set derived from Christensen's
!  FLTGEN program.

!                             INPUT
!   TERMS - maximum number of spherical degrees
!      IW - iw  angular frequency *(0.,1.)
!    RHOD - horizontal TX -> RX distance.
!     KER - stores kernel values from HSMD_KER    
!    
!

 USE FILT_COEF_Q

 IMPLICIT NONE
 INTEGER NTERMS,NLYR,I,J,SXLYR,RXLYR,NLAGS, M, INT
 REAL(KIND=QL) DEL_JN,Y,RHO_JN,LMBDA,RHO,ZRFD,ZTFD,THKD(NLYR-1),DPTHL(NLYR),ZS,ZR
 COMPLEX(KIND=QL) Eh(-2:2,3),Hh(-2:2,3), &
    RTM,RTE,P,YHAT(0:NLYR),ZHAT(0:NLYR), & 
    PSIA(NTERMS,-1:1), PSIF(NTERMS,-1:1), &
    FN(2,0:NLYR), AN(2,0:NLYR)
 LOGICAL JUMP

 DEL_JN = LOG (10.D0) / DBLE (NDEC_JN)
 RHO_JN = -LOG (RHO) - SHFTJN

 Eh = (0._QL, 0._QL)
 Hh = (0._QL, 0._QL)
 
 DO I = -50, JNHI             ! Start at I = -50 to pick up low values.
   Y = RHO_JN + DBLE(I) * DEL_JN
   LMBDA = EXP (Y)
   IF ( NLYR > 0 ) THEN
        CALL HSSPHERE_KER (LMBDA,SXLYR,NLYR,THKD,DPTHL,ZS,YHAT,ZHAT,AN,FN)
   END IF
   CALL MTSPHERE_JMP_SingleSource
   IF (JUMP .AND. I > -40) EXIT
 END DO

 JUMP = .FALSE.           ! Finish off the low end 
 DO I = -51, JNLO, -1
   Y = RHO_JN + DBLE(I) * DEL_JN
   LMBDA = EXP (Y)
   IF ( NLYR > 0 ) THEN
        CALL HSSPHERE_KER (LMBDA,SXLYR,NLYR,THKD,DPTHL,ZS,YHAT,ZHAT,AN,FN)
   END IF
   CALL MTSPHERE_JMP_SingleSource
   IF (JUMP .AND. I < -60) EXIT
 END DO
 
Eh = Eh / ( RHO * 4._QL * PI )
Hh = Hh / ( RHO * 4._QL * PI )

    CONTAINS
    
   SUBROUTINE MTSPHERE_JMP_SingleSource
!  ---------------------

!***  Called by: HSMD_HNK

!  Accumulates function calls for the Hankel transformation &
!  checks convergence.
   
     REAL(KIND=QL), PARAMETER :: TOL=1.D-6, TOL2=1.D-35, EXP_TOL=80.D0
     COMPLEX(KIND=QL), PARAMETER :: ZERO = CMPLX (0.D0, 0.D0, KIND=QL)
     INTEGER JINT,N,M,IC
     REAL(KIND=QL) QR,QI  
     COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL)
     REAL(KIND=QL) WJ(0:2)
     COMPLEX(KIND=QL) emphi, ESX, ERX
     COMPLEX(KIND=QL) AZ (2,-2:2), FZ(2,-2:2), EM(3,-2:2), HM(3,-2:2), EW(-2:2,3), HW(-2:2,3)
      
    ! Calcul des champs

    CALL MTSingleSourceCoefficients(NTERMS,PSIA,PSIF,LMBDA,YHAT(SXLYR),ZHAT(SXLYR),AZ,FZ)
    CALL MTSingleSourceFields(NLYR,DPTHL,RXLYR,ZR,LMBDA,yhat(RXLYR),zhat(RXLYR), &
         AZ,FZ,AN,FN,EM,HM)
    EW = (0._QL, 0._QL)
    HW = (0._QL, 0._QL)
    WJ = 0._QL
    WJ(0) = WJ0(I)
    WJ(1) = WJ1(I)
    WJ(2) = 2._QL * WJ(1) / ( RHO * LMBDA ) - WJ(0)
    DO M = -2, 2
        DO ic = 1, 3
            EW(M,ic) = EW(M,IC) + EM(IC,M) * WJ(ABS(M))
            HW(M,ic) = HW(M,IC) + HM(IC,M) * WJ(ABS(M))
        END DO
    END DO
    
    EW = EW * LMBDA 
    HW = HW * LMBDA  
        
    Eh = Eh + EW
    Hh = Hh + HW
    
    JUMP = .TRUE.
    DO M = -2, 2
        DO JINT = 1,3
            QR = ABS (REAL  (Eh(M,JINT), KIND=QL) )
            QI = ABS (AIMAG (Eh(M,JINT) ) )
            IF (QR > TOL2 .AND. ABS (REAL  (EW(M,JINT))) > TOL * QR) JUMP = .FALSE.
            IF (QI > TOL2 .AND. ABS (AIMAG (EW(M,JINT))) > TOL * QI) JUMP = .FALSE.
            QR = ABS (REAL  (Hh(M,JINT), KIND=QL) )
            QI = ABS (AIMAG (Hh(M,JINT) ) )
            IF (QR > TOL2 .AND. ABS (REAL  (HW(M,JINT))) > TOL * QR) JUMP = .FALSE.
            IF (QI > TOL2 .AND. ABS (AIMAG (HW(M,JINT))) > TOL * QI) JUMP = .FALSE.
        END DO
    END DO
       
END SUBROUTINE MTSPHERE_JMP_SingleSource

END SUBROUTINE MTSPHERE_HNK_SingleSource

SUBROUTINE MTSingleSourceCoefficients(NTERMS,PSIA,PSIF,LMBDA,YHAT,ZHAT,AZ,FZ)

    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NTERMS, N, M, NP, ID, IC
    REAL(KIND=QL) LMBDA, LMBDASQ
    COMPLEX(KIND=QL) YHAT,ZHAT,K,U, &
                     PSIA(NTERMS, -1:1), PSIF(NTERMS, -1:1), &
                     PSIAC(NTERMS, -NTERMS:NTERMS), PSIFC(NTERMS, -NTERMS:NTERMS), &
                     ZETA(2,0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                     EZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                     HZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1), & 
                     AZ(2,-2:2), &
                     FZ(2,-2:2)  
    
    K = SQRT ( - YHAT * ZHAT )    
    LMBDASQ = LMBDA * LMBDA
    U = SQRT ( LMBDASQ + YHAT * ZHAT )
    PSIAC = (0._QL,0._QL)
    PSIFC = (0._QL,0._QL)
    do M = -1, 1, 2
        PSIAC(:,M) = PSIA(:,M)
        PSIFC(:,M) = PSIF(:,M)
    end DO
    CALL VerticalFieldSphericalCoefficients(NTERMS,PSIAC,PSIFC,YHAT,ZHAT,EZNM,HZNM)    
    CALL SPHERICALCYLINDRICALCONVERSIONCoefficients(NTERMS+1,K,LMBDA,U,ZETA)
    AZ = (0._QL,0._QL)
    FZ = (0._QL,0._QL)
    DO ID = 1, 2
        do M = -1, 1, 2
            do N = 1, NTERMS+1
                AZ(ID,M) = AZ(ID,M) + ZETA(ID,N,M) * EZNM(N,M)
                FZ(ID,M) = FZ(ID,M) + ZETA(ID,N,M) * HZNM(N,M)
            end do
        end DO
    end do
  
    AZ = AZ * YHAT / LMBDASQ
    FZ = FZ * ZHAT / LMBDASQ 
    
    END SUBROUTINE MTSingleSourceCoefficients

SUBROUTINE MTSingleSourceFields(NLYR,DPTHL,RLYR,ZR,LMBDA,YHAT,ZHAT,AZ,FZ,AN,FN,EM,HM)

! Computation of the spherical coefficients for cylindrical fields

!           Input
!           -----
!   NTERMS - maximum spherical harmonic degree number
!   NLYR - number of layers
!   DEPTH - layer top depths
!   RXLYR - layer number of the receiver
!   ZS - depth of the source
!   ZR - depth of the receiver
!   LMBDA - cylindrical wavenumber
!   YHAT - layer admittances
!   ZHAT - layer impedances
!   AZ - Transverse magnetic vertical potential
!   FZ - Transverse electric vertical potential
!   AN      - TM layered earth propagation
!   FN      - TE layered earth propagation

!           Output
!           ------
!   ET      - Vertical electrical field Hankel transform spherical coefficients
!   EZNM    - Vertical magnetic field Hankel transform spherical coefficients

!   In the various vectors, ID is the propagation direction, 1 : Downward, 2 : Upward

    use PRECISION

    IMPLICIT NONE
    INTEGER NLYR,RLYR,ID,IC,N,M
    REAL(KIND=QL)    Il, LMBDA, lmbdasq, DPTHL(NLYR), ZR
    COMPLEX(KIND=QL) EM(3,-2:2), HM(3,-2:2), &
                     AZ(2,-2:2), FZ(2,-2:2), &
                     yhat, ZHAT, U, DZ(2), DECAY, RTM, RTE, &
                     AN(2,0:NLYR), FN(2,0:NLYR), &
                     AZID(-2:2), FZID(-2:2), &
                     AZDX(-2:2), AZDY(-2:2), FZDX(-2:2), FZDY(-2:2)
        
    LMBDASQ = LMBDA ** 2
    U = SQRT ( LMBDASQ + YHAT * ZHAT ) 
    DZ(1) = - U
    DZ(2) =   U
    EM = (0._QL,0._QL)
    HM = (0._QL,0._QL)
    DO ID = 1, 2
        if ( ( ID == 1 .AND. RLYR == 0 ) .OR. ( ID == 2 .AND. RLYR == NLYR ) ) cycle
        if ( ID == 1 ) THEN
            DECAY = EXP ( - U * ( ZR - DPTHL(RLYR) ) )
        ELSE IF ( RLYR == 0 ) THEN
            DECAY = EXP (   U * ( ZR - DPTHL(1) ) )
        ELSE
            DECAY = EXP (   U * ( ZR - DPTHL(RLYR) ) )
        end IF
        RTM = AN(ID,RLYR) * DECAY
        RTE = FN(ID,RLYR) * DECAY
        AZID = ZERO
        FZID = ZERO
        do M = -1, 1, 2
            AZID(M) = AZ(ID,M)
            FZID(M) = FZ(ID,M)
        end do
        call CylindricalDerivatives(1,1,LMBDA,AZID(-1:1),AZDX)
        call CylindricalDerivatives(1,2,LMBDA,AZID(-1:1),AZDY)
        call CylindricalDerivatives(1,1,LMBDA,FZID(-1:1),FZDX)
        call CylindricalDerivatives(1,2,LMBDA,FZID(-1:1),FZDY)
        EM(1,:) = EM(1,:) + RTM * AZDX     * DZ(ID)  / yhat
        EM(2,:) = EM(2,:) + RTM * AZDY     * DZ(ID)  / yhat
        EM(3,:) = EM(3,:) + RTM * AZID     * LMBDASQ / yhat
        HM(1,:) = HM(1,:) + RTM * AZDY
        HM(2,:) = HM(2,:) - RTM * AZDX
        EM(1,:) = EM(1,:) - RTE * FZDY
        EM(2,:) = EM(2,:) + RTE * FZDX
        HM(1,:) = HM(1,:) + RTE * FZDX     * DZ(ID)  / zhat
        HM(2,:) = HM(2,:) + RTE * FZDY     * DZ(ID)  / zhat
        HM(3,:) = HM(3,:) + RTE * FZID     * LMBDASQ / zhat
    END DO

END SUBROUTINE MTSingleSourceFields   



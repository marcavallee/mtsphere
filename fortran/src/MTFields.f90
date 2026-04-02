SUBROUTINE PlaneWaveImpedance(NW,NLYR,THK,RES,NF,FREQ,ZHAT,E)
    
! Subroutine to compute the coefficients of a plane wave propagating
! in a layered earth
    
!          Input
!          -----
! NW: printing unit
! NLYR: number of layers
! THK: layer thicknesses
! RES: layer resistivities
! NF: number of frequencies
! FREQ: list of frequencies

!          Output
!          -----
! ZHAT: surface impedances
! E : electric field coefficients
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    REAL(KIND=QL), PARAMETER :: MU0=12.56637D-7,EPS0=8.854156D-12
    REAL(KIND=QL), PARAMETER :: PI = 3.141592653589793    
    COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL)
    INTEGER JL,NLYR,JF,NF,NW
    REAL FREQ(NF),THK(NLYR-1),RES(NLYR), W
    COMPLEX(KIND=QL) K(0:NLYR),Z(0:NLYR),ZHAT(NF,NLYR),EI,E(NF,0:NLYR,2),A,B,DEN
    
    E = (0._QL,0._QL)
    DO JF = 1, NF
        W = 2._QL * PI * FREQ(JF) 
        DO JL = 0,NLYR
            if ( JL == 0 ) then
                K(JL) = w * SQRT ( mu0 * eps0 )
            else
                K(JL) = SQRT ( - CI * MU0 * W / RES(JL) )
            end if
            Z(JL) = W * MU0 / K(JL)
        END DO
        ZHAT(JF,NLYR) = Z(NLYR)
        DO JL = NLYR-1,1,-1
            A = ZHAT(JF,JL+1)+Z(JL)*TANH(CI*K(JL)*THK(JL))
            B = Z(JL)+ZHAT(JF,JL+1)*TANH(CI*K(JL)*THK(JL))
            ZHAT(JF,JL) = Z(JL)*A/B                              
        END DO
        DO JL = 0,NLYR
            if ( JL == 0 ) then
                E(JF,JL,1) = 1._QL
                E(JF,JL,2) = ( ZHAT(JF,1) - Z(0) ) / ( ZHAT(JF,1) + Z(0) )
            ELSE 
            	IF ( JL == 1 ) THEN
               		EI = E(JF,0,1) + E(JF,0,2)
               	ELSE
               		EI = E(JF,JL-1,1) * EXP ( - CI * K(JL-1) * THK(JL-1) ) + &
                         E(JF,JL-1,2) * EXP (   CI * K(JL-1) * THK(JL-1) ) 
                END IF
                if ( JL == NLYR ) then
                	E(JF,JL,1) = EI
                	E(JF,JL,2) = (0._QL,0._QL)
                ELSE
                	DEN = ( ZHAT(JF,JL+1) - Z(JL) ) * EXP ( - CI * 2._QL * K(JL) * THK(JL) ) + &
                		  ( ZHAT(JF,JL+1) + Z(JL) )
                	E(JF,JL,1) = ( ZHAT(JF,JL+1) + Z(JL) ) * EI / DEN
                	E(JF,JL,2) = ( ZHAT(JF,JL+1) - Z(JL) ) * EXP(- CI * 2._QL * K(JL)* THK(JL) ) * EI / DEN
                END IF
             END IF	
        END DO
    END DO
    
    WRITE(NW,'(/'' Layer impedance'')')
    DO JF = 1, NF
        DO JL = 1, NLYR
            WRITE(NW,'('' Frequency:'',G15.7,'', Layer:'',I5,'', Impedance:'',2G15.7)') &
                FREQ(JF),JL,ZHAT(JF,JL)
        END DO
    END DO
    
    WRITE(NW,'(/'' Propagation coefficients'')')
    DO JF = 1, NF
        DO JL = 0, NLYR
             WRITE(NW,100)FREQ(JF),JL,E(JF,JL,1),E(JF,JL,2)      
        END DO
    END DO
    
100 FORMAT(' Frequency:',G15.7,', Layer:',I5,', Earthward propagation:',2G15.7,', Sunward propagation:',2G15.7)
    
END SUBROUTINE PlaneWaveImpedance 
    
SUBROUTINE APPARENTRESISTIVITY(NF, NX, NY, FREQ, IMP, APPRES, PHASE)

! Subroutine to compute the apparent resistivity and phase

!          Input
!          -----
! NF: number of frequencies
! NX: number of x locations
! NY: number of y locations
! IMP: ground surface impedance

!          Output
!
! APPRES: Apparent resistivity
! PHASE: phase

IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    REAL(KIND=QL), PARAMETER :: MU0=12.56637D-7
    REAL(KIND=QL), PARAMETER :: PI = 3.141592653589793    
    INTEGER NF, NX, NY, JF, JX, JY, I, J
    REAL FREQ(NF)
    REAL(KIND=QL) APPRES(NF,NX,NY,2,2), PHASE(NF,NX,NY,2,2), DEN
    COMPLEX(KIND=QL) IMP(NF,NX,NY,3,2), TEST 
    
    DO JF = 1, NF
        DEN = 2._QL * PI * FREQ(JF) * MU0
        DO JX = 1, NX
            DO JY = 1, NY
                DO I = 1, 2
                    DO J = 1, 2
                        APPRES(JF,JX,JY,I,J) =  ABS ( IMP(JF,JX,JY,I,J) ) ** 2 / DEN
                        PHASE(JF,JX,JY,I,J) = ATAN2 ( AIMAG( IMP(JF,JX,JY,I,J) ), REAL( IMP(JF,JX,JY,I,J) )  ) * 180 / PI
                    END DO
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE APPARENTRESISTIVITY

SUBROUTINE PLANEWAVESPHERICALHARMONICS(NTERMS,YHAT,ZHAT,H,E,PSIA,PSIF)
    
! Subroutine to compute the spherical coefficients of the upward and downward 
! propagation of a sphere located in a layered earth

!          Input
!          -----
! NTERMS: number of spherical coefficients.
! YHAT: admittivity
! ZHAT: impeditivity
! H : elevation of the source
! E : source electric field

!          Output
!          -----
! PSIA: transverse magnetic radial coefficients
! PSIF: transverse electric radial coefficients

    use sphericalfunctions

    implicit none
    
    INTEGER NTERMS,IC,SIC,N,M,ID
    REAL(KIND=QL) FACTNM,H
    COMPLEX(KIND=QL) YHAT,ZHAT,K,E(2), &
                    PSIA(2,NTERMS,-1:1), &
                    PSIF(2,NTERMS,-1:1), E0, &
                    psiat(nterms,-1:1), &
                    psift(nterms,-1:1), IMPEDANCE
     
    PSIA = (0._QL,0._QL)
    PSIF = (0._QL,0._QL)
    K = SQRT ( - YHAT * ZHAT )
    DO IC = 1, 2 ! Downward and upward orientation
        IF ( E(IC) == 0 ) CYCLE
        IF ( IC == 1 ) THEN
            SIC = -1
        ELSE  
            SIC = 1
        END IF
        E0 = E(IC) * EXP ( SIC * CI * K * H )
        do ID = 1, 2
            call PlaneWaveFieldssphericalharmonics(e0,nterms,IC,id,K,yhat,zhat,psiat,psift)
            PSIA(Id,:,:) = PSIA(Id,:,:) + psiat
            PSIF(Id,:,:) = PSIF(Id,:,:) + psift
        end do
    END DO

END SUBROUTINE PLANEWAVESPHERICALHARMONICS
    
subroutine planewavefieldssphericalharmonics(E0,nterms,ic,id,k,yhat,zhat,psia,psif)

! Kernel of the representation of a plane wave as spherical expansion

!           Input
!           _____
!   E0 : Incident electric field.

    implicit none
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    REAL(KIND=QL), PARAMETER :: PI = 3.141592653589793    
    COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL)
    integer nterms,ic,ID,n,m,sic
    COMPLEX(KIND=QL) PSIA(NTERMS,-1:1), K, &
                     PSIF(NTERMS,-1:1), E0, AN, &
                     AFACT, FFACT, yhat, zhat, IMPEDANCE
    
    PSIA = (0._QL,0._QL)
    PSIF = (0._QL,0._QL)
    impedance = SQRT(zhat/yhat)
    IF ( IC == 1 ) THEN
        SIC = -1
    ELSE  
        SIC = 1
    END IF
    do n = 1, NTERMS
        do m = -1, 1, 2
            AN     = E0 * (SIC * CI ) ** N * SQRT(PI*(2*N+1)/N/(N+1)) / k
            IF ( ID == 1 ) THEN
                AFACT = - SIC * CI
                FFACT = M
            ELSE
                AFACT = - SIC * M
                FFACT = - CI
            END IF
            PSIA(N,M) = AFACT * AN * YHAT
            PSIF(N,M) = FFACT * AN * ZHAT / IMPEDANCE 
        END DO
    end do   

end subroutine planewavefieldssphericalharmonics    

SUBROUTINE MTSPHERE3D(NW, NF, NLYR, NTERMS, NX, NY, FREQ, PXY, THK, RES, DEPTH, &
                      RADIUS, SRES, ZHAT, E, ES, HS)

! Main subroutine to compute the sphere response

!          Input
!          -----
! NW: printing unit
! NF: number of frequencies
! NLYR: number of layers
! NTERMS: number of terms of the spherical computation
! NX: number of X locations
! NY: number of Y locations
! FREQ: vector of frequencies
! PXY: array of locations (NX,NY)
! THK: number of thicknesses
! RES: number of resistivities
! DEPTH: depth of the sphere
! RADIUS: radius of the sphere
! SRES: sphere resistivity
! ZHAT: layer surface admittances
! E: electric field coefficients

!          Output
!          -----
! IMPEDANCE: ground surface impedance
! ET: total electric field
! HT: total magnetic field

    USE PRECISION

    INTEGER NW, NF, NLYR, JF, JX, JY, JZ, SLYR, N, M, NP, MP, NTERMS, NX, NY, ID, I, TP, IP
    REAL THK(NLYR-1), RES(NLYR), DEPTH, FREQ(NF), SRES, RADIUS
    REAL(KIND=QL) THKD(NLYR-1),DPTHL(NLYR),ZS,DEPTHD, W, RADIUSD, PXY(NX,NY,3)
    COMPLEX(KIND=QL) E(NF,0:NLYR,2),EL(2),YHATL(0:NLYR),ZHATL(0:NLYR),YHATS,ZHATS, &
                     PSIAI(2,NTERMS,-1:1), & 
                     PSIFI(2,NTERMS,-1:1), &
                     PSIAC(NTERMS,-1:1), & 
                     PSIFC(NTERMS,-1:1), &
                     PSIAR(NTERMS,-1:1), &
                     PSIFR(NTERMS,-1:1), &
                     ES(NF,NX,NY,3,2), HS(NF,NX,NY,3,2), &
                     RTM(NTERMS), RTE(NTERMS), ZHAT(NF,NLYR)
    COMPLEX(KIND=QL), ALLOCATABLE :: AI(:,:,:)
    
    THKD = 0._QL
    THKD = REAL (THK, KIND=QL)
    RADIUSD = REAL (RADIUS, KIND=QL)
    DEPTHD = REAL ( DEPTH, KIND=QL)
    ALLOCATE ( AI (-1:1, 2*NTERMS, 2*NTERMS ) )
    
     DPTHL = 0._QL
     DO JZ = 1, NLYR-1
       DPTHL(JZ+1) = DPTHL(JZ) + THKD(JZ)
     END DO
     ZS = REAL (DEPTH,QL)
     SLYR = 0                ! Identify layer containing loop or GW source
     DO JZ = NLYR,1,-1
        IF (ZS > DPTHL(JZ)) THEN
            SLYR = JZ
            EXIT
        END IF
     END DO
    
    DO JF = 1, NF
        WRITE(NW,'(/''Frequency:'',G15.7)')FREQ(JF)
        W = 2._QL * PI * FREQ(JF)
        DO JZ = 0, NLYR
            IF ( JZ == 0 ) THEN
                YHATL(JZ) = CI * W * EPS0
            ELSE
                YHATL(JZ) = 1 / RES(JZ) + CI * W * EPS0
            END IF
        END DO
        ZHATL = CI * W * MU0
        YHATS = 1 / SRES + CI * W * EPS0
        ZHATS = CI * W * MU0
        CALL REFLECTIONCOEFFICIENTS(NTERMS,RADIUSD,YHATL(SLYR),ZHATL(SLYR),YHATS,ZHATS,RTM,RTE)
        call MTLayeredEarthCorrection(NW,NTERMS,NLYR,THKD,ZS,RADIUSD,SLYR,YHATL,ZHATL,YHATS,ZHATS,RTM,RTE,AI)
        WRITE(NW,'(/'' Reflection coefficients:'')')
        DO N = 1, NTERMS
            WRITE(NW,'(I5,100G15.7)')N,RTM(N),RTE(N)
        END DO
        EL = E(JF,SLYR,:)
        CALL PLANEWAVESPHERICALHARMONICS(NTERMS,YHATL(SLYR),ZHATL(SLYR),ZS-DPTHL(SLYR),EL,PSIAI,PSIFI)
        DO ID = 1, 2 ! X and Y directions
            WRITE(NW,'(/'' Plane wave spherical harmonics:'')')
            IF ( ID == 1 ) THEN
                WRITE(NW,'('' X directed propagation'')')
            ELSE
                WRITE(NW,'('' Y directed propagation'')')
            END IF
            DO N = 1, NTERMS
                DO M = -1, 1, 2
                    WRITE(NW,'(2I5,10G15.7)')N,M,PSIAI(ID,N,M),PSIFI(ID,N,M)
                END DO
            END DO
            CALL MTSphereReflectionCORRECTION(NTERMS,AI,PSIAI(ID,:,:),PSIFI(ID,:,:),PSIAC,PSIFC)
            WRITE(NW,'(/'' Corrected potentials'')')
            DO N = 1, NTERMS
                DO M = -1, 1, 2
                    WRITE(NW,'(2I5,10G15.7)')N,M,PSIAC(N,M),PSIFC(N,M)
               END DO
            END DO
            PSIAR = ZERO
            PSIFR = ZERO
           DO N =  1, NTERMS
                do M = -1, 1, 2
                    PSIAR(N,M) = RTM(N) * PSIAC(N,M)
                    PSIFR(N,M) = RTE(N) * PSIFC(N,M)
                end do
            END DO
            CALL MTSPHERE_SingleSource(NTERMS,PSIAR,PSIFR,DEPTHD,NX,NY,PXY,NLYR,THKD,YHATL,ZHATL,ES(JF,:,:,:,ID),HS(JF,:,:,:,ID))
            WRITE(NW,'(/''Induced fields'')')
            DO JX = 1, NX
                DO JY = 1, NY
                    WRITE(NW,'(G15.7,2I5,20G15.7)')FREQ(JF),JX,JY,ES(JF,JX,JY,:,ID),HS(JF,JX,JY,:,ID)
                END DO
            END DO
        END DO
    end do
    
    DEALLOCATE ( AI )

    END SUBROUTINE MTSPHERE3D
    
    SUBROUTINE GETIMPEDANCE(NF,NX,NY,ZHAT,E,ES,HS,IMPEDANCE,ET,HT)
    
! Subroutine to compute the ground impedance and total fields from the
! primary and secondary fields

!          Input
!          -----
! NX: number of X locations
! NY: number of Y locations
! E: primary electric field
! ES: secondary electric field
! HS: secondary magnetic field

!          Output
!          -----
! IMPEDANCE: surface impedance
! ET: total electric field
! HT: total magnetic field
    
        IMPLICIT NONE
        INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
        REAL(KIND=QL), PARAMETER :: EPS0=8.854156D-12, MU0=12.56637D-7
        INTEGER NX, NY, IX, IY, IC, ID, I, J, NF, JF
        COMPLEX(KIND=QL) IMPEDANCE(NF,NX,NY,3,2), ZHAT(NF), E(NF,2), &
                         ES(NF,NX,NY,3,2), HS(NF,NX,NY,3,2), &
                         ET(NF,NX,NY,3,2), HT(NF,NX,NY,3,2), &
                         A(2), B(2), C(2), DEN, EP, HP, EE, HH
        ET = ES
        HT = HS
        do JF = 1, NF
            EE = E(JF,1) + E(JF,2)
            HH = EE / ZHAT(JF)
            do IX = 1, nx
                do iy = 1, ny
                    ET(JF,ix,iy,1,1) = ET(JF,ix,iy,1,1) + EE
                    ET(JF,ix,iy,2,2) = ET(JF,ix,iy,2,2) + EE
                    HT(JF,ix,iy,1,2) = HT(JF,ix,iy,1,2) - HH
                    HT(JF,ix,iy,2,1) = HT(JF,ix,iy,2,1) + HH
                end do
            end do
            DO IX = 1, NX
                DO IY = 1, NY
                    DO ID = 1, 2 ! X and Y direction
                        A(ID) = HT(JF,IX,IY,1,ID)
                        B(ID) = HT(JF,IX,IY,2,ID)
                    END DO
                    DEN = A(1) * B(2) - B(1) * A(2)
                    DO IC = 1, 3 ! X, Y 
                        IF ( IC < 3 ) THEN
                            DO ID = 1, 2
                                C(ID) = ET(JF,IX,IY,IC,ID)
                            END DO
                        ELSE
                            DO ID = 1, 2
                                C(ID) = HT(JF,IX,IY,3,ID)
                            END DO
                        END IF
                        IMPEDANCE(JF,IX,IY,IC,1) = ( C(1) * B(2) - B(1) * C(2) ) / DEN !  Impedance x
                        IMPEDANCE(JF,IX,IY,IC,2) = ( A(1) * C(2) - C(1) * A(2) ) / DEN !  Impedance y
                    END DO
                END DO
            END DO
        end DO
        
    END SUBROUTINE GETIMPEDANCE
    
SUBROUTINE MTLayeredEarthCorrection(NW,NTERMS,NLYR,THK,DEPTH,RADIUS,SLYR,YHATL,ZHATL,YHATS,ZHATS,RTM,RTE,AI)

use precision

! This subroutine computes the layered earth correction for the reflection
! of the sphere field on the sphere

!          Input
!          -----

! NW: unit for printing
! NTERMS: degree number
! NLYR: number of layers
! THK: layer thicknesses
! DEPTH: sphere center depth
! RADIUS: sphere radius
! SLYR: Sphere layer
! YHATL: Layer admittivity
! ZHATL: Layer impedivity
! YHATS: Sphere admittivity
! ZHATS: Sphere impedivity

!          Output
!          -----
! C: Correction tensor

   IMPLICIT NONE
   INTEGER NTERMS, NLYR, N, M, ICR, IC, IH, IHR, NLAT, NLON, &
       SLYR, NS, IR, INFO, NP, NW, JLAt, jlon, lwork
   REAL(KIND=QL) THK(NLYR-1), DEPTH, RADIUS, centre(3), TOL
   COMPLEX(KIND=QL) YHATL(0:NLYR), ZHATL(0:NLYR), YHATS, ZHATS, &
                    AI(-1:1,2*NTERMS,2*NTERMS), &
                    RTM(NTERMS), RTE(NTERMS), KR
   integer, allocatable :: IPIV(:)
   REAL(KIND=QL), ALLOCATABLE :: PS(:,:,:)
   COMPLEX(KIND=QL), ALLOCATABLE :: E(:,:,:,:,:,:), H(:,:,:,:,:,:), &
                 psiA(:,:), psif(:,:), A(:,:), C(:,:), &
                 WORK(:), ER(:,:), HR(:,:)
   
   NLAT = 4 * NTERMS + 1
   NLON = 4 * NTERMS
   
   centre = 0._QL
   centre(3) = depth
   allocate ( A(2*NTERMS,2*NTERMS), C(2*NTERMS,2*NTERMS), &
       PS(NLAT,NLON,3), ER(NLAT,NLON), HR(NLAT,NLON), psia(nterms,-1:1), psif(nterms,-1:1), &
       E(2,NTERMS,-1:1,NLAT,NLON,3), H(2,NTERMS,-1:1,NLAT,NLON,3), &
       IPIV(2*NTERMS), WORK(2*NTERMS) ) 
   
   call SphericalGeometryInitialisation(centre,RADIUS,nlat,nlon,ps) 
   
   CALL MTSphere_MultipleSources(NW,NTERMS,NLAT,NLON,depth,depth,NLYR,THK,RADIUS,YHATL,ZHATL,RTM,RTE,E,H) 
   
   WRITE(NW,'(/''Spherical Harmonic Analysis'')') 

   do M = -1, 1, 2
        IF ( M < 0 ) THEN
            WRITE(NW,'(/''Sphere retroaction, negative order'')')
        ELSE
            WRITE(NW,'(/''Sphere retroaction, positive order'')')
        END IF
        ICR = 0
        do IC = 1, 2
            do n = 1, NTERMS
                iCR = icr + 1
                Call RadialFields(NLAT,NLON,PS,CENTRE,E(IC,n,m,:,:,:),H(IC,n,m,:,:,:),ER,HR)
                Call MTSphericalHarmonicAnalysis(nterms,nlat,nlon,radius,YHATL(SLYR),zhatL(SLYR),er,hr,.FALSE.,psia,psif)
                ihr = 0
                do Ih = 1, 2
                   do np = 1, nterms
                       ihr = ihr + 1
                       if ( ih == 1 ) then
                           C(ihr,icr) = psia(NP,M)
                       else
                           C(ihr,icr) = psif(NP,M)
                       end if
                       WRITE(NW,'(4I5,20G15.7)')IC,N,ih,np,C(ihr,icr)
                    end do
                end do
           end do
       end do           
       DO IR = 1, 2 * NTERMS
            DO IH = 1, 2 * NTERMS
                IF ( IR == IH ) THEN
                    A(IR,IH) = 1._QL - C(IR,IH)
                ELSE
                    A(IR,IH) = - C(IR,IH)
                END IF
            END DO
        END DO  
        CALL ZGETRF(2*NTERMS,2*NTERMS,A,2*NTERMS,IPIV,INFO)
        CALL ZGETRI(2*NTERMS,A,2*NTERMS,IPIV,WORK,2*NTERMS,INFO)    
        AI(M,:,:) = A
        IF ( M < 0 ) THEN
            WRITE(NW,'(/''Sphere correction, negative order'')')
        ELSE
            WRITE(NW,'(/''Sphere correction, positive order'')')
        END IF
        ICR = 0
        do IC = 1, 2
            do n = 1, NTERMS
                iCR = icr + 1
                ihr = 0
                do Ih = 1, 2
                   do np = 1, nterms
                       ihr = ihr + 1
                       WRITE(NW,'(4I5,20G15.7)')IC,n,ih,NP,A(ihr,icr)
                   end do
                end do
            END DO
        END DO
    END DO
    DEALLOCATE ( PS, ER, HR, psia, psif, E, H, A, C, IPIV, WORK )

END SUBROUTINE MTLayeredEarthCorrection    
       
SUBROUTINE MTSphereReflectionCORRECTION(NTERMS,AI,PSIAI,PSIFI,PSIAC,PSIFC)

! This subroutine applies the sphere reflection correction on the 
! spherical potential coefficients

!          Input
!          -----

! NTERMS: degree number
! C: Correction matrix
! PSIAI: original transverse magnetic spherical potential
! PSIFI: original transverse electric spherical potential

!          Output
!          -----
! PSIAC: corrected transverse magnetic spherical potential
! PSIFC: corrected transverse electric spherical potential
   
    IMPLICIT NONE
        INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
        INTEGER NTERMS,  IC, N, M, I, J, IA, NA, MA
        COMPLEX(KIND=QL) PSIAI(NTERMS,-1:1), & 
                    PSIFI(NTERMS,-1:1), &
                    PSIAC(NTERMS,-1:1), &
                    PSIFC(NTERMS,-1:1), &
                    AI(-1:1, 2 * NTERMS, 2 * NTERMS )
        
        complex(kind=QL), allocatable :: B(:), X(:)
        
        allocate ( b(2 * NTERMS ), x(2 * NTERMS ) )
    DO M = -1, 1, 2
        I = 0              
        DO IC = 1, 2
            DO N = 1, NTERMS
                I = I + 1
                IF( IC == 1 ) THEN
                    b(I) = PSIAI(N,M)
                ELSE
                    b(I) = PSIFI(N,M)
                END IF
            END DO
        END DO
        X = MATMUL(AI(M,:,:),B)
        I = 0              
        DO IC = 1, 2
            DO N = 1, NTERMS
                I = I + 1
                IF( IC == 1 ) THEN
                    PSIAC(N,M) = X(I)
                ELSE
                    PSIFC(N,M) = X(I)
                END IF
            END DO
        END DO
    END DO
    deallocate ( B, X )
    
    END SUBROUTINE MTSphereReflectionCORRECTION   
    
   
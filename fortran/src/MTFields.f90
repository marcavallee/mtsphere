SUBROUTINE PlaneWaveImpedance(NW,NLYR,THK,RES,NF,FREQ,K,Z,E)
    
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
    COMPLEX(KIND=QL) K(NF,0:NLYR),Z(NF,0:NLYR),ZHAT(NLYR),EI,E(NF,0:NLYR,2),A,B,DEN
    
    E = (0._QL,0._QL)
    DO JF = 1, NF
        W = 2._QL * PI * FREQ(JF) 
        DO JL = 0,NLYR
            if ( JL == 0 ) then
                K(JF,JL) = w * SQRT ( mu0 * eps0 )
            else
                K(JF,JL) = SQRT ( - CI * MU0 * W / RES(JL) )
            end if
            Z(JF,JL) = W * MU0 / K(JF,JL)
        END DO
        ZHAT(NLYR) = Z(JF,NLYR)
        DO JL = NLYR-1,1,-1
            A = ZHAT(JL+1)+Z(JF,JL)*TANH(CI*K(JF,JL)*THK(JL))
            B = Z(JF,JL)+ZHAT(JL+1)*TANH(CI*K(Jf,JL)*THK(JL))
            ZHAT(JL) = Z(JF,JL)*A/B                              
        END DO
        DO JL = 0,NLYR
            if ( JL == 0 ) then
                E(JF,JL,1) = 1._QL
                E(JF,JL,2) = ( ZHAT(1) - Z(JF,0) ) / ( ZHAT(1) + Z(JF,0) )
            ELSE 
            	IF ( JL == 1 ) THEN
               		EI = E(JF,0,1) + E(JF,0,2)
               	ELSE
               		EI = E(JF,JL-1,1) * EXP ( - CI * K(JF,JL-1) * THK(JL-1) ) + &
                         E(JF,JL-1,2) * EXP (   CI * K(JF,JL-1) * THK(JL-1) ) 
                END IF
                if ( JL == NLYR ) then
                	E(JF,JL,1) = EI
                	E(JF,JL,2) = (0._QL,0._QL)
                ELSE
                	DEN = ( ZHAT(JL+1) - Z(JF,JL) ) * EXP ( - CI * 2._QL * K(JF,JL) * THK(JL) ) + &
                		  ( ZHAT(JL+1) + Z(JF,JL) )
                	E(JF,JL,1) = ( ZHAT(JL+1) + Z(JF,JL) ) * EI / DEN
                	E(JF,JL,2) = ( ZHAT(JL+1) - Z(JF,JL) ) * EXP(- CI * 2._QL * K(JF,JL)* THK(JL) ) * EI / DEN
                END IF
             END IF	
        END DO
    END DO
    
    WRITE(NW,'(/'' Layer impedance'')')
    DO JF = 1, NF
        DO JL = 0, NLYR
            WRITE(NW,'('' Frequency:'',G15.7,'', Layer:'',I5,'', Impedance:'',2G15.7)') &
                FREQ(JF),JL,Z(JF,JL)
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
    
SUBROUTINE APPARENTRESISTIVITY(NF, NX, NY, ND, FREQ, IMP, APPRES, PHASE)

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
    INTEGER NF, NX, NY, ND, JF, JX, JY, JD, I, J
    REAL FREQ(NF)
    REAL(KIND=QL) APPRES(NF,NX,NY,ND,2,2), PHASE(NF,NX,NY,ND,2,2), DEN
    COMPLEX(KIND=QL) IMP(NF,NX,NY,ND,3,2), TEST 
    
    DO JF = 1, NF
        DEN = 2._QL * PI * FREQ(JF) * MU0
        DO JX = 1, NX
            DO JY = 1, NY
                do JD = 1, ND
                    DO I = 1, 2
                        DO J = 1, 2
                            APPRES(JF,JX,JY,JD,I,J) =  ABS ( IMP(JF,JX,JY,JD,I,J) ) ** 2 / DEN
                            PHASE(JF,JX,JY,JD,I,J) = ATAN2 ( AIMAG( IMP(JF,JX,JY,JD,I,J) ), &
                                                              REAL( IMP(JF,JX,JY,JD,I,J) )  ) * 180 / PI
                        END DO
                    END DO
                end DO
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

SUBROUTINE MTSPHERE3D(NW, NF, NLYR, NTERMS, NX, NY, ND, FREQ, PXYD, THK, RES, DEPTH, &
                      RADIUS, SRES, E, DPTHL, ES, HS)

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
! PXYZ: array of locations (NX,NY,NZ)
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

    INTEGER NW, NF, NLYR, JF, JX, JY, JZ, SXLYR, N, M, NP, MP, NTERMS, NX, NY, ND, JD, I, TP, IP
    REAL THK(NLYR-1), RES(NLYR), DEPTH, FREQ(NF), SRES, RADIUS
    REAL(KIND=QL) THKD(NLYR-1),DPTHL(NLYR),ZS,W, RADIUSD, PXYD(NX,NY,ND,3), R
    COMPLEX(KIND=QL) E(NF,0:NLYR,2),EL(2),YHATL(0:NLYR),ZHATL(0:NLYR),YHATS,ZHATS, &
                     PSIAI(2,NTERMS,-1:1), & 
                     PSIFI(2,NTERMS,-1:1), &
                     PSIAC(NTERMS,-1:1), & 
                     PSIFC(NTERMS,-1:1), &
                     PSIAR(NTERMS,-1:1), &
                     PSIFR(NTERMS,-1:1), &
                     PSIAT(NTERMS,-1:1), &
                     PSIFT(NTERMS,-1:1), &
                     ES(NF,NX,NY,ND,3,2), HS(NF,NX,NY,ND,3,2), &
                     RTM(NTERMS), RTE(NTERMS), TTM(NTERMS), TTE(NTERMS)
    COMPLEX(KIND=QL), ALLOCATABLE :: AI(:,:,:)
    
    THKD = 0._QL
    THKD = REAL (THK, KIND=QL)
    RADIUSD = REAL (RADIUS, KIND=QL)
    DEPTHD = REAL ( DEPTH, KIND=QL)
    ALLOCATE ( AI (-1:1, 2*NTERMS, 2*NTERMS ) )
    
     ! Identify sphere layer
     ZS = REAL (DEPTH,QL)
     SXLYR = 0                
     DO JZ = NLYR,1,-1
        IF (ZS > DPTHL(JZ)) THEN
            SXLYR = JZ
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
        CALL REFLECTIONCOEFFICIENTS(NTERMS,RADIUSD,YHATL(SXLYR),ZHATL(SXLYR), &
            YHATS,ZHATS,RTM,RTE,TTM,TTE)
        call MTLayeredEarthCorrection(NW,NTERMS,NLYR,THKD,DPTHL,ZS,RADIUSD, &
            SXLYR,YHATL,ZHATL,YHATS,ZHATS,RTM,RTE,AI)
        WRITE(NW,'(/'' Reflection coefficients:'')')
        DO N = 1, NTERMS
            WRITE(NW,'(I5,100G15.7)')N,RTM(N),RTE(N)
        END DO
        EL = E(JF,SXLYR,:)
        CALL PLANEWAVESPHERICALHARMONICS(NTERMS,YHATL(SXLYR),ZHATL(SXLYR),ZS-DPTHL(SXLYR),EL,PSIAI,PSIFI)
        DO IP = 1, 2 ! X and Y directions
            WRITE(NW,'(/'' Plane wave spherical harmonics:'')')
            IF ( IP == 1 ) THEN
                WRITE(NW,'('' X directed propagation'')')
            ELSE
                WRITE(NW,'('' Y directed propagation'')')
            END IF
            DO N = 1, NTERMS
                DO M = -1, 1, 2
                    WRITE(NW,'(2I5,10G15.7)')N,M,PSIAI(IP,N,M),PSIFI(IP,N,M)
                END DO
            END DO
            CALL MTSphereReflectionCORRECTION(NTERMS,AI,PSIAI(IP,:,:),PSIFI(IP,:,:),PSIAC,PSIFC)
            WRITE(NW,'(/'' Corrected potentials'')')
            DO N = 1, NTERMS
                DO M = -1, 1, 2
                    WRITE(NW,'(2I5,10G15.7)')N,M,PSIAC(N,M),PSIFC(N,M)
               END DO
            END DO
            PSIAR = ZERO
            PSIFR = ZERO
            PSIAT = ZERO
            PSIFT = ZERO
            DO N =  1, NTERMS
                do M = -1, 1, 2
                    PSIAR(N,M) = RTM(N) * PSIAC(N,M)
                    PSIFR(N,M) = RTE(N) * PSIFC(N,M)
                    PSIAT(N,M) = TTM(N) * PSIAC(N,M)
                    PSIFT(N,M) = TTE(N) * PSIFC(N,M)
                end do
            END DO
            CALL MTSPHERE_SingleSource(NTERMS,PSIAR,PSIFR,PSIAT,PSIFT,NX,NY,ND,PXYD,NLYR, &
                THKD,YHATL,ZHATL,YHATS,ZHATS,DPTHL,ZS,RADIUS,ES(JF,:,:,:,:,IP),HS(JF,:,:,:,:,IP))
            call dynamiccleaning(ES(JF,:,:,:,:,IP))
            call dynamiccleaning(HS(JF,:,:,:,:,IP))
            WRITE(NW,'(/''Induced fields'')')
            DO JX = 1, NX
                DO JY = 1, NY
                    DO JD = 1, ND
                        WRITE(NW,'(G15.7,3I5,20G15.7)')FREQ(JF),JX,JY,JD,&
                            ES(JF,JX,JY,JD,:,IP),HS(JF,JX,JY,JD,:,IP)
                    end DO
                END DO
            END DO
        END DO
    end do
    
    DEALLOCATE ( AI )
    
    contains

    subroutine dynamiccleaning(F)
    
    integer jx, jy, jd
    real(kind=QL) max_local, ceiling_local
    complex(kind=QL) F(NX,NY,ND,3), T(3)
    
    do jx = 1, nx
        do jy = 1, ny
            do jd = 1, nd
                T = f(jx,jy,jd,:)   
                max_local = MAX(ABS(real(t(1))),ABS(AIMAG(t(1))), &
                                ABS(real(t(2))),ABS(AIMAG(t(2))), &
                                ABS(real(t(2))),ABS(AIMAG(t(3))))
                ceiling_local = max_local * 1.D-12
                do IC = 1,3
                    if ( ABS(real(T(IC))) < ceiling_local) T(IC) = CMPLX(0._QL, AIMAG(T(IC)))
                    if ( ABS(aimag(T(IC))) < ceiling_local) T(IC) = CMPLX(real(T(IC)),0._QL)
                end do
                f(Jx,jy,jd,:) = T
            end do
        end do
    end do    
                 
    end subroutine dynamiccleaning
    
    END SUBROUTINE MTSPHERE3D
    
    
    SUBROUTINE GETIMPEDANCE(NF,NX,NY,ND,NLYR,NTERMS,PXYD,DPTHL,K,Z,DEPTH,RADIUS,E,ES,HS,IMPEDANCE,ET,HT)
    
! Subroutine to compute the point impedances and total fields from the
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
    
        use precision
    
        INTEGER NX, NY, ND, JX, JY, JD, IC, IP, J, NF, JF, NLYR, RXLYR
        real DEPTH, RADIUS
        real(KIND=QL) DR, PXYD(NX,NY,ND,3), DPTHL(NLYR), R
        COMPLEX(KIND=QL) IMPEDANCE(NF,NX,NY,ND,3,2), K(NF,0:NLYR), Z(NF,0:NLYR), E(NF,0:NLYR,2), &
                         ES(NF,NX,NY,ND,3,2), HS(NF,NX,NY,ND,3,2), &
                         ET(NF,NX,NY,ND,3,2), HT(NF,NX,NY,ND,3,2), &
                         A(2), B(2), C(2), DEN, EP, HP, EE, HH, &
                         YHATS, ZHATS, Hmat(2,2), Hplus(2,2), &
                         PSIAT(2,NTERMS,-1:1),PSIFT(2,NTERMS,-1:1)
        
        ! Compute the primary field at each vertical location
        ! Search
        
        ET = ES
        HT = HS
        do JF = 1, NF
            do Jd = 1, ND                        
                DR = PXYD(1,1,JD,3)
                RXLYR = 0
                do Jl = NLYR,1,-1
                    if ( DR > DPTHL(Jl) ) then
                        RXLYR = JL
                        exit
                    end if
                end do
                do JX = 1, nx
                    do Jy = 1, ny
                        R = SQRT ( PXYD(JX,Jy,JD,1) ** 2 + PXYD(JX,JY,JD,2) ** 2 + ( PXYD(JX,JY,JD,3) - DEPTH) ** 2 )
                        if ( R > RADIUS ) THEN
                            if ( rXLYR == NLYR ) then
                                EE = E(JF,NLYR,1) * EXP ( - CI * K(JF,RXLYR) * ( DR - DPTHL(RXLYR) ) )
                                HH = E(JF,NLYR,1) / Z(JF,RXLYR)
                            else if ( RXLYR == 0 ) THEN
                                EE = E(JF,0,1) * EXP ( - CI * K(JF,RXLYR) * DR ) &
                                   + E(JF,0,2) * EXP (   CI * K(JF,RXLYR) * DR )
                                HH = ( E(JF,0,1) * EXP ( - CI * K(JF,RXLYR) * DR ) &
                                     - E(JF,0,2) * EXP (   CI * K(JF,RXLYR) * DR ) ) / Z(JF,RXLYR)    
                            else
                                EE = E(JF,RXLYR,1) * EXP ( - CI * K(JF,RXLYR) * ( DR - DPTHL(RXLYR) ) ) &
                                   + E(JF,RXLYR,2) * EXP (   CI * K(JF,RXLYR) * ( DR - DPTHL(RXLYR) ) )
                                HH = ( E(JF,RXLYR,1) * EXP ( - CI * K(JF,RXLYR) * ( DR - DPTHL(RXLYR) ) ) &
                                     - E(JF,RXLYR,2) * EXP (   CI * K(JF,RXLYR) * ( DR - DPTHL(RXLYR) ) ) ) / Z(JF,RXLYR)    
                            end if
                            ET(JF,jx,jy,jd,1,1) = ET(JF,jx,jy,jd,1,1) + EE
                            ET(JF,jx,jy,jd,2,2) = ET(JF,jx,jy,jd,2,2) + EE
                            HT(JF,jx,jy,jd,1,2) = HT(JF,jx,jy,jd,1,2) - HH
                            HT(JF,jx,jy,jd,2,1) = HT(JF,jx,jy,jd,2,1) + HH
                        end IF
                        DO Ip = 1, 2 ! X and Y direction
                            A(Ip) = HT(JF,jx,jy,jd,1,Ip)
                            B(Ip) = HT(JF,jx,jy,jd,2,Ip)
                        END DO
                        
                        ! Construction de la matrice H 2x2
                        Hmat(1,1) = A(1)
                        Hmat(1,2) = A(2)
                        Hmat(2,1) = B(1)
                        Hmat(2,2) = B(2)

                        ! Pseudo-inverse (tolérance relative typique)
                        CALL ZPINV2(Hmat, Hplus, 1.0D-12)
                        
                        DO IC = 1, 3 ! X, Y 
                            IF ( IC < 3 ) THEN
                                DO Ip = 1, 2
                                    C(Ip) = ET(JF,jx,jy,jd,IC,Ip)
                                END DO
                            ELSE
                                DO Ip = 1, 2
                                    C(Ip) = HT(JF,jx,jy,jd,3,Ip)
                                END DO
                            END IF
                            IMPEDANCE(JF,jx,jy,jd,IC,1) = C(1) * Hplus(1,1) + C(2) * Hplus(2,1) !  Impedance x
                            IMPEDANCE(JF,jx,jy,jd,IC,2) = C(1) * Hplus(1,2) + C(2) * Hplus(2,2) !  Impedance y
                        END DO
                    END DO
                END DO
            end DO
        end DO
        
END SUBROUTINE GETIMPEDANCE
    
SUBROUTINE ZPINV2(A, Aplus, tol)

  ! Moore-Penrose for 2D complex matrix

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  COMPLEX(DP), INTENT(IN)  :: A(2,2)
  COMPLEX(DP), INTENT(OUT) :: Aplus(2,2)
  REAL(DP),    INTENT(IN)  :: tol

  COMPLEX(DP) :: U(2,2), V(2,2), VT(2,2), W(2,2)
  REAL(DP)    :: S(2), smax
  COMPLEX(DP) :: work(10)
  REAL(DP)    :: rwork(10)
  INTEGER     :: info, i

  ! Working copy
  
  W = A

  CALL ZGESVD('A','A', 2, 2, W, 2, S, U, 2, VT, 2, &
              work, 10, rwork, info)

  IF (info /= 0) THEN
     Aplus = (0.0_DP, 0.0_DP)
     RETURN
  END IF

  smax = MAX(S(1), S(2))
  DO i = 1, 2
     IF (S(i) > tol * smax .AND. S(i) > 0.0_DP) THEN
        S(i) = 1.0_DP / S(i)
     ELSE
        S(i) = 0.0_DP
     END IF
  END DO

  ! Sigma+
  W = (0.0_DP, 0.0_DP)
  W(1,1) = S(1)
  W(2,2) = S(2)

  ! A+ = V * Sigma+ * U^H
  ! VT contains V^H → V = CONJG(TRANSPOSE(VT))
  
  Aplus = MATMUL(W, CONJG(TRANSPOSE(U)))
  Aplus = MATMUL(CONJG(TRANSPOSE(VT)), Aplus)
  
END SUBROUTINE ZPINV2
    
SUBROUTINE MTLayeredEarthCorrection(NW,NTERMS,NLYR,THK,DPTHL,DEPTH,RADIUS, &
    SXLYR,YHATL,ZHATL,YHATS,ZHATS,RTM,RTE,AI)

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
! SXLYR: Sphere layer
! YHATL: Layer admittivity
! ZHATL: Layer impedivity
! YHATS: Sphere admittivity
! ZHATS: Sphere impedivity

!          Output
!          -----
! C: Correction tensor

   IMPLICIT NONE
   INTEGER NTERMS, NLYR, N, M, ICR, IC, IH, IHR, NLAT, NLON, &
       SXLYR, NS, IR, INFO, NP, NW, JLAt, jlon, lwork
   REAL(KIND=QL) THK(NLYR-1), DEPTH, RADIUS, centre(3), TOL, DPTHL(NLYR)
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
   
   CALL MTSphere_MultipleSources(NW,NTERMS,NLAT,NLON,depth,depth,NLYR,THK,&
       DPTHL,RADIUS,YHATL,ZHATL,RTM,RTE,E,H) 
   
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
                Call MTSphericalHarmonicAnalysis(nterms,nlat,nlon,radius,&
                    YHATL(SXLYR),zhatL(SXLYR),er,hr,.FALSE.,psia,psif)
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
    
   
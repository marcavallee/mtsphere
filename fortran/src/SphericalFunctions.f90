
module sphericalfunctions

    use precision 

    implicit none    

    contains 
    
COMPLEX(KIND=QL) FUNCTION YNM(N,M,theta,phi)

    IMPLICIT NONE
    INTEGER N,M
    REAL(KIND=QL) theta,phi,factnm,plmr !CX(1,0:N),costhe(1)
        
    FACTNM = (-1) ** M * SQRT(float(2*n+1)*FACTORIAL(N-ABS(M))/(4._QL*PI*FACTORIAL(N+ABS(M))))
!    costhe(1) = cos(theta)
!    call pm_polynomial_value (1, N, abs(m), costhe, cx)
!    YNM = factnm * cx(1,N) * EXP ( CI * M * PHI )
    plmr = plm(n,ABS(M),COS(THETA))
    YNM = factnm * plm(n,abs(M),COS(theta)) * EXP ( ci * M * phi )

END FUNCTION YNM    
    
function double_factorial(n)
    integer, intent(in) :: n
    integer :: i
    real(kind=QL) double_factorial
    double_factorial = 1.0_ql
    do i = 1, n, 2
       double_factorial = double_factorial * i
    end do
 end function double_factorial

real(kind=ql) function plm(n, m, x)
    integer, intent(in) :: n, m
    real(kind=ql), intent(in) :: x
    real(kind=ql) :: pm, pm1, pm2
    integer :: i
    if ( ABS(M) > N ) then
        plm = 0._QL
        return
    end if
    if ( M == 0 ) then
        pm = 1._QL
    else
        pm = (-1)**m * double_factorial(2 * m - 1) * (1.0_ql - x**2)**(m / 2.0_ql)
    end if
    if (n == m) then
       plm = pm
       return
    end if
    pm1 = x * (2.0_ql * m + 1) * pm
    if (n == m + 1) then
       plm = pm1
       return
    end if
    do i = m + 2, n
       pm2 = (x * (2.0_ql * i - 1) * pm1 - (i + m - 1) * pm) / (i - m)
       pm = pm1
       pm1 = pm2
    end do
    plm = pm1
end function plm

real(kind=ql) function factorial(n)
    integer, intent(in) :: n
    integer :: i
    real(kind=ql) :: fact
    if (n < 0) then
       write(*,*) 'Erreur : factorielle d un nombre négatif'
       stop
    end if
    fact = 1.0_ql
    do i = 1, n
       fact = fact * i
    end do
    factorial = fact
  end function factorial

    end module sphericalfunctions
  
SUBROUTINE SphericalFactorInitialisation(A,B,C,NTERMS)

! Function to compute factorial

!          Input
!          -----
! NTERMS: number of degrees

!          Output
!          -----
! A, B, C: spherical transformation coefficients

        IMPLICIT NONE
        INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
        INTEGER N,M,NTERMS
        REAL(KIND=QL) A(0:NTERMS,-NTERMS:NTERMS),B(0:NTERMS,-NTERMS:NTERMS),C(0:NTERMS,-NTERMS:NTERMS)
    
        A = 0._QL 
        B = 0._QL
        C = 0._QL
        DO N = 0, NTERMS
            DO M = -N,N
                A(N,M) = SQRT ( FLOAT ( N + 1 + M ) *  float ( N + 1 - M ) / ( float ( 2 * N + 1 ) * float ( 2 * N + 3 ) ) )
                IF ( M >= 0 ) THEN
                    b(N,M) = sqrt ( FLOAT ( n - m - 1) * float ( n - m ) / ( float ( 2 * n - 1 ) * float ( 2 * n + 1 ) ) )
                    c(N,M) = sqrt ( FLOAT ( n - m ) * float ( n + m + 1 ) )            
                ELSE
                    b(N,M) = - sqrt ( FLOAT ( n - m - 1) * float ( n - m ) / ( float ( 2 * n - 1 ) * float ( 2 * n + 1 ) ) )
                    c(N,M) = - sqrt ( FLOAT ( n - m ) * float ( n + m + 1 ) )
                END IF
            END DO
        END DO

    END SUBROUTINE SphericalFactorInitialisation 
    
 SUBROUTINE CylindricalDerivatives(NTERMS,IDER,LMBDA,F,DF)

 IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
     INTEGER NTERMS, IDER, M, MP
     REAL(KIND=QL) LMBDA
     COMPLEX(KIND=QL) F(-NTERMS:NTERMS), DF(-NTERMS-1:NTERMS+1),&
                      DXPIY(-NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                      DXMIY(-NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                      DX(   -NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                      DY(   -NTERMS-1:NTERMS+1,-NTERMS:NTERMS)

     CALL CylindricalDerivativeMatrixInitialisation(NTERMS,LMBDA,DXPIY,DXMIY,DX,DY)
     DF = (0._QL,0._QL)
     IF ( IDER == 1 ) THEN
         DO M = -NTERMS-1, NTERMS+1
            DO MP = -NTERMS, NTERMS
               DF(M) = DF(M) + DX(M,MP) * F(MP)
            END DO
        END DO
     ELSE IF ( IDER == 2 ) THEN
        DO M = -NTERMS-1, NTERMS+1
            DO MP = -NTERMS, NTERMS
               DF(M) = DF(M) + DY(M,MP) * F(MP)
            END DO
        END DO
    END IF

    END SUBROUTINE CylindricalDerivatives
   
SUBROUTINE REFLECTIONCOEFFICIENTS(NTERMS,RADIUS,YHAT,ZHAT,YHATS,ZHATS,RTM,RTE)

! Computation of the reflection coefficient of a conductive sphere

!          Input
!          -----
! NTERMS: number of degrees
! RADIUS: radius of the sphere
! YHAT: layer admittivity
! ZHAT: layer impedivity
! YHATS: sphere admittivity
! ZHATS: sphere impedivity

!          Output
!          -----
! RTM: radial transverse magnetic reflection coefficient
! RTE: radial transverse electric reflection coefficient

    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    COMPLEX(KIND=QL), PARAMETER :: CI=(0._ql,1._ql)
    INTEGER NTERMS,NM,N
    REAL(KIND=QL) RADIUS
    COMPLEX(KIND=QL) YHAT,ZHAT,YHATS,ZHATS,K,KS,RTM(NTERMS),RTE(NTERMS),DRJ,DRH,DRJS
    DOUBLE COMPLEX KA,KSA,JKA(0:NTERMS), DJKA(0:NTERMS), JKSA(0:NTERMS), DJKSA(0:NTERMS), &
                   YKA(0:NTERMS), DYKA(0:NTERMS), YKSA(0:NTERMS), DYKSA(0:NTERMS), & 
                   hkA(0:NTERMS), DHKa(0:NTERMS), hkSA(0:NTERMS), DHKSA(0:NTERMS),A,B,C,D
    logical inductive_limit
    
    K  = SQRT ( - YHAT  * ZHAT  )
    KS = SQRT ( - YHATS * ZHATS )
    KA =  K  * RADIUS
    KSA = KS * RADIUS

    CALL CSPHJY(NTERMS,KA,NM,JKA,DJKA,YKA,DYKA,HKA,DHKA)
    CALL CSPHJY(NTERMS,KSA,NM,JKSA,DJKSA,YKSA,DYKSA,HKSA,DHKSA)
        
    DO N = 1,NTERMS
                
        IF(ISNAN(REAL(JKSA(N))).OR.ISNAN(IMAG(JKSA(N))).OR. &
           ISNAN(REAL(YKSA(N))).OR.ISNAN(IMAG(YKSA(N))).OR.JKSA(N).EQ.0) THEN
            inductive_limit = .true.
        ELSE
            inductive_limit = .false.
        END IF
        if(inductive_limit) then
            RTM(n) = - ( jka(n) + ka * djka(n) ) / ( hka(n) + ka * dhka(n) )
            RTE(n) = - jka(n) / hka(n)
        else
            A = yhat  * jka(n)  * ( jksa(n) + ksa * djksa(n) )
            B = yhats * jksa(n) * ( jka(n)  + ka  * djka(n) )
            C = yhat  * hka(n)  * ( jksa(n) + ksa * djksa(n) )
            D = yhats * jksa(n) * ( hka(n)  + ka  * dhka(n) )
            RTM(n) = - ( A - B ) / ( C - D )
            
            A = zhat  * jka(n)  * ( jksa(n) + ksa * djksa(n) )
            B = zhats * jksa(n) * ( jka(n)  + ka  * djka(n) )
            C = zhat  * hka(n)  * ( jksa(n) + ksa * djksa(n) )
            D = zhats * jksa(n) * ( hka(n)  + ka  * dhka(n) )
            RTE(n) = - ( A - B ) / ( C - D )
       end if
        
    END DO
   
    END SUBROUTINE REFLECTIONCOEFFICIENTS
    
    SUBROUTINE SPHERICALCYLINDRICALCONVERSIONCoefficients(NTERMS,K,LMBDA,U,ZETA)

IMPLICIT NONE
INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)  
REAL(KIND=QL), PARAMETER :: PI = 3.141592653589793    
COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL) 
INTEGER N, M, MP, NTERMS, ID
REAL(KIND=QL) A(0:NTERMS,-NTERMS:NTERMS),B(0:NTERMS,-NTERMS:NTERMS),C(0:NTERMS,-NTERMS:NTERMS),LMBDA
COMPLEX(KIND=QL) ZETA(2,0:NTERMS,-NTERMS:NTERMS), K, U, TZ

! Development of the matrix converting spherical coefficients to cylindrical coefficients.
! Each element is a cylindrical function coefficient. The two first indices relate to spehrical functions
! and the last two to cylindrical functions.

!               Input
!
!       NTERMS: maximum spherical function degree
!       K: propagation constant
!       LAMBDA: cylindrical wavenumber
!       U: cylindrical propagation number

!               Output
!       CHI: conversion matrix

CALL SphericalFactorInitialisation(A,B,C,NTERMS)

ZETA = (0._QL,0._QL)
DO ID = 1, 2
    ZETA(ID,0,0) = CI * SQRT ( 4._QL * PI ) / ( K * U )
    DO N = 1, NTERMS
        do M = -n, n
            if ( M == n ) THEN
                ZETA(ID,N,M)   = &
                - LMBDA * ZETA(ID,N-1,M-1) / ( K * B(N,-N) )
            else if ( M == -n ) THEN
                ZETA(ID,N,M) = &
                - LMBDA * ZETA(ID,N-1,M+1) / ( K * B(N,-N) )
            else
                TZ = (0._QL,0._QL)
                IF ( N > ABS(M) + 1 ) THEN
                    TZ = ZETA(ID,N-2,M) * A(N-2,M) / A(N-1,M)
                END IF
                IF ( ID == 1 ) THEN
                    ZETA(ID,N,M) = TZ + U * ZETA(1,N-1,M) / ( K * A(N-1,M) )
                ELSE
                    ZETA(ID,N,M) = TZ - U * ZETA(2,N-1,M) / ( K * A(N-1,M) )
                END IF
            end if
        END DO
    end do
end DO
                
END SUBROUTINE SPHERICALCYLINDRICALCONVERSIONCoefficients    

subroutine CylindricalDerivativeMatrixInitialisation(NTERMS,LMBDA,DXPIY,DXMIY,DX,DY)

    implicit none
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NTERMS, M, MP
    COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL)
    REAL(KIND=QL) LMBDA
    COMPLEX(KIND=QL) DXPIY(-NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                     DXMIY(-NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                     DX   (-NTERMS-1:NTERMS+1,-NTERMS:NTERMS), &
                     DY   (-NTERMS-1:NTERMS+1,-NTERMS:NTERMS)
    
    DXPIY = (0._QL,0._QL)
    DXMIY = (0._QL,0._QL)
    DX    = (0._QL,0._QL)
    DY    = (0._QL,0._QL)
    DO M = -NTERMS, NTERMS
        IF ( M == 0 ) THEN
            DXPIY( 1, 0) = -1._QL
            DXMIY(-1, 0) = -1._QL
        ELSE IF ( M < 0 ) THEN    
            MP = M + 1
            DXPIY(MP,M) =   1._QL ! Zero si M == N, M = -N, MP = M + 1, NP = N - 1
            MP = M - 1
            DXMIY(MP,M) = - 1._QL ! Zero si M == N, M = -N, MP = M - 1, NP = N + 1
        ELSE IF ( M > 0 ) THEN
            MP = M + 1
            DXPIY(MP,M) = - 1._QL ! Zero si M == -N, M = N, MP = M + 1, NP = N + 1
            MP = M - 1
            DXMIY(MP,M) =   1._QL ! Zero si M == -N, M = N, MP = M - 1, NP = N - 1
        END IF
    END DO
    DXPIY = DXPIY * LMBDA 
    DXMIY = DXMIY * LMBDA
    DO M = -NTERMS-1, NTERMS+1
        DO MP = -NTERMS, NTERMS
            DX(M,MP) =        0.5_QL * ( DXPIY(M,MP) + DXMIY(M,MP) )
            DY(M,MP) = - CI * 0.5_QL * ( DXPIY(M,MP) - DXMIY(M,MP) )
        END DO
    END DO

    end subroutine CylindricalDerivativeMatrixInitialisation
    
    SUBROUTINE HSSPHERE_KER (LMBDA,SXLYR,NLYR,THKD,DPTHL,ZS,YHAT,ZHAT,AN,FN)
!-----------------------------------------------------------------------------

!  Kernel for the computation of the layer propagation potentials.
!
!          Input
!          -----
!   LMBDA - Hankel transform variable
!   SXLYR - source layer index
!   RXLYR - receiver layer index
!    NLYR - number of layers
!    THKD - layer thicknesses
!   DPTHL - depth to the top of each layer
!      ZP - depth of the source
!    YHAT - layer admittance
!    ZHAT - layer impedance

!          Output
!          ------
!      AN - TM layer propagation coefficients
!      FN - TE layer propagation coefficients


 IMPLICIT NONE 
 INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
 REAL(KIND=QL), PARAMETER :: EPS0=8.854156D-12, MU0=12.56637D-7, EXP_TOL=80.D0
 COMPLEX(KIND=QL), PARAMETER :: ONE=(1._QL,0._QL),ZERO=(0._QL,0._QL)
 INTEGER NLYR,J,SXLYR,INC
 REAL, DIMENSION(NLYR) :: CALF,CTAU,CFREQ,REPS
 REAL(KIND=QL) THKD(NLYR-1),LMBDA, LMBSQ, DPTHL(NLYR), ZS
 COMPLEX(KIND=QL)  EP,T(NLYR),EJ,EH,EMH,EPH,E2H,ES1,ES2,FN(2,0:NLYR),AN(2,0:NLYR), &
     FDEN,ADEN,RTE(2,1:NLYR-1),RTM(2,1:NLYR-1)
 COMPLEX(KIND=QL), DIMENSION (NLYR) :: FT,AT,ADMHAT,IMPHAT
 COMPLEX(KIND=QL), DIMENSION (0:NLYR) :: S,YHAT,ZHAT,ADM,IMP
 
 AN = (0._QL,0._QL)
 FN = (0._QL,0._QL)
 AT = (0._QL,0._QL)
 FT = (0._QL,0._QL)

 IF ( NLYR == 0 ) RETURN
 
 LMBSQ = LMBDA ** 2
 DO J = 0, NLYR
    S(J) = SQRT ( LMBSQ + YHAT(J) * ZHAT(J) )
    ADM(J) = S(J) / ZHAT(J)
    IMP(J) = S(J) / YHAT(J)
 END DO
 
 ! Admittance and impedance computation for layers higher than the source
 IF ( SXLYR < NLYR ) THEN 
    ADMHAT(NLYR) = ADM(NLYR)
    IMPHAT(NLYR) = IMP(NLYR)
    DO J = NLYR-1, SXLYR+1, -1
        EP = TANH ( S(J) * THKD(J) )
        ADMHAT(J) = ADM(J) * ( ADMHAT(J+1) + ADM(J) * EP ) / &
                             ( ADM(J) + ADMHAT(J+1) * EP )
        IMPHAT(J) = IMP(J) * ( IMPHAT(J+1) + IMP(J) * EP ) / &
                             ( IMP(J) + IMPHAT(J+1) * EP )
    END DO
 END IF 
 
 ! Admittance and impedance computation for layers lower than the source
 
IF ( SXLYR > 0 ) THEN
    ADMHAT(1) = - ADM(0)
    IMPHAT(1) = - IMP(0)
    DO J = 1, SXLYR-1
        EP = TANH ( S(J) * THKD(J) )
        ADMHAT(J+1) = ADM(J) * ( ADMHAT(J) - ADM(J)    * EP ) / &
                               ( ADM(J)    - ADMHAT(J) * EP )
        IMPHAT(J+1) = IMP(J) * ( IMPHAT(J) - IMP(J)    * EP ) / &
                               ( IMP(J)    - IMPHAT(J) * EP )
    END DO
END IF    

  IF ( SXLYR == 0 ) THEN
    EP = S(0) * ( DPTHL(1) - ZS )
    IF ( REAL (EP) < EXP_TOL) THEN
        EJ = EXP (-EP)
    ELSE
        EJ = ZERO
    END IF
    FN(2,0) = EJ * ( ADM(0) - ADMHAT(1) ) / ( ADM(0) + ADMHAT(1) )
    AN(2,0) = EJ * ( IMP(0) - IMPHAT(1) ) / ( IMP(0) + IMPHAT(1) )
    FT(1) = FN(2,0) + EJ
    AT(1) = AN(2,0) + EJ
 ELSE IF ( SXLYR == NLYR ) THEN
    EP = S(NLYR) * ( ZS - DPTHL(NLYR) )
    IF ( REAL (EP) < EXP_TOL) THEN
        EJ = EXP (-EP)
    ELSE
        EJ = ZERO
    END IF
    IF ( ADM(NLYR) .NE. ADMHAT(NLYR) ) THEN
        FN(1,NLYR) = EJ * ( ADM(NLYR) + ADMHAT(NLYR) ) / ( ADM(NLYR) - ADMHAT(NLYR) )
    END IF
    IF ( IMP(NLYR) .NE. IMPHAT(NLYR) ) THEN
        AN(1,NLYR) = EJ * ( IMP(NLYR) + IMPHAT(NLYR) ) / ( IMP(NLYR) - IMPHAT(NLYR) )
    END IF
    FT(NLYR) = FN(1,NLYR) + EJ
    AT(NLYR) = AN(1,NLYR) + EJ
ELSE
    EP = 2.0_QL * S(SXLYR) * THKD(SXLYR) 
    IF ( REAL(EP) < EXP_TOL ) THEN
        EMH = EXP(-EP/2.0_QL)
        EPH = EXP(EP/2.0_QL)
        E2H = EXP(-EP)
        ES1 = EXP( - S(SXLYR) * ( ZS - DPTHL(SXLYR) ) )
        ES2 = EXP( - S(SXLYR) * ( DPTHL(SXLYR+1) - ZS ) ) 
    ELSE
        EPH = ZERO
        EMH = ZERO
        ES1 = ZERO
        ES2 = ZERO
    END IF
    FDEN = ( ADMHAT(SXLYR) - ADM(SXLYR) ) * ( ADMHAT(SXLYR+1) + ADM(SXLYR) ) - &
           ( ADMHAT(SXLYR) + ADM(SXLYR) ) * ( ADMHAT(SXLYR+1) - ADM(SXLYR) ) * E2H
    ADEN = ( IMPHAT(SXLYR) - IMP(SXLYR) ) * ( IMPHAT(SXLYR+1) + IMP(SXLYR) )  - &
           ( IMPHAT(SXLYR) + IMP(SXLYR) ) * ( IMPHAT(SXLYR+1) - IMP(SXLYR) ) * E2H
    FN(1,SXLYR) = - ( ADMHAT(SXLYR) + ADM(SXLYR) ) * & 
        ( ( ADMHAT(SXLYR+1) + ADM(SXLYR) ) * ES1 - ( ADMHAT(SXLYR+1) - ADM(SXLYR) ) * EMH * ES2 ) / ADEN
    FN(1,SXLYR) = - ( IMPHAT(SXLYR) + IMP(SXLYR) ) * & 
        ( ( IMPHAT(SXLYR+1) + IMP(SXLYR) ) * ES1 - ( IMPHAT(SXLYR+1) - IMP(SXLYR) ) * EMH * ES2 ) / FDEN
    FN(2,SXLYR) = ( ADMHAT(SXLYR+1) - ADM(SXLYR) ) * & 
        ( ( ADMHAT(SXLYR) + ADM(SXLYR) ) * E2H * ES1 - ( ADMHAT(SXLYR) - ADM(SXLYR) ) * EMH * ES2 ) / ADEN
    AN(2,SXLYR) = ( IMPHAT(SXLYR+1) - IMP(SXLYR) ) * & 
        ( ( IMPHAT(SXLYR) + IMP(SXLYR) ) * E2H * ES1 - ( IMPHAT(SXLYR) - IMP(SXLYR) ) * EMH * ES2 ) / FDEN
    FT(SXLYR) = FN(1,SXLYR) + FN(2,SXLYR) + ES1
    AT(SXLYR) = AN(1,SXLYR) + AN(2,SXLYR) + ES1
    FT(SXLYR+1) = FN(1,SXLYR) * EMH + FN(2,SXLYR) * EPH + ES2
    AT(SXLYR+1) = AN(1,SXLYR) * EMH + AN(2,SXLYR) * EPH + ES2
END IF

do J = SXLYR-1, 0, -1
    if ( J == 0 ) then
        AN(2,0) = AT(1)
        FN(2,0) = FT(1)
    else
        EP = 2._QL * S(J) * THKD(J) 
        IF ( REAL(EP) < EXP_TOL ) THEN
            EMH = EXP(-EP/2.0_QL)
            E2H = EXP(-EP)
        ELSE
            EMH = ZERO
            E2H = ZERO
        END if
        AT(J) = AT(J+1) * 2 * imp(J) * EMH / ( imp(J) * ( 1.0_QL + E2H ) - imphat(J) * ( 1.0_QL - E2H ) )
        FT(J) = FT(J+1) * 2 * ADM(J) * EMH / ( ADM(J) * ( 1.0_QL + E2H ) - ADMhat(J) * ( 1.0_QL - E2H ) )
        AN(1,J) = AT(J) * ( imp(J) + imphat(J) ) / 2.0_QL / imp(J)
        AN(2,J) = AT(J) * ( imp(J) - imphat(J) ) / 2.0_QL / imp(J)
        FN(1,J) = FT(J) * ( adm(J) + admhat(J) ) / 2.0_QL / adm(J)
        FN(2,J) = FT(J) * ( adm(J) - admhat(J) ) / 2.0_QL / adm(J)
    end if
end do
do J = SXLYR+1, NLYR
    if ( J == NLYR ) then
        AN(1,NLYR) = AT(NLYR)
        FN(1,NLYR) = FT(NLYR)
    else
        EP = S(J) * THKD(J) 
        IF ( REAL(EP) < EXP_TOL ) THEN
            EMH = EXP(-EP)
            EPH = EXP(EP)
        else
            EMH = ZERO
            EPH = ZERO
        end if
        AN(1,J) = AT(J) * ( imp(J) + imphat(J) ) / 2.0_QL / imp(J)
        AN(2,J) = AT(J) * ( imp(J) - imphat(J) ) / 2.0_QL / imp(J)
        FN(1,J) = FT(J) * ( adm(J) + admhat(J) ) / 2.0_QL / adm(J)
        FN(2,J) = FT(J) * ( adm(J) - admhat(J) ) / 2.0_QL / adm(J)
        AT(J+1) = AN(1,J) * EMH + AN(2,J) * EPH
        FT(J+1) = FN(1,J) * EMH + FN(2,J) * EPH
    end if
end do

END SUBROUTINE HSSPHERE_KER
    
SUBROUTINE VerticalFieldSphericalCoefficients(NTERMS,PSIA,PSIF,YHAT,ZHAT,EZNM,HZNM)

! Computation of the spherical field coefficients from the potential spherical coefficients

!          Input
!          -----
! NTERMS: number of degrees
! PHIA: transverse magnetic spherical potential
! PSIF: transverse electric spherical potential
! YHAT: admittivity
! ZHAT: impedivity

!          Output
!          -----
! EZNM: electric field spherical coefficients
! HZNM: magnetic field spherical coefficients

    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL)
    INTEGER NTERMS, I, N, M 
    COMPLEX(KIND=QL) EZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1),HZNM(0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                     PSIA(NTERMS,-NTERMS:NTERMS), PSIF(NTERMS,-NTERMS:NTERMS),K,yhat,zhat
    REAL(KIND=QL) A(0:NTERMS+1,-NTERMS-1:NTERMS+1), B(0:NTERMS+1,-NTERMS-1:NTERMS+1), &
                  C(0:NTERMS+1,-NTERMS-1:NTERMS+1)

    CALL SphericalFactorInitialisation(A,B,C,NTERMS+1)
    K = sqrt ( - yhat * zhat )
    CALL ComputeField(PSIF,PSIA,EZNM,YHAT,1)
    CALL ComputeField(PSIA,PSIF,HZNM,ZHAT,-1)  
        
    CONTAINS
    
    SUBROUTINE ComputeField(PHI,PSI,FNM,PHAT,SIGNE)
    
    INTEGER SIGNE
    COMPLEX(KIND=QL) PHAT, PHI(NTERMS,-NTERMS:NTERMS), PSI(NTERMS,-NTERMS:NTERMS), &
        FNM(0:NTERMS+1,-NTERMS-1:NTERMS+1), TZ
    
        FNM = (0._QL,0._QL)
        DO N = 0, NTERMS + 1
            DO M = - N, N
                TZ    = (0._QL,0._QL) 
                IF ( N > 0 .AND. N .LE. NTERMS ) THEN
                    TZ    = TZ    + CI * M                      * PHI(N,M)   * SIGNE
                ENDIF
                IF ( N < NTERMS  ) THEN
                    TZ    = TZ    + K * ( N + 2 ) * A(N,M)      * PSI(N+1,M)   / PHAT
                END IF
                IF ( N > 1 .AND. ABS(M) < N )THEN                      
                    TZ    = TZ    + K * ( N - 1 ) * A(N-1,M)    * PSI(N-1,M)   / PHAT
                END IF
                FNM(N,M) = TZ
            END DO
        END DO    
    
    END SUBROUTINE ComputeField 

end subroutine VerticalFieldSphericalCoefficients   
    
    
    
    
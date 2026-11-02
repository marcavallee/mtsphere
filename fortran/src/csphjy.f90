        SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY,CSH,CDH)
!
!       ==========================================================
!       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!                and their derivatives for a complex argument
!       Input :  z --- Complex argument
!                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!       Output:  CSJ(n) --- jn(z)
!                CDJ(n) --- jn'(z)
!                CSY(n) --- yn(z)
!                CDY(n) --- yn'(z)
!                CSH(n) --- jn(z) - i y(n)
!                CDY(n) --- jn'(z) - i y'(n)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ==========================================================
!
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0
        DIMENSION CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N),CSH(0:N),CDH(0:N)
        CSJ = 0.0D0
        CDJ = 0.0D0
        CSY = 0.0D0
        CDY = 0.0D0
        CSH = 0.0D0
        CDH = 0.0D0
        CI = (0.0D0,1.0D0)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           CSJ(0)=(1.0D0,0.0D0)
           CDJ(0)=(0.0D0,0.0D0)
           CSY(0)=(-1.0D+300, 0.0d0 ) 
           CDY(0)=(1.0D+300, 0.0d0 )
           CSH(0)=CSJ(0) - CI * CSY(0)
           CDH(0)=CDJ(0) - CI * CDY(0)
           DO K=1,N
              CSJ(K)=0.0D0
              CDJ(K)=0.0D0
              CSY(K)=-1.0D+300
              CDY(K)=1.0D+300
           END DO
           RETURN
        ENDIF
    ! --- Étape 1 : Calcul de j_n(z) par récurrence descendante ---
    csj(0) = SIN(z) / z
    csj(1) = (csj(0) - COS(z)) / z
    
    mtop = msta2(a0, n, 15) 
    IF (mtop < n) THEN
        nm = mtop
    ELSE
        mtop = n
    END IF

    cf1 = (0.0d0, 0.0d0)
    cf0 = (1.0d-100, 0.0d0)
    cf = (0.0d0, 0.0d0)
    
    DO k = mtop, 0, -1
        cf = (2.0d0 * k + 3.0d0) * cf0 / z - cf1
        IF (k <= nm) csj(k) = cf0
        cf1 = cf0
        cf0 = cf
    END DO

    a0 = ABS(csj(0))
    IF (a0 >= ABS(csj(1))) THEN
        cs = SIN(z) / (z * csj(0))
    ELSE
        cs = (SIN(z) / z - COS(z)) / (z * csj(1))
    END IF
    
    DO k = 0, nm
        csj(k) = csj(k) * cs
    END DO

    ! --- Étape 2 : Calcul de y_n(z) par récurrence ascendante ---
    cy0 = -COS(z) / z
    cy1 = (cy0 - SIN(z)) / z
    csy(0) = cy0
    csy(1) = cy1
    
    DO k = 2, nm
        csy(k) = (2.0d0 * k - 1.0d0) * csy(k-1) / z - csy(k-2)
    END DO

    ! --- Étape 3 : Calcul des dérivées premières j'_n(z) et y'_n(z) ---
    cdj(0) = (COS(z) - SIN(z) / z) / z
    DO k = 1, nm
        cdj(k) = csj(k-1) - (k + 1.0d0) * csj(k) / z
    END DO
    
    cdy(0) = (SIN(z) + COS(z) / z) / z
    DO k = 1, nm
        cdy(k) = csy(k-1) - (k + 1.0d0) * csy(k) / z
    END DO

    ! --- Étape 4 : Calcul de la fonction de Hankel sphérique H2 et sa dérivée ---
    DO k = 0, nm
        csh(k) = csj(k) - ci * csy(k)
        cdh(k) = cdj(k) - ci * cdy(k)
    END DO
        RETURN
        END



!        PROGRAM MCSPHJY
!
!       ================================================================
!       Purpose: This program computes the spherical Bessel functions 
!                jn(z), yn(z), and their derivatives for a complex
!                argument using subroutine CSPHJY
!       Input :  z --- Complex argument
!                n --- Order of jn(z) & yn(z) ( 0 ó n ó 250 )
!       Output:  CSJ(n) --- jn(z)
!                CDJ(n) --- jn'(z)
!                CSY(n) --- yn(z)
!                !DY(n) --- yn'(z)
!       Example: z = 4.0+i 2.0
!
!     n     Re[jn(z)]       Im[jn(z)]       Re[jn'(z)]      Im[jn'(z)]
!   --------------------------------------------------------------------
!     0  -.80651523D+00  -.18941093D+00  -.37101203D-01   .75210758D+00
!     1   .37101203D-01  -.75210758D+00  -.67093420D+00   .11885235D+00
!     2   .60314368D+00  -.27298399D+00  -.24288981D+00  -.40737409D+00
!     3   .42955048D+00   .17755176D+00   .18848259D+00  -.24320520D+00
!     4   .12251323D+00   .22087111D+00   .19660170D+00   .17937264D-01
!     5  -.10242676D-01   .10975433D+00   .68951842D-01   .83020305D-01
!
!     n     Re[yn(z)]       Im[yn(z)]       Re[yn'(z)]      Im[yn'(z)]
!   --------------------------------------------------------------------
!     0   .21734534D+00  -.79487692D+00  -.77049661D+00  -.87010064D-02
!     1   .77049661D+00   .87010064D-02  -.92593503D-01  -.64425800D+00
!     2   .24756293D+00   .56894854D+00   .45127429D+00  -.25839924D+00
!     3  -.23845941D+00   .43646607D+00   .26374403D+00   .12439192D+00
!     4  -.27587985D+00   .20902555D+00  -.67092335D-01   .89500599D-01
!     5  -.70001327D-01   .18807178D+00  -.30472133D+00  -.58661384D-01
!       ================================================================
!
!        IMPLICIT COMPLEX*16 (C,Z)
!        REAL(KIND=QL) X,Y
!        DIMENSION CSJ(0:250),CDJ(0:250),CSY(0:250),CDY(0:250)
!        WRITE(*,*)'Please enter n,x,y (z=x+iy) '
!        READ(*,*)N,X,Y
!        WRITE(*,30)N,X,Y
!        Z=CMPLX(X,Y)
!        IF (N.LE.8) THEN
!           NS=1
!        ELSE
!           WRITE(*,*)'Please enter order step Ns '
!           READ(*,*)NS
!        ENDIF
!        CALL CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY,CSH,CDH)
!        WRITE(*,*)
!        WRITE(*,*)'  n      Re[jn(z)]        Im[jn(z)]',
!     &  '        Re[jn''(z)]       Im[jn''(z)]'
!        WRITE(*,*)'--------------------------------------------',
!     &  '----------------------------'
!        DO 10 K=0,NM,NS
!10         WRITE(*,20)K,CSJ(K),CDJ(K)
!        WRITE(*,*)
!        WRITE(*,*)'  n      Re[yn(z)]        Im[yn(z)]',
!     &  '        Re[yn''(z)]       Im[yn''(z)]'
!        WRITE(*,*)'--------------------------------------------',
!     &  '----------------------------'
!        DO 15 K=0,NM,NS
!15         WRITE(*,20)K,CSY(K),CDY(K)
!20      FORMAT(1X,I3,4D17.8)
!30      FORMAT(3X,6HNmaz =,I3,',     ','z = ',F8.1,'+ i',F8.1)
!        END
!
!
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
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ==========================================================
!
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0
        DIMENSION CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N),CSH(0:N),CDH(0:N)
        CSJ = 0
        CDJ = 0
        CSY = 0
        CDY = 0
        CSH = 0
        CDH = 0
        CI = (0.0D0,1.0D0)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO K=0,N
              CSJ(K)=0.0D0
              CDJ(K)=0.0D0
              CSY(K)=-1.0D+300
              CDY(K)=1.0D+300
           end DO
           CSJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(.333333333333333D0,0.0D0)
           RETURN
        ENDIF
        IF(AIMAG(Z).GT.-LOG(HUGE(A0)))THEN
            CSJ(0)=CDSIN(Z)/Z
            CSJ(1)=(CSJ(0)-CDCOS(Z))/Z
            IF (N.GE.2) THEN
               CSA=CSJ(0)
               CSB=CSJ(1)
               M=MSTA1(A0,200)
               IF (M.LT.N) THEN
                  NM=M
               ELSE
                  M=MSTA2(A0,N,15)
               ENDIF
               IF ( M > 0 ) THEN
                   CF0=0.0D0
                   CF1=1.0D0-100
                   DO K=M,0,-1
                      CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
                      IF (K.LE.NM) CSJ(K)=CF
                      CF0=CF1
                      CF1=CF
                   END DO
                   IF (CDABS(CSA).GT.CDABS(CSB)) CS=CSA/CF
                   IF (CDABS(CSA).LE.CDABS(CSB)) CS=CSB/CF0
                   DO K=0,NM
                      CSJ(K)=CS*CSJ(K)
                   END DO
               ELSE
                   DO K = 2,NM
                       CSJ(K) = (2*K-1)*CSJ(K-1)/Z-CSJ(K-2)
                   END DO
                END IF
            ENDIF
            CDJ(0)=(CDCOS(Z)-CDSIN(Z)/Z)/Z
        END IF
        DO K=1,NM
            CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
        end DO
        CSY(0)=-CDCOS(Z)/Z
     !   CSY(1)=(CSY(0)-CDSIN(Z))/Z
        CSY(1)=-CDCOS(Z)/(Z*Z)-CDSIN(Z)/Z
        DO K=2,NM
           IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
              CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
           ELSE
              CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
           ENDIF
          !    CSY(K) = (2*K-1)*CSY(K-1)/Z-CSY(K-2)
        END DO
        CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
        CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
        DO K=2,NM
            CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
        END DO
        CSH(0)=CI*CDEXP(-CI*Z)/Z
        CSH(1)=CDEXP(-CI*Z)*(CI/Z**2-1.0D0/Z)
        DO K=2,NM
         !  IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
         !     CSH(K)=(CSJ(K)*CSH(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
         !  ELSE
         !     CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
         !  ENDIF
            CSH(K) = (2*K-1)*CSH(K-1)/Z-CSH(K-2)
        END DO
        CDH(0)=CDEXP(-CI*Z)*(1.0D0/Z-CI/(Z**2))
        CDH(1)=CDEXP(-CI*Z)*(CI/Z+2.0D0/(Z**2)-2.0D0*CI/(Z**3))
        DO K=2,NM
            CDH(K)=CSH(K-1)-(K+1.0D0)*CSH(K)/Z
        END DO
        RETURN
        END


!        INTEGER FUNCTION MSTA1(X,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!       ===================================================
!
!        IMPLICIT REAL(KIND=QL) (A-H,O-Z)
!        A0=DABS(X)
!        N0=INT(1.1*A0)+1
!        F0=ENVJ(N0,A0)-MP
!        N1=N0+5
!        F1=ENVJ(N1,A0)-MP
!        DO 10 IT=1,20             
!           NTERMS=N1-(N1-N0)/(1.0D0-F0/F1)                  
!           F=ENVJ(NTERMS,A0)-MP
!           IF(ABS(NTERMS-N1).LT.1) GO TO 20
!           N0=N1
!           F0=F1
!           N1=NTERMS
! 10        F1=F
! 20     MSTA1=NTERMS
!        RETURN
!        END


!        INTEGER FUNCTION MSTA2(X,N,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
!        IMPLICIT REAL(KIND=QL) (A-H,O-Z)
!        A0=DABS(X)
!        HMP=0.5D0*MP
!        EJN=ENVJ(N,A0)
!        IF (EJN.LE.HMP) THEN
!           OBJ=MP
!           N0=INT(1.1*A0)
!        ELSE
!           OBJ=HMP+EJN
!           N0=N
!        ENDIF
!        F0=ENVJ(N0,A0)-OBJ
!        N1=N0+5
!        F1=ENVJ(N1,A0)-OBJ
!        DO 10 IT=1,20
!           NTERMS=N1-(N1-N0)/(1.0D0-F0/F1)
!           F=ENVJ(NTERMS,A0)-OBJ
!           IF (ABS(NTERMS-N1).LT.1) GO TO 20
!           N0=N1
!           F0=F1
!           N1=NTERMS
!10         F1=F
!20      MSTA2=NTERMS+10
!        RETURN
!        END

 !       REAL*8 FUNCTION ENVJ(N,X)
 !       REAL(KIND=QL) X
 !       ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
 !       RETURN
 !       END

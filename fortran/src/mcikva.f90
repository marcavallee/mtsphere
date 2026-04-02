
        SUBROUTINE CIKVA(V,Z,VM,CBI,CDI,CBK,CDK)

!       ============================================================
!       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
!                and their derivatives for an arbitrary order and
!                complex argument
!       Input :  z --- Complex argument
!                v --- Real order of Iv(z) and Kv(z)
!                      ( v = n+v0, n = 0,1,2,���, 0 � v0 < 1 )
!       Output:  CBI(n) --- In+v0(z)
!                CDI(n) --- In+v0'(z)
!                CBK(n) --- Kn+v0(z)
!                CDK(n) --- Kn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting 
!                point for backward recurrence
!       ============================================================

        IMPLICIT DOUBLE PRECISION (A,G,P,R,V,W)
        IMPLICIT DOUBLE COMPLEX (C,Z)
        DIMENSION CBI(0:*),CDI(0:*),CBK(0:*),CDK(0:*)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (N.EQ.0) N=1
        IF (A0.LT.1.0D-100) THEN
           DO K=0,N
              CBI(K)=0.0D0
              CDI(K)=0.0D0
              CBK(K)=-1.0D+300
              CDK(K)=1.0D+300
           END DO
           IF (V0.EQ.0.0) THEN
              CBI(0)=(1.0D0,0.0D0)
              CDI(1)=(0.5D0,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        K0=14
        IF (A0.GE.35.0) K0=10
        IF (A0.GE.50.0) K0=8
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LT.18.0) THEN
           IF (V0.EQ.0.0) THEN
              CA1=(1.0D0,0.0D0)
           ELSE
              V0P=1.0D0+V0
              CALL GAMMA(V0P,GAP)
              CA1=(0.5D0*Z1)**V0/GAP
           ENDIF
           CI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO K=1,50
              CR=0.25D0*CR*Z2/(K*(K+V0))
              CI0=CI0+CR
              IF (CDABS(CR).LT.CDABS(CI0)*1.0D-15)EXIT
           END DO
           CBI0=CI0*CA1
        ELSE
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO K=1,K0
              CR=-0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
              CS=CS+CR
           END DO
           CBI0=CA*CS
        ENDIF
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1+CF2
           IF (K.LE.N) CBI(K)=CF
           CF2=CF1
           CF1=CF
        END DO
        CS=CBI0/CF
        DO K=0,N
           CBI(K)=CS*CBI(K)
        END DO
        IF (A0.LE.9.0) THEN
           IF (V0.EQ.0.0) THEN
              CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
              CS=(0.0D0,0.0D0)
              W0=0.0D0
              CR=(1.0D0,0.0D0)
              DO K=1,50
                 W0=W0+1.0D0/K
                 CR=0.25D0*CR/(K*K)*Z2
                 CP=CR*(W0+CT)
                 CS=CS+CP
                 IF (K.GE.10.AND.CDABS(CP/CS).LT.1.0D-15)EXIT
              END DO
              CBK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA(V0N,GAN)
              CA2=1.0D0/(GAN*(0.5D0*Z1)**V0)
              CA1=(0.5D0*Z1)**V0/GAP
              CSU=CA2-CA1
              CR1=(1.0D0,0.0D0)
              CR2=(1.0D0,0.0D0)
              DO K=1,50
                 CR1=0.25D0*CR1*Z2/(K*(K-V0))
                 CR2=0.25D0*CR2*Z2/(K*(K+V0))
                 CSU=CSU+CA2*CR1-CA1*CR2
                 WS=CDABS(CSU)
                 IF (K.GE.10.AND.DABS(WS-WS0)/WS.LT.1.0D-15)EXIT
                 WS0=WS
              END DO
              CBK0=0.5D0*PI*CSU/DSIN(PIV)
           ENDIF
        ELSE
           CB=CDEXP(-Z1)*CDSQRT(0.5D0*PI/Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO K=1,K0
              CR=0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
              CS=CS+CR
           END DO
           CBK0=CB*CS
        ENDIF
        CBK1=(1.0D0/Z1-CBI(1)*CBK0)/CBI(0)
        CBK(0)=CBK0
        CBK(1)=CBK1
        CG0=CBK0
        CG1=CBK1
        DO K=2,N
           CGK=2.0D0*(V0+K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CGK
           CG0=CG1
           CG1=CGK
        END DO
        IF (REAL(Z).LT.0.0) THEN
           DO K=0,N
              CVK=CDEXP((K+V0)*PI*CI)
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBK(K)=CVK*CBK(K)+PI*CI*CBI(K)
                 CBI(K)=CBI(K)/CVK
              ELSE IF (DIMAG(Z).GT.0.0) THEN
                 CBK(K)=CBK(K)/CVK-PI*CI*CBI(K)
                 CBI(K)=CVK*CBI(K)
              ENDIF
           END DO
        ENDIF
        CDI(0)=V0/Z*CBI(0)+CBI(1)
        CDK(0)=V0/Z*CBK(0)-CBK(1)
        DO K=1,N
           CDI(K)=-(K+V0)/Z*CBI(K)+CBI(K-1)
           CDK(K)=-(K+V0)/Z*CBK(K)-CBK(K-1)
        END DO
        VM=N+V0
        RETURN
        END

        SUBROUTINE GAMMA(X,GA)

!       ==================================================
!       Purpose: Compute gamma function �(x)
!       Input :  x  --- Argument of �(x)
!                       ( x is not equal to 0,-1,-2,���)
!       Output:  GA --- �(x)
!       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        DATA G/1.0D0,0.5772156649015329D0, &
           -0.6558780715202538D0, -0.420026350340952D-1, &
           0.1665386113822915D0,-.421977345555443D-1, &
           -.96219715278770D-2, .72189432466630D-2, &
           -.11651675918591D-2, -.2152416741149D-3, &
           .1280502823882D-3, -.201348547807D-4, &
           -.12504934821D-5, .11330272320D-5, &
           -.2056338417D-6, .61160950D-8, &
           .50020075D-8, -.11812746D-8, &
           .1043427D-9, .77823D-11, &
           -.36968D-11, .51D-12, &
           -.206D-13, -.54D-14, .14D-14, .1D-15/
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=INT(X)-1
              DO K=2,M1
                 GA=GA*K
              END DO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO K=1,M
                 R=R*(Z-K)
              END DO
              Z=Z-M
           ELSE
              Z=X
           ENDIF
          GR=G(26)
           DO K=25,1,-1
              GR=GR*Z+G(K)
           END DO
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)

!       ===================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO IT=1,20
           NTERMS=N1-INT((N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NTERMS,A0)-MP
           IF(ABS(NTERMS-N1).LT.1)THEN
              MSTA1=NTERMS
              RETURN
           END IF
           N0=N1
           F0=F1
           N1=NTERMS
           F1=F
        END DO
        MSTA1=NTERMS
        RETURN
        END

        INTEGER FUNCTION MSTA2(X,N,MP)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO IT=1,20
           NTERMS=N1-INT((N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NTERMS,A0)-OBJ
           IF (ABS(NTERMS-N1).LT.1)EXIT
           N0=N1
           F0=F1
           N1=NTERMS
           F1=F
        END DO
        MSTA2=NTERMS+10
        RETURN
        END

        DOUBLE PRECISION FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END

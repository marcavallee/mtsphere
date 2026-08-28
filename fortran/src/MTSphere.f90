!  PROGRAM MTSPHERE Version 2.2.0 September 2026
!    
!    
!        Developed by: Marc A. Vallée
!                 For: Geo Data Solutions GDS Inc.
!
! 
!  Input file name: MTSphere.cfl
!  Output file names: MTSphere.out, MTSphere.mf1
!
!  MTSphere.cfl:
!    
!  1 : Comment line
!
!  2 :  NLYR, NTERMS, NX, NY, NZ
!       NLYR : Layer number
!       NTERMS : Number of spherical harmonics computed ( 6 recommended )
!       ITM : If ITM = 1, horizontal transverse magnetic solution is included.
!       NX : Number of X values
!       NY : Number of Y values
!       ND : Number of D values
!       IF NX and NY and ND > 0, XMIN, XMAX, YMIN, YMAX, DMIN and DMAX are provided.
!       Otherwise, all X, Y, and D locations must be provided.
!
!  3 : NF, MINFREQ, MAXFREQ, LOGARITHMIC
!      NF : Frequency number
!      MINFREQ : Minimum frequency
!      MAXFREQ : Maximum frequency
!      LOGARITHMIC: 0 (linear) or 1 (logarithmic) frequency spacing
!
!  4 : RADIUS, DEPTH, SPHRES
!    
!  5.1-5.NLYR-1 : RES, THK
!    
!  5.NLYR : RES(NLYR)
!
!  6 : XMIN XMAX OR X1, X2 ... XN
!
!  7 : YMIN YMAX OR Y1, Y2 ... YN
!    
!  8 : DMIN DMAX OR D1, D2 ... DN
!
!  Note: the sphere is located at X=0, Y=0, D=DEPTH.
! 
    
MODULE INPUT_DATA_FOR_MTSPHERE

IMPLICIT NONE

    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NR,NW,NW1,NW2,NLG,NSPH,LOGARITHMIC,NF,NLYR,QQDT(8),QQHMS(2),NTERMS, &
        NX, NY, ND, JX, JY, JD
    REAL MINFREQ, MAXFREQ, RADIUS, DEPTH, SPHRES, XMIN, XMAX, YMIN, YMAX, DMIN, DMAX, DX, DY, DD
    REAL, DIMENSION(:),ALLOCATABLE :: RES, THK, FREQ
    CHARACTER(LEN=10) TIME,DATE,ZONE    
    CHARACTER(LEN=3) MONTH(12)
    CHARACTER(LEN=60) PVC
    CHARACTER(LEN=120) INP,TITLE
    DATA MONTH /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
    DATA PVC /'MTSphere - Version 2.2.0 - September 2026'/
    REAL(KIND=QL), ALLOCATABLE :: X(:), Y(:), D(:), PXYD(:,:,:,:), APPRES(:,:,:,:,:,:), PHASE(:,:,:,:,:,:), DPTHL(:)
    COMPLEX(KIND=QL), ALLOCATABLE :: K(:,:), Z(:,:), E(:,:,:),IMP(:,:,:,:,:,:), &
                    ES(:,:,:,:,:,:),HS(:,:,:,:,:,:),ET(:,:,:,:,:,:),HT(:,:,:,:,:,:), &
                    PSIAT(:,:,:), PSIFT(:,:,:)
    LOGICAL POSITION
    
    CONTAINS
    
    SUBROUTINE READ_PARAMETERS
    
    ! Subroutine to read parameters from an ASCII file
    
    INTEGER JF,I
    
    NR =  3     !  Input unit number for MTSphere.cfl
    NW =  4     !  Output unit number for MTSphere.out
    NLG = 9     !  Log file unit number for MTSphere.log
    NW1 = 14    !  Output unit number for MTSphere.mf1 (apparent resistivity and phase)
    NW2 = 15    !  Output unit number for MTSphere.mf2 (fields)
    
    OPEN(NR,FILE = 'MTSphere.cfl',STATUS = 'OLD')
    OPEN(NW,FILE = 'MTSphere.out',STATUS = 'REPLACE') 
    OPEN(NW1,FILE = 'MTSphere.mf1',STATUS = 'REPLACE')
    OPEN(NW2,FILE = 'MTSphere.mf2',STATUS = 'REPLACE')
    
    WRITE(NW,1) PVC
    WRITE(*,1) PVC
    WRITE(NW,'(T11,A/T11,A/)') 'INPUT DATA', '----------'
    REFLECT_DATA: DO JF = 1,10000
    READ(NR,'(A)',END = 100) INP
    WRITE(NW,'(1X,A)') INP
    END DO REFLECT_DATA

    100 REWIND NR
        WRITE(NW,2)
        
    READ(NR,'(A)') TITLE
    WRITE(NW,'(/1X,A)') TRIM (TITLE)
        
    READ(NR,*) NLYR, NTERMS, NX, NY, ND
    if ( NX < 0 .OR. NY < 0 .OR. ND < 0 ) then
        POSITION = .TRUE.
        NX = ABS(NX)
        NY = ABS(NY)
        ND = ABS(ND)
    else
        POSITION = .FALSE.
    end IF
    READ(NR,*) NF, MINFREQ, MAXFREQ, LOGARITHMIC
    READ(NR,*) RADIUS, DEPTH, SPHRES
    
    ALLOCATE ( RES(NLYR), THK(NLYR-1), FREQ(NF), X(NX), Y(NY), D(ND), PXYD(NX,NY,ND,3) )
    ALLOCATE ( K(NF,0:NLYR), Z(NF,0:NLYR), E(NF,0:NLYR,2), APPRES(NF,NX,NY,ND,2,2), DPTHL(NLYR), & 
        PHASE(NF,NX,NY,ND,2,2), IMP(NF,NX,NY,ND,3,2), ES(NF,NX,NY,ND,3,2), HS(NF,NX,NY,ND,3,2), &
        PSIAT(2,NTERMS,-1:1),PSIFT(2,NTERMS,-1:1),ET(NF,NX,NY,ND,3,2), HT(NF,NX,NY,ND,3,2) )
    
    WRITE(NW,*)'Number of layers:',NLYR
    WRITE(NW,*)'Number of frequencies:',NF
    DO I = 1, NLYR - 1
        READ(NR,*)RES(I),THK(I)
        WRITE(NW,*)'Layer ',I,' : Resistivity:',RES(I),', Thickness:',THK(I)
    END DO
    READ(NR,*)RES(NLYR)
    WRITE(NW,*)'Layer ',NLYR,' : Resistivity:',RES(NLYR)
    IF ( LOGARITHMIC == 0 ) THEN
        DO I = 1, NF
            FREQ(I) = MINFREQ + ( I - 1 ) * ( MAXFREQ - MINFREQ ) / ( NF - 1 )
        END DO
    ELSE
        DO I = 1, NF
            FREQ(I) = MINFREQ * ( MAXFREQ / MINFREQ ) ** ( FLOAT ( I - 1 ) / FLOAT ( NF - 1 ) )
        END DO
    END IF
    WRITE(NW,'('' Frequencies:'',100G15.7)')(FREQ(I),I=1,NF)
    
    if ( POSITION ) THEN
        read(NR,*)X
        read(NR,*)Y
        read(NR,*)D
    else
        READ(NR,*)XMIN, XMAX
        READ(NR,*)YMIN, YMAX
        read(NR,*)DMIN, DMAX
        if ( NX > 1 ) THEN
            DX = ( XMAX - XMIN ) / ( NX - 1 )
        else
            DX = 0
        end if
        DO JX = 1, NX
            X(JX) = XMIN + ( JX - 1 ) * DX
        end DO
        if ( NY > 1 ) THEN
            DY = ( YMAX - YMIN ) / ( NY - 1 )
        else
            dy = 0
        end if
        do JY = 1, NY
            Y(JY) = YMIN + ( JY - 1 ) * DY
        end do
        if ( ND > 1 ) THEN
            DD = ( DMAX - DMIN ) / ( ND - 1 )
        else
            DD = 0
        end if
        do JD = 1, ND
            D(JD) = DMIN + ( JD - 1 ) * DD
        end do
    end IF
        
    PXYD = 0._QL
    DO JX = 1, NX
        DO JY = 1, NY
            do JD = 1, ND
                PXYD(JX,JY,JD,1) = X(JX)
                PXYD(JX,JY,JD,2) = Y(JY)
                PXYD(JX,JY,JD,3) = D(JD)
            end DO
        END DO
    END do
        
  1 FORMAT (T25,A/T25,'Developed by: Marc A. Vallee'/T25,'for: Geo Data Solutions GDS Inc.'///)   
  2 FORMAT (T1,79('-'))    
    
    END SUBROUTINE READ_PARAMETERS

    END MODULE INPUT_DATA_FOR_MTSPHERE
    
PROGRAM MAIN

! Main program that call all the required subroutines

USE INPUT_DATA_FOR_MTSPHERE

IMPLICIT NONE

    integer JZ
    REAL CMP_START, CMP_END, ELAPSED
    
     CALL CPU_TIME (CMP_START)
     CALL READ_PARAMETERS
 
    CALL PLANEWAVEIMPEDANCE(NW, NLYR, THK, RES, NF, FREQ, K, Z, E)
    
    ES = (0._QL,0._QL)
    HS = (0._QL,0._QL)
    DPTHL = 0._QL
    DO JZ = 1, NLYR-1
        DPTHL(JZ+1) = DPTHL(JZ) + THK(JZ)
    END DO
   
    CALL MTSPHERE3D(NW, NF, NLYR, NTERMS, NX, NY, ND, FREQ, PXYD, THK, RES, DEPTH, RADIUS, SPHRES, E, DPTHL, ES, HS)    
    CALL GETIMPEDANCE(NF,NX,NY,ND,NLYR,NTERMS,PXYD,DPTHL,K,Z,DEPTH,RADIUS,E,ES,HS,IMP,ET,HT)
    CALL APPARENTRESISTIVITY(NF, NX, NY, ND, FREQ, IMP, APPRES, PHASE )
    CALL WRITERESULTS(NW, NW1, NW2, NF, NX, NY, ND, FREQ, PXYD, APPRES, PHASE, IMP, ES, HS, ET, HT)
    
    deallocate ( RES, THK, FREQ, PXYD, K, Z, DPTHL, E, APPRES, PHASE, IMP, ES, HS, ET, HT )
   
    CALL DATE_AND_TIME (DATE, TIME, ZONE, QQDT)
    QQHMS(1:2) = QQDT(5:6)

    CALL CPU_TIME (CMP_END)
    ELAPSED = CMP_END - CMP_START   

   WRITE(NW,98) QQHMS(1:2),QQDT(3),MONTH(QQDT(2)),QQDT(1),ELAPSED
   WRITE(*,98)  QQHMS(1:2),QQDT(3),MONTH(QQDT(2)),QQDT(1),ELAPSED
   
   CLOSE(NW)
   STOP
   
 98 FORMAT(//T3,'MTSphere model completed at ',I2.2,':',I2.2,' on',I3.2,1X,A,I5 &
           //T3,'Computation time = ',F10.2,' seconds.'//)
   
    END PROGRAM MAIN
    
SUBROUTINE WRITERESULTS(NW, NW1, NW2, NF, NX, NY, ND, FREQ, PXYD, APPRES, PHASE, IMP, ES, HS, ET, HT)

! Subroutine to write all the results in an ASCII file

!          Input
!          -----
! NW: main printing unit
! NW1: impedance and phase printing unit
! NW2: field printing unit
! NF: number of frequencies
! NX: number of X locations
! NY: number of Y locations
! FREQ: frequency vector
! PXY: location array (NX,NY)
! APPRES: apparent resistivity
! PHASE: ground surface phase
! IMP: ground surface impedance
! ET: total electric field
! HT: total magnetic field

    IMPLICIT NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NW, NW1, NW2, NF, NX, NY, ND, JF, JX, JY, JD, I, J
    REAL FREQ(NF)
    REAL(KIND=QL) PXYD(NX,NY,ND,3), APPRES(NF,NX,NY,ND,2,2), PHASE(NF,NX,NY,ND,2,2)
    COMPLEX(KIND=QL) IMP(NF,NX,NY,ND,3,2), ES(NF,NX,NY,ND,3,2), HS(NF,NX,NY,ND,3,2), &
                                        ET(NF,NX,NY,ND,3,2), HT(NF,NX,NY,ND,3,2)
        
    DO I = 1, 2
        DO J = 1, 2
            SELECT CASE(I)
            CASE (1)
                SELECT CASE(J)
                CASE(1)
                    WRITE(NW,'(/'' XX Apparent resistivity, phase and impedance'')')
                    WRITE(NW1,'(/'' XX Apparent resistivity, phase and impedance'')')
               CASE(2)
                    WRITE(NW,'(/'' XY Apparent resistivity, phase and impedance'')')
                    WRITE(NW1,'(/'' XY Apparent resistivity, phase and impedance'')')
                END SELECT
            CASE (2)
                SELECT CASE(J)
                CASE(1)
                    WRITE(NW,'(/'' YX Apparent resistivity, phase and impedance'')')
                    WRITE(NW1,'(/'' YX Apparent resistivity, phase and impedance'')')
                CASE(2)
                    WRITE(NW,'(/'' YY Apparent resistivity, phase and impedance'')')
                    WRITE(NW1,'(/'' YY Apparent resistivity, phase and impedance'')')
                END SELECT
            END SELECT
            WRITE(NW,100)
            WRITE(NW1,100)
            DO JF = 1, NF
                do JD = 1, ND
                    DO JX = 1, NX
                        DO JY = 1, NY
                            WRITE(NW,'(G15.7,3F10.1,10G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),APPRES(JF,JX,JY,JD,I,J), &
                                                                      PHASE (JF,JX,JY,JD,I,J), &
                                                                      IMP(JF,JX,JY,JD,I,J)
                            WRITE(NW1,'(G15.7,3F10.1,10G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),APPRES(JF,JX,JY,JD,I,J), &
                                                                      PHASE (JF,JX,JY,JD,I,J), &
                                                                      IMP(JF,JX,JY,JD,I,J)
                        END DO
                    END DO
                END DO
            end do
        END DO
    END DO

    WRITE(NW,'(/'' Tipper'')')
    WRITE(NW1,'(/'' Tipper'')')
    WRITE(NW,120)
    WRITE(NW1,120)
    DO JF = 1, NF
        do JD = 1, ND
            DO JX = 1, NX
                DO JY = 1, NY
                    WRITE(NW,'(G15.7,3F10.1,10G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),IMP(JF,JX,JY,JD,3,:)
                    WRITE(NW1,'(G15.7,3F10.1,10G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),IMP(JF,JX,JY,JD,3,:)
                END DO
            END DO
        END DO
    end do
    
      DO I = 1, 2
        SELECT CASE(I)
        CASE (1)
            WRITE(NW,'(/'' X polarization electric and magnetic fields'')')
            WRITE(NW2,'(/'' X polarization electric and magnetic fields'')')
        CASE(2)
            WRITE(NW,'(/'' Y polarization electric and magnetic fields'')')
            WRITE(NW2,'(/'' Y polarization electric and magnetic fields'')')
        END SELECT
        WRITE(NW,140)
        WRITE(NW2,140)
        write(NW,'(/''Secondary'')')
        write(nw2,'(/''Secondary'')')
        DO JF = 1, NF
            do JD = 1, ND
                DO JX = 1, NX
                    DO JY = 1, NY
                        WRITE(NW,'(G15.7,3F10.1,22G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),&
                            (ES(JF,JX,JY,JD,J,I),J=1,3),(HS(JF,JX,JY,JD,J,I),J=1,3)
                        WRITE(NW2,'(G15.7,3F10.1,22G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),&
                            (ES(JF,JX,JY,JD,J,I),J=1,3),(HS(JF,JX,JY,JD,J,I),J=1,3)
                    END DO
                END DO
            END DO
        end do
        write(NW,'(/''Total'')')
        write(Nw2,'(/''Total'')')
        DO JF = 1, NF
            do JD = 1, ND
                DO JX = 1, NX
                    DO JY = 1, NY
                        WRITE(NW,'(G15.7,3F10.1,22G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),&
                            (ET(JF,JX,JY,JD,J,I),J=1,3),(HT(JF,JX,JY,JD,J,I),J=1,3)
                        WRITE(NW2,'(G15.7,3F10.1,22G15.7)')FREQ(JF),PXYD(JX,JY,JD,:),&
                            (ET(JF,JX,JY,JD,J,I),J=1,3),(HT(JF,JX,JY,JD,J,I),J=1,3)
                    END DO
                END DO
            END DO
        end do
    END DO  
    
100 FORMAT(T7,'Frequency',T25,'X',T35,'Y',T45,'D',T48,'Apparent_res.',T67,'Phase',T81,'Imp.(real)',T96,'Imp.(imag)')
120 FORMAT(T7,'Frequency',T25,'X',T35,'Y',T45,'D',T52,'Kzx(real)',T67,'Kzx(imag)',T82,'Kzy(real)',T97,'Kzy(imag)')
    140 FORMAT(T7,'Frequency',T25,'X',T35,'Y',T45,'D',T52,'Ex(real)',T67,'Ex(imag)',T82,'Ey(real)',T97,'Ey(imag)', &
    T112,'Ez(real)',T127,'Ez(imag)', &
    T142,'Hx(real)',T157,'Hx(imag)',T172,'Hy(real)',T187,'Hy(imag)',T202,'Hz(real)',T217,'Hz(imag)')

END SUBROUTINE WRITERESULTS
  
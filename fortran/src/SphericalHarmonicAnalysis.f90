  SUBROUTINE SphericalGeometryInitialisation(Centre,Radius,nlat,nlon,PS)
    
! This subroutine initializes the position matrix for a given sphere geometry

!          Input
!          -----
! Centre: location of the center of the sphere in cartesian coordinates
! Radius: radius of the sphere
! nlat: number of latitudes
! nlon: number of longitudes

!          Output
!          -----
! PS: matrix with the location of regularly distributed points on the
!     surface of the sphere

    use precision
    
    implicit none
    integer nlat,nlon,ilat,ilon
    REAL(KIND=QL) centre(3),ps(nlat,nlon,3),radius,theta,phi,xgauss(nlat),wgauss(nlat)
   
    IF ( NLAT == 1 ) THEN
        ps(1,1,:) = centre
    else    
        call gauss_legendre(nlat,xgauss,wgauss)
        do ilat = 1,nlat
            theta = ACOS(xgauss(ilat))
            do ilon = 1,nlon
                phi = (ilon-1) * 2._ql * pi / nlon
                ps(ilat,ilon,1) = centre(1) + radius * sin(theta) * cos(phi)
                ps(ilat,ilon,2) = centre(2) + radius * sin(theta) * sin(phi)
                ps(ilat,ilon,3) = centre(3) + radius * cos(theta)
            end do 
        end do
    end if
    
END SUBROUTINE SphericalGeometryInitialisation
    
  SUBROUTINE RadialFields(NLAT,NLON,PS,CENTRE,E,H,ER,HR)
  
! This subroutine computes the radial fields from the cartesian fields.

!          Input
!          -----
! nlat: number of latitudes
! nlon: number of longitudes
! PS: matrix with the location of regularly distributed points on the
!     surface of the sphere
! Centre: location of the center of the sphere in cartesian coordinates
! E: Cartesian components of the electric field at the PS points
! H: Cartesian components of the magnetic field at the PS points

!          Output
!          -----
! ER: Radial component of the electric field
! HR: Radial component of the magnetic field
    
    implicit NONE
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    INTEGER NLAT,NLON,JLAT,JLON
    REAL(KIND=QL) X,Y,Z,R,THETA,PHI,CENTRE(3),PS(NLAT,NLON,3)
    COMPLEX(KIND=QL) E(NLAT,NLON,3),ER(NLAT,NLON),H(NLAT,NLON,3),HR(NLAT,NLON)
       
    DO JLAT = 1,NLAT
        DO JLON = 1,NLON
            X = PS(JLAT,JLON,1) - CENTRE(1)
            Y = PS(JLAT,JLON,2) - CENTRE(2)
            Z = PS(JLAT,JLON,3) - CENTRE(3)
            R = SQRT(X*X+Y*Y+Z*Z)
            THETA = ACOS ( Z / R )
            PHI = ATAN2(Y,X)
            ER(JLAT,JLON) = SIN(THETA)*COS(PHI)*E(JLAT,JLON,1) + &
                            SIN(THETA)*SIN(PHI)*E(JLAT,JLON,2) + &
                            COS(THETA)         *E(JLAT,JLON,3)
            HR(JLAT,JLON) = SIN(THETA)*COS(PHI)*H(JLAT,JLON,1) + &
                            SIN(THETA)*SIN(PHI)*H(JLAT,JLON,2) + &
                            COS(THETA)         *H(JLAT,JLON,3)
            
        END DO
    END DO

END SUBROUTINE RadialFields
    
SUBROUTINE MTSphericalHarmonicAnalysis(nterms,nlat,nlon,radius,yhat,zhat,er,hr,singular,psia,psif)

! This subroutine computes the spherical coefficients of an electromagnetic field from the 
! radial components using spherical harmonic analysis

!          Input
!          -----
! nterms: degree number
! nlat: number of latitudes
! nlon: number of longitudes
! radius: radius of the sphere
! yhat: background admittivity
! zhat: background impedivity
! er: radial electric field matrix
! hr: radial magnetic field matrix
! singular: select singular or regular solution

!          Output
!          -----
! psia: transverse magnetic potential coefficients
! psif: transverse electric potential coefficients

    use precision

    implicit none 
    integer nterms,nlat,nlon,jlat,jlon,n,m,nm
    real(kind=ql) radius
    complex(kind=ql) yhat,zhat,fact,er(nlat,nlon),hr(nlat,nlon), &
        psia(NTERMS,-1:1),psif(NTERMS,-1:1)
    DOUBLE COMPLEX KR,JKR(0:NTERMS),DJKR(0:NTERMS), &
        YKR(0:NTERMS),DYKR(0:NTERMS),HKR(0:NTERMS),DHKR(0:NTERMS)
    LOGICAL singular
    
    call MTSphericalAnalysis(nterms,nlat,nlon,ER,psia)
    call MTSphericalAnalysis(nterms,nlat,nlon,Hr,psif)
    
    do n = 1, nterms
        psia(n,:) = psia(n,:) * radius * yhat / ( n * ( n + 1 ) )
        psif(n,:) = psif(n,:) * radius * Zhat / ( n * ( n + 1 ) )
    end do

    kr = sqrt ( - yhat * zhat ) * radius
    call csphjy(NTERMS,KR,NM,JKR,DJKR,YKR,DYKR,HKR,DHKR)
    do n = 1,NTERMS
        IF (singular) THEN
            PSIA(N,:) = PSIA(N,:) / HKR(N)
            PSIF(N,:) = PSIF(N,:) / HKR(N)
        ELSE
            PSIA(N,:) = PSIA(N,:) / JKR(N)
            PSIF(N,:) = PSIF(N,:) / JKR(N)
        END IF
    end do
    
    end subroutine MTSPHERICALHARMONICANALYSIS
    
    subroutine MTSphericalAnalysis(nterms,nlat,nlon,f,psi)
    
    use precision
    use sphericalfunctions 
    
    implicit none
    integer nterms,nlat,nlon,n,m,jlat,jlon
    real(kind=QL) theta, phi, deltalon, deltalat, xgauss(nlat), wgauss(nlat)
    complex(kind=QL) f(nlat,nlon),psi(nterms,-1:1), integral, term, f1, f2
    
    call gauss_legendre(nlat, xgauss, wgauss)
    psi = (0._QL,0._QL)
    deltalon = 2._QL * PI / NLON
    do n = 1, nterms
        do M = -1, 1, 2
            integral = (0.0d0, 0.0d0)
            do jlon = 1, nlon
                phi = ( jlon - 1 ) * 2._QL * PI / NLON
                do jlat = 1, nlat
                    deltalat = wgauss(jlat)
                    theta = ACOS(xgauss(jlat))
                    integral = integral + f(jlat,jlon) * YNM(N,-M,theta,phi) * deltalon * deltalat
                end do
            end do
            psi(n,M) = integral
        end do
    end do

    end subroutine MTSphericalAnalysis
    
    subroutine gauss_legendre(n, x, w)
  
   ! Subroutine to calcule nodes and weight of Gauss-Legendre integration
  
  !         Input
  !         -----
  !     x1: minimum
  !     x2: maximum
  !     n : number of nodes
  !         Output
  !         ------
  !     x : node locations
  !     w : node weights.
  
    implicit none
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80) 
    real(kind=ql), parameter :: pi = 3.141592653589793_ql    
    integer, intent(in) :: n
    real(8), intent(out) :: x(n), w(n)
    integer :: i, j, iter
    real(8) :: z, p1, p2, p3, pp, z1
    real(8), parameter :: eps = 1.0d-15
    do i = 1, (n + 1) / 2
       z = cos(pi * (i - 0.25d0) / (n + 0.5d0))
       do iter = 1, 100
          p1 = 1.0d0
          p2 = 0.0d0
          do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2.0d0 * j - 1.0d0) * z * p2 - (j - 1.0d0) * p3) / j
          end do
          pp = n * (z * p1 - p2) / (z**2 - 1.0d0)
          z1 = z
          z = z1 - p1 / pp
          if (abs(z - z1) < eps) exit
       end do
       x(i) = -z
       x(n + 1 - i) = z
       w(i) = 2.0d0 / ((1.0d0 - z**2) * pp**2)
       w(n + 1 - i) = w(i)
    end do
    
end subroutine gauss_legendre
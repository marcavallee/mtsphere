module precision
    implicit none
    integer, parameter :: QL=SELECTED_REAL_KIND(12,80)    
    real(kind=ql), parameter :: pi = 3.141592653589793_ql    
    real(kind=QL), parameter :: mu0 = 4._QL * PI * 1.0d-7
    real(kind=QL), parameter :: c_light = 299792458.0d0
    real(kind=QL), parameter :: eps0 = 1._QL / ( mu0 * C_light * C_light )
    COMPLEX(KIND=QL), PARAMETER :: CI = CMPLX (0.D0, 1.D0, KIND=QL) 
    COMPLEX(KIND=QL), PARAMETER :: ONE=(1._QL,0._QL)
    COMPLEX(KIND=QL), PARAMETER :: ZERO=(0._QL,0._QL)
end module precision
FUNCTION Gamma_CR(zeta,y_a,y_m,y_He)
    !     calculates heating rate (ergs/s/g) due to CR
    IMPLICIT REAL*8(a-h,o-z)
    !double precision :: zeta,y_a,y_m,y_He
    Gamma_CR=3.26d12*(0.46d0*y_a+0.50d0*y_He+0.94d0*y_m)*zeta/(1.d0+4.d0*y_He)
    return
END FUNCTION Gamma_CR
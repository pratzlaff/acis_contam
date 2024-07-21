      SUBROUTINE runpow(ear, ne, param, ifl, photar)
      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne)

      integer nemax
      parameter (nemax=100000)

C---
      double precision gamma, beta, e, alpha
      integer ie

      gamma = param(1)
      beta = param(2)

      DO ie = 1, ne
        e = 0.5*(ear(ie-1)+ear(ie))
        alpha = gamma + beta*log(e)
        photar(ie) = (ear(ie)-ear(ie-1))*e**(alpha)
      ENDDO

      RETURN
      END

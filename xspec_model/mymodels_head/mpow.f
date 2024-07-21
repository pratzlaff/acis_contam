      SUBROUTINE mpow(ear, ne, param, ifl, photar)
      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne)

      integer nemax
      parameter (nemax=100000)

C---
      double precision gamma,e
      integer ie

      gamma = param(1)


      DO ie = 1, ne
        e = 0.5*(ear(ie-1)+ear(ie))
        photar(ie) = e**(-gamma)
      ENDDO

      RETURN
      END

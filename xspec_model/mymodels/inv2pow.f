      SUBROUTINE inv2pow (ear, ne, param, ifl, photar)
      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne)

      integer nemax
      parameter (nemax=100000)

C---
      double precision gamma1,gamma2,norm1,norm2,e, de
      integer ie

      gamma1 = param(1)
      gamma2 = param(2)
      norm2  = param(3)

      DO ie = 1, ne
        e = 0.5*(ear(ie-1)+ear(ie))
        de = ear(ie)-ear(ie-1)
        photar(ie) = de/(e**gamma1+norm2*e**gamma2)
      ENDDO

      RETURN
      END

      SUBROUTINE exptail(ear, ne, param, ifl, photar)
      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne)

      integer nemax
      parameter (nemax=100000)

C---
      double precision e0, etau, norm, e
      integer ie

      e0 = param(1)
      etau = param(2)

c Trap the case of etau=0

      IF ( param(2) .LE. 0.0 ) THEN
        DO ie = 1, ne
          if (param(1).gt.ear(ie-1).and.param(1).le.ear(ie)) then
            photar(ie) = 1.
          else
            photar(ie) = 0.
          endif
        ENDDO
        RETURN
      ENDIF

      norm = etau*(1-exp(-e0/etau))
      norm = 1/norm

      DO ie = 1, ne
        e = (0.5*(ear(ie-1)+ear(ie)) - e0)/etau
        if (e.lt.0.0) then
          photar(ie) = norm*(ear(ie)-ear(ie-1))*exp(e)
        else
          photar(ie) = 0
        endif
      ENDDO

      RETURN
      END

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      SUBROUTINE exptai2 (ear, ne, param, ifl, photar)

      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(3), photar(ne)

      integer nemax
      parameter (nemax=100000)

C---
      double precision e0, etau, norm, e, slope
      integer ie

      e0 = param(1)
      etau = param(2)
      slope = param(3)

c Trap the case of etau=0

      IF ( param(2) .LE. 0.0 ) THEN
        DO ie = 1, ne
          if (param(1).gt.ear(ie-1).and.param(1).le.ear(ie)) then
            photar(ie) = 1.
          else
            photar(ie) = 0.
          endif
        ENDDO
        RETURN
      ENDIF

      norm = 0
      DO ie = 1, ne
        e = (0.5*(ear(ie-1)+ear(ie)) - e0)/etau
        if (e.lt.0.0) then
          photar(ie) = (ear(ie)-ear(ie-1))*exp(-abs(e)**slope)
        else
          photar(ie) = 0
        endif
        norm = norm + photar(ie)
      ENDDO

      do ie=1,ne
        photar(ie) = photar(ie)/norm
      enddo

      RETURN
      END

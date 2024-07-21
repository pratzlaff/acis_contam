      subroutine locateAV(xx,n,x,j)
c      
c Given an array xx(n) and given value of x returns j such as x 
c between xx(j) and xx(j+1). xx(n) must be monotonic, either increasing 
c or decreasing. j=0 or j=n return to indicate that x out of range. 
c
c--------
      implicit none
      integer n,j
      real xx(n),x
      integer jl, ju, jm
      
      jl=0
      ju=n+1
      do while (ju-jl.gt.1)
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      enddo
      j=jl
      return
      end


      subroutine locatei4(xx,n,x,j)
c      
c Given an array xx(n) and given value of x returns j such as x 
c between xx(j) and xx(j+1). xx(n) must be monotonic, either increasing 
c or decreasing. j=0 or j=n return to indicate that x out of range. 
c
c--------
      implicit none
      integer n,j
      integer xx(n),x
      integer jl, ju, jm
      
      jl=0
      ju=n+1
      do while (ju-jl.gt.1)
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      enddo
      j=jl
      return
      end




cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine locate8(xx,n,x,j)
c      
c Given an array xx(n) and given value of x returns j such as x 
c between xx(j) and xx(j+1). xx(n) must be monotonic, either increasing 
c or decreasing. j=0 or j=n return to indicate that x out of range. 
c
c--------
      implicit none
      integer n,j
      real*8 xx(n),x
      integer jl, ju, jm
      
      jl=0
      ju=n+1
      do while (ju-jl.gt.1)
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      enddo
      j=jl
      return
      end


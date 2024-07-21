      subroutine huntAV(xx,n,x,jlo)
c
c      
c Given an array x(n) and given value of x returns jlo such as x 
c between xx(jlo) and xx(jlo+1). xx(n) must be monotonic, either increasing 
c or decreasing. jlo=0 or jlo=n return to indicate that x out of range. 
c jlo on input is taken as initial guess for jlo on output.
c
c
c
      implicit none
      integer n,jlo
      real xx(n),x
      logical ascnd
      integer jhi,inc,jm
      
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1      jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
 2      jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine huntAV8(xx,n,x,jlo)
c
c      
c Given an array x(n) and given value of x returns jlo such as x 
c between xx(jlo) and xx(jlo+1). xx(n) must be monotonic, either increasing 
c or decreasing. jlo=0 or jlo=n return to indicate that x out of range. 
c jlo on input is taken as initial guess for jlo on output.
c
c
c
      implicit none
      integer n,jlo
      real*8 xx(n),x
      logical ascnd
      integer jhi,inc,jm
      
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1      jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
 2      jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end

      
      subroutine hunti(xx,n,x,jlo)
c
c      
c Given an array x(n) and given value of x returns jlo such as x 
c between xx(jlo) and xx(jlo+1). xx(n) must be monotonic, either increasing 
c or decreasing. jlo=0 or jlo=n return to indicate that x out of range. 
c jlo on input is taken as initial guess for jlo on output.
c
c
c
      implicit none
      integer n,jlo
      integer xx(n),x
      logical ascnd
      integer jhi,inc,jm
      
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1      jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
 2      jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end








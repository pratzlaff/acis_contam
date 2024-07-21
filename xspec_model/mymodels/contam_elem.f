      subroutine contamCOF (ear,ne,param,ifl,photar)
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl
      real contam_elem_func

      real tauC, tauF, tauO, e, eold, eref, tauref, tau0
      integer i, ibin
      real w1, w2

      tauref = param(1)
      tauC = 1
      tauO = param(2) ! OtoC
      tauF = param(3) ! FtoC
      eref = param(4)

      ibin = -1
      eold = -1
      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        if (eref.gt.eold.and.eref.le.e) then
          ibin = i
          call get_weights(eref,eold,e,w1,w2)
        endif
        photar(i) = contam_elem_func (e,tauC, tauO, tauF)
        eold = e
      enddo

      if (ibin.le.1) then
        write (0,*) 
     ~       'ERROR: REFERENCE ENERGY OUTSIDE ARF BOUNDS in CONTAM_ELEM'
        call exit(1)
      else
        tau0 = w1*photar(ibin-1)+w2*photar(ibin)
      endif
      
      do i=1,ne
        photar(i) = exp(-photar(i)*tauref/tau0)
      enddo
      
      return
      end

c -----------------------------------------------------------------------
      subroutine contamCOFL (ear,ne,param,ifl,photar)
c just like contam-C-O-F but normalized to the ECS-L line complex
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl
      real contam_elem_func

      real tauC, tauF, tauO, e, eref, tauref, tau0
      integer i
      real e1(100000)

      tauref = param(1)
      tauC = 1
      tauO = param(2) ! OtoC
      tauF = param(3) ! FtoC

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        e1(i)=e
        photar(i) = contam_elem_func (e,tauC, tauO, tauF)
      enddo

      call get_tau_L (e1,photar,ne,tau0)
      
      do i=1,ne
        photar(i) = exp(-photar(i)*tauref/tau0)
      enddo
      
      return
      end

c---------------------------------------------------------------
      real function contam_elem_func (e, tauC, tauO, tauF)
      implicit none
      real e, tauC, tauO, tauF

      logical firstcall
      data firstcall /.true./
      real abs0C(20000), abs0O(20000), abs0F(20000)

      save firstcall, abs0C, abs0O, abs0F

      real ee(10000), tt(10000)
      integer nn
      character*1024 filename, dir
      integer status

      integer ibin

      if (firstcall) then
        firstcall = .false.

        call get_environment_variable ('CONTAM_COF_DIR',dir,
     ~       status=status)
        if (status.ne.0) then
        !  dir='/data/alexey/cal/contam/C-O-F'
          dir='/data/legs/rpete/contam/xspec_model/dat/C-O-F'
        endif

        ! load and process C
        filename = trim(dir)//'/C.dat'
        call contam_elem_load1file (filename,ee,tt,nn)
        call contam_elem_prepare_array (ee,tt,nn,abs0C,20000)

        ! ... and O
        filename = trim(dir)//'/O.dat'
        call contam_elem_load1file (filename,ee,tt,nn)
        call contam_elem_prepare_array (ee,tt,nn,abs0O,20000)

        ! ... and F
        filename = trim(dir)//'/F.dat'
        call contam_elem_load1file (filename,ee,tt,nn)
        call contam_elem_prepare_array (ee,tt,nn,abs0F,20000)

      endif
      
      ibin = nint(e*1000)
      if (ibin.lt.1) then
        contam_elem_func = -10
      elseif (ibin.gt.20000) then
        contam_elem_func = 0
      else
        contam_elem_func = tauC*abs0C(ibin) + 
     ~       tauO*abs0O(ibin) +
     ~       tauF*abs0F(ibin)
      endif

      return
      end

c---------------------------------------------------------------------
      subroutine contam_elem_load1file (filename,ee,tt,nn)
      implicit none
      character*(*) filename
      real ee(*), tt(*)
      integer nn

      integer unit, newunit, status

      unit = newunit()
      open (unit,file=filename,status='old')
      read (unit,*)
      nn = 0
      status = 0
      do while (status.eq.0)
        read (unit,*,iostat=status) ee(nn+1), tt(nn+1)
        if (status.eq.0) then
          nn = nn + 1
        endif
      enddo
      close (unit)
      
      return
      end

c---------------------------------------------------------------------
      subroutine contam_elem_prepare_array (ee,tt,nn,abs0C,n)
      implicit none
      integer nn, n
      real ee(nn), tt(nn), abs0C(n)
      
      integer i, ibin 
      real eee, w1, w2

      do i=1,n
        eee = i*0.001
        if (eee.lt.ee(1)) then
          ibin = 0
        elseif (eee.gt.ee(nn)) then
          ibin = nn
        else
          ibin = min(nn-1,max(1,ibin))
          call locateAV(ee,nn,eee,ibin)
        endif
        if (ibin.le.0) then
          ibin = 1
          w1 = 1
          w2 = 0
        else if (ibin.ge.nn) then
          ibin = nn-1
          w1 = 0
          w2 = 1
        else
          call get_weights (eee,ee(ibin),ee(ibin+1),w1,w2)
        endif
        
        abs0C(i) = -log(w1*tt(ibin)+w2*tt(ibin))
c        print*,eee,abs0C(i)
      enddo
      
      return
      end


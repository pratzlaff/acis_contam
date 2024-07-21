      subroutine contam2func (ear,ne,param,ifl,photar)
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl

      logical firstcall
      data firstcall /.true./
      save firstcall

      integer ne0, ne0max
      parameter (ne0max=1000000)
      real e0(ne0max), tau_elem(ne0max), tau_fluf(ne0max)
      save ne0, e0, tau_elem, tau_fluf

      integer ibin, i
      
      real tau_K, tau_F, w1, w2, e

      if (firstcall) then
        call load_cont_data (e0,tau_elem,tau_fluf,ne0)
        firstcall = .false.
      endif

      tau_K = param(1)
      tau_F = param(2)

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        if (i.eq.1) then
          call locateAV (e0,ne0,e,ibin)
        else
          if (e.lt.e0(1)) then
            ibin = 0
          elseif (e.gt.e0(ne0)) then
            ibin = ne0
          else
            ibin = min(ne0-1,max(1,ibin))
            call huntAV (e0,ne0,e,ibin)
          endif
        endif
        if (ibin.le.0) then
          ibin = 1
          w1 = 1
          w2 = 0
        else if (ibin.ge.ne0) then
          ibin = ne0-1
          w1 = 0
          w2 = 1
        else
          call get_weights (e,e0(ibin),e0(ibin+1),w1,w2)
        endif

        photar(i) = exp(
     ~      -(w1*tau_elem(ibin)+w2*tau_elem(ibin+1))*tau_K
     ~      -(w1*tau_fluf(ibin)+w2*tau_fluf(ibin+1))*tau_F
     ~      )

      enddo

      return
      end

      
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       subroutine load_cont_data (e,tau_elem,tau_fluf,ne)
       implicit none
       integer ne
       real e(*), tau_elem(*), tau_fluf(*)
       character file*200
       integer unit, newunit
       integer status

       integer ibin
       real w1,w2
       real f_F,E_F, tau_F_e, tau_F_f
       real f_Mn,E_Mn, tau_Mn_e, tau_Mn_f
       real tau_L
       parameter (f_F  = 0.6, E_F=0.679)
       parameter (f_Mn = 0.4, E_Mn=0.638)
       integer i

       unit = newunit()

       file = '/data/alexey/cal/contam/tau_K_fluf.dat'
       
       open (unit,file=file,status='old')
       ne = 0
       status = 0
       do while (status.eq.0)
         read (unit,*,iostat=status) e(ne+1), tau_elem(ne+1),
     ~       tau_fluf(ne+1)
         if (status.eq.0) then
           ne = ne + 1
         endif
       enddo

       close(unit)

c Renorm to unit tau for the L-complex
       call locateAV (e,ne,E_F,ibin)
       call get_weights (E_F,e(ibin),e(ibin+1),w1,w2)
       tau_F_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
       tau_F_f = w1*tau_fluf(ibin)+w2*tau_fluf(ibin+1)

       call locateAV (e,ne,E_Mn,ibin)
       call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
       tau_Mn_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
       tau_Mn_f = w1*tau_fluf(ibin)+w2*tau_fluf(ibin+1)

       tau_L = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
       print*,tau_L
       do i=1,ne
         tau_elem(i) = tau_elem(i)/tau_L
       enddo

       tau_L = -log(f_F*exp(-tau_F_f)+f_Mn*exp(-tau_Mn_f))
       do i=1,ne
         tau_fluf(i) = tau_fluf(i)/tau_L
       enddo


       


       return
       end


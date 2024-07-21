      subroutine contam32lev (ear,ne,param,ifl,photar)
c
c Just like contam2, but F is proportional to tau_K, 
c and C and O tied to produce tauK; and extra parameter to specify suppression of 
c the C-L EXAFS structure
c
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl
      real areafrac
      integer i
      real tparam(4)

      areafrac=param(4)
*      print*,param(1),areafrac,(exp(-param(1))-(1-areafrac))/areafrac
      tparam(1)=param(1)
      tparam(2)=param(2)
      tparam(3)=param(3)
      tparam(4)=param(4)

      call contam3func(ear,ne,tparam,ifl,photar)

      do i=1,ne
         photar(i)=areafrac*photar(i)+(1-areafrac)
!         exp(-tauL) = a*exp(-tauL2)+(1-a) 
!         (exp(-tauL)-(1-a))/a = exp(-tauL2)
!         -log((exp(-tauL)-(1-a))/a) = tauL2
      enddo
      return
      end

      subroutine contam3func (ear,ne,param,ifl,photar)
c
c Just like contam2, but F is proportional to tau_K, 
c and C and O tied to produce tauK; and extra parameter to specify suppression of 
c the C-L EXAFS structure
c
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
      real e0(ne0max), tau_elem(ne0max), tau_C(ne0max), tau_O(ne0max), 
     ~ tau_F(ne0max), exafs_C(ne0max)
      real  tau_C0(ne0max), tau_O0(ne0max), 
     ~ tau_F0(ne0max), exafs_C0(ne0max)
      real tau_L_C, tau_L_O
      save ne0, e0, tau_elem, tau_C, tau_O, tau_F, exafs_C
      save tau_C0, tau_O0, tau_F0, exafs_C0, tau_L_C, tau_L_O

      integer ibin, i
      
      real tau_K, w1, w2, e
      
      real ex, OtoC
      
      real ee

      if (firstcall) then
        call load_cont_data_3 (e0,tau_elem,tau_C,tau_O,tau_F,exafs_C,
     ~      ne0,
     ~      tau_L_C,tau_L_O)
        firstcall = .false.

        do i=1,20000
           ee = i/1000.0
           if (ee.lt.e0(1)) then
              ibin = 0
           elseif (ee.gt.e0(ne0)) then
              ibin = ne0
           else
              ibin = min(ne0-1,max(1,ibin))
              call locateAV (e0,ne0,ee,ibin)
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
              call get_weights (ee,e0(ibin),e0(ibin+1),w1,w2)
           endif

           tau_C0(i) = w1*tau_C(ibin)+w2*tau_C(ibin+1)
           tau_O0(i) = w1*tau_O(ibin)+w2*tau_O(ibin+1)
           tau_F0(i) = w1*tau_F(ibin)+w2*tau_F(ibin+1)
           exafs_C0(i) = w1*exafs_C(ibin)+w2*exafs_C(ibin+1)
        enddo
      endif

      tau_K = param(1)
      OtoC  = param(2)
      ex    = param(3)

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        ibin = nint(e*1000)
        if (ibin.lt.1.or.ibin.gt.20000) then
           photar(i)=0
        else
           photar(i)=tau_K/(tau_L_C+OtoC*tau_L_O)*
     ~       (tau_C0(ibin)*((exafs_C0(ibin)-1)*ex+1) + 
     ~         tau_O0(ibin)*OtoC + 
     ~         tau_F0(ibin))
        endif
        photar(i)=exp(-photar(i))
      enddo

      return
      end

      
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine load_cont_data_3 (e,tau_elem,tau_C,tau_O,tau_F,
     ~    exafs_C,ne,
     ~    tau_L_C, tau_L_O)
      implicit none
      integer ne
      real e(*), tau_elem(*), tau_C(*), tau_O(*), tau_F(*), exafs_C(*)
      real tau_L_C, tau_L_O
      character file*200
      integer unit, newunit
      integer status
      
      integer ibin
      real w1,w2
      real f_F,E_F, tau_F_e, tau_F_C, tau_F_O
      real f_Mn,E_Mn, tau_Mn_e, tau_Mn_C, tau_Mn_O
      real tau_L
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      integer i
      
      unit = newunit()
      
      file = '/data/alexey/cal/contam/tau_components.dat'
      
      open (unit,file=file,status='old')
      ne = 0
      status = 0
      do while (status.eq.0)
         read (unit,*,iostat=status) e(ne+1), tau_elem(ne+1),
     ~        tau_C(ne+1), tau_O(ne+1), tau_F(ne+1), exafs_C(ne+1)
         if (status.eq.0) then
            ne = ne + 1
         endif
      enddo
      
      close(unit)
      
c     Renorm to unit tau for the L-complex
      call locateAV (e,ne,E_F,ibin)
      call get_weights (E_F,e(ibin),e(ibin+1),w1,w2)
      tau_F_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
      tau_F_C = w1*tau_C(ibin)+w2*tau_C(ibin+1)
      tau_F_O = w1*tau_O(ibin)+w2*tau_O(ibin+1)
      

      call locateAV (e,ne,E_Mn,ibin)
      call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
      tau_Mn_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
      tau_Mn_C = w1*tau_C(ibin)+w2*tau_C(ibin+1)
      tau_Mn_O = w1*tau_O(ibin)+w2*tau_O(ibin+1)
      
      tau_L = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
      tau_L_C = -log(f_F*exp(-tau_F_C)+f_Mn*exp(-tau_Mn_C))
      tau_L_O = -log(f_F*exp(-tau_F_O)+f_Mn*exp(-tau_Mn_O))
      print*,tau_L, tau_L_C, tau_L_O
      do i=1,ne
         tau_elem(i) = tau_elem(i)/tau_L
         tau_C(i) = tau_C(i)/tau_L/exafs_C(i)
         tau_O(i) = tau_O(i)/tau_L
         tau_F(i) = tau_F(i)/tau_L
      enddo
      tau_L_C = tau_L_C / tau_L
      tau_L_O = tau_L_O / tau_L
      
      return
      end
      

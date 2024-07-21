      subroutine contam542func  (ear,ne,param,ifl,photar)
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl

      real param1(100)
      real e1(100000), tau_L, tau_L_1
      integer i

! param1 = target L absorption
! param2 = sigma

      param1(1)=1 ! C
      param1(2)=param(5) ! exsC
      param1(3)=param(3) ! O
      param1(4)=0 ! Cl
      param1(5)=param(4) ! F
      
      call contam54func(ear,ne,param1,ifl,photar)
      
c Convert back to tau
      do i=1,ne
        e1(i)=0.5*(ear(i-1)+ear(i))
        photar(i)=-log(photar(i))
      enddo
      
c Compute tau_L for C=1 absorption
      call get_tau_L (e1,photar,ne,tau_L)

c Find a target tau_L absorption from equation
c
c t-1/2*sigma**2*t**2=TauL
c
      tau_L_1 = (1-sqrt(1-2*param(2)**2*param(1)))/param(2)**2
!      print*,tau_L, tau_L_1, param(1), param(2)
      do i=1,ne
        photar(i)=photar(i)*tau_L_1/tau_L
        photar(i)=exp(-photar(i)+0.5*param(2)**2*photar(i)**2)
      enddo
      return
      end

      
      
      subroutine contam54func (ear,ne,param,ifl,photar)
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
      parameter (ne0max=100000)
      real e00(0:ne0max), www(ne0max)
      real e0(ne0max), tau_elem(ne0max), tau_C(ne0max), tau_O(ne0max), 
     ~ tau_F(ne0max), exafs_C(ne0max), tau_Cl(ne0max)
      real  tau_C0(ne0max), tau_O0(ne0max), 
     ~ tau_F0(ne0max), exafs_C0(ne0max), tau_Cl0(ne0max)
      real tau_L_C, tau_L_O, tau_L_Cl, tau_L_elem
      save ne0, e0, tau_elem, tau_C, tau_O, tau_F, exafs_C
      save tau_C0, tau_O0, tau_F0, exafs_C0, tau_Cl0
      save tau_L_C, tau_L_O, tau_L_Cl, tau_L_elem

      integer ibin, i
      
      real w1, w2, e
      
      real exC, tauC, tauO, tauCl, Ffact
      real param0(30)
      integer ifl0
      
      real ee
      integer unit, newunit

      if (firstcall) then
        call load_cont_data_4 (e0,tau_elem,tau_C,tau_O,tau_F,exafs_C,
     ~      tau_Cl,ne0,
     ~      tau_L_C,tau_L_O, tau_L_Cl, tau_L_elem)
        firstcall = .false.

*c       
*        e00(0)=e0(1)
*        do i=1,ne0
*          e00(i)=e0(i)
*        enddo
*        do i=1,18
*          param0(i)=0
*        enddo
*        param0(12)=50
*        ifl0=1
*        call xsabsv(e00, ne0, param0, ifl0, tau_Cl0, www)
*        do i=1,20000
*          tau_Cl0(i)=-log(tau_Cl0(i))
*        enddo
*        unit = newunit()
*        open (unit,file='/data/alexey/cal/contam/cl.dat')
*        do i=1,ne0
*          write (unit,*) e0(i),tau_cl0(i)
*        enddo
*        close(unit)
*c

        e00(0)=0.001
        do i=1,20000
           ee = i/1000.0
           e00(i)=ee
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
*           tau_N0(i) = w1*tau_N(ibin)+w2*tau_N(ibin+1)
        enddo
        
        do i=1,18
          param0(i)=0
        enddo
        param0(12)=50
        ifl0=1
        call xsabsv(e00, 20000, param0, ifl0, tau_Cl0, www)
        do i=1,20000
           tau_Cl0(i)=-log(tau_Cl0(i))
        enddo
        call get_tau_L (e00(1),tau_Cl0,20000,tau_L_Cl)
        print*,'tau_L_cl=',tau_L_Cl
      endif

      tauC = param(1)
      exC   = param(2)
      tauO  = param(3)
      tauCl = param(4)
      Ffact = param(5)

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        ibin = nint(e*1000)
        if (ibin.lt.1.or.ibin.gt.20000) then
           photar(i)=0
        else
           photar(i)=tauC/tau_L_C*
     ~       tau_C0(ibin)*((exafs_C0(ibin)-1)*exC+1) +
     ~       tauO/tau_L_O*tau_O0(ibin) +
     ~       tauCl/tau_L_Cl*tau_Cl0(ibin) +
     ~          (tauC+tauO)*Ffact*tau_F0(ibin)/tau_L_elem
        endif
        photar(i)=exp(-photar(i))
      enddo

      return
      end

      
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       
      subroutine contam4func (ear,ne,param,ifl,photar)
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
     ~ tau_F(ne0max), exafs_C(ne0max), tau_N(ne0max)
      real  tau_C0(ne0max), tau_O0(ne0max), 
     ~ tau_F0(ne0max), exafs_C0(ne0max), tau_N0(ne0max)
      real tau_L_C, tau_L_O, tau_L_N, tau_L_elem
      save ne0, e0, tau_elem, tau_C, tau_O, tau_F, exafs_C
      save tau_C0, tau_O0, tau_F0, exafs_C0, tau_N0
      save tau_L_C, tau_L_O, tau_L_N, tau_L_elem

      integer ibin, i
      
      real tau_K, w1, w2, e
      
      real ex, OtoC, NtoC
      
      real ee

      if (firstcall) then
        call load_cont_data_4 (e0,tau_elem,tau_C,tau_O,tau_F,exafs_C,
     ~      tau_N,ne0,
     ~      tau_L_C,tau_L_O, tau_L_N, tau_L_elem)
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
           tau_N0(i) = w1*tau_N(ibin)+w2*tau_N(ibin+1)
           exafs_C0(i) = w1*exafs_C(ibin)+w2*exafs_C(ibin+1)
        enddo
      endif

      tau_K = param(1)
      OtoC  = param(2)
      NtoC  = param(3)
      ex    = param(4)
      

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        ibin = nint(e*1000)
        if (ibin.lt.1.or.ibin.gt.20000) then
           photar(i)=0
        else
           photar(i)=(tau_K-OtoC-NtoC)/tau_L_C*
     ~       tau_C0(ibin)*((exafs_C0(ibin)-1)*ex+1) +
     ~       (tau_K-NtoC)*tau_F0(ibin)/tau_L_elem  +
     ~       OtoC/tau_L_O*tau_O0(ibin) +
     ~       NtoC/tau_L_N*tau_N0(ibin)
c$$$           photar(i)=tau_K/(tau_L_C+OtoC*tau_L_O+NtoC*tau_L_N)*
c$$$     ~       (tau_C0(ibin)*((exafs_C0(ibin)-1)*ex+1) + 
c$$$     ~         tau_O0(ibin)*OtoC + 
c$$$     ~         tau_N0(ibin)*NtoC + 
c$$$     ~         tau_F0(ibin))
        endif
        photar(i)=exp(-photar(i))
      enddo

      return
      end

      
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine load_cont_data_4 (e,tau_elem,tau_C,tau_O,tau_F,
     ~    exafs_C,tau_N,ne,
     ~    tau_L_C, tau_L_O, tau_L_N, tau_L)
      implicit none
      integer ne
      real e(*), tau_elem(*), tau_C(*), tau_O(*), tau_F(*), exafs_C(*)
      real tau_N(*)
      real tau_L_C, tau_L_O, tau_L_N
      character file*200
      integer unit, newunit
      integer status
      
      integer ibin
      real w1,w2
      real f_F,E_F, tau_F_e, tau_F_C, tau_F_O, tau_F_N
      real f_Mn,E_Mn, tau_Mn_e, tau_Mn_C, tau_Mn_O, tau_Mn_N
      real tau_L
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      integer i
      
      unit = newunit()
      
      file = '/data/alexey/cal/contam/tau_components_caldbgrid.dat'
      
      open (unit,file=file,status='old')
      ne = 0
      status = 0
      do while (status.eq.0)
         read (unit,*,iostat=status) e(ne+1), tau_elem(ne+1),
     ~        tau_C(ne+1), tau_O(ne+1), tau_F(ne+1), exafs_C(ne+1),
     ~        tau_N(ne+1)
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
      tau_F_N = w1*tau_N(ibin)+w2*tau_N(ibin+1)
      

      call locateAV (e,ne,E_Mn,ibin)
      call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
      tau_Mn_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
      tau_Mn_C = w1*tau_C(ibin)+w2*tau_C(ibin+1)
      tau_Mn_O = w1*tau_O(ibin)+w2*tau_O(ibin+1)
      tau_Mn_N = w1*tau_N(ibin)+w2*tau_N(ibin+1)
      
      tau_L = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
      tau_L_C = -log(f_F*exp(-tau_F_C)+f_Mn*exp(-tau_Mn_C))
      tau_L_O = -log(f_F*exp(-tau_F_O)+f_Mn*exp(-tau_Mn_O))
      tau_L_N = -log(f_F*exp(-tau_F_N)+f_Mn*exp(-tau_Mn_N))
      print*,tau_L, tau_L_C, tau_L_O, tau_L_N
      do i=1,ne
         tau_C(i) = tau_C(i)/exafs_C(i)
      enddo
      
      return
      end
      

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       
      subroutine get_tau_L (e,tau,n,tau_L)
      implicit none
      integer n
      real e(n), tau(n), tau_L
      real f_F, e_F, f_Mn, E_Mn
      real tau_F_e, tau_Mn_e
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      integer ibin
      real w1,w2

c     Renorm to unit tau for the L-complex
      call locateAV (e,n,E_F,ibin)
      call get_weights (E_F,e(ibin),e(ibin+1),w1,w2)
      tau_F_e = w1*tau(ibin)+w2*tau(ibin+1)
      

      call locateAV (e,n,E_Mn,ibin)
      call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
      tau_Mn_e = w1*tau(ibin)+w2*tau(ibin+1)
      
      tau_L = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))

      return
      end

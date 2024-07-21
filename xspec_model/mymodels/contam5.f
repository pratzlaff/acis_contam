      subroutine contam53func (ear,ne,param,ifl,photar)
c
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl
      
      real param1(30)
      integer i

      param1(1)=param(1)
      param1(2)=param(2)
      param1(3)=param(3)
      do i=4,19
        param1(i)=0
      enddo
      param1(5)=param(4)
      param1(12)=param(5)
      param1(19)=param(6)
      call contam5func (ear,ne,param1,ifl,photar) 
      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       
      subroutine contam5func (ear,ne,param,ifl,photar)
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
      parameter (ne0max=300000)
      real e0(0:ne0max), tau_elem(ne0max,20), www(ne0max)
      real tau_L(20)
      real exafs_C(ne0max)
      real param0(30)

      save ne0, e0, tau_elem, tau_L, exafs_C

      real abund_angers_TTT(19)
      data abund_angers_TTT /
     ~     1.000,               ! H
     ~    9.77e-2,              ! He
     ~    3.63e-4,              ! C
     ~    1.12e-4,              ! N
     ~    8.51e-4,              ! O
     ~    1.23e-4,              ! Ne
     ~    2.14e-6,              ! Na
     ~    3.80e-5,              ! Mg
     ~    2.95e-6,              ! Al
     ~    3.55e-5,              ! Si
     ~    1.62e-5,              ! S 
     ~    1.88e-7,              ! Cl
     ~    3.63e-6,              ! Ar
     ~    2.29e-6,              ! Ca
     ~    4.87e-7,              ! Cr
     ~    4.68e-5,              ! Fe
     ~    8.6e-8,               ! Co
     ~    1.78e-6,              ! Ni
     ~    3.63e-8               ! F
     ~    /


      real f_F,E_F, tau_F_e
      real f_Mn,E_Mn, tau_Mn_e
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)

      real exC, tauLtot, tauLtot1, sum, ee
      integer ibin, ielem, i, ifl0

      ifl0 = 1

      if (firstcall) then
         firstcall = .false.
         ne0 = 20000
         e0(0)=0.001
         do i=1,ne0
            e0(i)=0.001+(i-1)/2000.0
         enddo
         do ielem=1,16
            do i=1,18
               param0(i)=0
            enddo
            param0(ielem+2)=1e-5/abund_angers_TTT(ielem+2)
            call xsabsv(e0, ne0, param0, ifl0, tau_elem(1,ielem),
     ~        www)
            do i=1,ne0
               tau_elem(i,ielem)=-log(tau_elem(i,ielem))
            enddo
            print*,tau_elem(800,ielem)

c           compute norm at the L-shell

            ibin=nint(E_F*2000)
            tau_F_e = tau_elem(ibin,ielem)

            ibin=nint(E_Mn*2000)
            tau_Mn_e = tau_elem(ibin,ielem)
            
            tau_L(ielem) = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
            print*,'ini elem ',ielem,' tau_L_i=',tau_L(ielem)
         enddo

c        Load F data separately
        call load_cont_data_F (e0,ne0,tau_elem(1,17),exafs_C)

        tau_L(17)=0
        print*,'ini elem ',ielem,' tau_L_i=',tau_L(ielem)

        do i=1,17
           print*,'tau_4.5=',tau_elem(9000,i),' for ',i
        enddo
      endif
            
c Load parameters in a convenient way
      do i=1,17
         param0(i)=param(i+2)
      enddo
      tauLtot = param(1)
      exC = param(2)

      tauLtot1 = 0
      do ielem=1,17
         tauLtot1=tauLtot1+param0(ielem)*tau_L(ielem)
      enddo
      do i=1,ne
         ee = 0.5*(ear(i-1)+ear(i))
         ibin = nint(ee*2000)
         ibin = max(1,min(ibin,ne0))
         sum = 0
         do ielem=1,16
            if (ielem.eq.1) then
               sum = sum + param0(ielem)*tau_elem(ibin,ielem)*
     ~              ((exafs_C(ibin)-1)*exC+1)
            else
               sum = sum + param0(ielem)*tau_elem(ibin,ielem)
            endif
         enddo
         sum = sum*tauLtot/tauLtot1
         sum = sum + tauLtot*param0(17)*tau_elem(ibin,17) ! F
         photar(i) = exp(-sum)
      enddo

      return
      end

      
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine load_cont_data_F (e0,ne,tau_F,exafs_C)
      implicit none
      integer ne
      real e0(0:ne), tau_F(ne), exafs_C(ne)
      real e(100000), tau_F1(100000), tau_elem(100000), exafs_C1(100000)
      character file*200
      integer unit, newunit
      integer status
      
      integer ibin
      real w1,w2
      real f_F,E_F, tau_F_e
      real f_Mn,E_Mn, tau_Mn_e
      real tau_L
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      integer i, ne1
      real g1, g2, ee
      
      unit = newunit()
      
      file = '/data/alexey/cal/contam/tau_components_caldbgrid.dat'
      
      open (unit,file=file,status='old')
      ne1 = 0
      status = 0
      do while (status.eq.0)
         read (unit,*,iostat=status) e(ne1+1), tau_elem(ne1+1),
     ~        g1, g2, tau_F1(ne1+1), exafs_C1(ne1+1)
         if (status.eq.0) then
            ne1 = ne1 + 1
         endif
      enddo
      close(unit)
      
c     Renorm to unit tau for the L-complex
      call locateAV (e,ne1,E_F,ibin)
      call get_weights (E_F,e(ibin),e(ibin+1),w1,w2)
      print*,ne1,ibin,tau_elem(ibin),w1,w2
      tau_F_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)

      call locateAV (e,ne1,E_Mn,ibin)
      call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
      tau_Mn_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
      
      tau_L = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
      print*,'tau_L = ',tau_L, f_F, tau_F_e, f_Mn, tau_Mn_e

      do i=1,ne
         ee = 0.5*(e0(i-1)+e0(i))
         if (ee.lt.e(1)) then
            ibin = 1
            w1 = 1
            w2 = 0
         else if (ee.ge.e(ne1)) then
            ibin = ne1-1
            w1 = 0
            w2 = 1
         else
            call locateAV(e,ne1,ee,ibin)
            call get_weights (ee,e(ibin),e(ibin+1),w1,w2)
         endif
         tau_F(i) = (w1*tau_F1(ibin)+w2*tau_F1(ibin+1))/tau_L
         exafs_C(i) = (w1*exafs_C1(ibin)+w2*exafs_C1(ibin+1))
      enddo
      
      return
      end
      

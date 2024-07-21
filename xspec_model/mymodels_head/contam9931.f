      subroutine contam9931func (ear,ne,param,ifl,photar)
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
      real e00(0:ne0max), tau0(ne0max)
      real e0(ne0max), tau_elem(ne0max), tau_C(ne0max), tau_O(ne0max), 
     ~ tau_F(ne0max), exafs_C(ne0max), tau_Cl(ne0max)
      real  tau_C0(ne0max), tau_O0(ne0max), 
     ~ tau_F0(ne0max), exafs_C0(ne0max), tau_Cl0(ne0max)
      real tau_L_C, tau_L_O, tau_L_Cl, tau_L_elem
      save ne0, e0, tau_elem, tau_C, tau_O, tau_F, exafs_C
      save tau_C0, tau_O0, tau_F0, exafs_C0, tau_Cl0
      save tau_L_C, tau_L_O, tau_L_Cl, tau_L_elem
      save e00,tau0

      integer ibin, i
      
      real w1, w2, e
      
      real exC, tC, tE, tO, tF, eC
      real param0(30)
      integer ifl0
      
      real ee
      integer unit, newunit

      real tauL

      real f_F,E_F, tau_F_e, tau_F_f
      real f_Mn,E_Mn, tau_Mn_e, tau_Mn_f
      real tau_L
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      real efold, tauLexp

      if (firstcall) then
        call load_cont_data_5 (e0,tau_elem,tau_C,tau_O,tau_F,exafs_C,
     ~      tau_Cl,ne0,
     ~      tau_L_C,tau_L_O, tau_L_Cl, tau_L_elem)

        firstcall = .false.

        efold = 0.63
        tau_F_e=exp(-(E_F/efold)**2)
        tau_Mn_e=exp(-(E_Mn/efold)**2)
        tauLexp = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))


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

           tE = exp(-(e00(i)/efold)**2)/tauLexp
           tC = (w1*tau_C(ibin)+w2*tau_C(ibin+1))/tau_L_C
           eC = w1*exafs_C(ibin)+w2*exafs_C(ibin+1)
           tC = tC*((eC-1)*1.0+1.0)
           tO = 0.18*(w1*tau_O(ibin)+w2*tau_O(ibin+1))/tau_L_O
           tF = (w1*tau_F(ibin)+w2*tau_F(ibin+1)) /
     ~         tau_L_elem*(1+0.288256)
           tau0(i) = (tC+tO+0.75*tE+tF)/(1+0.75+0.18)
        enddo
        
      endif

      tauL = param(1)

      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        ibin = nint(e*1000)
        if (ibin.lt.1.or.ibin.gt.20000) then
           photar(i)=0
        else
          photar(i)=tauL*tau0(ibin)
        endif
        photar(i)=exp(-photar(i))
      enddo

      return
      end

      

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine load_cont_data_5 (e,tau_elem,tau_C,tau_O,tau_F,
     ~    exafs_C,tau_Cl,ne,
     ~    tau_L_C, tau_L_O, tau_L_Cl, tau_L_elem)
      implicit none
      integer ne
      real e(*), tau_elem(*), tau_C(*), tau_O(*), tau_F(*), exafs_C(*)
      real tau_Cl(*)
      real tau_L_C, tau_L_O, tau_L_Cl
      character file*200
      integer unit, newunit
      integer status
      
      integer ibin
      real w1,w2
      real f_F,E_F, tau_F_e, tau_F_C, tau_F_O, tau_F_Cl
      real f_Mn,E_Mn, tau_Mn_e, tau_Mn_C, tau_Mn_O, tau_Mn_Cl
      real tau_L_elem
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)
      integer i
      real gov
      
      unit = newunit()
      
      file = '/data/alexey/cal/contam/tau_components_caldbgrid.dat'
      
      open (unit,file=file,status='old')
      ne = 0
      status = 0
      do while (status.eq.0)
         read (unit,*,iostat=status) e(ne+1), tau_elem(ne+1),
     ~        tau_C(ne+1), tau_O(ne+1), tau_F(ne+1), exafs_C(ne+1),
     ~        gov, tau_Cl(ne+1)
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
      tau_F_Cl = w1*tau_Cl(ibin)+w2*tau_Cl(ibin+1)
      

      call locateAV (e,ne,E_Mn,ibin)
      call get_weights (E_Mn,e(ibin),e(ibin+1),w1,w2)
      tau_Mn_e = w1*tau_elem(ibin)+w2*tau_elem(ibin+1)
      tau_Mn_C = w1*tau_C(ibin)+w2*tau_C(ibin+1)
      tau_Mn_O = w1*tau_O(ibin)+w2*tau_O(ibin+1)
      tau_Mn_Cl = w1*tau_Cl(ibin)+w2*tau_Cl(ibin+1)
      
      tau_L_elem = -log(f_F*exp(-tau_F_e)+f_Mn*exp(-tau_Mn_e))
      tau_L_C = -log(f_F*exp(-tau_F_C)+f_Mn*exp(-tau_Mn_C))
      tau_L_O = -log(f_F*exp(-tau_F_O)+f_Mn*exp(-tau_Mn_O))
      tau_L_Cl = -log(f_F*exp(-tau_F_Cl)+f_Mn*exp(-tau_Mn_Cl))
      print*,tau_L_elem, tau_L_C, tau_L_O, tau_L_Cl
      do i=1,ne
         tau_C(i) = tau_C(i)/exafs_C(i)
      enddo
      
      return
      end

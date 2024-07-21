      subroutine exptau (ear,ne,param,ifl,photar)
      implicit none
      integer ne
      real ear(0:ne), photar(ne)
      real param(*)
      integer ifl
      
      real f_F, e_F, f_Mn, E_Mn
      real tau_F, tau_Mn, tau_L1
      parameter (f_F  = 0.6, E_F=0.679)
      parameter (f_Mn = 0.4, E_Mn=0.638)

      real tau_L, efold, e
      integer i


      tau_L = param(1)
      efold = param(2)

      
      tau_F=exp(-(E_F/efold)**2)
      tau_Mn=exp(-(E_Mn/efold)**2)
      tau_L1 = -log(f_F*exp(-tau_F)+f_Mn*exp(-tau_Mn))


      
      do i=1,ne
        e = 0.5*(ear(i-1)+ear(i))
        photar(i) = tau_L/tau_L1*exp(-(e/efold)**2)
        photar(i)=exp(-photar(i))
      enddo

      return
      end

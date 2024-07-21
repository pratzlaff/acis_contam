      subroutine get_weights (x,x1,x2,w1,w2)
      implicit none
      real x,x1,x2,w1,w2

      w1=(x2-x)/(x2-x1)
      w2=(x-x1)/(x2-x1)
      
      return
      end

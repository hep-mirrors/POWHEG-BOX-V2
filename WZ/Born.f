      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born,colcf
      integer j,k,mu,nu
c                    id1      iad1     id2      idw
      call setzcoupl(bflav(5),bflav(6),bflav(7),bflav(3))


      call compborn(p,bflav,born,bmunu,bornjk)

      end


      subroutine compborn(pin,bflav,born,bmunu,bornjk)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 pin(0:3,nlegborn)
      integer bflav(nlegborn)
      real * 8 born,bmunu(0:3,0:3,nlegs),bornjk(nlegs,nlegs)
      real *8 amp2,colcf
      integer ferm_type(nlegborn)
      integer i,j,k,mu,nu
      real * 8 ferm_charge(nlegborn)
      real * 8 p(12,1:4)  ! 12 = mxpart in MCFM
      double precision msq(-5:5,-5:5)
      real * 8 suppfact2e
c     vector boson id and decay
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode

      do i=1,2
         p(i,4) = pin(0,i)
         p(i,1:3) = pin(1:3,i)
      enddo

      do i=3,6
         p(i,4) = pin(0,i+2)
         p(i,1:3) = pin(1:3,i+2)
      enddo

      p(1,:) = -p(1,:) 
      p(2,:) = -p(2,:) 

      call qqb_wz(p,msq)

      born = msq(bflav(1),bflav(2))
      born = born * normbr


C     -- no gluons, so no spin correlated Born  
      do i=0,3
         do j=0,3
            bmunu(i,j,:)=0d0
         enddo
      enddo
      
      bornjk=0d0
C     -- here the bornjk is very simple B_i j = B C_i
      do j=1,nlegs
         if(abs(bflav(j)).le.6) then
            do k=j+1,nlegs
               if(abs(bflav(k)).le.6) then
                  colcf=cf
                  bornjk(j,k)=born*colcf
                  bornjk(k,j)=bornjk(j,k)
               endif
            enddo
         endif
      enddo
      


      end


      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      logical isquark
c     neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      icolup(1,5)=0
      icolup(2,5)=0
      icolup(1,6)=0
      icolup(2,6)=0
c     colored particles
      if((idup(1).gt.0).and.(idup(2).lt.0)) then
         icolup(1,1)=501
         icolup(2,1)=0
         icolup(1,2)=0
         icolup(2,2)=501
      elseif((idup(1).lt.0).and.(idup(2).gt.0)) then
         icolup(1,1)=0
         icolup(2,1)=501
         icolup(1,2)=501
         icolup(2,2)=0
      else
         write(*,*) ' invalid flavour'
         stop
      endif
      if(isquark(5)) then
         icolup(1,5)=502
         icolup(2,5)=0
         icolup(1,6)=0
         icolup(2,6)=502
      else
         icolup(:,5)=0
         icolup(:,6)=0
      endif
      if(isquark(7)) then
         icolup(1,7)=503
         icolup(2,7)=0
         icolup(1,8)=0
         icolup(2,8)=503
      else
         icolup(:,7)=0
         icolup(:,8)=0
      endif      
      end

      subroutine finalize_lh
      implicit none
      include 'LesHouches.h'
      integer i1,i2
      real * 8 v(5),v1
      integer itmp2(2)
      real * 8 random
      external random
      call lhefinitemasses
      end

      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'zcouple.f'  !TM now set the z-coupling parameters here
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'cvecbos.h'
      logical verbose
      parameter(verbose=.true.)
      physpar_ml(1)=0.511d-3
      physpar_ml(2)=0.1057d0
      physpar_ml(3)=1.777d0
      physpar_mq(1)=0.33d0     ! up
      physpar_mq(2)=0.33d0     ! down
      physpar_mq(3)=0.50d0     ! strange
      physpar_mq(4)=1.50d0     ! charm
      physpar_mq(5)=4.80d0     ! bottom

      call smcouplings

c     number of light flavors
      st_nlight = 5

      !TM added QCD couplings

      write(*,*)'alpha',st_alpha
      gsq = st_alpha*4d0*pi
      as  = st_alpha
      ason2pi = st_alpha/2d0/pi
      ason4pi = st_alpha/4d0/pi


      
      !TM added z couplings
      Q(+1)=-0.333333333333333d0
      Q(+2)=+0.666666666666667d0
      Q(+3)=-0.333333333333333d0
      Q(+4)=+0.666666666666667d0
      Q(+5)=-0.333333333333333d0
      tau=(/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/)
      esq = ph_unit_e**2
      zmass = ph_Zmass
      zwidth = ph_Zwidth

      call couplz(ph_sthw2)

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'e**2  = ',ph_unit_e**2
      write(*,*) '*************************************'
      endif

      end





      subroutine setzcoupl(id1,id2)
      implicit none
      integer id1,id2
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'zcouple.f'  !TM now set the z-coupling parameters here
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'cvecbos.h'
      include 'vvsettings.f'
      logical withinterference,ini
      real * 8 mllmin,m1,m2
      data ini/.true./
      save withinterference,ini,mllmin
      real * 8 powheginput
      logical isquark,islepton,isnu
      external powheginput,isquark,islepton,isnu
      if(ini) then
         if(powheginput("#withinterference").eq.0) then
            withinterference = .false.
         else
            withinterference = .true.
         endif
         ini=.false.
         mllmin=powheginput("#mllmin")
         if(mllmin.le.0) mllmin=0.1d0
      endif

      normbr = 1

      if(isquark(id1)) then
         q1=q(id1)
         l1=l(id1)
         r1=r(id1)
         normbr=normbr*(1d0+ph_deltas)*3
      elseif(islepton(id1)) then
         q1=-1
         l1=le
         r1=re
      elseif(isnu(id1)) then
         q1=0
         l1=ln
         r1=rn
      else
         write(*,*) 'setzcoupl: invalid Z decay product' ,id1
         call pwhg_exit(-1)
      endif

      if(isquark(id2)) then
         q2=q(id2)
         l2=l(id2)
         r2=r(id2)
         normbr=normbr*(1d0+ph_deltas)*3
      elseif(islepton(id2)) then
         q2=-1
         l2=le
         r2=re
      elseif(isnu(id2)) then
         q2=0
         l2=ln
         r2=rn
      else
         write(*,*) 'setzcoupl: invalid Z decay product' ,id2
         call pwhg_exit(-1)
      endif

      if(.not.isnu(id1)) then
         m1=sqrt(kn_cmpborn(0,3)**2-kn_cmpborn(1,3)**2
     1        -kn_cmpborn(2,3)**2-kn_cmpborn(3,3)**2)
         if(m1.lt.mllmin) then
            kn_jacborn = 0
         endif
      endif

      if(.not.isnu(id2)) then
         m2=sqrt(kn_cmpborn(0,4)**2-kn_cmpborn(1,4)**2
     1        -kn_cmpborn(2,4)**2-kn_cmpborn(3,4)**2)
         if(m2.lt.mllmin) then
            kn_jacborn = 0
         endif
      endif

      if(id1.eq.id2) then
         vsymfact=0.5d0
         if(withinterference) then
            interference=.true.
         else
            interference=.false.
         endif
         if(.not.isnu(id1)) then
c impose mllmin cut also upon the crossed pairs
            m1=sqrt(
     1           (kn_cmpborn(0,5)+kn_cmpborn(0,8))**2-
     2           (kn_cmpborn(1,5)+kn_cmpborn(1,8))**2-
     3           (kn_cmpborn(2,5)+kn_cmpborn(2,8))**2-
     4           (kn_cmpborn(3,5)+kn_cmpborn(3,8))**2  )
            m2=sqrt(
     1           (kn_cmpborn(0,6)+kn_cmpborn(0,7))**2-
     2           (kn_cmpborn(1,6)+kn_cmpborn(1,7))**2-
     3           (kn_cmpborn(2,6)+kn_cmpborn(2,7))**2-
     4           (kn_cmpborn(3,6)+kn_cmpborn(3,7))**2  )
            if(min(m1,m2).lt.mllmin) then
               kn_jacborn = 0
            endif
         endif
      else
         vsymfact=1
         interference=.false.
      endif

      end

      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      integer i1,i2,i9,id1,id2
      logical debug
      parameter (debug=.true.)
      integer j,k
      real * 8 powheginput
      external powheginput
c     vector boson id and decay
c     lepton masses
      include 'cvecbos.h'
      logical condition,alloweddec
      external alloweddec
      include 'vvsettings.f'

      z1_to_ch=.false.
      z2_to_ch=.false.

      if(powheginput('#cutallpairs').eq.0) then
         cutallpairs = .false.
      else
         cutallpairs = .true.
      endif

      if (powheginput("#zerowidth").eq.1) then 
         zerowidth = .true. 
      else
         zerowidth = .false. 
      endif

      if (powheginput("#dronly").eq.1) then 
         dronly = .true. 
         write(*,*) 'Double resonant diagrams only'
c  cant have srdiags if zerowidth is true
      elseif (zerowidth) then
         write(*,*) 'Single resonant diagrams require off-shell Z'
         write(*,*) 'Double resonant diagrams only'
         dronly = .true.
      else
         dronly=.false.
         write(*,*) 'Including single resonant diagrams'
      endif

      if(powheginput("#withdamp").ne.0) then
         flg_withdamp = .true.
         flg_bornzerodamp = .true.
      else
         flg_withdamp = .false.
         flg_bornzerodamp = .false.
      endif

      interference = (powheginput('#withinterference').ne.0)
      interference=(interference).and.(.not. zerowidth)

      idvecbos1 = 23
      idvecbos2 = 23

      lprup(1)=10000


c     index of the first coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=9

*********************************************************************
***********            BORN SUBPROCESSES              ***************
*********************************************************************
      flst_nborn=0
      do i1=-5,5
         do i2=-5,5
            if(i1.ne.0.and.i1+i2.eq.0) then
c     q qbar
               do id1=1,16
                  do id2=1,16
                     if(.not.alloweddec(id1,id2)) cycle
                     flst_nborn=flst_nborn+1
                     if(flst_nborn.gt.maxprocborn) goto 999
                     flst_born(1,flst_nborn)=i1
                     flst_born(2,flst_nborn)=i2
                     flst_born(3,flst_nborn)=23
                     flst_born(4,flst_nborn)=23
                     flst_born(5,flst_nborn)=id1
                     flst_born(6,flst_nborn)=-id1
                     flst_bornres(5,flst_nborn)=3
                     flst_bornres(6,flst_nborn)=3
                     flst_born(7,flst_nborn)=id2
                     flst_born(8,flst_nborn)=-id2
                     flst_bornres(7,flst_nborn)=4
                     flst_bornres(8,flst_nborn)=4
                  enddo
               enddo
            endif
         enddo
      enddo
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

*********************************************************************
***********            REAL SUBPROCESSES              ***************
*********************************************************************
      flst_nreal=0
      do i1=-5,5
         do i2=-5,5
            do i9=-5,5
               condition=.false.
               if(.not.(i1.eq.0.and.i2.eq.0)) then
c     exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.0).and.(i9.eq.0)) then
c     q qbar -> g
                     condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i9)) then
c     g q
                     condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i9)) then
c     q g
                     condition=.true.
                  endif
               endif
               if(condition) then
                  do id1=1,16
                     do id2=1,16
                        if(.not.alloweddec(id1,id2)) cycle
                        flst_nreal=flst_nreal+1
                        if(flst_nreal.gt.maxprocreal) goto 998
                        flst_real(1,flst_nreal)=i1
                        flst_real(2,flst_nreal)=i2
                        flst_real(3,flst_nreal)=23
                        flst_real(4,flst_nreal)=23
                        flst_real(5,flst_nreal)=id1
                        flst_real(6,flst_nreal)=-id1
                        flst_realres(5,flst_nreal)=3
                        flst_realres(6,flst_nreal)=3
                        flst_real(7,flst_nreal)=id2
                        flst_real(8,flst_nreal)=-id2
                        flst_realres(7,flst_nreal)=4
                        flst_realres(8,flst_nreal)=4
                        flst_real(9,flst_nreal)=i9
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif
      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      call exit(-1)
 999  write(*,*) 'init_processes: increase maxprocborn'
      call exit(-1)
      end

      function alloweddec(id1,id2)
      implicit none
      logical alloweddec
      integer id1,id2
      include 'cvecbos.h'
      logical isquark,isutype,isdtype,islepton,isnu
      real * 8 powheginput
      external powheginput,isquark,isutype,isdtype,islepton,isnu
      alloweddec = .false.
      if(.not.(isquark(id1).or.islepton(id1).or.isnu(id1))) return
      if(.not.(isquark(id2).or.islepton(id2).or.isnu(id2))) return
c if the two id are different, then order them
c as follows
      if(id1.ne.id2) then
c in semileptonic decays, id1 is always the lepton
         if(islepton(id2).and.isquark(id1)) return
c in lepton - neutrino, first is the lepton
         if(islepton(id2).and.isnu(id1)) return
c if up and down types, ups are first
         if(isutype(id2).and.isdtype(id1)) return
c if they are both of the same type, order by increasing id
         if((islepton(id1).and.islepton(id2).or.
     1        isutype(id1).and.isutype(id2).or.
     2        isdtype(id1).and.isdtype(id2).or.
     3        isnu(id1).and.isnu(id2)).and.id1.gt.id2) return
      endif

c Now decay types are ordered: leptons before hadrons, ups before downs,
c leptons before nus, and everything being equal, by increasing id.
c Check which decay configurations we want to keep
      if(powheginput("#semileptonic").eq.1.and.
     1     .not.(islepton(id1).and.isquark(id2))) return
      if(powheginput("#e-mu").eq.1.and..not.(id1.eq.11.and.id2.eq.13))
     1     return
      if(powheginput("#e-tau").eq.1.and..not.(id1.eq.11.and.id2.eq.15))
     1     return
      if(powheginput("#mu-tau").eq.1.and..not.(id1.eq.13.and.id2.eq.15))
     1     return
      if(powheginput("#leptonic").eq.1.and.
     1     .not.(islepton(id1).and.islepton(id2))) return
      if(powheginput("#leptons-nu").eq.1.and.
     1     .not.(islepton(id1).and.isnu(id2))) return
      if(powheginput("#hadrons-nu").eq.1.and.
     1     .not.(isquark(id1).and.isnu(id2))) return
      if(powheginput("#only-e").eq.1.and. (id1.ne.11.or.id2.ne.11))
     1 return
      if(powheginput("#only-mu").eq.1.and.(id1.ne.13.or.id2.ne.13))
     2     return
      if(powheginput("#only-tau").eq.1.and.(id1.ne.15.or.id2.ne.15))
     2     return

c Add here other options at will, if you like;


      alloweddec = .true.

c Here we want to know if there can be charged particle decays in Z1
      if(.not.isnu(id1)) then
         Z1_to_ch=.true.
      endif
c Here we want to know if there can be charged particle decays in Z2
      if(.not.isnu(id2)) then
         Z2_to_ch=.true.
      endif

      end


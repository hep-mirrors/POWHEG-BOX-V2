      subroutine pwhgnewweight(iunin,iunrwgt)
      implicit none
      integer iunin,iunrwgt
      character * 500 stringin
      integer iret
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
      integer maxev
      integer k
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      real * 8 newweight
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 3 whichpdfpk
      character * 200 string
      call lhefreadevnew(iunin,iunrwgt,iret,stringin)
      if(iret.lt.0) then
         write(*,*) ' End of event file! Aborting ...'
         call exit(-1)
      endif
      call setrandom(gen_seed,gen_n1,gen_n2)
      if(rad_type.eq.1) then
         call gen_btilderw
         newweight=rad_btilde_arr(rad_ubornidx)*
     1        rad_btilde_sign(rad_ubornidx)
      elseif(rad_type.eq.2) then
         call gen_sigremnantrw
         newweight=rad_damp_rem_arr(rad_realalr)
      elseif(rad_type.eq.3) then
         call gen_sigremnantrw
         newweight=rad_reg_arr(rad_realreg)
      else
         write(*,*) 'Error in pwhgnewweight, invalid rad_type: ',
     $        rad_type
         call exit(-1)
      endif

      if(.not.pwhg_isfinite(newweight)) newweight=0d0
      write(string,*) '#new weight,renfact,facfact,pdf1,pdf2',
     1        xwgtup*newweight/rad_currentweight,st_renfact,
     2        st_facfact,pdf_ndns1,pdf_ndns2,' ',whichpdfpk()
      write(iunrwgt,'(a)') trim(adjustl(string))

      if(adjustl(stringin).eq.'<rwgt>') then
         write(iunrwgt,'(a)') trim(stringin)
         do k=1,1000000
            read(unit=iunin,fmt='(a)',end=999,err=999) stringin
            if(stringin.ne.'</rwgt>') then
               write(iunrwgt,'(a)') trim(stringin)
            else
               call lhrwgt_writeweight
     1              (iunrwgt,xwgtup*newweight/rad_currentweight)
               write(iunrwgt,'(a)') trim(stringin)
               return
            endif
         enddo
      else
         write(iunrwgt,'(a)') trim(stringin)
      endif
      return
 999  continue
      write(*,*) 'Did not find </rwgt> after <rwgt>'
      write(*,*) 'in Les Houches file. Exiting ...'
      call exit(-1)
      end


      subroutine gen_btilderw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer mcalls,icalls
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)      
      real * 8 btilde
      external btilde
      mcalls=0
      icalls=0
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,2,
     #    mcalls,icalls,xx)
      end

      subroutine gen_sigremnantrw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      real * 8 sigremnant
      external sigremnant
      mcalls=0
      icalls=0
      call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1 xmmmrm,ifoldrm,2,mcalls,icalls,xx)
      end



      subroutine openoutputrw(iunrwgt)
      implicit none
      integer iunrwgt
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      call newunit(iunrwgt)
      if(rnd_cwhichseed.ne.'none') then
         open(unit=iunrwgt,file=pwgprefix(1:lprefix)//'events-rwgt-'
     1        //rnd_cwhichseed//'.lhe'
     2     ,status='unknown')
      else
         open(unit=iunrwgt,file=pwgprefix(1:lprefix)//'events-rwgt.lhe'
     1     ,status='unknown')
      endif
      end

c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadevnew(nlf,nuo,iret,string)
      implicit none
      integer nlf,nuo,iret
      character * 500 string
      include 'LesHouches.h'
      integer i,j
      iret=0
 1    continue
      string=' '
      read(nlf,fmt='(a)',err=777,end=666) string
      write(nuo,'(a)') trim(string)
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         read(nlf,fmt='(a)',err=777,end=666) string
         write(nuo,'(a)') trim(string)
         read(string,fmt=*,end=998,err=1)
     1        nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,'(a)',err=777,end=666) string
            write(nuo,'(a)') trim(string)
            read(string,fmt=*,end=998,err=1)
     1           idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadextrarw(nlf,nuo,iret,string)
c Go on reading, find the end of the weights section
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      write(*,*) "Error in reading"
      write(*,*) string
      call exit(-1)
 666  continue
      iret=-1
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      iret=-1
      nup=0      
 999  end


      subroutine lhefwritetrailernw(iunin,iunout)
      implicit none
      integer iunin,iunout
      character * 500 string
 1    continue
      string=' '
      read(iunin,fmt='(a)',err=998,end=998) string
      write(iunout,'(a)') trim(string)
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      goto 1
 998  continue
      end



      subroutine lhefreadextrarw(nlf,nou,iret,string)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      character * 500 string
      integer nlf,nou,iret
      logical readrw
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      readrw = .false.
 1    continue
      read(unit=nlf,fmt='(a)',end=998) string
      if(string.eq.'<rwgt>'.or.
     1     string.eq.'</event>') then
c Don't write the end event record; first we must output the new weight
         return
      endif
      if(string.eq.'<event>') then
         if(.not.readrw) then
            write(*,*) 
     $ 'Error in lhefreadextra, while reading rwg informations'
            write(*,*)'Abort run'
            call exit(-1)
         endif
         backspace nlf
         return
      endif
      if(string.eq.'# Start extra-info-previous-event') then
         write(nou,'(a)') trim(string)
         read(nlf,'(a)') string
         write(nou,'(a)') trim(string)
         read(string(3:),*) rad_kinreg
         read(nlf,'(a)') string
         read(string(3:),*) rad_type
      endif
      write(nou,'(a)') trim(string)
      if(flg_newweight) then
c read a string; if it starts with #rwgt, read first rad_type from the
c string, then all other information, depending upon rad_type.
c set readrw to true
         string=adjustl(string)
         if(string(1:5).eq.'#rwgt') then
            string(1:5)=' '
c     do things
c            print*, 'FOUND'
            read(string,*) rad_type
            if(rad_type.eq.1) then
c     btilde
               read(string,*)rad_type,
     $              rad_ubornidx,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.2) then
c     remnant
               read(string,*)rad_type,
     $              rad_realalr,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.3) then
c     regular
               read(string,*)rad_type,
     $              rad_realreg,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            else
               write(*,*) 'Invalid rad_type in lhefwriteevrw: ',rad_type
               call exit(-1)
            endif
c     if all went ok, set readrw to true
            readrw=.true.
         endif
      endif
      goto 1
 998  continue
      iret=-1
      end


c The following are routines for the implementation of the
c Les Houches standard for reweighting


      subroutine lhrwgt_writeweight(iun,weight)
      implicit none
      integer iun
      include 'pwhg_lhrwgt.h'
      double precision weight
      write(iun,'(a,e16.9,a)')
     1     "<wgt id='"//trim(adjustl(lhrwgt_id))//"'> ",
     2     weight,' </wgt>'
      end

      subroutine lhrwgt_writeheader(iun)
      implicit none
      integer iun
      include 'pwhg_lhrwgt.h'
      character * 200 tmpstring
      integer j,iret
      logical written,in_group
c This flag is set to true when the current weight has
c been written out
      written = .false.
      write(iun,'(a)') '<initrwgt>'
      if(lhrwgt_group_name.eq.' ') then
c Write all lines firs
         do j=1,lhrwgt_nheader
            write(iun,'(a)') trim(lhrwgt_header(j))
         enddo
      else
c this flag is set to true when we are writing out weights in
c the current group
         in_group = .false.
         do j=1,lhrwgt_nheader
            if(.not.written.and.in_group
     1           .and.adjustl(lhrwgt_header(j)).eq.'</weightgroup>')
     2           then
               write(iun,'(a)') "<weight id='"//
     1              trim(adjustl(lhrwgt_id))//
     2              "'>"//' '//trim(adjustl(lhrwgt_descr))
     3              //' </weight>'
               in_group = .false.
               written = .true.
            endif
            write(iun,'(a)') trim(lhrwgt_header(j))
            tmpstring = adjustl(lhrwgt_header(j))
            if(tmpstring(1:18).eq.'<weightgroup name=') then
               tmpstring = tmpstring(19:)
               call getquotedstring(tmpstring,tmpstring,iret)
               if(tmpstring.eq.lhrwgt_group_name) in_group = .true.
            endif
         enddo
      endif
      if(.not.written) then
         if(lhrwgt_group_name.ne.' ') then
            if(lhrwgt_group_combine.ne.' ') then
               write(iun,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"' combine='"
     1              //trim(lhrwgt_group_combine)//"'>"
            else
               write(iun,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"'>"
            endif
         endif
         write(iun,'(a)') "<weight id='"//
     1        trim(adjustl(lhrwgt_id))//
     1        "'>"//' '//trim(adjustl(lhrwgt_descr))//' </weight>'
         if(lhrwgt_group_name.ne.' ') then
            write(iun,'(a)') '</weightgroup>'
         endif
      endif
      write(iun,'(a)') '</initrwgt>'
      end

      subroutine lhrwgt_checkheader
      implicit none
      include 'pwhg_lhrwgt.h'
      character * (lhrwgt_max_header_columns) string
      integer j,iret
      do j=1,lhrwgt_nheader
         string = adjustl(lhrwgt_header(j))
         if(string(1:11).eq.'<weight id=') then
            call getquotedstring(string(12:),string,iret)
            if(string.eq.lhrwgt_id) then
               write(*,*) ' lhrwgt_checkheader: ERROR'
               write(*,*) ' lhrwgt_id ',trim(string),' already in use'
               call exit(-1)
            endif
         endif
      enddo
      end

      subroutine lhrwgt_readheader(iun)
      implicit none
      include 'pwhg_lhrwgt.h'
      integer iun
      character * (lhrwgt_max_header_columns+200) string,string0
      integer j,k
      do j=1,1000000
         read(iun,'(a)') string0
         string=adjustl(string0)
         if(string.eq.'<initrwgt>') then
            lhrwgt_nheader = 0
            do k=1,1000000
               read(iun,'(a)') string0
               string=adjustl(string0)
               if(string.ne.'</initrwgt>') then
                  lhrwgt_nheader = lhrwgt_nheader + 1
                  if(lhrwgt_nheader.gt.lhrwgt_maxnheader) then
                     write(*,*) '******** ERROR ******'
                     write(*,*) ' lhrwgt_readheader: '
                     write(*,*) 
     1                    ' found too many strings in header'
                     write(*,*) ' increase lhrwgt_maxnheader'
     1                    //' in include/pwg_lhrwgt.h'
                     call exit(-1)
                  endif                    
c check that the string fits
                  if(len(trim(string0)).gt.
     1                 lhrwgt_max_header_columns) then
                     write(*,*) '******** ERROR ******'
                     write(*,*) ' lhrwgt_readheader: '
                     write(*,*) 
     1                    ' found too long a string in the header'
                     write(*,*) 
     1                    ' increase lhrwgt_max_header_columns'
     1                    //' in include/pwg_lhrwgt.h'
                     call exit(-1)
                  endif
                  lhrwgt_header(lhrwgt_nheader) = trim(string0)
               else
                  return
               endif
            enddo
            write(*,*) '******** ERROR ******'
            write(*,*) ' lhrwgt_readheader: '
            write(*,*) ' didn;t find the end of the header'
            call exit(-1)
         endif
      enddo
      write(*,*) '******** ERROR ******'
      write(*,*) ' lhrwgt_readheader: '
      write(*,*) ' didn;t find the header'
      call exit(-1)
      end

      subroutine lhrwgt_copyheader(iunin,iunout)
      implicit none
      include 'pwhg_lhrwgt.h'
      integer iunin,iunout
      character * (400) string,string0
      integer j,k
      logical written,in_group
      do j=1,1000000
         read(unit=iunin,fmt='(a)',end=999,err=999) string0
         string=adjustl(string0)
         if(string.eq.'<initrwgt>') then
            backspace(iunin)
            call lhrwgt_readheader(iunin)
            call lhrwgt_checkheader
            call lhrwgt_writeheader(iunout)
            return
         else
            write(iunout,'(a)') trim(string0)
         endif
      enddo
 999  continue
      write(*,*) ' ERROR: did not find <initrwgt> in lhe file'
      write(*,*) ' Most likely, the initial run did not specify'
      write(*,*) ' a rwgt_id line'
      call exit(-1)
      end


      subroutine getreweightinput
      implicit none
      include 'pwhg_lhrwgt.h'
      call powheginputstring("#lhrwgt_group_name",lhrwgt_group_name)
      call powheginputstring("#lhrwgt_group_combine",
     1     lhrwgt_group_combine)
      call powheginputstring("#lhrwgt_id",lhrwgt_id)
      call powheginputstring("#lhrwgt_descr",lhrwgt_descr)
      end

      subroutine printrwghthdr(iun)
      implicit none
      include 'pwhg_lhrwgt.h'
      INTEGER IUN
      if(lhrwgt_id.ne.' ') then
         write(iun,'(a)') '<initrwgt>'
         if(lhrwgt_group_name.ne.' ') then
            if(lhrwgt_group_combine.ne.' ') then
               write(iun,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"' combine='"
     1              //trim(lhrwgt_group_combine)//"'>"
            else
               write(iun,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"'>"
            endif
         endif
         write(iun,'(a)') "<weight id='"//trim(lhrwgt_id)
     1       //"'> "//trim(lhrwgt_descr)//" </weight>"
         if(lhrwgt_group_name.ne.' ') then
             write(iun,'(a)') "</weightgroup>"
          endif
         write(iun,'(a)') "</initrwgt>"
      endif
      end

      subroutine printrwgtev(nlf,weight)
      implicit none
      integer nlf
      double precision weight
      include 'pwhg_lhrwgt.h'
      write(nlf,'(a,e16.9,a)')"<wgt id='"//trim(lhrwgt_id)//"'> ",
     1     weight,' </wgt>' 
      end

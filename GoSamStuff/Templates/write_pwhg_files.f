      program fix_everything
      implicit none
      call write_virtual
      call write_init_couplings
      end


      subroutine write_init_couplings
      implicit none
      integer maxnumlines
      integer maxnum_physpar
      parameter (maxnum_physpar=37)  ! To be increased if further parameters are added
      parameter (maxnumlines=100000)
      integer incoup, inproc, outcoup, outproc, i
      character * 150 stringa
      logical found
      character * 15 phys_par_names(2*maxnum_physpar)
      data phys_par_names/
     $     'Nf','st_nlight',
     $     'alpha','ph_alphaem',
     $     'GF','ph_gfermi',
     $     'mZ','ph_Zmass',
     $     'mW','ph_Wmass',
     $     'mH','ph_Hmass',

     $     'mC','ph_cmass',
     $     'mB','ph_bmass',
     $     'mBMS','ph_bmass',   ! b-quark Yukawa coupling. Set to the b mass here
     $     'mT','ph_tmass',

     $     'me','ph_emass',
     $     'mmu','ph_mumass',
     $     'mtau','ph_taumass',

     $     'wZ','ph_Zwidth',
     $     'wW','ph_Wwidth',
     $     'wH','ph_hwidth',
     $     'wB','ph_bwidth',
     $     'wT','ph_twidth',
     $     'wtau','ph_tauwidth',

     $     'VUD','ph_CKM(1,1)',
     $     'CVDU','ph_CKM(1,1)',
     $     'VUS','ph_CKM(1,2)',
     $     'CVSU','ph_CKM(1,2)',
     $     'VUB','ph_CKM(1,3)',
     $     'CVBU','ph_CKM(1,3)',
     $     'VCD','ph_CKM(2,1)',
     $     'CVDC','ph_CKM(2,1)',
     $     'VCS','ph_CKM(2,2)',
     $     'CVSC','ph_CKM(2,2)',
     $     'VCB','ph_CKM(2,3)',
     $     'CVBC','ph_CKM(2,3)',
     $     'VTD','ph_CKM(3,1)',
     $     'CVDT','ph_CKM(3,1)',
     $     'VTS','ph_CKM(3,2)',
     $     'CVST','ph_CKM(3,2)',
     $     'VTB','ph_CKM(3,3)',
     $     'CVBT','ph_CKM(3,3)'/
   
      found=.false.
      incoup=10
      outcoup=11
      inproc=12
      outproc=13
      open(unit=incoup,status='old',file='../init_couplings.f')
      open(unit=outcoup,file='init_couplings_new.f')
      do i=1,maxnumlines
         read(incoup,'(a)',err=666,end=777) stringa
         if (stringa(1:31).eq.'      subroutine init_couplings') 
     $        found=.true.   
         if (found.and.stringa(1:10).eq.'      end') then            
            write(outcoup,1) '      call golem_initialize'
            found=.false.
         endif
         write(outcoup,1) stringa(1:len_trim(stringa))
      enddo
 666  continue
 777  continue
c     writing golem_initialize subroutine
      write(outcoup,4) ""
      write(outcoup,6) "initializes all the couplings in the "//
     $     " virtual, code"
      write(outcoup,6) "generated by GoSam and sets them equal to the"
      write(outcoup,6) "ones defined in the POWHEG BOX.              "
      write(outcoup,4) "subroutine golem_initialize"
      write(outcoup,4) "implicit none"
      write(outcoup,4) "include 'PhysPars.h'"
      write(outcoup,4) "include 'pwhg_st.h'"
      write(outcoup,4) "include 'pwhg_rnd.h'"
      write(outcoup,4) "integer ierr"
      write(outcoup,4) "integer ioerr"
      write(outcoup,4) "character * 20 param"
      write(outcoup,4) "character * 20 value"
      write(outcoup,4) "character * 50 line"
      write(outcoup,4) "character * 29 path"
      write(outcoup,4) "real * 8 powheginput"
      write(outcoup,4) "external powheginput"
      write(outcoup,4) "integer parallelstage,rndiwhichseed"
      write(outcoup,4) "common/cpwhg_info/parallelstage,rndiwhichseed"
      write(outcoup,4) ""
      write(outcoup,4) "rndiwhichseed=rnd_iwhichseed"
      write(outcoup,4) "parallelstage=powheginput('#parallelstage')"
      write(outcoup,4) ""
      write(outcoup,6) "Parameter definition"
      write(outcoup,4) ""
      do i=1,2*maxnum_physpar-1,2
         write(outcoup,4) "param = '"//
     $        phys_par_names(i)(1:len_trim(phys_par_names(i)))//"='"
         if (i .eq. 1) then
            write(outcoup,4) "write(value,'(I1)') "//
     $           phys_par_names(i+1)(1:len_trim(phys_par_names(i+1)))
         else
            write(outcoup,4) "write(value,'(F20.10)') "//
     $           phys_par_names(i+1)(1:len_trim(phys_par_names(i+1)))
         endif
         write(outcoup,4) "line = trim(param)//trim(adjustl(value))"
         write(outcoup,4) "call OLP_Option(line,ierr)"
         write(outcoup,4) "call check_gosam_err(param,ierr)"
         write(outcoup,4) ""
      enddo
      write(outcoup,6) "Initialize virtual code"
      write(outcoup,4) ""
      write(outcoup,4) "path = '../GoSam_POWHEG/orderfile.olc'"
      write(outcoup,4) ""
      write(outcoup,41)"call OLP_Start(path,ioerr,parallelstage,",
     $     "rndiwhichseed)"
      write(outcoup,4) "call check_gosam_err('olp_start routine',ierr)"
      write(outcoup,4) "end"

      write(outcoup,4) ""
      write(outcoup,4) ""
      write(outcoup,4) "subroutine check_gosam_err(param,ierr)"
      write(outcoup,4) "implicit none"
      write(outcoup,4) "character *(*) param"
      write(outcoup,4) "integer ierr"
      write(outcoup,4) "if (ierr.lt.0) then"
      write(outcoup,4) "   write(*,*)"
      write(outcoup,5) "    'GoSam '//param(1:len_trim(param))//"//
     $     " ' reports an error.'"
      write(outcoup,4) "   write(*,*) 'The POWHEG BOX aborts'"
      write(outcoup,4) "   call exit(1)"
      write(outcoup,4) "endif"
      write(outcoup,4) "end"
      close(incoup)
      close(outcoup)
      
c$$$      open(unit=inproc,status='old',file='../init_processes.f')
c$$$      open(unit=outproc,file='init_processes_new.f')
c$$$      found=.false.
c$$$      do i=1,maxnumlines
c$$$         read(inproc,'(a)',err=888,end=999) stringa
c$$$         if (stringa(1:31).eq.'      subroutine init_processes') 
c$$$     $   then
c$$$            found=.true.
c$$$         endif
c$$$         if (found.and.stringa(1:26).eq.'      call gosam_flst_born')
c$$$     $   then
c$$$            stringa(1:26)='c     call gosam_flst_born'
c$$$c            found=.false.
c$$$         endif
c$$$         write(outproc,1) stringa(1:len_trim(stringa))
c$$$      enddo
c$$$ 888  continue
c$$$ 999  continue
c$$$      if (.not.found) then
c$$$         write(*,*) "************************************************"
c$$$         write(*,*) "Have NOT found the string 'call gosam_flst_born'"
c$$$         write(*,*) "Check that this string is commented in the code!"
c$$$         write(*,*) "************************************************"
c$$$      endif
c$$$      close(inproc)
c$$$      close(outproc)
 1    format(A)
 4    format('      ',A)
 41   format('      ',A,A)
 5    format('     $',4x,A)
 6    format('C     ',A)
 7    format('C     ',A,A)
      end


      subroutine write_virtual
      implicit none
      include '../nlegborn.h'
      integer maxnumlines
      parameter (maxnumlines=100000)
      integer maxprocvirt
      parameter (maxprocvirt=10000)

      integer in, out,i,j,k
      integer vflav(nlegborn),code, proc
      integer coup,apower,aspower
      character * 11 coupling
      character * 2 ok
      character * 2 arrow
      character * 1 pipe     
      integer vflav_gosam(1:nlegborn,0:maxprocvirt-1)
 
      in=10
      out=11
      open(unit=in,status='old',file='orderfile.olc', err=200)
      open(unit=out,file='virtual_new.f')

c     START WRITING      
      write(out,6) "This file is generated AUTOMATICALLY by the"
      write(out,6) "write_pwhg_files.f program."
      write(out,6) "Do NOT edit by hand"
      write(out,4) ""
      write(out,6) "returns 2 Re(M_B * M_V)/(as/(2pi)), "
      write(out,6) "where M_B is the Born amplitude and "
      write(out,6) "M_V is the finite parte of the virtual amplitude."
      write(out,6) "The as/(2pi) factor is attached at a later point."
      write(out,6) "The virtual amplitude is generated using GoSam. "
      write(out,4) "subroutine setvirtual(p,vflav,virtual)"
      write(out,4) "implicit none"
      write(out,4) "include 'nlegborn.h'"
      write(out,4) "include 'PhysPars.h'"
      write(out,4) "include 'pwhg_st.h'"
      write(out,4) "include 'pwhg_flst.h'"
      write(out,4) "include 'pwhg_math.h'"
      write(out,4) "real * 8 p(0:3,nlegborn)"
      write(out,4) "integer vflav(nlegborn)"
      write(out,4) "real * 8 virtual"
      write(out,4) ""
      write(out,4) "integer proc, i"
      write(out,4) "integer vflav_gosam(1:nlegborn,0:maxprocborn-1)"
      write(out,4) "logical equalintlists"
      write(out,4) "external equalintlists"
      write(out,4) "integer dim_mom_array"
      write(out,4) "parameter (dim_mom_array=50)"
      write(out,4) "real * 8 pgosam(dim_mom_array)"     
      write(out,6) "real * 8 pgosam(5*nlegborn)"
      write(out,4) "real * 8 params(10),muren,res(4)"

c READ CONTRACT FILE
      do i=1,maxnumlines
         read(in,*,err=666,end=777) coupling, coup, pipe, ok  
         if (coupling .eq. 'AlphaPower') then
            apower = coup
         endif
         if (coupling .eq. 'AlphasPower') then
            aspower = coup
         endif
c         write(*,*) coupling, ' ==>', coup
 666     continue
      enddo
      write(*,*) 'The orderfile.olc contains more than ',maxnumlines
      write(*,*) ' lines. Increase it and run again'
      call exit(1)

 777  continue
      rewind(in)
      do i=1,maxnumlines
         read(in,*,err=888,end=999) vflav(1), vflav(2), arrow, 
     $        (vflav(k),k=3,nlegborn), pipe, code, proc
c         write(*,*) (vflav(k),k=1,nlegborn), ' ==>',proc
         write(out,1) proc
         do k=1,nlegborn
            if (vflav(k).eq.21) vflav(k)=0
            if (k.ne.nlegborn) then 
               write(out,2) vflav(k)
            else
               write(out,3) vflav(k)
            endif
         enddo
 888     continue
      enddo

c CONTRACT FILE READ, CONTINUE WRITING virtual.f FILE:
 999  continue
      write(out,4) ""
      write(out,4) "do i=0,flst_nborn-1"
      write(out,7) "if (equalintlists(nlegborn,vflav,vflav_gosam",
     $"(1,i))) then"
      write(out,4) "      proc=i"
      write(out,4) "      goto 222"
      write(out,4) "   endif"
      write(out,4) "enddo"
      write(out,*) "111  write(*,*) 'NO matching flavour string "//
     $     "between POWHEG and GoSam'"
      write(out,4) "write(*,*) 'PROGRAM ABORT'"
      write(out,4) "call exit(1)"
      write(out,4) ""
      write(out,*) "222  call gosam_momenta(p,pgosam)"
      write(out,4) ""
      write(out,4) "muren=sqrt(st_muren2)"
      write(out,4) "params(1)=1d0"
      write(out,4) ""
      write(out,5) "call OLP_EvalSubProcess(proc,pgosam,",
     $"muren,params,res)"
      write(out,4) "virtual=res(3)"
      write(out,4) ""
      write(out,11)"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
     $"CCCCCCCCCCCCCCCCCCCCCCCCC"
      write(out,6) "GOSAM returns Virtual with NO gs factor ==>"
      write(out,6) "virt_gosam ->  virt_gosam * (gs^2)^AlphasPower =  "
      write(out,6) "virt_gosam * (4*pi*st_alpha)^AlphasPower"
      write(out,11)"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
     $"CCCCCCCCCCCCCCCCCCCCCCCCC"
      write(out,8) "virtual=virtual * (4*pi*ph_alphaem)**",apower
      write(out,9) "* (4*pi*st_alpha)**",aspower
      write(out,6) 'The as/(2pi) factor is attached at a later point'
      write(out,6) 'We have then to multiply for 2*pi'
      write(out,10)'* (2*pi)' 
      write(out,4) "end"
      close(in)
      close(out)

 1    format('      data(vflav_gosam(i,',i5,'),i=1,nlegborn)/')
 2    format('     $   ',i5,',') 
 3    format('     $   ',i5,'/')
 4    format('      ',A)
 5    format('      ',A,A)
 6    format('C     ',A)
 7    format('         ',A,A)
 8    format('      ',A,I2)
 9    format('     $                ',A,I2)
 10   format('     $                ',A)
 11   format('CCCCCC',A,A)
      return
 200  write(*,*) 'orderfile.olc NOT found'
      call exit(1)
      end
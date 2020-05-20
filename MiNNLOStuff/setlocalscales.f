      subroutine setlocalscales(iuborn,imode,rescfac)
c     returns the rescaling factor including sudakov form factors and
c     coupling rescaling, for Born (imode=1) and NLO corrections (imode=2)
c     (imode=3) is as imode=1, but with the computation of the d3 term. 
c      use libnnlops
      use rad_tools
      use pdfs_tools
      use sudakov_radiators
      use internal_parameters, pi_hoppet => pi, cf_hoppet => cf, ca_hoppet => ca, tf_hoppet => tf
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'minnlo_flg.h'
      integer iuborn,imode
      real * 8 rescfac,expsud,sudakov,pwhg_alphas,resc_coupl
      real * 8 ptb2,mb2,mu2,alphas,alphas_ptb2,b0,optb2,omb2,orescfac,omuf2
      real * 8 pb(0:3), ptmax, yh, sudakov_form_factor
      integer oimode,i,flav
      common/flav_initial_state/ flav
      save optb2,omb2,orescfac,oimode,omuf2
      data optb2/-1d0/
      logical ini,bmass_in_minlo_flg
      common/c_bmass_in_minlo_flg/bmass_in_minlo_flg
      data ini/.true./
      save ini
      real * 8 powheginput,factsc2min,frensc2min,as,y,b1,tmp,bfact
      save factsc2min,frensc2min,b0,b1
      integer imax
      double precision D1, D2, D3, xx1, xx2
      integer nflav, pdf_set
      parameter(nflav=5) ! BEWARE OF THE FLAVOUR: nflav must match LHAPDF's set
      double precision msqB(-nflav:nflav,-nflav:nflav), msqV1(-nflav:nflav,-nflav:nflav), msqV2(-nflav:nflav,-nflav:nflav)
      double precision pborn(1:nlegborn-1,1:4), L, pt, jacob
      character*100 pdf_name
      double precision d3terms
      common/d3terms/d3terms
      double precision pborn_UUB(0:3,nlegborn-1),e2,e,ebeam
      double precision Delta1,Delta2,modlog
      double precision s3, kappar, kappaf
      double precision mb2_sav, kappar_sav, kappaf_sav
      save mb2_sav, kappar_sav, kappaf_sav
      double precision, parameter :: accuracy=1d-2
      double precision alphas_cutoff_fact, lambda
      save alphas_cutoff_fact
      logical use_analytic_alphas
      parameter (use_analytic_alphas = .false.)
      double precision, save :: eps
      logical flg_profiledscales
      save flg_profiledscales
      
      
c     the only purpose of this is to invalidate the currently stored results
c     of setlocalscales (see the line  if(imode.eq.oimode ...))
      if(imode == -1) then
         oimode = -1
         return
      endif
      

c     All the process dependence in this subroutine is and must be encoded through 
c     flg_minnloproc

      if (flg_minnloproc == 'H') then 
c     gg -> H production
c     Sudakov for a gluon 
         flav=0
      elseif (flg_minnloproc == 'Z') then 
c     Sudakov for a quark
         flav=1                 ! any value different from zero
      else
         write(*,*) ' setlocalscales: flg_minnloproc is ',flg_minnloproc
         write(*,*) ' not handled now, exiting ...'
         call exit(-1)
      endif

      if(ini) then

c     set the number of flavour used in the Dterms
         call set_nflav(nflav)
         
c     this to make sure that at the line   if (.. abs(kappar-kappar_sav) > 1d-7 )
c     the if block is entered on the first call
         mb2_sav=-1d0
         kappar_sav=-1d0
         kappaf_sav=-1d0

c     (process dependent) precision on the colourless system mass difference
c     to decide if events are identical
         if(flav == 0) then
            eps = 1d-3
         elseif(flav == 1) then
            eps = 1d-7
         endif

         bmass_in_minlo_flg = powheginput("#bmass_in_minlo").eq.1

         b0=(33d0-2d0*st_nlight)/(12*pi)
         b1=(153d0-19d0*st_nlight)/(24*pi**2)
         ini = .false.

c     read the modified logarithms parameter p, if it's not found, set it to its default (p=6)
         modlog_p = powheginput("#modlog_p")
         if(modlog_p < 0) modlog_p=6
         if(modlog_p.gt.0) then
            flg_modlog=.true.
            write(*,*) "modified log, exponent p = ",modlog_p
         else
            flg_modlog = .false.
            write(*,*) "ERROR: modlog_p must be > 0"
            write(*,*) "minnlops without modified log ",
     c           "not fully tested yet"
            call exit(-1)
            write(*,*) "sharp cutoff at hard scale, no mod log"
         endif


c     alphas_cutoff_fact**2 * pdf_q2min is the scale (in GeV^2) at which we freeze the running of alphas
c     pdf_q2min is the LHAPDF cutoff (e.g. for NNPDF30 it's 1GeV^2)
            
         if (.not. flg_hoppet_initialized) then
            pdf_cutoff_fact = powheginput('#pdf_cutoff_fact')

c     The conditions below should ultimately be replaced by a condition on pt/M
            if (flg_minnloproc == 'Z') then 
               if(pdf_cutoff_fact.lt.0d0) pdf_cutoff_fact=1.8d0
            elseif (flg_minnloproc == 'H') then 
               if(pdf_cutoff_fact.lt.0d0) pdf_cutoff_fact=2.5d0
            else
               write(*,*) ' setlocalscales: flg_minnloproc is ',flg_minnloproc
               call exit(-1)
            endif
         endif

         alphas_cutoff_fact = powheginput('#alphas_cutoff_fact')
         if(alphas_cutoff_fact.lt.0d0) alphas_cutoff_fact=pdf_cutoff_fact


         write(*,*) 'freezing of alphas in setlocalscales at [GeV] '
     $        ,sqrt(alphas_cutoff_fact**2 *pdf_q2min)

         write(*,*)'freezing of PDF evolution in '
     $        ,'setlocalscales at [GeV]'
     $        ,sqrt(pdf_cutoff_fact**2 *pdf_q2min)

         
         if (flg_minnlo) then
            write(*,*) '========================='
            write(*,*) '*MINNLO* ACTIVATED'
            write(*,*) '========================='
!     set the MiNNLO flags
            flg_include_delta_terms = .false.
            if(powheginput('#inc_delta_terms').eq.1)
     $           flg_include_delta_terms=.true.
            write(*,*) "include_delta_terms is ",flg_include_delta_terms

            flg_distribute_by_ub = .true.
            if(powheginput('#distribute_by_ub').eq.0)
     $           flg_distribute_by_ub = .false.
            write(*,*) "distribute_by_ub is ",flg_distribute_by_ub

            flg_distribute_by_ub_AP = .true.
            if(powheginput('#distribute_by_ub_AP').eq.0)
     $           flg_distribute_by_ub_AP = .false.
            write(*,*) "distribute_by_ub_AP is ",flg_distribute_by_ub_AP

            
            flg_dtermsallorders = .true.
            if(powheginput('#d3allorders').eq.0)
     $           flg_dtermsallorders = .false.
            write(*,*) '============================================='
            write(*,*) "flg_dtermsallorders = ",flg_dtermsallorders
            write(*,*) '============================================='
         else
c$$$            write(*,*) "ERROR: running at nlo (flg_minnlo=false)",
c$$$     c           " not fully tested yet"
c$$$            call exit(-1)
            write(*,*) '====================================='
            write(*,*) 'MINLO ACTIVATED (BUT *MINNLO* IS OFF)'
            write(*,*) '====================================='
         endif
c This also calls init_anom_dim
         call init_Dterms(flav,flg_minnlo)

c Setup profiled scales (main switch and parameters)    
         flg_profiledscales=powheginput('#profiledscales').eq.1
         Q0 = powheginput('#Q0')
         if(Q0.lt.0d0) Q0 = 2d0
         npow = powheginput('#npow')
         if(npow.lt.0d0) npow = 1d0              
         call init_profiled_scales_parameters(flg_profiledscales,Q0,npow)

         if(flg_profiledscales) then
            write(*,*) '============================================='
            write(*,*) 'Using profiled scales: (Q0,npow) = ',Q0,npow
            write(*,*) '============================================='
         endif
      endif

      rescfac = 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   These lines are process dependent!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     pb(0:3) is the colourless "boson" momentum
      pb(:)=0d0
c     sum over colourless particles (they must all come from a single boson decay)
      do i=3,nlegborn
c     the sequence of colourless particles is unchanged in the Born and in the real flavour list
         if (abs(flst_born(i,iuborn)).gt.6) then
            if(flg_minlo_real) then
               pb(:)=pb(:) + kn_cmpreal(:,i)
!               pb(:)=pb(:) + kn_cmpborn(:,i) ! alternative (consistent) option to define pT in reals
            else
               pb(:)=pb(:) + kn_cmpborn(:,i)
            endif
         endif
      enddo

      ptb2 = pb(1)**2 + pb(2)**2 ! transverse momentum squared
      mb2  = pb(0)**2 - pb(3)**2 - ptb2 ! invariant mass squared

      if(imode.eq.oimode.and.ptb2.eq.optb2.and.mb2.eq.omb2) then
         rescfac=orescfac
         st_mufact2=omuf2
         return
      else
         optb2=ptb2
         omb2=mb2
         oimode=imode
      endif


c     if we are running the original MiNLO, the code will do what it used to do.
c     If we run MiNNLO *with* modified logs, we can't have a theta function here, the Sudakov is
c     automatically switched off by the modified logs (see sudakov, expsudakov and expsudakov2).
      if(ptb2.gt.mb2.and.(.not.flg_modlog)) then
c     OLD MiNLO: this might be slightly inconsistent, with the upper limit of the sudakov integral (see subroutine sudakov)
c     Here there's a slight issue (only if running without
c     modlogs and both with MiNLO and MiNNLO): we are switching off the Sudakov and the D3 at pt = M.
c     This will cause a small step in the differential pT spectrum at this scale, that technically is
c     not smooth. We expect in practice to be small, but we haven't tested it so far.
         rescfac = 1d0
         expsud  = 0d0
         d3terms = 0d0
         return
      end if 


c     first define the logarithms
      pt = sqrt(ptb2)


c     Here we assume that B(2) is fixed (not dependent upon the kinematics).
c     In this case, all the following is needed only when the D3 contribution
c     is evaluated.
      if(imode == 3) then
         call invISRmap(kn_pborn,nlegborn,pborn_UUB)
         call get_B_V1_V2(pborn_UUB,msqB,msqV1,msqV2)
         e2 =  kn_sbeams
         e = sqrt(e2)
         ebeam = e/2d0        
         xx1 = pborn_UUB(0,1)/ebeam
         xx2 = pborn_UUB(0,2)/ebeam   
      endif
C     this call arranges for the A and B resummation coefficients to be set
C     as well as the virtual corrections
      kappar = st_renfact
      kappaf = st_facfact
      if (abs(mb2-mb2_sav) > 0 .or. abs(kappar-kappar_sav) > 0
     1     .or. abs(kappaf-kappaf_sav) > 0) then
         call reset_dynamical_parameters(dsqrt(mb2), kappar, kappaf)
         call init_anom_dim
         mb2_sav = mb2
         kappar_sav = kappar
         kappaf_sav = kappaf
      endif

c     Here we compute the modified log with KR=1 as the KR dependence is
c     explicitly accounted for in the computation of the D terms
      if(flg_modlog) then
         L = modlog(pt,cs%Q,1d0) ! modified logarithm
         jacob = 1d0/pt*(cs%Q/pt)**modlog_p/(1d0 + (cs%Q/pt)**modlog_p) ! dL/dpt
      else
         L = 0d0
         jacob = 0d0
         if (pt <= cs%Q) then
            L = log(cs%Q/pt)
            jacob = 1d0/pt      ! dL/dpt
         endif
      endif

      if(imode == 3) then
         call virtual_scale_dep(msqB, msqV1, msqV2)
      endif
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     alpha_s reweighting
      if(flg_profiledscales) then
         mu2 = st_renfact**2 * (sqrt(mb2) * exp(-L) + Q0 / (1d0 +
     $        (sqrt(mb2)/Q0*exp(-L))**npow))**2
         st_mufact2 = st_facfact**2 * (sqrt(mb2) * exp(-L) + Q0 / (1d0 +
     $        (sqrt(mb2)/Q0*exp(-L))**npow))**2
      else
         mu2        = st_renfact**2*mb2*exp(-2d0*L)
         st_mufact2 = st_facfact**2*mb2*exp(-2d0*L)
      endif

      
      omuf2=st_mufact2
      if(pdf_alphas_from_pdf) then
         if (use_analytic_alphas) then           
            if(flg_profiledscales) then
               lambda=st_alpha*b0*log(sqrt(mb2) / (sqrt(mb2)*exp(-L) +
     $              Q0 / (1d0 + (sqrt(mb2)/Q0*exp(-L))**npow)))
            else
               lambda = st_alpha*b0*L
            endif
            alphas = st_alpha / (1-2d0*lambda) * (1d0 
     &           - st_alpha / (1-2d0*lambda) * b1/b0 * log(1d0-2d0*lambda))
         else
            alphas = pwhg_alphas(mu2,st_lambda5MSB,st_nlight)
            if((mu2.lt.alphas_cutoff_fact**2 *pdf_q2min) .and. (.not.flg_profiledscales)) then
                 alphas = pwhg_alphas(alphas_cutoff_fact**2 *pdf_q2min,st_lambda5MSB,st_nlight)
            endif
         end if  
      else
         write(*,*) "ERROR: pdf_alphas_from_pdf=false ",
     c        "not supported for minnlops"
         call exit(-1)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      !>> ER+PM: use numerical integration for Sudakov
      sudakov_form_factor = Sudakov_pt_exact(L, alphas_cutoff_fact*sqrt(pdf_q2min), use_analytic_alphas)
      Delta1  = expsudakov_pt(L)  * alphas
      Delta2  = exp2sudakov_pt(L) * alphas**2
      rescfac = sudakov_form_factor

      if(imode.eq.2) then
         rescfac = rescfac * (alphas/st_alpha)**2
         if (flg_include_delta_terms.and.flg_minnlo) rescfac = rescfac * (1d0 + Delta1)
      elseif(imode == 1) then
         rescfac = rescfac * (alphas/st_alpha)
         if(.not.flg_bornonly) then
            if (.not. flg_minnlo) then 
               rescfac = rescfac * 
     1              (1 + Delta1 + alphas*b0*log(mu2/st_muren2))
            else
               if (flg_include_delta_terms) then
                  rescfac = rescfac *
     1                 (1d0 + (Delta1+alphas*b0*log(mu2/st_muren2) * (1d0 + Delta1)) + (Delta1**2/2d0+Delta2))
               else
                  rescfac = rescfac * (1d0 + (Delta1+alphas*b0*log(mu2/st_muren2)))
               end if
            endif
         endif
      endif

      if (flg_minnlo .and. imode .eq. 3) then
         
         call reset_dynamical_C_matxs()

         if(flg_dtermsallorders) then
            call DtermsAllOrders(D1, D2, D3, exp(-L), xx1, xx2, msqB, msqV1, msqV2, alphas)
         else
            call Dterms(D1, D2, D3, exp(-L), xx1, xx2, msqB, msqV1, msqV2, alphas)
         endif
        
         
         if (flg_include_delta_terms) then
            d3terms = D3 - D2 * Delta1 + D1 * Delta1**2/2d0 - D1 * Delta2
         else if (flg_uubornonly) then
            d3terms = D1+D2+D3
         else
            d3terms = D3
         end if
         d3terms = d3terms * jacob * sudakov_form_factor
         return
      endif

      if(bmass_in_minlo_flg) then
         call bmass_in_minlo(bfact,alphas)
         rescfac = rescfac * bfact
      endif
      orescfac=rescfac

      end


      subroutine init_Dterms(flav,minnlo)
      use rad_tools
      use pdfs_tools
      use internal_parameters, pi_hoppet => pi, cf_hoppet => cf, ca_hoppet => ca, tf_hoppet => tf
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'minnlo_flg.h'
      double precision D1, D2, D3, xx1, xx2, get_M_for_init_Dterms
      integer pdf_set, flav
      character*100 pdf_name
      character*100 string
      integer stringlength
      real *8 kappar,kappaf
      logical ini,minnlo
      data ini/.true./
      save ini
      real *8 powheginput
      logical check_pdf
      parameter (check_pdf = .true.)
      


      !************************************************************************
      ! initialise pdf (already done in hoppetif.f, see there)
      string = " "
      call lhapdfname(pdf_ndns1,string,pdf_set)
      pdf_name = trim(string(1:stringlength(string)-1))//trim(".LHgrid")
      ! set the global mass of the colour singlet (this has to happen before the pdf initialisation)
      cs%M = get_M_for_init_Dterms()
      cs%Q = get_M_for_init_Dterms()
      call init_masses_Dterms() ! initializes hard-coded mass values in NNLOPS_plugin
      if(.not. flg_hoppet_initialized) then
         if(flg_use_NNLOPS_pdfs) then
            call init_pdfs_NNLOPS(pdf_name, pdf_set, cutoff_high_fact=pdf_cutoff_fact)
         else
            call init_pdfs_from_LHAPDF(pdf_name, pdf_set, pdf_cutoff_fact)
         endif         
      endif
      
      if(check_pdf) call checkpdf


      !************************************************************************

      kappar=st_renfact
      kappaf=st_facfact

      if (flav.eq.0) then ! gg-initiated
         call set_process_and_parameters('pp', 'gg', sqrt(kn_sbeams), get_M_for_init_Dterms(), kappar, kappaf) ! M, KR, KF
      else ! qqbar-initiated
         call set_process_and_parameters('pp', 'qq', sqrt(kn_sbeams), get_M_for_init_Dterms(), kappar, kappaf) ! M, KR, KF
      endif

c     initialise anomalous dimensions
      call init_anom_dim

      if(minnlo)then
         ! initialise coefficient functions
         call init_C_matxs()
      endif

      end

      double precision function modlog(pt,m_singlet,KR)
      include 'minnlo_flg.h'
      double precision pt, m_singlet, KR
      modlog = 1d0/modlog_p*log((m_singlet/pt)**modlog_p + 1d0) - log(KR)
      end

      subroutine sudakov_exponent
      write(*,*) ' sudakov_exponent should not be called here'
      call exit(-1)
      end
      

      subroutine checkpdf
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      real *8 xxx(1:4),muF,pdf1(-pdf_nparton:pdf_nparton)
      integer i

      write(*,*) '==============================================='
      write(*,*) 'LHAPDF'
      write(*,*) '==============================================='

      xxx(1)=0.001d0
      xxx(2)=0.01d0
      xxx(3)=0.05d0
      xxx(4)=0.1d0

      muF=0.75d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo

      muF=1.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo

      muF=2.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      muF=5.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      muF=91.d0
      st_mufact2 = muF**2
      write(*,*) 'muF = ',muF
      do i=1,4
         call pdfcall(1,xxx(i),pdf1)
         write(*,*) 'x,g,d,u,c,b -> ',xxx(i),pdf1(0),pdf1(1),pdf1(2),pdf1(4),pdf1(5)
      enddo
      
      write(*,*) '==============================='
      end
      
      


!== see "readme" at bottom of file for explanation of code

! External Functions  
! CHARACTER*8        INT2STRING
! COMPLEX*16         DETERMINANTC
! INTEGER            MATMUL
! REAL*8             DETERMINANT
! REAL               ETIME,SUM,TRANSPOSE
! subroutines:
!CLOSEALLFILES
!DISTRIBUTIONSCRITICALGAIN,
!FLUSH
!FREEMATRIXCREATION, SCATMATRIXCREATION
!INPUTVARIABLES, OUTPUTFILES,
!INVERSE_C_MATRIX,       
!JWDUNITFLUX
!MEDIUMCREATION,   

program transfer_matrix_2d

  implicit none !== all of the variables should be of an explicit type

  logical, save :: writeAccuracy
  logical, save :: writeReflectionMatrix
  logical, save :: writeRunningStatus
  logical, save :: writeTmRmDet
  logical, save :: writeTransmissionMatrix
  logical, save :: writeMedium
  logical, save :: readMedium
  logical, save :: finddx
  logical, save :: writeDeterminats
  logical, save :: writeSelfEmbedTmRm
  logical, save :: writeScatFreeMatrix
  logical, save :: writeFrequency
  logical, save :: writeGain
  logical, save :: readGain
  logical, save :: varyStrength
  logical, save :: writeFinish
  logical, save :: writeElecField
  logical, save :: writeElecFieldError
  logical, save :: findElecField
  logical, save :: writeABcoefficientsUF
  logical, save :: writeABcoefficientsPW
  logical, save :: findJWDplanewave
  logical, save :: findJWDunitflux
  logical, save :: writeEachJW
  logical, save :: writeSummedJW
  logical, save :: doefse
  logical, save :: pausable
  logical, save :: writeTEuf
  logical, save :: writegEpw
  logical, save :: writegEuf
  logical, save :: writeTEufdist
  logical, save :: findcriticalgain
  logical, save :: repeatgainloop
  logical, save :: writecriticalgain
  logical, save :: writedistcriticalgain
  logical, save :: writeconductionpwdist
  logical, save :: writeconductionufdist
  logical, save :: incrementedgainlastcycle
  logical, save :: writeaveTab
  logical, save :: writeufJWDchan
  logical, save :: writepwJWDchan
  logical, save :: findeignvalues
  logical, save :: firstpropagation
  logical, save :: writecriticalgainduetovfwderror
  logical, save :: writeEzyPW
  logical, save :: writeEzyUF
  logical, save :: find_dzy
  logical, save :: dzy_gain_resolved
  logical :: finished_status_exists  !== boolean to test whether "out_finished.status" file exists
  logical :: inquire ! function to see if files exist
  
  !===================Definitions=========================
                          !== random number generator function by Jeff. platform/compiler independent
                          !== Dr Yamilov says:
  real            :: ran2 !== RAN2 SHOULD BE DEFINED AS "REAL", NOT "REAL*8"
                          !== IF YOU DO THIS SAME MISTAKE ONE MORE TIME, I WILL BE REALLLLY UPSET
                          !== note the number of "L"s.  
  real,pointer    :: Medium(:,:)     !== matrix that contains scatterers, is saved to a file called medium.dat
  real*8,pointer  :: FreeMatrix(:,:),& !== Transfer matrix the space between scatterers
                     DetFree(:),&  !== Vector used to store determinants of Free matrices
                     Kperp(:),&      !== k_perpendicular
                     Kpara(:),&      !== k_parallel
                     kparaMat(:,:),&  !== for matrix multiplication
                     freq(:),&       !== frequencies
                     gain(:),&       !== gain or absorption values
                     dxFixed(:,:),&
                     expon_temp(:,:),&
                     Jflux_l_sum(:,:),& !== J, W are indexed by position,gain
                     Jflux_r_sum(:,:),&
!                     Wdist_sum(:,:),&
                     Wdist_r_sum(:,:),&
                     Wdist_l_sum(:,:),&
                     energypassive(:),& !== passive waveguide energy for each channel
                     aveT(:,:,:),&
                     aveR(:,:,:),&
                     Wdist_r_n_sum(:,:,:,:),&
                     Wdist_l_n_sum(:,:,:,:),&
                     Jflux_r_n_sum(:,:,:,:),&
                     Jflux_l_n_sum(:,:,:,:),&
                     Wdist_pw_r_n_sum(:,:,:),&
                     Wdist_pw_l_n_sum(:,:,:),&
                     Jflux_pw_r_n_sum(:,:,:),&
                     Jflux_pw_l_n_sum(:,:,:),&
                     correlationE_sum(:,:),&
                     jflux_z(:,:),&
                     jflux_y(:,:),&
                     wdist_zy(:,:),&
                     wdist_zy_incohrnt(:,:),&
                     jflux_z_gr(:,:,:),&
                     jflux_y_gr(:,:,:),&
                     wdist_zy_gr(:,:,:),&
                     wdist_zy_gr_incohrnt(:,:,:),&
                     chi_y(:,:),&
                     dchi_y(:,:),&
                     uzy2(:,:),&
                     partialuwrty2(:,:),&
                     partialuwrtz2(:,:),&
                     y(:)

   real*8,pointer :: eignvalues(:),rwork(:) !== for eigenvalues zheev
   complex*16,pointer :: tthc(:,:),work(:) !== for eigenvalues zheev
   
  integer,pointer :: scatr_indx(:),&
                     tebinsuf(:,:,:,:),& !== T,E,T/E,log() distributions for each channel (input channel resolved)
                     gainbin(:),& !== critical gain
                     condbinspw(:,:,:),& !== g,E summed over input channels
                     condbinsuf(:,:,:),& !== g,E summed over input channels
                     alphasign(:),& !== for varying alpha
                     validgaincount(:)

  complex*16,pointer :: Rm (:,:),&   !== Self-embedding matrix (Reflection) (S(0,J+1) in notes)
                     Tm   (:,:),&    !== Self-embedding matrix (Transmission) (S(J+1,J+1) in notes)
                     Zm   (:,:),&    !== Sigma(J) in notes (part of Self-embedding calculation): total medium matrix
                     Zmt  (:,:),&    !== Used as temp matrix for the entire medium result
                     Gm   (:,:),&    !== Used for initial Tm
                     Hm   (:,:),&    !== Used for initial Tm
                     Im   (:,:),&    !== Identity Matrix
                     Ref  (:,:),&    !== Final Reflection Matrix
                     Trans(:,:),&    !== Final Transmission Matrix
                     Tempse(:,:),&     !== Dr Yamilov update for electric-field self-embedding, 20090705
                     SnN(:,:,:),&    !== Dr Yamilov update for electric-field self-embedding, 20090705
                     ChunkMatrix(:,:),& !== Transfer matrix for a self-embedding segment of waveguide
                     ScatteringMatrix(:,:),& !== Transfer matrix across a scatterer
                     InvFreeMatrix(:,:),& !== for finding electric field
                     InvScatMatrix(:,:),& !== for finding electric field
                     DetScat(:),&  !== Vector used to store determinants of Scattering matrices
                     DetChunk(:),& !== Vector used to store determinants of Chunk Matrix
                     Vr   (:),& !== Vector with reflected E-field (and derivatives) info for each channel
                     Vt   (:),& !== Vector with transmitted E-field(and derivatives) info for each channel
                     Vbc  (:),&         !== Vector used for boundary conditions
                     expon(:,:),& !== exponent for J,W (finding D)
                     Apwscatsum(:,:),&
                     Bpwscatsum(:,:),&
                     uzy(:,:),&
                     partialuwrty(:,:),&
                     partialuwrtz(:,:),&
                     AB_inside_trans(:,:) !== find_dzy

!== jwp, array
  complex*16,pointer  :: vector(:,:)
  complex*16,pointer  :: A_raw_jpw(:,:),&
                 B_raw_jpw(:,:),&
                 A_inside_jpw(:,:),&
                 B_inside_jpw(:,:)
  real*8,pointer      :: Apw2_jpw(:,:),&
                 Bpw2_jpw(:,:),&
                 Apwscat2_jpw(:,:),&
                 Bpwscat2_jpw(:,:),&
                 energy_integrand_jpw(:)

  !== jwp, scalar
  real*8  :: pwchannorm_jpw,k_jpw,integrand_jpw
  real*8  :: condpw_jpw,energypw_jpw,emptypwwg_jpw,energy2pw_jpw
  real :: DeltaX_jpw, DeltaX
  integer :: a_jpw,c_jpw,d_jpw,h_jpw,f_jpw,nz_jpw
  integer :: numscatterersplustwo!,inputchan  

  complex*16      :: determinantc    !== function to find determinant of a complex matrix
  integer,pointer :: Indx(:)    !== used in multiple subroutines
  integer         :: a,&        !== index for many small loops. Also used in medium creation
                     i,&        !== 
                     b,&        !== used in medium creation
                     c,&        !== used in medium creation
                     d,&        !== ?
                     ntot,&     !== size of all matrices, =2*(total number of channels)
                     g,&        !== FreeSpace matrix creation
                     p,&        !== FreeSpace matrix creation
                     q,&        !== FreeSpace matrix creation
                     M,&        !== total number of scatterers; in the .input
                     seed,&     !== random number generator seed [must be negative]
                     seed_interrupt,& !== used to read in from "seed.input" 
                                          !== file to see if seed was set to positive. 
                                          !== If yes, the pause [dump parameters to file and quit]
                     nmax,&     !== total number of channels, both open and closed
                     remainder,&!== for Free space calculation
                     se_step,&  !== number of scatterers in a self-embedding step. Set in parameters.input
                     M_se,&     !== number of chunks (how many self embedding steps)
                     z,&        !== index for self-embedding interval
                     zz,&       !== index for Dr Yamilov's addon, 20090705
                     fli,&      !== frequency loop index
                     gli,&      !== gain loop index
                     sli,&      !== scatterer loop index, for electric field
                     n_omg,&    !== number of frequencies
                     n_alphaa,& !== number of gain/absorption steps
                     n_closed,& !== number of closed channels
                     n_open,&   !== number open channels
                     rlz,& !== which realization is being propagated. Loop index
                     n_rlz,& !== number of realizations. Set in parameters.input
                     seed_counter,& !== used for "pause" functionality. 
                                          !== counts how many times the "ran2" function is called inside medium loop
                     start_rlz,&!==  used for "pause" functionality. which realization to start on
                     inputchan,&
                     binsize,&    !== T,E,T/E,log() distributions
                     n_decgs    !== number of decrements of gain step (for finding critical gain)
  
  real*8          :: L,&     !== system length; in the .input
                     W,&     !== width of the waveguide; in the .input
                     alphas,&   !== scattering strength; in the .input
                     alphaa,&   !== active scattering
                     alphaai,&  !== initial gain for loop; in the .input
                     alphaaf,&  !== final gain for loop; in the .input
                     omg,&      !== current frequency being calculated
                     omgi,&     !== initial frequency
                     omgf,&     !== final frequency
                     pi,&       !== pi = 3.14 (constant)
                     ABstep,&   !== size(A_inside2) = 0:lambda/ABstep:L
                     lambda,&   !== wavelength
                     v_fwd_error,&
                     Lmin,&     !== minimum allowable distance between any two scatterers
                     T,&        !== 
                     determinant,&!== function to find determinant of a real matrix
                     K2,&       !== k^2 = (2*pi/lambda)^2
                     !!! WARNING: do not use rand(), since the output is compiler-dependent !!! 
                     !rand,&    !== "built-in" random number generator, 
                        !== but compiler-dependent. different 
                        !== output for ifort and gfortran.  
                     deviation_count,& !== count the number of times the 
                             !== determinant exceeds 1 by more than 9E-9
                     percntwriteJWD,&
                     initialgainstep,&  !==finding critical gain
                     condpw,& !== conductance (for finding distribution) g = sum_a T_a
                     conduf,&
                     total_energypw,& !== energy summed over all changels (for finding distribution)
                     total_energyuf,&
                     total_passive_energypw,&
                     total_passive_energyuf,&
                     conductvty,& !== is conductivity not defined in active media?
                     oneminusconductvty,&
                     Emin,Emax,Tmin,Tmax,&  !== unit flux distributions
                     Epwmin,Epwmax,gmin,gmax,& !== plane wave distributions
                     Eamin,Eamax,Tamin,Tamax !== channel resolved unit flux distributions
  real            :: FirstNumber,& !== used for seed manipulation
                     etime,&    !== etime, elapsed, and totaltime are 
                        !== three variable for timing the script
                     elapsed(2),& !== For receiving user and system time
                     totaltime  !== For receiving total time
  double precision starttime,endtime

  !== finding critical gain
  integer :: decgs !== decremented gain step index
  real*8 :: deltagain,T_oldest,T_old,T_prev,T_current,criticalgain  

  !== for output file names
  character(len=50) :: file_name
  integer*4       :: ch0,ch1,ch2,ch3,ch4,ch5,ch6,ch7
  character*8     :: ch_n
  character*8              :: int2String
  !== for electric field
  complex*16,pointer :: v_fwd(:,:),v_bck(:,:),& !== forward and backward prop vectors
                        v_fwd_check(:)

  integer*4 rank,nproc,ierr,k  !== MPI 

  parameter (pi=3.14159265358979323846d+0) !== set the constant
  
  !== initialize deviation count (accuracy of determinant)
  deviation_count=0 !== used for file 126, "out_accuracy.dat"

  rank = 0 !== initialization of these MPI-specific variables is unnecessary, but simplifies serial f90 synchronization
  ierr = 0
  ! always gets written, no matter what [no "if-then-endif" needed]
  open(125,file='out_screen.dat',POSITION='APPEND') 
  if ((n_alphaa.gt.1).AND.findcriticalgain) then
     open(163,file='out_how_critical_gain_was_found.dat',POSITION='APPEND')
  endif
  open(165,file='out_where_jw_pw_was_written_from.dat',POSITION='APPEND')

  seed_counter=0 !== "FirstNumber always gets called, regardless of whether program is "fresh" or "resume" state

  !====== INPUT FILES ======
!  if (rank.eq.0) then 
  call inputvariables(seed,rank,omgi,omgf,n_omg,percntwriteJWD,&
                alphas,varyStrength,alphaai,alphaaf,n_alphaa,&
                L,W,M,n_closed,start_rlz,n_rlz,lambda,se_step,&
                nmax,ntot,n_open,M_se,Lmin,writeElecField,findeignvalues,&
                writeRunningStatus,writeFrequency,writeReflectionMatrix,&
                writeTransmissionMatrix,writeTmRmDet, writeAccuracy,& 
                writeScatFreeMatrix,writeMedium,writeDeterminats,& 
                writeSelfEmbedTmRm,writeGain,readGain,writeFinish,readMedium,&
                findElecField,writeABCoefficientsPW,writeABCoefficientsUF,writeElecFieldError,&
                findJWDplanewave,findjwdunitflux,doefse,ABstep,writeSummedJW,writeEachJW,pausable,&
                writeTEuf,writegEuf,writegEpw,binsize,initialgainstep,n_decgs,findcriticalgain,&
                writeTEufdist,writecriticalgain,writedistcriticalgain,&
                writeconductionufdist,writeconductionpwdist,writeaveTab,&
                writeufJWDchan,writepwJWDchan,writeEzyPW,writeEzyUF,find_dzy,dzy_gain_resolved)
!  endif

  !== distribution limits
     !== unit flux, input channel resolved
  Tamax = 4 !== tranmission never should exceed 1 (physical)
  Tamin = 0
  Eamax = 14  !== Emax probably scales as  W*L*M*n_open
  Eamin = 0  
     !== unit flux
  Tmax = 15  !== in diffusive regime, Tmax = n_open*10
  Tmin = 0
  Emax = 15  !== in diffusive regime, Emax ~ n_open
  Emin = 0
     !== planewave
  gmax = 2  !== in diffusive regime, gmax = 1
  gmin = 0
  Epwmax = 15
  Epwmin = 0

  !== now that all scalars have been set by input file, allocate arrays

  a = int(L*ABstep)+1
!== jwp, array
  allocate( vector(1:ntot,1:(M+2)))
  allocate( A_raw_jpw(1:nmax,1:(M+2)))
  allocate( B_raw_jpw(1:nmax,1:(M+2)))
  allocate( A_inside_jpw(1:nmax,1:a))
  allocate( B_inside_jpw(1:nmax,1:a))
  allocate( Apw2_jpw(1:nmax,1:a))
  allocate( Bpw2_jpw(1:nmax,1:a))
  allocate( Apwscat2_jpw(1:nmax,1:(M+2)))
  allocate( Bpwscat2_jpw(1:nmax,1:(M+2)))
  allocate( energy_integrand_jpw(1:a))

  !===== transfer-matrices =======================
  allocate(ChunkMatrix     (1:ntot,1:ntot))     !==
  allocate(DetFree         (1:2*M+2)      )     !==
  allocate(DetScat         (1:2*M+2)      )     !==
  allocate(DetChunk        (1:2*M+2)      )     !==
  allocate(Indx            (1:nmax)       )     !== nmax=n_open+n_closed
  allocate(ScatteringMatrix(1:ntot,1:ntot))     !== ntot=2*nmax
  allocate(FreeMatrix      (1:ntot,1:ntot))     !== free space propagation
  allocate(Im              (1:ntot,1:ntot))     !==

  allocate(Kperp(nmax)           )      !== 
  allocate(Kpara(nmax)           )      !== 
  allocate(freq (1:n_omg)        )      !== frequency
  allocate(gain (1:n_alphaa)     )      !== gain
  allocate(Ref  (1:nmax,1:n_open))      !== nmax=n_open+n_closed
  allocate(Trans(1:nmax,1:n_open))
  allocate(aveT (1:nmax,1:n_open,1:n_alphaa))
  allocate(aveR (1:nmax,1:n_open,1:n_alphaa))
  allocate(Vbc  (1:ntot)         )      !== ntot=2*nmax
  allocate(Vr   (1:ntot)         )
  allocate(Vt   (1:ntot)         )

  allocate(InvFreeMatrix(1:ntot,1:ntot))        !== 
  allocate(InvScatMatrix(1:ntot,1:ntot))        !== 

  !====== matrices to be used in self-embedding ====
  allocate(Gm    (1:ntot,1:ntot))       !== boundary condition
  allocate(Hm    (1:ntot,1:ntot))       !== boundary condition
  allocate(Tm    (1:ntot,1:ntot))       !== 
  allocate(Zm    (1:ntot,1:ntot))       !== "Z matrix"
  allocate(Zmt   (1:ntot,1:ntot))       !== "Z matrix" temporary
  allocate(Rm    (1:ntot,1:ntot))       !== 
  allocate(Medium(M,2)          )       !== M=total number of scatterers
  allocate(alphasign(M))

  !== electric field ============
  allocate(Tempse  (1:ntot,1:ntot))     !== self-embedding electric field
  allocate(SnN   (1:ntot,1:ntot,1:M_se))        !== self-embedding electric field
  allocate(v_fwd(1:ntot,1:M+2)) !== electric field, just after every scatterer
  allocate(v_bck(1:ntot,1:M+2)) !== same, but backward propagating (uses inverses)
  allocate(v_fwd_check(1:ntot))
     
  allocate(validgaincount(1:n_alphaa))

  allocate(tthc(1:nmax,1:n_open)) !== same size as transmission matrix
  allocate(eignvalues(nmax))
  allocate(rwork(3*nmax-2))
  allocate(work(2*nmax-1))

  !== converting from A,B to J,W,D =====
  !== note: there is an error catch to see if int(L*ABstep).NE.L
  allocate( dxFixed(1,1:a)        )       !== distance after each scatterer for propagating A,B [per rlz]
  allocate( expon(1:nmax,1:a) )  !== exponential propagation for A,B inside medium  [per rlz, per frequency]
  allocate( expon_temp(1:nmax,1:a)    )
  allocate( scatr_indx(1:a)          )  !== which scatterer is A,B supposed to use? [per rlz]
  allocate( kparaMat(1:nmax,1))   !== version of kpara that can be multiplied by dxFixed to get a matrix

  allocate( tebinsuf(1:8,1:binsize,1:n_open,1:n_alphaa))  !== input channel resolved
  allocate( gainbin(1:binsize) )
  allocate( condbinsuf(1:8,1:binsize,1:n_alphaa))
  allocate( condbinspw(1:8,1:binsize,1:n_alphaa))
  allocate( energyPassive(1:n_open) )


  if (findjwdunitflux) then
!    allocate( Jflux_l(1:a,1:n_alphaa) )
!    allocate( Jflux_r(1:a,1:n_alphaa) )
!    allocate( Wdist  (1:a,1:n_alphaa) )
    allocate( Jflux_l_sum(1:a,1:n_alphaa) )
    allocate( Jflux_r_sum(1:a,1:n_alphaa) )
!    allocate( Wdist_sum  (1:a,1:n_alphaa) )
    allocate( Wdist_r_sum  (1:a,1:n_alphaa) )
    allocate( Wdist_l_sum  (1:a,1:n_alphaa) )

!    allocate( Jflux_l_n    (1:n_open,1:n_open,1:a,1:n_alphaa) )
!    allocate( Jflux_r_n    (1:n_open,1:n_open,1:a,1:n_alphaa) )
!    allocate( Wdist_l_n    (1:n_open,1:nmax  ,1:a,1:n_alphaa) )
!    allocate( Wdist_r_n    (1:n_open,1:nmax  ,1:a,1:n_alphaa) )
    allocate( Jflux_l_n_sum(1:n_open,1:n_open,1:a,1:n_alphaa) )
    allocate( Jflux_r_n_sum(1:n_open,1:n_open,1:a,1:n_alphaa) )
    allocate( Wdist_l_n_sum(1:n_open,1:nmax  ,1:a,1:n_alphaa) )
    allocate( Wdist_r_n_sum(1:n_open,1:nmax  ,1:a,1:n_alphaa) )
  endif
  if (findJWDplanewave) then
    allocate( Apwscatsum(1:nmax,1:(M+2)) )
    allocate( Bpwscatsum(1:nmax,1:(M+2)) )
!     allocate( Wdist_pw_r_n    (1:nmax,1:a,1:n_alphaa) )
    allocate( Jflux_pw_l_n_sum(1:n_open,1:a,1:n_alphaa) )
    allocate( Jflux_pw_r_n_sum(1:n_open,1:a,1:n_alphaa) )
    allocate( Wdist_pw_l_n_sum(1:nmax,1:a,1:n_alphaa) )
    allocate( Wdist_pw_r_n_sum(1:nmax,1:a,1:n_alphaa) )
    allocate( correlationE_sum(1:int(L*ABstep),1:n_alphaa) )
    if (find_dzy.AND.(.NOT.dzy_gain_resolved)) then
      allocate( jflux_z(1:a,1:int(W*ABstep)+1))
      allocate( jflux_y(1:a,1:int(W*ABstep)+1))
      allocate( wdist_zy(1:a,1:int(W*ABstep)+1))
      allocate( wdist_zy_incohrnt(1:a,1:int(W*ABstep)+1))
      allocate( uzy(1:a,1:int(W*ABstep)+1))
      allocate( uzy2(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrty(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrtz(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrty2(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrtz(1:a,1:int(W*ABstep)+1))
      allocate( AB_inside_trans(1:a,1:nmax))
      allocate( chi_y(1:nmax,1:a))
      allocate( dchi_y(1:nmax,1:int(W*ABstep)+1))
      allocate( y(1:int(W*ABstep)+1))
    elseif (find_dzy.AND.dzy_gain_resolved) then
      allocate( jflux_z_gr(1:a,1:int(W*ABstep)+1,1:n_alphaa))
      allocate( jflux_y_gr(1:a,1:int(W*ABstep)+1,1:n_alphaa))
      allocate( wdist_zy_gr(1:a,1:int(W*ABstep)+1,1:n_alphaa))
      allocate( wdist_zy_gr_incohrnt(1:a,1:int(W*ABstep)+1,1:n_alphaa))
      allocate( uzy(1:a,1:int(W*ABstep)+1))
      allocate( uzy2(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrty(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrtz(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrty2(1:a,1:int(W*ABstep)+1))
      allocate( partialuwrtz2(1:a,1:int(W*ABstep)+1))
      allocate( AB_inside_trans(1:a,1:nmax))
      allocate( chi_y(1:nmax,1:a))
      allocate( dchi_y(1:nmax,1:int(W*ABstep)+1))
      allocate( y(1:int(W*ABstep)+1))
    endif
  endif

  if (readGain) then 
    open(167,file='quasi1d_rect_gain_array.input')
    do a=1,n_alphaa
      read(167,*) gain(a)
    enddo
  endif

  validgaincount(:) =0
  tebinsuf(:,:,:,:) = 0
  gainbin(:) = 0
  condbinsuf(:,:,:) = 0
  condbinspw(:,:,:) = 0
  if (findjwdunitflux) then
    Jflux_l_sum(:,:) = 0.0d+0
    Jflux_r_sum(:,:) = 0.0d+0
!    Wdist_sum  (:,:) = 0.0d+0
    Wdist_r_sum  (:,:) = 0.0d+0
    Wdist_l_sum  (:,:) = 0.0d+0
    Jflux_l_n_sum(:,:,:,:) = 0.0d+0
    Jflux_r_n_sum(:,:,:,:) = 0.0d+0
    Wdist_l_n_sum(:,:,:,:) = 0.0d+0
    Wdist_r_n_sum(:,:,:,:) = 0.0d+0
  endif
  if (findjwdplanewave) then
    Apwscatsum(:,:) = 0.0d+0
    Bpwscatsum(:,:) = 0.0d+0
    Jflux_pw_l_n_sum(:,:,:) = 0.0d+0
    Jflux_pw_r_n_sum(:,:,:) = 0.0d+0
    Wdist_pw_l_n_sum(:,:,:) = 0.0d+0
    Wdist_pw_r_n_sum(:,:,:) = 0.0d+0
    correlationE_sum(:,:)   = 0.0d+0
  endif

  Im(:,:) = 0.0d+0 !== Identity Matrix. Used for initialization and self-embed
  do a=1, ntot 
     Im(a,a) = (1.0d+0,0.0d+0)
  enddo

  do a=1, nmax !== Define Kpara and Kperp for all channels =======
     Kperp(a) = a*pi/W
  enddo
  
  aveR(:,:,:) = (0.0d+0)     !== Zero the matrices and vectors               
  aveT(:,:,:) = (0.0d+0)
  !====== Initializing self-embedding matrices "G" and "H" =====
  Gm(:,:) = (0.0d+0,0.0d+0)
  Hm(:,:) = (0.0d+0,0.0d+0)
  do a=1, ntot                               !creates Gm and Hm
     p = a + nmax
     q = a - nmax
     if (a.le.n_open) then
        Gm(a,a) = ( 1.0d+0, 0.0d+0)
        Gm(a,p) = ( 0.0d+0,-1.0d+0)
     else if (a.le.nmax) then   
        Gm(a,a) = ( 1.0d+0, 0.0d+0)
        Gm(a,p) = (-1.0d+0, 0.0d+0)
     else if (a.le.nmax+n_open) then
        Hm(a,q) = ( 0.0d+0, 1.0d+0)
        Hm(a,a) = (-1.0d+0, 0.0d+0)        
     else
        Hm(a,q) = ( 1.0d+0, 0.0d+0)
        Hm(a,a) = ( 1.0d+0, 0.0d+0)    
     endif
  enddo !== a loop, boundary condition matrices for self-embedding

  do rlz=start_rlz,n_rlz
     !== chop up the realization index
     ch0=int(    int(1.0*(rlz+(n_rlz*rank*1.0))+0.1)/10000000)
     ch1=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),10000000)/1000000)
     ch2=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),1000000)/100000)
     ch3=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),100000)/10000)
     ch4=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),10000)/1000)
     ch5=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),1000 )/100)
     ch6=int(mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),100  )/10)
     ch7=    mod(int(1.0*(rlz+(n_rlz*rank*1.0))+0.1),10   )

     !== append all four characters to a string
  ch_n=char(ch0+48)//char(ch1+48)//char(ch2+48)//char(ch3+48)//char(ch4+48)//char(ch5+48)//char(ch6+48)//char(ch7+48)

     !how to do collapse multiple strings: file_name="e-dft"//ch_n//ch_f//".dat" 
     !open(87,file=file_name,status='old',POSITION='append')

     !== for a specific realization, create output files
     call outputfiles(rank,ch_n,rlz,n_rlz,writeRunningStatus,writeFrequency,&
                         writeReflectionMatrix,writeTransmissionMatrix,&
                         writeTmRmDet,writeAccuracy, writeScatFreeMatrix, &
                         writeMedium,writeDeterminats,writeElecField, &   
                         writeSelfEmbedTmRm, writeGain, writeFinish,&
                         findElecField,writeABcoefficientsPW,writeABCoefficientsUF,writeElecFieldError,&
                         findJWDplanewave,findjwdunitflux,writeSummedJW,writeEachJW,pausable,&
                         percntwriteJWD,writeTEuf,writegEpw,&
                         writeTEufdist,writecriticalgain,writedistcriticalgain,&
                         writeconductionufdist,writeconductionpwdist,writeaveTab,&
                         writeufJWDchan,writepwJWDchan,findcriticalgain,writegEuf,find_dzy,dzy_gain_resolved)

     if(readMedium) then !== read the medium in from file "out_%rlz%_medium.dat"
        file_name = "out_"//ch_n//"_Medium.dat"
!OLD:   open (100, file='out_Medium.dat',POSITION='APPEND')
        open (100, file=file_name)
        do a=1,M
           read(100,*) Medium(a,1),Medium(a,2)
        enddo
     else   !== create new Medium
        Medium(:,:) = (0.0d+0) !== initialize

        if (findJWDplanewave.OR.findjwdunitflux.OR.writeABcoefficientsPW.OR.writeABcoefficientsUF.OR.writeTEuf &
           .OR.writeconductionpwdist.OR.writeconductionufdist.OR.writeEzyPW.OR.writeEzyUF) then
           finddx = .true.
        else
           finddx = .false.
        endif

        call mediumcreation(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
             finddx,dxFixed,scatr_indx,ABstep,rank) 
        !call mediumcreationBAD(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
        !      findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
        !call mediumcreationNOCOLLISIONCHECK(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
        !      findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
        !call mediumcreationSTACKING(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
        !     findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
        !call mediumcreationOLD(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
        !     findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
     endif

     !=== frequency loop ================
     do fli=1,n_omg

        !== set the current frequency
        if(n_omg.eq.1) then
           omg=omgi
        else
           omg=omgi+(omgf-omgi)*((fli-(1.0d+0))/(n_omg-(1.0d+0))) !== the 1.0 is to change the variable type to real
        endif
        !== record the frequencies to an array
        freq(fli)=omg
        if (writeFrequency.AND.(rlz.eq.start_rlz)) then        
           write(135,*) freq(fli)
           call flush(135)
        endif

        !== k^2
        K2 = (2*pi*omg/lambda)**2 !== omg is the dimensionless frequency

        !======= Define Kpara and Kperp for all channels =======
        Kpara(1:nmax) = sqrt(abs(K2 - Kperp(1:nmax)**2))   !== Kpara can be complex if Kperp^2 > K2
        ! Note Kpara is defined as |Kpara| so it is Kpara for open   channels
        !                                  and Kappa_para for closed channels

        expon(:,:) = (0.0d+0,0.0d+0)  !== complex precise zero
        expon_temp(:,:) = 0.0d+0

        if (finddx) then !== if "dxFixed" was found, it indicates you will be using it for expon()
           kparaMat(1:nmax,1) = kpara(1:nmax)
           expon_temp = matmul(KparaMat,dxFixed) !== per rlz, per frequency. This matrix is real-valued
        
!           do a= 1,(int(L*ABstep)+1)
!             write(*,*) expon_temp(1:nmax,a)
!           enddo

           do a = 1,nmax
              do b = 1,(int(L*ABstep)+1)
                 expon(a,b) = exp((0.0d+0,1.0d+0)*omg*expon_temp(a,b))
              enddo
           enddo
           !== original matlab version:
           !expon(:,:,fli) = exp((0,1)*omg*  matmul(Kpara(1:nmax),dxFixed)) !== per rlz, per frequency
        endif

        !=== gain loop ================

        gli = 0
        repeatgainloop = .true.
!        write(163,*) 'starting gain loop. fli,rlz=,',fli,rlz
        if (findcriticalgain) then
           write(111,*) 'starting gain loop. fli,rlz=,',fli,rlz
        endif
        do while (repeatgainloop)   !== previously, do gli=1,n_alphaa was used. Replaced to find critical gain
           gli=gli+1
           !== set the gain
           if (.NOT.findcriticalgain) then !== when not finding critical gain
              if((n_alphaa.eq.1).AND.(.NOT.readGain)) then !== only doing one gain value, use the initial gain value
                 alphaa=alphaai
              elseif (.NOT.readGain) then !== blindly increment gain a fixed amount
                 alphaa=alphaai+(alphaaf-alphaai)*((gli-(1.0d+0))/(n_alphaa-(1.0d+0))) !== the 1.0 is to change the variable type to real
              elseif (readGain) then
                 alphaa = gain(gli) ! read in array from file
              endif
              if (gli.eq.1) then
                 T_prev =0
              else
                 T_prev = sum((cdabs(Trans))**2)
              endif
           elseif (findcriticalgain.AND.(gli.eq.1)) then !== find passive conductance (g=0)
              !write(111,*) 'gli is 1'
              decgs = 0 !== decremented gain step index
              deltagain = initialgainstep !== set in input file
              alphaa = 0 !== passive
              repeatgainloop = .true.
              incrementedgainlastcycle = .false.
              T_oldest = 0 !== the following values allow the program to compare and then reset gain
              T_old    = 0 
              T_prev = 0
              T_current = 0
           elseif (findcriticalgain.AND.(gli.eq.2)) then !== first active system
              !write(111,*) 'gli is 2, fli is ',fli
              T_prev = sum((cdabs(Trans))**2)
              write(111,*) 'conductance is ',T_prev,' when gain is 0. gli,fli,rlz is ',gli,fli,rlz
              alphaa = alphaa + deltagain
              repeatgainloop = .true.
              !write(111,*) 'gain is ',alphaa
           elseif (findcriticalgain.AND.(gli.ge.3)) then
              T_current = sum((cdabs(Trans))**2)
              write(111,*) 'conductance is ',T_current,' when gain is ',alphaa,' and T_prev ',T_prev
              if (T_current.gt.1000) then   !== maximum conductance, far above diffusive maximum
                 criticalgain = alphaa !== since T>>T_diffusive_max, then call this the critical gain
                 write(111,*) 'TIME TO GO TO NEXT RLZ (T>1000). gli,fli,rlz',gli,fli,rlz
                 totaltime = etime(elapsed)
                 write(111,*) 'run time so far:', totaltime/60.0,' minutes'
                 write(148,*) rlz,fli,gli, criticalgain, decgs, 2 !'due to T>1000'
                 gli = n_alphaa !== to get J,W to write to file, tell them that gain is done
                 call distributionscriticalgain(criticalgain,gainbin,binsize)
                 call flush(148)
                 incrementedgainlastcycle = .false.
                 repeatgainloop = .false. !== go to next rlz
              elseif (incrementedgainlastcycle.AND.(T_current.ge.T_oldest)) then !== keep incrementing gain
                    T_old    = T_oldest
                    T_prev   = T_current
                    alphaa = alphaa + deltagain
                    incrementedgainlastcycle = .false.
                    !write(111,*) 'but this time, Tc>Toldest'
                    repeatgainloop = .true.
              elseif (incrementedgainlastcycle.AND.(T_current.lt.T_oldest)) then !== too far, go back in gain
                    alphaa = alphaa - 2*deltagain
                    deltagain = deltagain/10.0
                    decgs = decgs+1
                    incrementedgainlastcycle = .true.
                    repeatgainloop = .true.
              elseif ((.NOT.incrementedgainlastcycle).AND.(gli.eq.3).AND.(T_current.le.T_prev)) then !== not able to decrement gain because we have not taken enough steps
                 !== reset gain to zero, decrease step size, start over
                 deltagain = deltagain/10.0
                 decgs = decgs+1
                 alphaa = 0
                 T_prev = 0
                 repeatgainloop = .true.
              elseif ((.NOT.incrementedgainlastcycle).AND.(gli.ge.3).AND.(T_current.ge.T_prev)) then !== keep incrementing gain 
                 T_oldest = T_old
                 T_old    = T_prev
                 T_prev   = T_current
                 alphaa = alphaa + deltagain
                 incrementedgainlastcycle = .false.
                 repeatgainloop = .true.
              elseif ((.NOT.incrementedgainlastcycle).AND.(gli.gt.3).AND.(T_current.lt.T_prev)) then !== too far
                 write(111,*) 'too far: T_current= ',T_current,' < Tprev= ',T_prev
                 criticalgain = alphaa - deltagain
                 alphaa = alphaa - 2*deltagain
                 deltagain = deltagain/10.0
                 decgs = decgs+1
                 incrementedgainlastcycle = .true.
                 repeatgainloop = .true.
                 if (decgs.eq.n_decgs) then !== decremented gain step maximum number of times
                    write(111,*) 'TIME TO GO TO NEXT RLZ (max decgs). gli,fli,rlz=',gli,fli,rlz
                    totaltime = etime(elapsed)
                    write(111,*) 'run time so far:', totaltime/60.0,' minutes'
                    write(148,*) rlz,fli,gli, criticalgain, decgs, 1 !'critical gain'
                    gli = n_alphaa !== to get J,W to write to file, tell them that gain is done
                    call distributionscriticalgain(criticalgain,gainbin,binsize)
                    call flush(148)
                    incrementedgainlastcycle = .false.
                    repeatgainloop = .false.
                 endif

              else !== gli ~ 3
                 write(111,*) 'ERROR in gain loop 1'
              endif
           else !== find critical gain
              write(111,*) 'ERROR in gain loop 2'
           endif
           !== record the gains to an array
           if (.NOT.readGain) then
             gain(gli)=alphaa
           endif

           !== zeroize the matrices and arrays (while within inner-most loop)
           Tm(:,:)    = (0.0d+0,0.0d+0)
           DetFree(:) = (0.0d+0,0.0d+0)
           DetScat(:) = (0.0d+0,0.0d+0)
           DetChunk(:)= (0.0d+0,0.0d+0) 
           Ref  (:,:) = (0.0d+0,0.0d+0)   !== Zero the matrices and vectors               
           Trans(:,:) = (0.0d+0,0.0d+0)
           Vr   (:  ) = (0.0d+0,0.0d+0)
           Vt   (:  ) = (0.0d+0,0.0d+0)

           Zmt = (Gm + Hm) !== Zm_temporary,Tm will be the result of combining the chunks from self-embedding
           call inverse_c_matrix(Zmt,ntot,Tm)   !== ntot is the size of Zmt 
           !== returns the inversed Zmt matrix which is Tm
           Rm = Tm !== initializes first Reflection matrix to be first Transmission matrix
!           if (writeScatFreeMatrix) then
!              write(131,*) 'detc(Tm,ntot)                       = ',determinantc(Tm,ntot)
!              write(131,*) '====== Loop over self-embedding intervals ============='
!           endif
           if (writeTmRmDet) then  !== note: "cdabs(A)" is absolute value, SQRT(dreal(A)**2+dimag(A)**2)
              write(124,*) rlz,0,&
                 cdabs(determinantc(Tm,ntot)),cdabs(determinantc(Tm(1:n_open,1:n_open),n_open)),&
                 cdabs(determinantc(Rm,ntot)),cdabs(determinantc(Rm(1:n_open,1:n_open),n_open))
           endif
           !====== Loop over self-embedding intervals =============   
           do z=1,M_se
              ChunkMatrix = Im
              !====== Finding transfer matrix within se-interval by simple multiplication
              do a=(z-1)*se_step+1,z*se_step            
                 !=== Generate transfer matrix for free space between scatterers
                 if (a.eq.1) then
                    DeltaX = Medium(1,1)
                 else
                    DeltaX = Medium(a,1) - Medium(a-1,1)
                 endif

                 call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)

                 !== call the "determinant" function (with two arguments passed) and set the result equal to a vector
                 if (writeDeterminats) then
                    DetFree(a)  = determinant(FreeMatrix,ntot)
                 endif
                 ChunkMatrix      = matmul(FreeMatrix,ChunkMatrix)
                 if (writeDeterminats) then
                    DetChunk(a) = determinantc(ChunkMatrix,ntot)
                 endif
           
                 if (writeScatFreeMatrix) then
!                    write(138,*) 'free       matrix dx=',DeltaX,' z=',z,' a=',a
!                    write(138,*) 'fm=[...'
                    do i=1,ntot  !== formatting: 20=columns, ES=scientific notation, 20=characters, 13=right of decimal
                                 !== see page 181 of "Fortran 90" by Chapman
                       !write(138,'(20ES21.13)') real(FreeMatrix(i,:)) !,';...'
                       write(138,*) FreeMatrix(i,:)
                    enddo
!                    write(138,*) ']'
                 endif
                 !begin Scattering Matrix Subroutine    
                 ScatteringMatrix = Im
                 if ((fli.eq.1).AND.(gli.eq.1)) then
                    firstpropagation=.true.
                 else
                    firstpropagation=.false.
                 endif
                 call scatmatrixcreation(ScatteringMatrix,varyStrength,seed,&
                                         alphasign,seed_counter,nmax,ntot,alphas,&
                                         alphaa,omg,W,Kpara,Kperp,Medium,M,a,firstpropagation)
                 if (writeScatFreeMatrix) then
!                    write(131,*) 'scattering matrix dx=',DeltaX,' z=',z,' a=',a   
!                    write(131,*) 'sm=[...'
                    do i=1,ntot
                       !write(131,'(40ES21.13)') real(ScatteringMatrix(i,:)),imag(ScatteringMatrix(i,:))
                       write(131,*) dreal(ScatteringMatrix(i,:)),dimag(ScatteringMatrix(i,:))
                    enddo
!                    write(131,*) ']'
                 endif
                 if (writeDeterminats) then
                    DetScat(a+M+1)  = determinantc(ScatteringMatrix,ntot)
                 endif
                 ChunkMatrix          = matmul(ScatteringMatrix,ChunkMatrix)
                 if (writeDeterminats) then
                    DetChunk(a+M+1) = determinantc(ChunkMatrix,ntot)
                 endif

              enddo     !== end of the "a" (chunk) indexed loop

              if (writeScatFreeMatrix) then
!                 write(139,*) 'chunk      matrix dx=',DeltaX,' z=',z,' a=',a
!                 write(139,*) 'chunk=[...'
                 do i=1,ntot
                    write(139,*) dreal(ChunkMatrix(i,:)),dimag(ChunkMatrix(i,:))
                 enddo
!                 write(139,*) ']'
              endif
              !==== Self-embedding procedure
!== NEW: [Dr Yamilov, 20090705]
             Zmt=Im - matmul(Tm,matmul(Hm,(Im-ChunkMatrix)))   !====Sigma(J) in notes
             call inverse_c_matrix(Zmt,ntot,Zm)                !====returns Zm as Sigma(J) from notes
             Tempse=Im+matmul(Hm,matmul(Im-ChunkMatrix,matmul(Zm,Tm))) 
             Rm = matmul(Rm,Tempse)
             Tm = matmul(ChunkMatrix,matmul(Zm,Tm))
             if (doefse) then
                do zz=1,z
                    if (zz.eq.z) then    !== at the chunk boundary
                         SnN(:,:,zz)=Tm             !==== Initialize S(n,N) as S(n,n)
                    else
                         SnN(:,:,zz)=matmul(SnN(:,:,zz),Tempse)     !==== Update S(n,N+1) based on S(n,N)
                    endif
                enddo
             endif
!== OLD:
!              Zmt=Im - matmul(Tm,matmul(Hm,(Im-ChunkMatrix)))   !====Sigma(J) in notes
!              call inverse_c_matrix(Zmt,ntot,Zm)                 !====returns Zm as Sigma(J) from notes
!              Rm = matmul(Rm,Im+matmul(Hm,matmul(Im-ChunkMatrix,matmul(Zm,Tm))))
!              Tm = matmul(ChunkMatrix,matmul(Zm,Tm))

              if (writeTmRmDet) then
                 write(124,*) rlz,a,&
                      cdabs(determinantc(Tm,ntot)),cdabs(determinantc(Tm(1:n_open,1:n_open),n_open)),&
                      cdabs(determinantc(Rm,ntot)),cdabs(determinantc(Rm(1:n_open,1:n_open),n_open))
              endif
              if (writeRunningStatus) then
                 write(127,*) 'self embedding loop step number (z) = ',z
                 write(127,*) 'free matrix number (a)              = ',a
                 write(127,*) 'det(ChunkMatrix,ntot)-1             = ',determinantc(ChunkMatrix,ntot) - (1.0d+0)
              endif
              if ((abs(determinantc(ChunkMatrix,ntot) - (1.0d+0))).gt.(1E-10)) then
                 deviation_count = deviation_count+1
                 if (writeAccuracy) then
                    write(126,*) 'deviation count is currently', deviation_count
                    write(126,*) ' '
                 endif
                 if (deviation_count.ge.(.05*n_rlz*M/se_step)) then
                    open(133,file='out_ERROR.log',POSITION='APPEND')
                    write(133,*) 'WARNING: number of deviations exceeds 5% of the Self-Embedding Steps'
                    close(133)
                    write(125,*) 'WARNING: number of deviations exceeds 5% of the Self-Embedding Steps'
                    close(125)
                    write(*,*) 'WARNING: number of deviations exceeds 5% of the Self-Embedding Steps'
                    !if (rlz.lt.(int(.5*n_rlz))) then
!*                    call mpi_finalize(ierr)
                    !call mpi_abort(MPI_COMM_WORLD,ierr)
!*                    call closeAllFiles(1)
!*                       stop ! die
!*                       write(*,*) char(7)
!*                       call sleep(1)
!*                       write(*,*) char(7)
                    !endif
                 endif
              endif
              !write(126,*) 'detc(Tm,ntot)                       = ',determinantc(Tm,ntot)
           enddo                !== end the "z" indexed loop for self-embedding

           !==== tack on one more free matrix ===============
           DeltaX = L - medium(M,1)
           call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)

           if (writeDeterminats) then
              DetFree(M+1)  = determinant(FreeMatrix,ntot)
           endif
           ChunkMatrix        = FreeMatrix        !== new chunk matrix for the last free space   
           if (writeDeterminats) then
              DetChunk(M+1) = determinantc(ChunkMatrix,ntot)
           endif

           if (writeScatFreeMatrix) then
!              write(138,*) 'free       matrix dx=',DeltaX,' z=',z,' a=',a
!              write(138,*) 'fm=[...'
              do i=1,ntot
                 !write(138,'(20ES21.13)') real(FreeMatrix(i,:)) !,';...'
                 write(138,*) FreeMatrix(i,:) !,';...'
              enddo
!              write(138,*) ']'
           endif

!== Dr Yamilov's correction, 20090707
           Zmt=Im - matmul(Tm,matmul(Hm,(Im-ChunkMatrix)))   !====Sigma(J) in notes
           call inverse_c_matrix(Zmt,ntot,Zm)                !====returns Zm as Sigma(J) from notes
           Tempse=Im+matmul(Hm,matmul(Im-ChunkMatrix,matmul(Zm,Tm))) 
           Rm = matmul(Rm,Tempse)
           Tm = matmul(ChunkMatrix,matmul(Zm,Tm))
           if (doefse) then
              do zz=1,M_se
                  SnN(:,:,zz)=matmul(SnN(:,:,zz),Tempse)        !==== Update S(n,N+1) based on S(n,N)
              enddo
           endif

!== OLD:
           !Zmt = Im - matmul(Tm,matmul(Hm,(Im-ChunkMatrix)))
           !call inverse_c_matrix(Zmt,ntot,Zm)
           !Rm = matmul(Rm,Im+matmul(Hm,matmul(Im-ChunkMatrix,matmul(Zm,Tm))))                 !==  S(0,J+1)
           !Tm = matmul(ChunkMatrix,matmul(Zm,Tm)) !== now Tm includes the last free space     !==  S(J+1,J+1)

           if (writeTmRmDet) then
              write(124,*) rlz,a,&
                      cdabs(determinantc(Tm,ntot)),cdabs(determinantc(Tm(1:n_open,1:n_open),n_open)),&
                      cdabs(determinantc(Rm,ntot)),cdabs(determinantc(Rm(1:n_open,1:n_open),n_open))
           endif
           if (writeRunningStatus) then
              !write(127,*) 'dx                                  = ',DeltaX
              write(127,*) 'self embedding loop step number (z) = ',z
              write(127,*) 'free matrix number (a)              = ',a
              write(127,*) 'det(ChunkMatrix,ntot)-1             = ',determinantc(ChunkMatrix,ntot) - (1.0d+0)
              !write(127,*) 'detc(Tm,ntot)                       = ',determinantc(Tm,ntot)
              write(127,*) 'rlz = ', rlz,'freq = ',fli,' gain = ',gli
              write(127,*) ' '   
           endif
           !=== At the opposite side of the incident beam
    
           !===== Saving Determinants of free and scattering matrices ==============
           if (writeDeterminats) then
              do b=1,M+1
                 if (b.eq.1) then
                    DeltaX = Medium(1,1)
                 else if (b.eq.M+1) then 
                    DeltaX = L - Medium(M,1)
                 else
                    DeltaX = Medium(b,1) - Medium(b-1,1)
                 endif
                 write (101,*) 'rlz=',rlz,'Freespace index=',b, &
                    'Det(Free)-1=',DetFree(b)    -(1.0d+0,0.0d+0),' ',&
                    'Det(Scat)-1=',DetScat(b+M+1)-(1.0d+0,0.0d+0),' ',&
                    'dx=',DeltaX
                 write (102,*) 'rlz',rlz,'Chunkmatrix index=',b,&
                    'Det(ChunkMatrix)-1 after Free matrix=',DetChunk(b)    -(1.0d+0,0.0d+0),&
                    'det(ChunkMatrix)-1 after Scat matrix=',DetChunk(b+M+1)-(1.0d+0,0.0d+0),&
                    'dx=',DeltaX
              enddo   
           endif

           !===== Saving self-embedding reflection and transmission matrices =======
           !write(125,*) 'FINAL: detc(Tm,ntot)         = ',cdabs(determinantc(Tm,ntot))
           if (writeSelfEmbedTmRm) then
              do a=1, ntot
                 write(103,*) dreal(Tm(a,:)),dimag(Tm(a,:))
              enddo
              !write(125,*) 'FINAL: detc(Rm,ntot)         = ',cdabs(determinantc(Rm,ntot))
              do a=1, ntot
                 write(104,*) dreal(Rm(a,:)),dimag(Rm(a,:))
              enddo
           endif
           !===== Loop to determine Trans and Ref (1:ntot,1:n_open) matrices =======
           conductvty =0 !== is conductivity not defined in active media?
           oneminusconductvty =0
           do a=1, n_open          !== Loop to make final Ref and Trans
              !== create V_{boundary condition} vector
              Vbc(:) = 0.0d+0
              Vbc(a) = 2.0d+0      !== INPUT CHANNEL. MULTIPLY  Initial Boundary 
                !== Condition(also alters it every step to make the next boundary cond.)
              Vt (:) = 0.0d+0
              Vr (:) = 0.0d+0
        
              Vt = matmul(Tm,Vbc)!== Makes a column of Trans (transmission) matrix
              do b=1, nmax         !== Loop to copy the column into Trans matrix
                 Trans(b,a) = Vt(b)*(Kpara(b)/Kpara(a))**0.5
                 conductvty = conductvty + cdabs(Trans(b,a))**2
              enddo
!              if (a.eq.2) then
!                 write(*,*) ' '
!                 write(*,*) 'rlz = ', rlz
!                 write(*,*) 'vt(1), from self-embed, input chan2 = ',dreal(Vt(1)),dimag(Vt(1))
!              endif
        
              Vr = matmul(Rm,Vbc)   !== Makes a column of Ref (reflection) matrix
              Vr(a) = Vr(a)-1      !== subtracts 1 from the current open channel 
                        !== step(makes it actually the Ref matrix column)
              do b=1, nmax         !== Loop to copy the column into Ref matrix
                 Ref  (b,a) = Vr(b)*(Kpara(b)/Kpara(a))**0.5
                 oneminusconductvty = oneminusconductvty + cdabs(Ref(b,a))**2
              enddo
!              if (a.eq.2) then
!                 write(*,*) 'rlz = ', rlz
!                 write(*,*) 'vr(1), from self-embed, input chan2 = ',dreal(Vr(1)),dimag(Vr(1))
!              endif

           enddo  !== a
           if ((sum((cdabs(Trans))**2).gt.1000).AND.(.NOT.findcriticalgain).AND.(alphaaf.lt.0)) then
              repeatgainloop = .false.  !== very close to critical gain, or exceeding it
!~~              write(163,*) gli,fli,rlz,'T>1000'
           endif
           if ((sum((cdabs(Trans))**2).lt.T_prev).AND.(.NOT.findcriticalgain).AND.(alphaaf.lt.0)) then
              repeatgainloop = .false.  !== very close to critical gain, or exceeding it
!~~              write(163,*) gli,fli,rlz,'T_current < T_prev'
           endif

          if (findeignvalues.AND.repeatgainloop) then
            != eignvalues() is the target vector for storing eigen values
            tthc = matmul(Trans,transpose(conjg(Trans)))
            call zheev('N','U',nmax,tthc,nmax,eignvalues(1:nmax),work,2*nmax-1,rwork,c)
            write(166,"("//int2String(nmax+3)//"ES24.15)") eignvalues(1:nmax),real(gli),real(fli),real(rlz)
          endif

           if (findElecField) then
              writecriticalgainduetovfwderror = .true.
              do b = 1,n_open !== input channel
                 !== after the last free space, we know t, r, S(N,N), and S(0,N)
                 Vbc(:) = 0.0d+0
                 Vbc(b) = 2.0d+0
    
 !               v_bck(:,M+2) = matmul(Tm,Vbc)!== Makes a column of Trans (transmission) matrix
                 v_fwd(:,1) = matmul(Rm,Vbc) !== tracking indicies 1 and M+2 are Vr and Vt
                  
                 if (repeatgainloop) then
                    !== forward propagation   
    !== NEW: [Dr Yamilov, 20090705]
                    do sli = 1,M
                      if (sli.eq.1) then
                          DeltaX = Medium(1,1)
                      else
                          DeltaX = Medium(sli,1) - Medium(sli-1,1)
                      endif
                      if (doefse) then
                          if ((mod(sli,se_step).eq.0).and.(sli.NE.M)) then !== sli is at a chunk boundary, but not the end
                            call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)
                            call scatmatrixcreation(ScatteringMatrix,varyStrength,seed,&
                                                    alphasign,seed_counter,nmax,ntot,alphas,&
                                                    alphaa,omg,W,Kpara,Kperp,Medium,M,sli,.false.)
                            v_fwd_check = matmul(ScatteringMatrix,matmul(FreeMatrix,v_fwd(:,sli))) 
                            v_fwd(:,sli+1) = matmul(SnN(:,:,sli/se_step),Vbc)
                            v_fwd_error=sum(cdabs(v_fwd_check-v_fwd(:,sli+1)))                     
                            if (v_fwd_error.gt.(1E-10).AND.(gli.eq.1)) then
                            open(133,file='out_ERROR.log',POSITION='APPEND')
                                write(133,*) 'WARNING: electric field error exceeds 1E-10, yields invalid A,B',gli,fli,rlz,b
                                close(133)
                                write(125,*) 'WARNING: electric field error exceeds 1E-10, yields invalid A,B',gli,fli,rlz,b
                                close(125)
                                write(*,*) 'WARNING: electric field error exceeds 1E-10, yields invalid A,B',gli,fli,rlz,b
                                write(*,*) 'diff is ',v_fwd_error
    !*                            call closeAllFiles(1)
                                !if (rlz.lt.(int(.5*n_rlz))) then
    !*                            call mpi_finalize(ierr)
                                !call mpi_abort(MPI_COMM_WORLD,ierr)
    !*                               stop ! die
    !*                              write(*,*) char(7) != beep
    !*                               call sleep(1)
    !*                               write(*,*) char(7)
                                !endif
                            elseif (v_fwd_error.gt.(1E-10)) then !== electric field error for active media
                                !== this realization errored out  
                                repeatgainloop = .false.
!~~                                write(163,*) gli,fli,rlz,'v_fwd error',b
                                if (writecriticalgainduetovfwderror.AND.findcriticalgain) then
                                  write(148,*) rlz,fli,gli, alphaa,decgs, 3 !'v_fwd error'
                                  call distributionscriticalgain(alphaa,gainbin,binsize)
                                  call flush(148)
                                  write(111,*) 'v_fwd error',rlz, 0, gli, decgs,b
                                  gli = n_alphaa !== to get J,W to write to file, tell them that gain is done
                                  writecriticalgainduetovfwderror = .false. !== only write to file once per input channel
                                elseif (.NOT.findcriticalgain) then
                                  !write(*,*) 'entered non-physical regime v_fwd error',rlz,fli,gli,b
                                  repeatgainloop = .false.
                                endif
                            endif
                          else !== between chunk boundaries, propagate forward as normal
                            call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)
                            call scatmatrixcreation(ScatteringMatrix,varyStrength,seed,&
                                                    alphasign,seed_counter,nmax,ntot,alphas,&
                                                    alphaa,omg,W,Kpara,Kperp,Medium,M,sli,.false.)
                            v_fwd(:,sli+1) = matmul(ScatteringMatrix,matmul(FreeMatrix,v_fwd(:,sli)))
                          endif
                      else
                          call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)
                          call scatmatrixcreation(ScatteringMatrix,varyStrength,seed,&
                                                  alphasign,seed_counter,nmax,ntot,alphas,&
                                                  alphaa,omg,W,Kpara,Kperp,Medium,M,sli,.false.)
                          v_fwd(:,sli+1) = matmul(ScatteringMatrix,matmul(FreeMatrix,v_fwd(:,sli)))
                      endif
                    enddo !== sli forwards, to the last scatterer
   
                    DeltaX = L - Medium(M,1)  !== tack on last free space to find Vt
                    call freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)
                    v_fwd(:,M+2) = matmul(FreeMatrix,v_fwd(:,M+1))
                 endif

                 if (findJWDplanewave.AND.(mod(b,2).eq.1).AND.repeatgainloop) then
!   vector = v_fwd
!   numscatterersplustwo=M+2
!   inputchan=b
!   k_jpw = sqrt(k2)
!   nz_jpw = int(L*ABstep)+1
! 
!   !== see notes, 20090827
!   pwchannorm_jpw = sqrt(8.0)/(k_jpw*pi*inputchan) !== plane wave input channel normalization
!   !write(*,*) 'inputchan=',inputchan,' k=',k,' pi=',pi,'pwchannorm = ',pwchannorm
! 
!   A_raw_jpw(1:nmax,1:numscatterersplustwo) = (0.5)*( vector(1:nmax,1:numscatterersplustwo)-&
!                       (0.0d+0,1.0d+0)*vector(nmax+1:2*nmax,1:numscatterersplustwo) )*(pwchannorm_jpw*(1.0d+0,0.0d+0))
!   B_raw_jpw(1:nmax,1:numscatterersplustwo) = (0.5)*( vector(1:nmax,1:numscatterersplustwo)+&
!                       (0.0d+0,1.0d+0)*vector(nmax+1:2*nmax,1:numscatterersplustwo) )*(pwchannorm_jpw*(1.0d+0,0.0d+0))
!   !== NOTE: including the (1,0) for pwchan norm does not affect A,B_raw. Both are still complex
!   
!   if (inputchan.eq.1) then !== first call of this subroutine, reset matrices
!      Apwscatsum(1:nmax,1:numscatterersplustwo) = (0.0d+0,0.0d+0)
!      Bpwscatsum(1:nmax,1:numscatterersplustwo) = (0.0d+0,0.0d+0)
!   endif
! 
!   !== summing over input channels 
!   Apwscatsum = Apwscatsum + A_raw_jpw
!   Bpwscatsum = Bpwscatsum + B_raw_jpw
!   if ((writepwJWDchan.OR.writegEpw.OR.writeconductionpwdist).AND.((inputchan.eq.n_open).OR.(inputchan.eq.n_open-1))&
!          .AND.(pi.lt.(5.0))) then
!      !== J,W pw are recorded after summing over input channels
!      !== begin finding distribution of g,E
!      if (writeconductionpwdist.OR.writegEpw) then
!         Apwscat2_jpw = cdabs(Apwscatsum)**2
!         Bpwscat2_jpw = cdabs(Bpwscatsum)**2
! 
!          do c_jpw = 1,(numscatterersplustwo-1) !== for every scatterer
!             !== note: there is no energyDensity associated with c=numscatterersplustwo
! 
!             integrand_jpw = k2*sum(Apwscat2_jpw(1:nmax,c_jpw)+Bpwscat2_jpw(1:nmax,c_jpw))
!             if (c_jpw.eq.1) then
!                energypw_jpw=0   
!                DeltaX_jpw = Medium(1,1)
!             elseif (c_jpw.lt.numscatterersplustwo-1) then
!                DeltaX_jpw = Medium(c_jpw,1) - Medium(c_jpw-1,1)
!               if (DeltaX_jpw.lt.0) then
!                   write(133,*) 'negative deltax',rlz,gli,fli,DeltaX_jpw,c_jpw
!               endif
!             elseif (c_jpw.eq.numscatterersplustwo-1) then
!                DeltaX_jpw = L-Medium(M,1)
!               if (DeltaX_jpw.lt.0) then
!                   write(133,*) 'negative deltax',rlz,gli,fli,DeltaX_jpw,c_jpw
!               endif
!             endif  
!             energypw_jpw = energypw_jpw + integrand_jpw*W*(DeltaX_jpw*(1.0d+0))
!             energy2pw_jpw = energy2pw_jpw + (integrand_jpw**2)*W*(DeltaX_jpw*(1.0d+0))
!          enddo 
!          !== note: energy2pw is for the inverse participation ratio criterion distribution
!  
!          !== NOTE: T,E normalizations are left over from uf calculations! [probably need to be changed for pw]
!          condpw_jpw = k_jpw*sum(Apwscat2_jpw(1:n_open,numscatterersplustwo)*kparaMat(1:n_open,1))
! 
!          !check: 
!          if (sum(Bpwscat2_jpw(1:n_open,numscatterersplustwo)).gt.(1E-10)) then
!             open(133,file='out_ERROR.log',POSITION='APPEND')
!             write(133,*) 'WARNING: incident flux on right side !?'
!             close(133)
!             write(125,*) 'WARNING: incident flux on right side !?'
!             close(125)
!             write(*,*) 'WARNING: incident flux on right side !?'
!            !if (rlz.lt.(int(.5*n_rlz))) then
! !*            call closeAllFiles(1)
! !*            stop ! die
!            !endif
!          endif
!  
!          !== reality check
!          if ((condpw_jpw.lt.0).OR.(energypw_jpw.lt.0)) then 
!             open(133,file='out_ERROR.log',POSITION='APPEND')
!             write(133,*) 'WARNING: energy or transmission less then 0',condpw_jpw,energypw_jpw
!             close(133)
!             write(125,*) 'WARNING: energy or transmission less then 0',condpw_jpw,energypw_jpw
!             close(125)
!             write(*,*) 'WARNING: energy or transmission less then 0',condpw_jpw,energypw_jpw
!            !if (rlz.lt.(int(.5*n_rlz))) then
!  !*           call closeAllFiles(1)
!  !*          stop ! die
!            !endif
!          endif
! 
!         !== normalize the energy per channel by empty waveguide
!            !== see email, 20090913
!            !== energy in empty waveguide = (L*W) * sum_{a=1:nopen} (pwchannorm_a)^2
!         emptypwwg_jpw = (0.0d+0)
!         do a_jpw = 1,n_open,2
!            emptypwwg_jpw = emptypwwg_jpw + (1.0d+0)/(a_jpw*a_jpw*(1.0d+0))
!         enddo
!         energypw_jpw = energypw_jpw*pi*pi/(L*W*8*emptypwwg_jpw)  
!         energy2pw_jpw = energy2pw_jpw*(pi**4)/(L*W*64*emptypwwg_jpw*emptypwwg_jpw)  !== ?? emptypwwg^2, or sum(a^4)?
!         if (writegEpw) then !.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
!            write(141,'(5ES24.15)') condpw_jpw,energypw_jpw,energy2pw_jpw,real(gli),real(fli)  !== per input channel
!            !== reminder: may need to eliminate "fli" to avoid wrapping columns
!         endif
! !        if (alphaa.eq.0) then   !== for the distributions, no longer in use
! !          energyPassivepw = energypw_jpw
! !        endif
!      endif !== writeconductionpwdist. Done with distributions
! 
!      !== given A, B after each scatterer, convert to static step size
!      A_inside_jpw(1:nmax,1:nz_jpw)=&
!        Apwscatsum(1:nmax,scatr_indx(1:nz_jpw))*expon(1:nmax,1:nz_jpw)
!      B_inside_jpw(1:nmax,1:nz_jpw)=&
!        Bpwscatsum(1:nmax,scatr_indx(1:nz_jpw))*(1/expon(1:nmax,1:nz_jpw))  !== this should be element-wise inversion. Tested & works
! 
!      if (find_dzy) then
! 
!         !== find u, du/dy, du/dz. These are the components for J,W
!         do a=1,int(W*ABstep)+1
!            y(a) = (a-1)*W/(int(W*ABstep)*(1.0d+0))
!         enddo
!         do a=1,nmax
!            chi_y(a,1:int(W*ABstep)+1)=sqrt(2/W)*sin(a*pi*y(:)/W)
!            dchi_y(a,1:int(W*ABstep)+1)=sqrt(2/W)*(a*pi/W)*cos(a*pi*y(:)/W)
!         enddo
! 
!         AB_inside_trans(:,:) = (0.0d+0,0.0d+0)
!         AB_inside_trans(1:nz_jpw,1:nmax) = transpose(A_inside_jpw+B_inside_jpw)
!         uzy(1:nz_jpw,1:int(W*ABstep)+1) = matmul(AB_inside_trans,chi_y(1:nmax,1:int(W*ABstep)+1))
! 
!         uzy2(1:nz_jpw,1:int(W*ABstep)+1) = &
!         matmul(abs(AB_inside_trans*AB_inside_trans),(chi_y(1:nmax,1:int(W*ABstep)+1)*chi_y(1:nmax,1:int(W*ABstep)+1)))
! 
!         partialuwrty(1:nz_jpw,1:int(W*ABstep)+1) = matmul(AB_inside_trans,dchi_y)
!         partialuwrty2(1:nz_jpw,1:int(W*ABstep)+1) = matmul(abs(AB_inside_trans*AB_inside_trans),(dchi_y*dchi_y))
!         AB_inside_trans(:,:) = (0.0d+0,0.0d+0)
!         do a=1,nmax 
!             AB_inside_trans(1:nz_jpw,a) = kpara(a)*(A_inside_jpw(a,1:nz_jpw)-B_inside_jpw(a,1:nz_jpw))
!         enddo
!         partialuwrtz(1:nz_jpw,1:int(W*ABstep)+1) = (0.0d+0,1.0d+0)*matmul(AB_inside_trans,chi_y(1:nmax,1:int(W*ABstep)+1))
!         partialuwrtz2(1:nz_jpw,1:int(W*ABstep)+1) = &
!      matmul(abs(AB_inside_trans*AB_inside_trans),(chi_y(1:nmax,1:int(W*ABstep)+1)*chi_y(1:nmax,1:int(W*ABstep)+1)))
! 
!         if (.NOT.dzy_gain_resolved) then
!           if (rlz.eq.1) then !== reset summed j,w
!             wdist_zy(:,:) = (0.0d+0)
!             wdist_zy_incohrnt(:,:) = (0.0d+0)
!             jflux_y(:,:)=(0.0d+0)
!             jflux_z(:,:)=(0.0d+0)
!           endif
! 
!           wdist_zy(1:nz_jpw,1:int(W*ABstep)+1) = wdist_zy(:,:)+.5*(k_jpw*k_jpw*abs(uzy)*abs(uzy)+&
!               abs(partialuwrty)*abs(partialuwrty)+abs(partialuwrtz)*abs(partialuwrtz))
!           wdist_zy_incohrnt(1:nz_jpw,1:int(W*ABstep)+1) = wdist_zy(:,:)+.5*(k_jpw*k_jpw*uzy2 + partialuwrty2+partialuwrtz2)
!           jflux_y(1:nz_jpw,1:int(W*ABstep)+1) = jflux_y(:,:)+k_jpw*dimag(conjg(uzy)*partialuwrty)
!           jflux_z(1:nz_jpw,1:int(W*ABstep)+1) = jflux_z(:,:)+k_jpw*dimag(conjg(uzy)*partialuwrtz)
!         elseif (dzy_gain_resolved) then
!           if (rlz.eq.1) then !== reset summed j,w
!             wdist_zy_gr(:,:,gli) = (0.0d+0)
!             wdist_zy_gr_incohrnt(:,:,gli) = (0.0d+0)
!             jflux_y_gr(:,:,gli)=(0.0d+0)
!             jflux_z_gr(:,:,gli)=(0.0d+0)
!           endif
! 
!           wdist_zy_gr(1:nz_jpw,1:int(W*ABstep)+1,gli) = wdist_zy_gr(:,:,gli)+.5*(k_jpw*k_jpw*abs(uzy)*abs(uzy)+&
!               abs(partialuwrty)*abs(partialuwrty)+abs(partialuwrtz)*abs(partialuwrtz))
!           wdist_zy_gr_incohrnt(1:nz_jpw,1:int(W*ABstep)+1,gli) = wdist_zy_gr(:,:,gli)+&
!                 .5*(k_jpw*k_jpw*uzy2+partialuwrty2+partialuwrtz2)
!           jflux_y_gr(1:nz_jpw,1:int(W*ABstep)+1,gli) = jflux_y_gr(:,:,gli)+k_jpw*dimag(conjg(uzy)*partialuwrty)
!           jflux_z_gr(1:nz_jpw,1:int(W*ABstep)+1,gli) = jflux_z_gr(:,:,gli)+k_jpw*dimag(conjg(uzy)*partialuwrtz)
!         endif
! 
!      endif !== find_dzy
! 
!      Apw2_jpw=cdabs(A_inside_jpw)**2
!      Bpw2_jpw=cdabs(B_inside_jpw)**2
! 
!      do a_jpw=1,n_open  !== interior channels
!         Jflux_pw_r_n_sum(a_jpw,1:nz_jpw,gli)=Jflux_pw_r_n_sum(a_jpw,1:nz_jpw,gli)+k_jpw*Apw2_jpw(a_jpw,1:nz_jpw)*kparaMat(a_jpw,1)
!         Jflux_pw_l_n_sum(a_jpw,1:nz_jpw,gli)=Jflux_pw_l_n_sum(a_jpw,1:nz_jpw,gli)+k_jpw*Bpw2_jpw(a_jpw,1:nz_jpw)*kparaMat(a_jpw,1)
!         !write(*,*) 'Jflux_pw_l_n_sum(5,5,1,1) = ',Jflux_pw_l_n_sum(5,5,1,1)
!       enddo
!       energy_integrand_jpw(1:nz_jpw) = (0.0d+0)
!       do a_jpw=1,nmax  !== interior channels
!          Wdist_pw_r_n_sum(a_jpw,1:nz_jpw,gli)=Wdist_pw_r_n_sum(a_jpw,1:nz_jpw,gli)+k2*Apw2_jpw(a_jpw,1:nz_jpw)
!          Wdist_pw_l_n_sum(a_jpw,1:nz_jpw,gli)=Wdist_pw_l_n_sum(a_jpw,1:nz_jpw,gli)+k2*Bpw2_jpw(a_jpw,1:nz_jpw)
!          energy_integrand_jpw(1:nz_jpw) = energy_integrand_jpw(1:nz_jpw) +k2* (Apw2_jpw(a_jpw,1:nz_jpw)+Bpw2_jpw(a_jpw,1:nz_jpw))
!       enddo
!       do a_jpw=1,(nz_jpw-1)
!          correlationE_sum(a_jpw,gli) = correlationE_sum(a_jpw,gli) + (sum( energy_integrand_jpw(1:(nz_jpw-a_jpw))*&
!                                   energy_integrand_jpw((a_jpw+1):nz_jpw) )/((nz_jpw-a_jpw)*(1.0d+0)))
!       enddo
! 
!       !== write J_+, J_-, W_+,W_- for each channel to file
!       if ((fli.eq.n_omg).AND.(gli.eq.n_alphaa).AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz))).AND.&
!          ((inputchan.eq.n_open).OR.(inputchan.eq.n_open-1))) then 
!          write(165,*) gli,fli,rlz,'wrote in subroutine'
!      
!          if (find_dzy.AND.(.NOT.dzy_gain_resolved)) then
!             do c_jpw=1,int(W*ABstep)+1
!                do a_jpw=1,nz_jpw
!      write(164,'(4ES24.15)') wdist_zy(a_jpw,c_jpw),wdist_zy_incohrnt(a_jpw,c_jpw),jflux_y(a_jpw,c_jpw),jflux_z(a_jpw,c_jpw)
!                enddo
!             enddo
!          elseif (find_dzy.AND.dzy_gain_resolved) then
!             do h_jpw=1,n_alphaa
!               do c_jpw=1,int(W*ABstep)+1
!                 do a_jpw=1,nz_jpw
!      write(164,'(4ES24.15)') wdist_zy_gr(a_jpw,c_jpw,h_jpw),wdist_zy_gr_incohrnt(a_jpw,c_jpw,h_jpw), &
!                    jflux_y_gr(a_jpw,c_jpw,h_jpw), jflux_z_gr(a_jpw,c_jpw,h_jpw)
!                 enddo
!               enddo
!             enddo
!          endif
! 
!          if (writepwJWDchan) then
!          do c_jpw = 1,n_alphaa
!             do h_jpw = 1,n_open  !== interior channel
!                do a_jpw = 1,nz_jpw
!                   write(156,*) Jflux_pw_l_n_sum(h_jpw,a_jpw,c_jpw),Jflux_pw_r_n_sum(h_jpw,a_jpw,c_jpw)  !== sum
!                enddo
!             enddo
!             do h_jpw = 1,nmax  !== interior channel
!                do a_jpw = 1,nz_jpw
!                   write(158,*) Wdist_pw_l_n_sum(h_jpw,a_jpw,c_jpw),Wdist_pw_r_n_sum(h_jpw,a_jpw,c_jpw) !== sum
! !!                     write(159,*) Wdist_pw_r_n_sum(h,a,c)/(n_rlz*(1.0d+0)),Wdist_pw_l_n_sum(h,a,c)/(n_rlz*(1.0d+0))   !== average
!                enddo  
!             enddo  
!             do a_jpw=1,(nz_jpw-1)
!                write(162,*) correlationE_sum(a_jpw,c_jpw)
!             enddo
!          enddo
!          endif
!       endif
!    endif
! 
!   call flush(156)
!   call flush(157)
!   call flush(158)
!   call flush(159)
!   call flush(165)

                 endif
                 if (findJWDunitflux.AND.repeatgainloop) then
                     call JWDunitflux(rank,ierr,v_fwd,M+2,ntot,nmax,kpara,b,pi,k2,W,&
                        scatr_indx,expon,L,ABstep,n_omg,fli,gli,n_open,writeufJWDchan,&
                        kparaMat,n_alphaa,alphaa,rlz,n_rlz,percntwriteJWD,M,Medium,&
                        writeEachJW,writeSummedJW,Jflux_r_sum,Jflux_l_sum,Wdist_r_sum,Wdist_l_sum,&
                        Wdist_r_n_sum,Wdist_l_n_sum,Jflux_r_n_sum,Jflux_l_n_sum,&
                        tebinsuf,condbinsuf,binsize,writeTEuf,writeTEufdist,gainbin,findcriticalgain,&
                        conduf,total_energyuf,total_passive_energyuf,writeconductionufdist)
                 endif

                  if (writeGain.AND.repeatgainloop.AND.(.NOT.findcriticalgain).AND.(b.eq.n_open)) then
                      write(137,*) gain(gli),gli,fli,rlz
                      validgaincount(gli) = validgaincount(gli)+1
                      call flush(137)
                  endif

                  ! ========= Saving Reflection and Transmission matrices ================
                  if (writeReflectionMatrix.AND.repeatgainloop.AND.(b.eq.n_open)) then
!                      write(*,*) 'R',gli,fli,rlz,n_alphaa
                      do c=1, nmax
                        !write(121,*) dreal(Ref(b,:)),dimag(Ref(b,:))  !== note: gfortran likes this
                        !write(121,"(20ES)") real(Ref(b,:)),imag(Ref(b,:))  
                         write(121,"("//int2String(n_open*2)//"ES24.15)") real(Ref(c,:)),imag(Ref(c,:))  
                        !== note: NIC uses this way, otherwise it wraps @ 3 numbers
                        !== "ES"=type:real, scientific notation [page 223, Fortran90 by Chapman]
                      enddo
                  endif
!                   if (writeaveTab.AND.repeatgainloop.AND.(b.eq.n_open)) then
!                       do a=1,n_open
!                         do c=1, nmax      !== Loop to copy the column into Trans matrix
!                             aveT(c,a,gli) = aveT(c,a,gli) + cdabs(Trans(c,a))**2
!                             aveR(c,a,gli) = aveR(c,a,gli) + cdabs(Ref(c,a))**2
!                         enddo
!                       enddo
!                   endif
                  if (writeTransmissionMatrix.AND.repeatgainloop.AND.(b.eq.n_open)) then
                      do c=1, nmax
                        !write(122,*) dreal(Trans(b,:)),dimag(Trans(b,:))  !== note: gfortran likes this
                        !write(122,"(20ES)") real(Trans(b,:)),imag(Trans(b,:))
                          write(122,"("//int2String(n_open*2)//"ES24.15)") real(Trans(c,:)),imag(Trans(c,:))    
                        !== note: NIC uses this way, otherwise it wraps @ 3 numbers
                      enddo
                  endif

                 if ((b.eq.n_open).AND.(.NOT.findcriticalgain).AND.(.NOT.repeatgainloop).AND.(n_alphaa.gt.1)) then
                    Ref(:,:) = (0.0d+0,0.0d+0)
                    if (writeTransmissionMatrix) then
                      !== need to reuse the initialized R because T is used later
                      do a=1,(1+n_alphaa-gli)
                          do c=1,nmax
                            write(122,"("//int2String(n_open*2)//"ES24.15)") real(Ref(c,:)),imag(Ref(c,:))  
                          enddo
                      enddo
                    endif
                    if (writeReflectionMatrix) then
                      !== initialize Reflection matrix to zero
                      do a=1,(1+n_alphaa-gli)
!                         write(*,*) '0',gli,fli,rlz,n_alphaa,a+gli-1
                          do c=1,nmax
                            write(121,"("//int2String(n_open*2)//"ES24.15)") real(Ref(c,:)),imag(Ref(c,:))  
                          enddo
                      enddo
                    endif
                 endif
 
                 if ((fli.eq.n_omg).AND.(.NOT.repeatgainloop).AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))&
                    .AND.(b.eq.n_open).AND.(n_alphaa.gt.1)) then
                    write(165,*) 'entered write statement in main body'
                    !== reasoning: 
                    !== fli.eq.n_omg: end of frequency loop
                    !== .not.repeatgainloop: one of the "critical gain" conditions was triggered
                    !== mod: write to file every percntwriteJWD for crash recovery
                    !== b.eq.n_open: end of input channels
                    !== n_alphaa.gt.1: not doing a passive-only media

                    !== the following was removed, because it prevented files from being written at all
                    !== gli.ne.n_alphaa: NOT at end of gain loop. [if you were, then these files would be written in subroutines]

                    !========== UNIT FLUX OUTPUT, AVERAGED OVER FREQUENCY AND REALIZATIONS ===============

                    if (writeSummedJW) then  !== write J_+, J_-, W to file
                      d = n_rlz*n_open*n_omg !== input channel, realization
                      do c = 1,n_alphaa  
                         !Jflux(:)=Jflux_r(:,b,c)-Jflux_l(:,b,c)  !== do not write to file, save disk space and write time
                          do a = 1,(int(L*ABstep)+1)
                    write(144,'(4ES24.15)') Wdist_l_sum(a,c),Wdist_r_sum(a,c),Jflux_l_sum(a,c),Jflux_r_sum(a,c)   !== sum !== 20090901: Ben removed /(n_open*(1.0d+0))
!!                    write(145,*) Wdist_sum(a,c)/(d*(1.0d+0)),     Jflux_r_sum(a,c)/(d*(1.0d+0)),     Jflux_l_sum(a,c)/(d*(1.0d+0))  !== average
                          enddo  
                      enddo  
                    endif

                    if (writeufJWDchan) then !== write J_+, J_-, W_+,W_- for each channel to file
                      do c = 1,n_alphaa
                         do q = 1,n_open !== input
                             do p = 1,n_open  !== interior
                                do a = 1,(int(L*ABstep)+1)
                                   write(152,*) Jflux_l_n_sum(q,p,a,c),Jflux_r_n_sum(q,p,a,c)  !== sum
!!                  write(153,*) Jflux_r_n_sum(q,p,a,c)/(n_rlz*n_omg*(1.0d+0)),Jflux_l_n_sum(q,p,a,c)/(n_rlz*n_omg*(1.0d+0))  !== average
                                enddo  
                            enddo  
                            do d = 1,nmax !== interior
                                do a = 1,(int(L*ABstep)+1)
                                   write(154,*) Wdist_l_n_sum(q,d,a,c),Wdist_r_n_sum(q,d,a,c)  !== sum
!!                  write(155,*) Wdist_r_n_sum(q,d,a,c)/(n_rlz*n_omg*(1.0d+0)),Wdist_l_n_sum(q,d,a,c)/(n_rlz*n_omg*(1.0d+0))  !== average
                                enddo  
                            enddo  
                          enddo
                      enddo
                    endif
                    call flush(143)
                    call flush(144)
                    call flush(145)
                    call flush(152)
                    call flush(153)
                    call flush(154)
                    call flush(155)

                    !========== PLANE WAVE OUTPUT, AVERAGED OVER FREQUENCY AND REALIZATIONS ===============

                    if (find_dzy.AND.(.NOT.dzy_gain_resolved)) then
                        do c_jpw=1,int(W*ABstep)+1
                          do a_jpw=1,nz_jpw
    write(164,'(4ES24.15)') wdist_zy(a_jpw,c_jpw),wdist_zy_incohrnt(a_jpw,c_jpw),jflux_y(a_jpw,c_jpw), jflux_z(a_jpw,c_jpw)
                          enddo
                        enddo
                    elseif (find_dzy.AND.dzy_gain_resolved) then
                        do h_jpw=1,n_alphaa
                          do c_jpw=1,int(W*ABstep)+1
                            do a_jpw=1,nz_jpw
    write(164,'(4ES24.15)') wdist_zy_gr(a_jpw,c_jpw,h_jpw),wdist_zy_gr_incohrnt(a_jpw,c_jpw,h_jpw), &
                              jflux_y_gr(a_jpw,c_jpw,h_jpw), jflux_z_gr(a_jpw,c_jpw,h_jpw)
                            enddo
                          enddo
                        enddo
                    endif

                    if (writepwJWDchan) then
                    write(165,*) 'calling jwpw from main body'
                      do c = 1,n_alphaa
                         do d = 1,n_open  !== interior channel
                            do a = 1,(int(L*ABstep)+1)
                                write(156,*) Jflux_pw_l_n_sum(d,a,c),Jflux_pw_r_n_sum(d,a,c)  !== sum
                            enddo
                         enddo
                          do d = 1,nmax  !== interior channel
                            do a = 1,(int(L*ABstep)+1)
                                write(158,*) Wdist_pw_l_n_sum(d,a,c),Wdist_pw_r_n_sum(d,a,c) !== sum
                            enddo  
                         enddo  
                         do a = 1,(int(L*ABstep))
                            write(162,*) correlationE_sum(a,c)
                         enddo
                      enddo
                    endif
                    call flush(156)
                    call flush(157)
                    call flush(158)
                    call flush(159)
                 endif

                 !== error calculation
                 !write(*,*) 'vt(1), from elec field = ',dreal(v_fwd(1,M+2)),dimag(v_fwd(1,M+2))
                 !== removed un-necessary backwards propagation
                 Vr = matmul(Rm,Vbc)
                 Vt = matmul(Tm,Vbc)
                 if (writeElecFieldError) then
                    write(136,*) (sum((cdabs(v_fwd(:,M+2)-Vt))**2)/sum((cdabs(Vt))**2)) !, &
!                                (sum((cdabs(v_bck(:,1)  -Vr))**2)/sum((cdabs(Vr))**2))
                 endif

                 if ((writeABCoefficientsPW.OR.writeABcoefficientsUF) .AND. &
                     ((sum((cdabs(v_fwd(:,M+2)-Vt))**2)/sum((cdabs(Vt))**2)) .gt. 0.1)) then
                     open(133,file='out_ERROR.log',POSITION='APPEND')
                     write(133,*) 'WARNING: electric field error yields invalid A,B',gli,fli,rlz,b
                     !write(133,*) 'program exited, DO NOT USE output data'
                     close(133)
                     write(125,*) 'WARNING: electric field error yields invalid A,B',gli,fli,rlz,b
                     !write(125,*) 'program exited, DO NOT USE output data'
                     !close(125)
                     write(*,*) 'WARNING: electric field error yields invalid A,B',gli,fli,rlz,b
                     write(*,*) 'ABcoefficients bad'
                     !write(*,*) 'program exited, DO NOT USE output data'
                     !if (rlz.lt.(int(.5*n_rlz))) then
!*                     call closeAllFiles(1)
                     !call mpi_finalize(ierr)
                     !call mpi_abort(MPI_COMM_WORLD,ierr)
!*                        stop
!*                        write(*,*) char(7)
!*                        call sleep(1)
!*                        write(*,*) char(7)
                     !endif
                 endif

                 if (writeElecField) then
                    do sli = 1,M+2
                       do a = 1,nmax !ntot  !== only write field (not field and deriv)
                          write(129,*) dreal(v_fwd(a,sli)),dimag(v_fwd(a,sli))
!                          write(141,*) dreal(v_bck(a,sli)),dimag(v_bck(a,sli))
                       enddo
                    enddo
                 endif !== writeElecField
              enddo !== b, open channels only
           endif !== findElecField
           call flush(125)
           !== allow for interrupt by changing seed to positive. Saves position
           if (mod(gli,1).eq.0) then
              open(123,file='seed.input')
              read(123,*) seed_interrupt 
              close(123)
              if (seed_interrupt.gt.0) then
                 open(132,file='quasi1d_rect_parameters.input') ! input what parameters to use next time
                 write(132,*) omgi,omgf,n_omg
                 write(132,*) 'omgi    omgf   n_omg' !== input variables documentation
                 write(132,*) ' ' !== second and third lines skipped when reading input
                 write(132,*) alphas,varyStrength,alphaai,alphaaf,n_alphaa
                 write(132,*) 'alphas  varyStrength   alphaai       alphaaf          n_alphaa'
                 write(132,*) ' '
                 write(132,*) L,W,M,n_closed,rlz,n_rlz,lambda,se_step
                 write(132,*) 'L    W      M   n_closed start_rlz n_rlz lambda se_step'
                 write(132,*) ' '
                 write(132,*) percntwriteJWD,binsize,initialgainstep,n_decgs
                 close(132)

                 !== once variables are saved, close all files before exiting
!*                 call closeAllFiles(1) !== value of integer does not matter
!*                 stop  !== exit the program since interrupt=true
!*                 write(*,*) char(7)
!*                 call sleep(1)
!*                 write(*,*) char(7)
              endif !== seed interrupt
           endif !== mod(rlz,1)
     
           !== keep track of how many time the ran2(seed) was called; for the 
           !== purpose of "pausing" and resuming
           if (pausable.AND.(gli.eq.1).AND.(fli.eq.1)) then  !== only need to write once per frequency, gain
             write(130,*) seed_counter !,fli,gli,rlz
             seed_counter=0 !== after writing, reset
             call flush(130)
           endif

           if (repeatgainloop.AND.(gli.eq.n_alphaa)) then
!~~              write(163,*) gli,fli,rlz,'did not exceed critical gain for this freq, rlz'  !j,w write in the subroutine
           endif

           if (gli.eq.n_alphaa) then !== maximum number of searches is set in the input file by n_alphaa
              !== n_alphaa serves as an upper limit on number of loops when finding critical gain
              repeatgainloop = .false.
           endif
        enddo !== gain loop

        !== Shows the percentage of the calculation done        for frequencies
        !== every 5 frequencies
        if (mod(fli,50).eq.0) then  !== if number of frequencies=1, then no screen output
!           write(*,*) (fli/(n_omg+0.0))*100,"% of frequencies are done for rlz",rlz
!           totaltime = etime(elapsed)
!           write(*,*) 'run time so far:', totaltime/60.0,' minutes'
!           write(*,*) 'projected end time for rlz', rlz,'is',(n_omg/(fli+0.0))*totaltime/60.0,'minutes'
!                      ' seconds, or ',(n_omg/(fli+0.0))*totaltime/60, &
!                      ' minutes, or ',((n_omg/(fli+0.0))*totaltime/60)/60,' hours'

           write(125,*) (fli/(n_omg+(0.0d+0)))*100,"% of frequencies are done for rlz",rlz
           totaltime = etime(elapsed)
           write(125,*) 'run time so far:', totaltime/60.0,' minutes'
!           write(125,*) 'projected end time for rlz', rlz,'is',(n_omg/(fli+0.0))*totaltime/60.0,'minutes'
!                        ' seconds, or ',(n_omg/(fli+0.0))*totaltime/60, &
!                        ' minutes, or ',((n_omg/(fli+0.0))*totaltime/60)/60,' hours'
        endif
     enddo !== frequency loop

     !== Shows the percentage of the calculation done   
     !== every 1 realization
     if (mod(rlz,100).eq.0) then  !== if the mod of rlz loop with 100 is zero, then write to the terminal
        totaltime = etime(elapsed)
        write(125,*) rlz/(n_rlz+(0.0d+0))*100,"% of realizations are done"
        write(125,*) 'run time so far:', totaltime/60.0,' minutes'
        write(*,*) rlz/(n_rlz+(0.0d+0))*100,"% of realizations are done"
        write(*,*) 'run time so far:', totaltime/60.0,' minutes'
     endif
  enddo !==for the realization loop

!== moved into gain loop, so that if program crashes, there's still a gain value in the file
!     if (writeGain) then
!        do gli=1,n_alphaa  
!           write(137,*) gain(gli)
!        enddo
!     endif 
 
!== moved into freq loop, so that if program crashes, there's still a frequency in the file    
!     if (writeFrequency) then
!        do fli=1,n_omg 
!           write(135,*) freq(fli)
!        enddo
!     endif 

  if (writeaveTab) then
     do c=1,n_alphaa
        do b=1,nmax
            do a=1, n_open
              aveT(b,a,gli) = aveT(b,a,gli)/(sum(condbinsuf(1,:,gli))*(1.0d+0))
              aveR(b,a,gli) = aveR(b,a,gli)/(sum(condbinsuf(1,:,gli))*(1.0d+0))
            enddo
            write(151,*) aveT(b,:,gli) !== note: there is no imaginary part at this point, 
                                            !== but dimag(0) gets written if nothing is specified
        enddo
     enddo
     !do a=1,n_open
     !   write(151,*) dreal(aveR(b,:))
     !enddo
  endif

  if (writedistcriticalgain.AND.findcriticalgain) then
    write(149,*) gainbin
  endif
 
  !== place a marker with the data that shows the program finished. Useful 
  !== for no-email NIC sessions and necessary for the Windows cluster crash recovery
  if (writeFinish) then
     open(134,file='out_finished.status',POSITION='APPEND')  !== note: should not 
     !== need append, as this is hypothetically the first time the file is 
     !== being written to. However, we would want to observe any deviations
     !== from that expectation
     write(134,*) 'finished run'
     close(134)
  endif
  totaltime = etime(elapsed)
  write(125,*) 'total run time= ', totaltime/60.0,' minutes'
  !== not as useful: , ' user=', elapsed(1), ' system=', elapsed(2)
  call closeAllFiles(1)  !== value of integer does not matter

  !== begin experimental section: using MPI_REDUCE(MPI_SUM) to do some output variable consolidation for J,W sums

  !== free up memory that is no longer needed
  deallocate(ChunkMatrix,DetFree, DetScat,DetChunk)  
  deallocate(Indx,ScatteringMatrix,FreeMatrix,Im)  
  deallocate(Kperp,Kpara,freq,gain,Ref,Trans,Vbc,Vr,Vt)  
  deallocate(InvFreeMatrix,InvScatMatrix)  
  deallocate(Gm,Hm,Tm,Zm,Zmt,Rm,Medium)  
  deallocate(Tempse,SnN,v_fwd,v_bck)  
  deallocate(dxFixed,expon,expon_temp,scatr_indx,kparaMat)

  close(205)
  close(206)
  close(207)

  deallocate(tebinsuf,gainbin,condbinspw,condbinsuf)
  if (findjwdunitflux) then
     deallocate( Jflux_l_sum,Jflux_r_sum,Wdist_r_sum,Wdist_l_sum,&
                 Jflux_l_n_sum,Jflux_r_n_sum,Wdist_l_n_sum,Wdist_r_n_sum)
  endif
  if (findJWDplanewave) then
     deallocate(Apwscatsum,Bpwscatsum,Jflux_pw_l_n_sum,Jflux_pw_r_n_sum,&
                Wdist_pw_l_n_sum,Wdist_pw_r_n_sum,correlationE_sum)
  endif
end program transfer_matrix_2d

!== begin functions and subroutines section

!======== Medium creation ====================
!== method 2: if M>100, then divide the waveguide into smaller
!== "sub-blocks" and attempt random placement within each sub-block
!== CALLED BY main body
!== DEPENDS ON ran2
subroutine mediumcreation(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
                   finddx,dxFixed,scatr_indx,ABstep,rank)
  implicit none !== all of the variables should be of an explicit type
  !== the following variable types need to match type from main body
  !== arguments, scalar
  integer :: M,&        !== number of scatterers
             seed,&
             rlz
  integer*4 :: rank
  real ::    ran2       !== function needs to be declared type
  real*8 :: W,L         !== width, length of waveguide
  real*8 :: ABstep
  !== arguments, array
  real  :: Medium(M,2) !== array for holding scatterer positions
  real*8 :: dxFixed(1,1:(int(L*ABstep)+1))
  !== arguments, logical
  logical :: finddx
  logical :: writeMedium
  integer :: scatr_indx(1:(int(L*ABstep)+1))
  !== local variables, scalar
  real*8  :: dL,dW,&    !== sub-block length, width
             Temp,&
             Lmin
  integer :: nM,nMmax,& !== maximum number of scatterers to put in a sub-block
             nL,nW,&    !== number of sub-blocks along length,width
             a1,a2,&
             a01,a02,&
             a11,a12,&  !== used for medium blocks
             a61,a62,&
             a71,a72,&
             a81,a82,&  !== disorder speedup parameters
             a,b,c,&    !== loop indices
             seed_rlz,&
             n_rlz,&
             seed_counter
  integer :: summer
  !== local variables, array
  real*8 :: scatpos(1:(M+2)),&
            dx(1:(int(L*ABstep)+1))

  !=== How to divide the random medium into small parts === 
  nMmax=100 ! maximum number of scatterers per sub-block
  if (M.gt.nMmax) then
     if (W.lt.L)then
        nW=int(W/(L*W/(int(M/nMmax)+(0.0d+0)))**0.5)
        if(nW.lt.1)then
           nW=1
        endif
        nL=int(M/(nW+(0.0d+0))/(nMmax+(0.0d+0)))
        if(nL.lt.1)then
           nL=1
        endif
     else
        nL=int(L/(L*W/(int(M/nMmax)+(0.0d+0)))**0.5)
        if(nL.lt.1)then
           nL=1
        endif
        nW=int(M/(nL+(0.0d+0))/(nMmax+(0.0d+0)))
        if(nW.lt.1)then
           nW=1
        endif
     endif
     nM=int(M/(nL+(0.0d+0))/(nW+(0.0d+0)))
  else
     nW=1 ! Number of divisions in W
     nL=1 ! Number of divisions in L
     nM=M ! Number of scatterers in a division
  endif
  dL=L/(nL+(0.0d+0))
  dW=W/(nW+(0.0d+0))
!  write(125,*) 'nW = ',nW
!  write(125,*) 'nL = ',nL
!  write(125,*) 'nM = ',nM
!  write(125,*) 'dW = ',dW
!  write(125,*) 'dL = ',dL
  !=== End of division ====================================

  !== open realization-specific output files
  
 ! write(125,*) 'Generating random positions of scatterers for rlz=',rlz,'  ...'
  !==================Create Medium========================
  ! generate a matrix of (x,y) position of point scatterers
  Lmin  = 0.5 * sqrt(L*W/M)
  do b=1,nW
        do c=1,nL
           !write(125,*) 'b=',b,' c=',c
           a1=((b-1)*nL+c-1)*nM+1
           a2=((b-1)*nL+c  )*nM
           a01=1
           a02=nM
           a11=a1-nM        !----------------------------
           a12=a2-nM        !_b1:c1!______!______!__->__!
           a61=a1-nM*(nL-1) !______!______!______!______!
           a62=a2-nM*(nL-1) !______!__8___!__7___!__6___!
           a71=a1-nM*(nL  ) !______!__1___!_X9X__!__5___!
           a72=a2-nM*(nL  ) !______!__2___!__3___!__4___!
           a81=a1-nM*(nL+1) !__I___!______!______!______!
           a82=a2-nM*(nL+1) !__V___!______!______!______!
           a=a1

           do while (a.le.a2) 
              Medium(a,1) = (c-1+ran2(seed))*dL
              Medium(a,2) = (b-1+ran2(seed))*dW
              seed_counter = seed_counter +2

        !== Check that no scatterers are closer then Lmin
              if    ((b.eq.1).and.(c.eq.1))then
                 if(&
minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#1
minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#2
minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#3
minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#4
minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#5
minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#6
minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#7
minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#8
minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
                    a = a - 1
                 endif
              elseif((b.gt.1).and.(c.eq.1))then
                 if(&
minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#1
minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
minval((Medium(a,1)      -Medium(a61:a62,1))**2+(Medium(a,2)      -Medium(a61:a62,2))**2).lt.Lmin**2 .or.& !#6
minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
minval((Medium(a,1)-(c-2)*dL-Medium(1:Nm,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#8
minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
                    a = a - 1
                 endif
              elseif((b.eq.1).and.(c.gt.1))then
                 if(&
minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#6
minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#7
minval((Medium(a,1)-(c-2)*dL-Medium(1:Nm,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#8
minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
                    a = a - 1
                 endif
              elseif((b.gt.1).and.(c.gt.1).and.(c.lt.nL))then
                if(&
minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
minval((Medium(a,1)      -Medium(a61:a62,1))**2+(Medium(a,2)      -Medium(a61:a62,2))**2).lt.Lmin**2 .or.& !#6
minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
minval((Medium(a,1)      -Medium(a81:a82,1))**2+(Medium(a,2)      -Medium(a81:a82,2))**2).lt.Lmin**2 .or.& !#8
minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
                    a = a - 1
                 endif
           elseif((b.gt.1).and.(c.eq.nL))then
                 if(&
minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#6
minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
minval((Medium(a,1)      -Medium(a81:a82,1))**2+(Medium(a,2)      -Medium(a81:a82,2))**2).lt.Lmin**2 .or.& !#8
minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
                    a = a - 1
                 endif
               endif
              a=a+1
           enddo !== while
        enddo  !== c
  enddo   !== b  (rlz)
     
     !== Sort positions of the scatterers in order of increasing x's
     Temp = 0
     do b=1, M
        do a=1, M-1
           if (Medium(a,1).gt.Medium(a+1,1)) then       
              Temp = Medium(a,1)                        !Move the x-values
              Medium(a,1) = Medium(a+1,1)               
              Medium(a+1,1) = Temp
              Temp = Medium(a,2)                   !Move the y-values
              Medium(a,2) = Medium(a+1,2)
              Medium(a+1,2) = Temp                                              
           endif
        enddo !== a
     enddo !== b
     !== Save generated positions of the scatterers into a file
     if (writeMedium) then
       do a=1, M
          write (100,*) Medium(a,:)
       enddo
       close(100)
     endif
     !=======================================================

  !== dx = 0:lambda/10:L;
  dx(1) = 0
  do a = 2,(int(L*ABstep)+1)
     dx(a) = dx(a-1)+ 0.1
  enddo
  !== scatpos = [0; Medium(1:M,1); L];
  scatpos(1) =0
  scatpos(2:(M+1)) = Medium(1:M,1)
  scatpos(M+2) = L

  if (finddx) then
     do a=1,(int(L*ABstep)+1) !== size of dx

        !== scatr_indx(a) = sum(dx(a).gt.Medium(:,1))+1;  !== used for A, B
        summer=0
        do b =1,M
           if (dx(a).gt.Medium(b,1)) then
              summer = summer+1
           endif
        enddo
        scatr_indx(a) = summer+1

        dxFixed(1,a) = dx(a) - scatpos(scatr_indx(a));

     enddo
  endif
!     write(125,*) '... done with generating scatterers, now onto the matrix math.'

end subroutine mediumcreation

!== CALLED BY 
!== DEPENDS ON 
! subroutine mediumcreationBAD(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
!             findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
!   implicit none !== all of the variables should be of an explicit type
!   !== the following variable types need to match type from main body
!   integer :: M,&        !== number of scatterers
!              seed,&
!              rlz
!   integer*4 :: rank
!   real ::    ran2       !== function needs to be declared type
!   real*8 :: W,L         !== width, length of waveguide
!   real*8 :: ABstep
!   real  :: Medium(M,2) !== array for holding scatterer positions
!   real*8 :: dxFixed(1,1:(int(L*ABstep)+1))
!   logical :: findJWDplanewave,findjwdunitflux
!   logical :: writeMedium
!   integer :: scatr_indx(1:(int(L*ABstep)+1))
! 
!   !== local variables
!   real*8  :: dL,dW,&    !== sub-block length, width
!              Temp,&
!              Lmin
!   real*8 :: scatpos(1:(M+2)),&
!             dx(1:(int(L*ABstep)+1))
!   integer :: nM,nMmax,& !== maximum number of scatterers to put in a sub-block
!              nL,nW,&    !== number of sub-blocks along length,width
!              a1,a2,&
!              a01,a02,&
!              a11,a12,&  !== used for medium blocks
!              a61,a62,&
!              a71,a72,&
!              a81,a82,&  !== disorder speedup parameters
!              a,b,c,&    !== loop indices
!              seed_rlz,&
!              n_rlz,&
!              seed_counter
!   integer :: summer
! 
!   !=== How to divide the random medium into small parts === 
!   nMmax=100 ! maxium number of scatterers per sub-block
!   ! REMOVED IF STATEMENT, SO NO SUB-BLOCKS!!!!
!   !if
!      nW=1 ! Number of divisions in W
!      nL=1 ! Number of divisions in L
!      nM=M ! Number of scatterers in a division
!   !endif
!   dL=L/(nL+(0.0d+0))
!   dW=W/(nW+(0.0d+0))
! !  write(125,*) 'nW = ',nW
! !  write(125,*) 'nL = ',nL
! !  write(125,*) 'nM = ',nM
! !  write(125,*) 'dW = ',dW
! !  write(125,*) 'dL = ',dL
!   !=== End of division ====================================
! 
!   !== open realization-specific output files
!   
!   write(125,*) 'Generating random positions of scatterers for rlz=',rlz,'  ...'
!   !==================Create Medium========================
!   ! generate a matrix of (x,y) position of point scatterers
!   Lmin  = 0.5 * sqrt(L*W/M)
!   do b=1,nW
!         do c=1,nL
!            !write(125,*) 'b=',b,' c=',c
!            a1=((b-1)*nL+c-1)*nM+1
!            a2=((b-1)*nL+c  )*nM
!            a01=1
!            a02=nM
!            a11=a1-nM        !----------------------------
!            a12=a2-nM        !_b1:c1!______!______!__->__!
!            a61=a1-nM*(nL-1) !______!______!______!______!
!            a62=a2-nM*(nL-1) !______!__8___!__7___!__6___!
!            a71=a1-nM*(nL  ) !______!__1___!_X9X__!__5___!
!            a72=a2-nM*(nL  ) !______!__2___!__3___!__4___!
!            a81=a1-nM*(nL+1) !__I___!______!______!______!
!            a82=a2-nM*(nL+1) !__V___!______!______!______!
!            a=a1
! 
!            do while (a.le.a2) 
!               Medium(a,1) = (c-1+ran2(seed))*dL
!               Medium(a,2) = (b-1+ran2(seed))*dW
!               seed_counter = seed_counter +2
! 
!         !== Check that no scatterers are closer then Lmin
!               if    ((b.eq.1).and.(c.eq.1))then
!                  if(&
! minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#1
! minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#2
! minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#3
! minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)+dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#4
! minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#5
! minval((Medium(a,1)+dL-Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#6
! minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#7
! minval((Medium(a,1)-dL-Medium(a1:a-1,1))**2+(Medium(a,2)-dW-Medium(a1:a-1,2))**2).lt.Lmin**2 .or.& !#8
! minval((Medium(a,1)   -Medium(a1:a-1,1))**2+(Medium(a,2)   -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
!                     a = a - 1
!                  endif
!               elseif((b.gt.1).and.(c.eq.1))then
!                  if(&
! minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#1
! minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
! minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
! minval((Medium(a,1)      -Medium(a61:a62,1))**2+(Medium(a,2)      -Medium(a61:a62,2))**2).lt.Lmin**2 .or.& !#6
! minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
! minval((Medium(a,1)-(c-2)*dL-Medium(1:Nm,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#8
! minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
!                     a = a - 1
!                  endif
!               elseif((b.eq.1).and.(c.gt.1))then
!                  if(&
! minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
! minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
! minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#6
! minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#7
! minval((Medium(a,1)-(c-2)*dL-Medium(1:Nm,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#8
! minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
!                     a = a - 1
!                  endif
!               elseif((b.gt.1).and.(c.gt.1).and.(c.lt.nL))then
!                 if(&
! minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
! minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
! minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
! minval((Medium(a,1)      -Medium(a61:a62,1))**2+(Medium(a,2)      -Medium(a61:a62,2))**2).lt.Lmin**2 .or.& !#6
! minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
! minval((Medium(a,1)      -Medium(a81:a82,1))**2+(Medium(a,2)      -Medium(a81:a82,2))**2).lt.Lmin**2 .or.& !#8
! minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
!                     a = a - 1
!                  endif
!            elseif((b.gt.1).and.(c.eq.nL))then
!                  if(&
! minval((Medium(a,1)      -Medium(a11:a12,1))**2+(Medium(a,2)      -Medium(a11:a12,2))**2).lt.Lmin**2 .or.& !#1
! minval((Medium(a,1)-(c-2)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#2
! minval((Medium(a,1)-(c-1)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#3
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-0)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#4
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-1)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#5
! minval((Medium(a,1)-(c-0)*dL-Medium(1:nM,1))**2+(Medium(a,2)-(b-2)*dW-Medium(1:nM,2))**2).lt.Lmin**2 .or.& !#6
! minval((Medium(a,1)      -Medium(a71:a72,1))**2+(Medium(a,2)      -Medium(a71:a72,2))**2).lt.Lmin**2 .or.& !#7
! minval((Medium(a,1)      -Medium(a81:a82,1))**2+(Medium(a,2)      -Medium(a81:a82,2))**2).lt.Lmin**2 .or.& !#8
! minval((Medium(a,1)       -Medium(a1:a-1,1))**2+(Medium(a,2)       -Medium(a1:a-1,2))**2).lt.Lmin**2 )then !#9
!                     a = a - 1
!                  endif
!                endif
!               a=a+1
!            enddo !== while
!         enddo  !== c
!   enddo   !== b  (rlz)
!      
!      !== Sort positions of the scatterers in order of increasing x's
!      Temp = 0
!      do b=1, M
!         do a=1, M-1
!            if (Medium(a,1).gt.Medium(a+1,1)) then       
!               Temp = Medium(a,1)                        !Move the x-values
!               Medium(a,1) = Medium(a+1,1)               
!               Medium(a+1,1) = Temp
!               Temp = Medium(a,2)                   !Move the y-values
!               Medium(a,2) = Medium(a+1,2)
!               Medium(a+1,2) = Temp                                              
!            endif
!         enddo !== a
!      enddo !== b
!      !== Save generated positions of the scatterers into a file
!      if (writeMedium) then
!        do a=1, M
!           write (100,*) Medium(a,:)
!        enddo
!        close(100)
!      endif
!      !=======================================================
! 
! 
!   !== dx = 0:lambda/10:L;
!   dx(1) = 0
!   do a = 2,(int(L*ABstep)+1)
!      dx(a) = dx(a-1)+ 0.1
!   enddo
!   !== scatpos = [0; Medium(1:M,1); L];
!   scatpos(1) =0
!   scatpos(2:(M+1)) = Medium(1:M,1)
!   scatpos(M+2) = L
! 
!   if (findJWDplanewave.OR.findjwdunitflux) then
!      do a=1,(int(L*ABstep)+1) !== size of dx
! 
!         !== scatr_indx(a) = sum(dx(a).gt.Medium(:,1))+1;  !== used for A, B
!         summer=0
!         do b =1,M
!            if (dx(a).gt.Medium(b,1)) then
!               summer = summer+1
!            endif
!         enddo
!         scatr_indx(a) = summer+1
! 
!         dxFixed(1,a) = dx(a) - scatpos(scatr_indx(a));
! 
!      enddo
!   endif
! !     write(125,*) '... done with generating scatterers, now onto the matrix math.'
! 
! end subroutine mediumcreationBAD
! 
! !==============================================================================
! !== Medium creation alternative: stacking scatterers
! !== The following method "stacks" the scatterers according to a desired 
! !== distribution, instead of random placement. Should be faster. 
! !== == == NEEDS TO BE TESTED == == ==!
! !== -verify that it creates a medium with the same characteristics as the 
! !== previous methods, and is faster
! !== -issue: the ability to specify number of scatterers may be compromised. ie, 
! !== the last scatterer(s) may be improperly placed, or not able to be placed at all.
! !== -once $\Delta x$ is chosen, choose a uniformly distributed y. Then check
! !== absolute separation with previous scatterer. If too close, regenerate y.
! !== -note: do not need to sort scatterers if using this method
! !== ==
!== CALLED BY 
!== DEPENDS ON 
! subroutine mediumcreationSTACKING(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
!             findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
!   implicit none !== all of the variables should be of an explicit type
!   !== the following variable types need to match type from main body
!   integer :: M,&        !== number of scatterers
!              seed,&
!              rlz
!   integer*4 :: rank
!   real ::    ran2       !== function needs to be declared type
!   real*8 :: W,L         !== width, length of waveguide
!   real*8 :: ABstep
!   real  :: Medium(M,2) !== array for holding scatterer positions
!   real*8 :: dxFixed(1,1:(int(L*ABstep)+1))
!   logical :: findJWDplanewave,findjwdunitflux
!   logical :: writeMedium
!   integer :: scatr_indx(1:(int(L*ABstep)+1))
! 
!   !== local variables
!   real*8 :: scatpos(1:(M+2)),&
!             dx(1:(int(L*ABstep)+1))
!   integer :: a,b,&
!              seed_counter
!   real*8 :: Lmin
!   logical :: tryAgain,&
!              badLength
!   integer :: summer
! 
!   Medium(:,:) = (0.0d+0)
!   write(125,*) 'Generating random positions of scatterers for rlz=',rlz,'  ...'
! 
!   Lmin = 0.5 * sqrt(L*W/M)
! 
!   badLength = .true.
!   do while (badLength)
!      do a=1,M 
!         !== according to Dr Yamilov, P(dx)=1/a*exp[-dx/a]
!         !== where a is the average value of dx.
!         !== a=[system length]/[number of scatterers]
!         if (a.eq.1) then
!            Medium(a,1) = abs(((L/M)*log(1-ran2(seed))))
!         else
!            Medium(a,1) = abs(((L/M)*log(1-ran2(seed))))+Medium(a-1,1)
!         endif
!         Medium(a,2) = (ran2(seed))*W
!         seed_counter = seed_counter +2
!         !==if separation with previous is small, then regenerate Medium(a,2) and check again
!         tryagain=.true.
!         do while (tryAgain)
!            if(((Medium(a,1)-Medium(a-1,1))**2+(Medium(a,2)-Medium(a-1,2))**2).lt.Lmin**2)then
!               Medium(a,2) = (ran2(seed))*W !== another random y
!               seed_counter = seed_counter +1
!               tryAgain=.true.
!            else
!               tryAgain=.false.
!            endif
!         enddo
!      enddo
!      if (Medium(M,1).gt.L) then
!         badLength = .true. !== too long (>L)
!      elseif ( (abs(((L/M)*log(1-ran2(seed))))+Medium(a-1,1)).lt.L ) then
!         badLength = .true. !== too short (?may not be needed)
!         seed_counter = seed_counter +1
!      else
!         badLength = .false.
!         seed_counter = seed_counter +1 !== if this was reached, then the elseif statement was used
!      endif
!   enddo
! 
!   !== Save generated positions of the scatterers into a file
!   if (writeMedium) then
!      do a=1, M
!         write (100,*) Medium(a,:)
!      enddo
!      close(100)
!   endif
! 
!   !== dx = 0:lambda/10:L;
!   dx(1) = 0
!   do a = 2,(int(L*ABstep)+1)
!      dx(a) = dx(a-1)+ 0.1
!   enddo
!   !== scatpos = [0; Medium(1:M,1); L];
!   scatpos(1) =0
!   scatpos(2:(M+1)) = Medium(1:M,1)
!   scatpos(M+2) = L
! 
!   if (findJWDplanewave.OR.findjwdunitflux) then
!      do a=1,(int(L*ABstep)+1) !== size of dx
! 
!         !== scatr_indx(a) = sum(dx(a).gt.Medium(:,1))+1;  !== used for A, B
!         summer=0
!         do b =1,M
!            if (dx(a).gt.Medium(b,1)) then
!               summer = summer+1
!            endif
!         enddo
!         scatr_indx(a) = summer+1
! 
!         dxFixed(1,a) = dx(a) - scatpos(scatr_indx(a));
! 
!      enddo
!   endif
! 
!   write(125,*) '... done with generating scatterers, now onto the matrix math.'
! end subroutine mediumcreationSTACKING
! !==============================================================================
! 
! !==== OLD MEDIUM CREATION METHOD, ONE BIG LUMP (not optimized) ====
!== CALLED BY 
!== DEPENDS ON 
! subroutine mediumcreationOLD(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
!          findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
!   implicit none !== all of the variables should be of an explicit type
!   !== the following variable types need to match type from main body
!   integer :: M,&        !== number of scatterers
!              seed,&
!              rlz
!   integer*4 :: rank
!   real ::    ran2       !== function needs to be declared type
!   real*8 :: W,L         !== width, length of waveguide
!   real*8 :: ABstep
!   real  :: Medium(M,2) !== array for holding scatterer positions
!   real*8 :: dxFixed(1,1:(int(L*ABstep)+1))
!   logical :: findJWDplanewave,findjwdunitflux
!   logical :: writeMedium
!   integer :: scatr_indx(1:(int(L*ABstep)+1))
! 
!   !== local variables
!   real*8 :: scatpos(1:(M+2)),&
!             dx(1:(int(L*ABstep)+1))
!   integer :: a,b,&
!              seed_counter
!   real*8 :: Lmin
!   integer :: summer
! 
!   Lmin = 0.5 * sqrt(L*W/M)
!   a = 1
!   do while (a.le.M) 
!      Medium(a,1) = (ran2(seed)*L)
!      Medium(a,2) = (ran2(seed)*W)
!      !     Medium(a,1) = rand(0)*L
!      !     Medium(a,2) = rand(0)*W
!      seed_counter = seed_counter+2
!      ! Check that no scatterers are closer then Lmin
!      if(minval((Medium(a,1)  -Medium(1:a-1,1))**2+(Medium(a,2)  -Medium(1:a-1,2))**2).lt.Lmin**2 .or.&
!           minval((Medium(a,1)  -Medium(1:a-1,1))**2+(Medium(a,2)-W-Medium(1:a-1,2))**2).lt.Lmin**2 .or.&
!           minval((Medium(a,1)  -Medium(1:a-1,1))**2+(Medium(a,2)+W-Medium(1:a-1,2))**2).lt.Lmin**2 .or.&
!           minval((Medium(a,1)-L-Medium(1:a-1,1))**2+(Medium(a,2)  -Medium(1:a-1,2))**2).lt.Lmin**2 .or.&
!           minval((Medium(a,1)+L-Medium(1:a-1,1))**2+(Medium(a,2)  -Medium(1:a-1,2))**2).lt.Lmin**2 )then 
!         a = a - 1
!      endif
!      a=a+1
!      !if (mod(a*10,M).eq.0)then
!      !   write(125,*) int(10*a/M)*10,'% done'
!      !endif
!   enddo !== while
! 
!   !== Save generated positions of the scatterers into a file
!   if (writeMedium) then
!      do a=1, M
!         write (100,*) Medium(a,:)
!      enddo
!      close(100)
!   endif
! 
!   !== dx = 0:lambda/10:L;
!   dx(1) = 0
!   do a = 2,(int(L*ABstep)+1)
!      dx(a) = dx(a-1)+ 0.1
!   enddo
!   !== scatpos = [0; Medium(1:M,1); L];
!   scatpos(1) =0
!   scatpos(2:(M+1)) = Medium(1:M,1)
!   scatpos(M+2) = L
! 
!   if (findJWDplanewave.OR.findjwdunitflux) then
!      do a=1,(int(L*ABstep)+1) !== size of dx
! 
!         !== scatr_indx(a) = sum(dx(a).gt.Medium(:,1))+1;  !== used for A, B
!         summer=0
!         do b =1,M
!            if (dx(a).gt.Medium(b,1)) then
!               summer = summer+1
!            endif
!         enddo
!         scatr_indx(a) = summer+1
! 
!         dxFixed(1,a) = dx(a) - scatpos(scatr_indx(a));
! 
!      enddo
!   endif
! 
! !  write(125,*) '... done with generating scatterers, now onto the matrix math.'
! end subroutine mediumcreationOLD
! 
! subroutine mediumcreationNOCOLLISIONCHECK(rlz,M,W,L,Medium,writeMedium,seed,seed_counter,&
!        findJWDplanewave,findjwdunitflux,dxFixed,scatr_indx,ABstep,rank)
!   implicit none !== all of the variables should be of an explicit type
!   !== the following variable types need to match type from main body
!   !== arguments, scalar
!   integer :: M,&        !== number of scatterers
!              seed,&
!              rlz
!   integer*4 :: rank
!   real ::    ran2       !== function needs to be declared type
!   real*8 :: W,L         !== width, length of waveguide
!   real*8 :: ABstep
!   !== arguments, array
!   real  :: Medium(M,2) !== array for holding scatterer positions
!   real*8 :: dxFixed(1,1:(int(L*ABstep)+1))
!   integer :: scatr_indx(1:(int(L*ABstep)+1))
!   !== arguments, logical
!   logical :: findJWDplanewave,findjwdunitflux
!   logical :: writeMedium
! 
!   !== local variables
!   real*8 :: scatpos(1:(M+2)),&
!             dx(1:(int(L*ABstep)+1))
!   integer :: a,b,&
!              seed_counter
!   integer :: summer
! 
!   do a=1,M 
!      Medium(a,1) = (ran2(seed)*L)
!      Medium(a,2) = (ran2(seed)*W)
!      seed_counter = seed_counter+2
!   enddo
! 
!   !== Save generated positions of the scatterers into a file
!   if (writeMedium) then
!      do a=1, M
!         write (100,*) Medium(a,:)
!      enddo
!      close(100)
!   endif
! 
!   !== dx = 0:lambda/10:L;
!   dx(1) = 0
!   do a = 2,(int(L*ABstep)+1)
!      dx(a) = dx(a-1)+ 0.1
!   enddo
!   !== scatpos = [0; Medium(1:M,1); L];
!   scatpos(1) =0
!   scatpos(2:(M+1)) = Medium(1:M,1)
!   scatpos(M+2) = L
! 
!   if (findJWDplanewave.OR.findjwdunitflux) then
!      do a=1,(int(L*ABstep)+1) !== size of dx
! 
!         !== scatr_indx(a) = sum(dx(a).gt.Medium(:,1))+1;  !== used for A, B
!         summer=0
!         do b =1,M
!            if (dx(a).gt.Medium(b,1)) then
!               summer = summer+1
!            endif
!         enddo
!         scatr_indx(a) = summer+1
! 
!         dxFixed(1,a) = dx(a) - scatpos(scatr_indx(a));
! 
!      enddo
!   endif
! 
! !  write(125,*) '... done with generating scatterers, now onto the matrix math.'
! 
! end subroutine mediumcreationNOCOLLISIONCHECK

!== ran2 generates a uniform random number between 0,1
!== the built-in fortran function rand(0) has the same function (and is
!== probably faster) but has different (repeatable) output depending on 
!== whether the ifort or gfortran compilers are used.  
!===
! long period (>2x10^18) random number generator of l'ecuyer with bays-durham shffle
! and added safeguards. returns a uniform random deviate between 0.0 and 1.0 (exclusive
! of the endpoint values). call with idum a negative integer to initialize; thereafter
! , do not alter idum between successive deviates in a sequence. rnmx should
! approximate the largest floating values that is less than 1
!
!== see http://books.google.com/books?id=gn_4mpdN9WkC&pg=PA271
!== CALLED BY main body, medium creation, scatmatrixcreation, inputvariables
!== DEPENDS ON NONE
function ran2(idum)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer*4 idum
  !== local, scalar
  integer*4 im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
  real ran2,am,eps,rnmx
  parameter(im1=2147483563, im2=2147483399,am=1./im1, imm1=im1-1,  &
       ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211,   &
       ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
  integer*4 idum2, j, k, iv(ntab), iy
  save iv,iy,idum2
  data idum2/123456789/,iv/ntab*0/, iy/0/
  if (idum .le. 0) then
     idum=max(-idum,1)
     idum2=idum
     do 11 j=ntab+8,1,-1
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if(idum .lt. 0) idum=idum+im1
        if (j .le. ntab) iv(j)=idum
     11 continue
     iy=iv(1)
  endif
  k=idum/iq1
  idum=ia1*(idum-k*iq1)-k*ir1
  if (idum .lt. 0) idum=idum+im1
  k=idum2/iq2
  idum2=ia2*(idum2-k*iq2)-k*ir2
  if (idum2 .lt. 0) idum2=idum2+im2
  j=1+iy/ndiv
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy .lt. 1) iy=iy+imm1
  ran2=min(am*iy,rnmx)
  return
end function ran2

!*****************************************************************************
!** input variables from files =================
!== called once per run
!== CALLED BY main body
!== DEPENDS ON ran2, int2String
subroutine inputvariables(seed,rank,omgi,omgf,n_omg,percntwriteJWD,&
                alphas,varyStrength,alphaai,alphaaf,n_alphaa,&
                L,W,M,n_closed,start_rlz,n_rlz,lambda,se_step,&
                nmax,ntot,n_open,M_se,Lmin,writeElecField,findeignvalues,&
                writeRunningStatus,writeFrequency,writeReflectionMatrix,&
                writeTransmissionMatrix,writeTmRmDet, writeAccuracy,& 
                writeScatFreeMatrix,writeMedium,writeDeterminats,& 
                writeSelfEmbedTmRm,writeGain,readGain,writeFinish,readMedium,&
                findElecField,writeABCoefficientsPW,writeABCoefficientsUF,writeElecFieldError,&
                findJWDplanewave,findjwdunitflux,doefse,ABstep,writeSummedJW,writeEachJW,pausable,&
                writeTEuf,writegEuf,writegEpw,binsize,initialgainstep,n_decgs,findcriticalgain,&
                writeTEufdist,writecriticalgain,writedistcriticalgain,&
                writeconductionufdist,writeconductionpwdist,writeaveTab,&
                writeufJWDchan,writepwJWDchan,writeEzyPW,writeEzyUF,find_dzy,dzy_gain_resolved)
  implicit none !== all of the variables should be of an explicit type
  !== external functions
  real    :: ran2
!  include 'mpif.h' !== MPIF90
  !== arguments, scalar
  integer :: seed,&
             seed_input,&
             n_alphaa,&
             n_omg,&
             M,&
             n_closed,&
             n_decgs,&
             start_rlz,&
             n_rlz,&
             se_step,&
             nmax,ntot,n_open,&
             binsize,&
             M_se,&
             fli, gli,rlz,& !== pause function
             tmp_indx,&  !== pause function
             tmp_fli_indx,&
             tmp_gli_indx,&
             num_times_call_seed,& !== used for "pause" functionality. 
                 !== Read from file to loop index correct number of times
             seed_rlz !== used for "pause" functionality. 
                 !== Read from file. the above "num_times_call_seed" 
                 !== is per realization. Seed_rlz specifies rlz
             integer*4 :: rank !== MPI
!  character(len=50)        :: file_name
!  integer*4                :: ch0,ch1,ch2,ch3,ch4,ch5,ch6,ch7
  character*8              :: ch_n
  character*8              :: int2String
  real*8  :: omgi,omgf,&
             alphas,alphaai,alphaaf,&
             L,W,&
             ABstep,&
             initialgainstep,&
             lambda,&
             Lmin,&
             percntwriteJWD
  !== arguments, logical
  logical :: writeRunningStatus,&
             writeFrequency,&
             writeReflectionMatrix,& 
             writeTransmissionMatrix,& 
             writeTmRmDet,& 
             writeAccuracy,& 
             writeScatFreeMatrix,& 
             writeMedium,& 
             writeDeterminats,& 
             writeSelfEmbedTmRm,& 
             writeGain,& 
             readGain,&
             writeFinish,&
             readMedium,&
             varyStrength,&
             writeElecField,&
             writeElecFieldError,&
             findeignvalues,&
             finished_status_exists,&
             findElecField,&
             writeABCoefficientsUF,&
             writeABCoefficientsPW,&
             findJWDplanewave,findjwdunitflux,&
             doefse,&
             writeSummedJW,&
             writeEachJW,&
             pausable,&
             writeTEuf,&
             writegEuf,&
             writegEpw,&
             writeTEufdist,&
             findcriticalgain,&
             writecriticalgain,&
             writedistcriticalgain,&
             writeconductionufdist,&
             writeconductionpwdist,&
             writeaveTab,&
             writeufJWDchan,&
             writepwJWDchan,&
             writeEzyPW,writeEzyUF,&
             find_dzy,&
             dzy_gain_resolved
  !== local variables, scalar
  real    :: FirstNumber
  real*8 :: tmp !== pause function
  !== local variables, logical
  logical ::     ran2_file_exists

  !== initialize all variables. Formerly, this was in the main body prior to 
  !== calling this subroutine, but that was found to be unnecessary. Further,
  !== initialization with the subroutine is also unnecessary, but it seems 
  !== safer to include it
  seed_input = 0
  seed = 0
  FirstNumber = 0.0
  omgi = (0.0d+0)
  omgf = (0.0d+0)
  n_omg = 0
  alphas = (0.0d+0)
  varyStrength = .false.
  alphaai = (0.0d+0)
  alphaaf = (0.0d+0)
  n_alphaa = 0
  L = (0.0d+0)
  W = (0.0d+0)
  M = 0
  n_closed = 0
  start_rlz = 0
  n_rlz = 0
  lambda = (0.0d+0)
  se_step = 0
  nmax = 0
  ntot = 0
  initialgainstep = (0.0d+0)
  n_decgs = 0
  n_open = 0
  writeRunningStatus = .true.
  writeFrequency = .true.
  writeReflectionMatrix = .true.
  writeTransmissionMatrix = .true.
  writeTmRmDet = .true.
  writeAccuracy = .true.
  writeScatFreeMatrix = .true.
  writeMedium = .true.
  writeDeterminats = .true.
  writeSelfEmbedTmRm = .true.
  writeGain = .true.
  readGain = .false.
  writeElecField = .true.
  findeignvalues = .true.  
  writeFinish = .true.
  writeABCoefficientsUF = .true.
  writeABCoefficientsPW = .true.
  writeSummedJW = .true.
  writeEachJW = .true.
  writeTEuf = .true.
  writegEuf = .true.
  writegEpw = .true.
  writeTEufdist = .true.
  writecriticalgain = .true.
  writedistcriticalgain = .true.
  writeconductionufdist = .true.
  writeconductionpwdist = .true.
  writeaveTab = .true.
  writeufJWDchan = .true.
  writepwJWDchan = .true.
  writeEzyPW = .true.
  writeEzyUF = .true.

  pausable = .true.
  readMedium = .true.

  find_dzy = .true.
  dzy_gain_resolved = .true.
  findcriticalgain = .true.
  findElecfield = .true.
  findJWDplanewave = .true.
  findjwdunitflux = .true.
  doefse = .true.

  open(123,file='seed.input')                   ! input seed
  read(123,*) seed_input
  if (seed_input.ge.0) then  !== positive seed error detection
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: seed was detected to be positive: ', seed_input
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: seed was detected to be positive: ', seed_input
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: seed was detected to be positive: ', seed_input
     write(*,*) 'program exited without making changes to existing files'
     !== warning: do not just set the seed to negative and continue on, as the 
     !== windows cluster uses the seed sign as an indicator of crash status [ask Ben]
     stop
  endif
  close(123)

  seed = seed_input !== the seed is specific to the node
  !== call ran2 the first time, which upsets the seed number. In place to 
  !== avoid using the first ran2 output, because otherwise seed gets screwed up (?!)
!  write(*,*) 'seed is ',seed
  !== ran2 needs to have a variable passed to it, because it passes that variable back
  FirstNumber=ran2(seed)  
!  write(*,*) 'seed is now ',seed
!  write(*,*) 'and FirstNumber is ',FirstNumber
  
  open(132,file='quasi1d_rect_parameters.input')        ! input what parameters to use
  read(132,*) omgi,omgf,n_omg
  read(132,*) ! skip second
  read(132,*) ! and third lines when reading input
  read(132,*) alphas,varyStrength,alphaai,alphaaf,n_alphaa
  read(132,*) ! skip fourth 
  read(132,*) ! and fifth lines when reading input (allows room for comments in file)
  read(132,*) L,W,M,n_closed,start_rlz,n_rlz,lambda,se_step
  read(132,*) !skip
  read(132,*) !skip
  read(132,*) percntwriteJWD,binsize,initialgainstep,n_decgs
  close(132)

  !== note: the following variable is not in the input file, but could be moved there:
  ABstep = 10  !== for dx=0:lambda/ABstep:L

  n_open = int(2.0*W/lambda)
  !== Total number of channels, both open and closed
  nmax = n_open+n_closed 
  !old: nmax =(W/pi) * sqrt(K2+((Lmin**(-2))*(log(error)**2)))
  ! note: removed "error" from input
  !== Size of all matrices 
  ntot = 2*nmax   !== size of matrices

  if ((omgi.lt.0).OR.(omgf.lt.0).OR.(omgf.lt.omgi).OR.(n_omg.lt.1).OR.&
     (n_alphaa.lt.1).OR.&
     (L.lt.0).OR.(W.lt.0).OR.(M.lt.0).OR.(n_closed.lt.0).OR.&
     (start_rlz.lt.1).OR.(n_rlz.lt.1).OR.(lambda.lt.0).OR.(se_step.lt.0).OR.&
     (n_open.lt.2).OR.(ntot.lt.0)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: one of the input variables is non-physical'
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: one of the input variables is non-physical'
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: one of the input variables is non-physical'
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if ((omgi.lt.(n_open+.25)/(n_open+.5)) .OR. (omgf.gt.(n_open+.75)/(n_open+.5))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: frequency range too wide. Creates a new open channel. '
     write(133,*) 'input = ',omgi,'lower limit = ',((n_open+.25)/(n_open+.5))
     write(133,*) 'input = ',omgf,'upper limit = ',((n_open+.75)/(n_open+.5))
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: frequency range too wide. Creates a new open channel. '
     write(125,*) 'input = ',omgi,'lower limit = ',((n_open+.25)/(n_open+.5))
     write(125,*) 'input = ',omgf,'upper limit = ',((n_open+.75)/(n_open+.5))
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: frequency range too wide. Creates a new open channel. '
     write(*,*) 'input = ',omgi,'lower limit = ',((n_open+.25)/(n_open+.5))
     write(*,*) 'input = ',omgf,'upper limit = ',((n_open+.75)/(n_open+.5))
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if (int(percntwriteJWD*n_rlz).eq.0) then 
     !== floating point exception, since mod(rlz,%writeJWD*nrlz) is dividing by zero
     write(125,*) 'resetting percntwriteJWD to 1 because rounding down causes divide by 0'
     write(*,*) 'resetting percntwriteJWD to 1 because rounding down causes divide by 0'
     percntwriteJWD=1
  endif

!== need to implement this correctly:
!  if ((1.0*int(percntwriteJWD*n_rlz)).NE.(percntwriteJWD*n_rlz)) then
  if (.false.) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: percntwriteJWD of n_rlz is not an integer, which screws up JWD output file',&
                n_rlz,(percntwriteJWD*n_rlz),real(int(percntwriteJWD*n_rlz))
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: percntwriteJWD of n_rlz is not an integer, which screws up JWD output file',&
                n_rlz,(percntwriteJWD*n_rlz),real(int(percntwriteJWD*n_rlz))
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: percntwriteJWD of n_rlz is not an integer, which screws up JWD output file',&
                n_rlz,(percntwriteJWD*n_rlz),real(int(percntwriteJWD*n_rlz))
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if ((findJWDplanewave.OR.findjwdunitflux).AND.((int(L*ABstep)).NE.L*ABstep)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: L*10 is not an integer (necessary for J,W,D calculation)',L,ABstep
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: L*10 is not an integer (necessary for J,W,D calculation)',L,ABstep
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: L*10 is not an integer (necessary for J,W,D calculation)',L,ABstep
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if (mod(M,se_step).NE.0) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: M is not divisible by se_step ',M,se_step
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: M is not divisible by se_step ',M,se_step
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: M is not divisible by se_step ',M,se_step
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  !======= IF start_rlz does not equal 1, then 
  INQUIRE(FILE="out_num_ran2_calls.dat", EXIST=ran2_file_exists)
  if ((start_rlz .NE. 1) .AND. (.NOT. ran2_file_exists )) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: Need the prior number of ran2 calls to start. '
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: Need the prior number of ran2 calls to start. '
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: Need the prior number of ran2 calls to start. '
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  open(128,file='quasi1d_rect_files.input')             ! input what files to write
  read(128,*) !== header line
  read(128,*) writeAccuracy
  read(128,*) writeRunningStatus
  read(128,*) writeReflectionMatrix, writeTransmissionMatrix
  read(128,*) writeMedium, readMedium
  read(128,*) writeDeterminats
  read(128,*) writeTmRmDet
  read(128,*) writeSelfEmbedTmRm
  read(128,*) writeScatFreeMatrix
  read(128,*) writeFinish
  read(128,*) writeFrequency
  read(128,*) writeGain, readGain
  read(128,*) findElecField, doefse
  read(128,*) findeignvalues
  read(128,*) writeElecField, writeElecFieldError
  read(128,*) writeABCoefficientsUF, writeABCoefficientsPW
  read(128,*) findJWDplanewave, findJWDunitflux 
  read(128,*) writeSummedJW, writeEachJW
  read(128,*) writeufJWDchan, writepwJWDchan
  read(128,*) pausable
  read(128,*) writeTEuf, writegEuf, writegEpw
  read(128,*) findcriticalgain, writecriticalgain, writedistcriticalgain
  read(128,*) writeconductionufdist, writeconductionpwdist, writeTEufdist
  read(128,*) writeaveTab
  read(128,*) writeEzyPW,  writeEzyUF
  read(128,*) find_dzy, dzy_gain_resolved
  close(128)  !== now we know which files to write

  !============ "DUH" CHECKS =====================
  !== check for logical inconstancies supplied by user

  !== "writeFinish" needs to come first, since it prevents the program 
  !== from continuing if the program finished previously
  if (writeFinish) then
     !== check to make sure a "finished" marker does not already exist. If it does, exit
     INQUIRE(FILE="out_finished.status", EXIST=finished_status_exists)
     if (finished_status_exists) then
        open(133,file='out_ERROR.log',POSITION='APPEND')
        write(133,*) 'WARNING: finish status marker already exists.',writeFinish
        write(133,*) 'program exited without making changes to existing files'
        close(133)
        write(125,*) 'WARNING: finish status marker already exists.',writeFinish
        write(125,*) 'program exited without making changes to existing files'
        close(125)     
        write(*,*)   'WARNING: finish status marker already exists.',writeFinish
        write(*,*)   'program exited without making changes to existing files'
        stop
     endif
  endif

  if (readGain) then
     !== check to make sure a "finished" marker does not already exist. If it does, exit
     INQUIRE(FILE="quasi1d_rect_gain_array.input", EXIST=finished_status_exists)
     if (.NOT.finished_status_exists) then
        open(133,file='out_ERROR.log',POSITION='APPEND')
        write(133,*) 'WARNING: gain array needs to exist',readGain
        write(133,*) 'program exited without making changes to existing files'
        close(133)
        write(125,*) 'WARNING: gain array needs to exist',readGain
        write(125,*) 'program exited without making changes to existing files'
        close(125)     
        write(*,*)   'WARNING: gain array needs to exist',readGain
        write(*,*)   'program exited without making changes to existing files'
        stop
     endif
  endif

  if ((start_rlz.ne.1).AND.(.NOT.(pausable))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: Pause functionality is turned off, but start_rlz .ne.1 '
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: Pause functionality is turned off, but start_rlz .ne.1 '
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: Pause functionality is turned off, but start_rlz .ne.1 '
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  !========== PAUSE FUNCTION: RESUME ===============
  tmp=(0.0d+0) !== tmp is only used in this loop, but reset it anyways
  if ((start_rlz.ne.1).AND.pausable) then
     open(130,file='out_num_ran2_calls.dat')
     do rlz=1,(start_rlz-1) 
        read(130,*) num_times_call_seed
        !== ran2 is used in for FirstNumber, scattering matrix, and medium creation.
        !== medium is called once per rlz, scattering matrix is called every rlz, fli, gli
        do tmp_indx=1,num_times_call_seed
                 tmp=ran2(seed)
        enddo
     enddo !== rlz
     close(130)
  endif
  tmp=(0.0d+0)
  readMedium = .false. !== account for the re-setting of a paused state

  if (find_dzy.AND.(.NOT.(findJWDplanewave))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: can not find D(z,y) without finding plane wave J,W '
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: can not find D(z,y) without finding plane wave J,W '
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: can not find D(z,y) without finding plane wave J,W '
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if ((.NOT.find_dzy).AND.dzy_gain_resolved) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: can not find D(z,y) gain resolved without finding D(z,y) '
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: can not find D(z,y) gain resolved without finding D(z,y) '
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: can not find D(z,y) gain resolved without finding D(z,y) '
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  if (find_dzy.AND.(.NOT.dzy_gain_resolved).AND.(n_alphaa.gt.1)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: finding D(z,y) but not gain resolved when varying gain give useless data ',n_alphaa
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: finding D(z,y) but not gain resolved when varying gain give useless data ',n_alphaa
     write(125,*) 'program exited without making changes to existing files'
     close(125)
     write(*,*) 'WARNING: finding D(z,y) but not gain resolved when varying gain give useless data ',n_alphaa
     write(*,*) 'program exited without making changes to existing files'
     stop
  endif

  !== critical gain booleans:
  
  if (findcriticalgain.AND.(.NOT.writecriticalgain).AND.(.NOT.writedistcriticalgain)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: unnecessarily finding critical gain ',findcriticalgain,writecriticalgain,writedistcriticalgain
     write(133,*) 'probably need to turn on critical gains or distributions'
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: unnecessarily finding critical gain ',findcriticalgain,writecriticalgain,writedistcriticalgain
     write(125,*) 'probably need to turn on critical gains or distributions'
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: unnecessarily finding critical gain ',findcriticalgain,writecriticalgain,writedistcriticalgain
     write(*,*)   'probably need to turn on critical gains or distributions'
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

  if ((.NOT.findcriticalgain).AND.(writecriticalgain.OR.writedistcriticalgain)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: cannot write cg or dist cg without finding critical gain',&
                 findcriticalgain,writecriticalgain,writedistcriticalgain
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: cannot write cg or dist cg without finding critical gain',&
                 findcriticalgain,writecriticalgain,writedistcriticalgain
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: cannot write cg or dist cg without finding critical gain',&
                 findcriticalgain,writecriticalgain,writedistcriticalgain
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

  if (findcriticalgain.AND.(.NOT.writeGain)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'warning: when finding critical gain, one should record what gains were. Reset writeGain to .true.'
     close(133)
     write(133,*) 'warning: when finding critical gain, one should record what gains were. Reset writeGain to .true.'
     close(125) 
     write(133,*) 'warning: when finding critical gain, one should record what gains were. Reset writeGain to .true.'
     writeGain = .true.
  endif

  if (findcriticalgain.AND.(writeSummedJW.OR.writeEachJW.OR.writeufJWDchan.OR.writepwJWDchan.OR.&
        findJWDplanewave.OR.findJWDunitflux)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: writing J,W from non-uniform gain step results in non-sense'
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: writing J,W from non-uniform gain step results in non-sense'
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: writing J,W from non-uniform gain step results in non-sense'
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

  if ((writeEzyPW.OR.writeEzyUF).AND.((.NOT.findElecField).OR.(.NOT.doefse))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: cannot find E(z,y) without finding electric field and doing self-embed'
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: cannot find E(z,y) without finding electric field and doing self-embed'
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: cannot find E(z,y) without finding electric field and doing self-embed'
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

!  if (writeMedium .AND. readMedium) then
!     open(133,file='out_ERROR.log',POSITION='APPEND')
!     write(133,*) 'WARNING: medium logic variables conflict: ',writeMedium,readMedium
!     write(133,*) 'program exited without making changes to existing files'
!     close(133)
!     write(125,*) 'WARNING: medium logic variables conflict: ',writeMedium,readMedium
!     write(125,*) 'program exited without making changes to existing files'
!     close(125)     
!     write(*,*)   'WARNING: medium logic variables conflict: ',writeMedium,readMedium
!     write(*,*)   'program exited without making changes to existing files'
!     stop
!  endif

  if (((.NOT. findJWDplanewave).AND.(.NOT.findjwdunitflux) .OR. (.NOT. findElecField)) .AND. doefse) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'warning: booleans indicate you may be performing unnecessary electric field self-embedding'
     !write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'warning: booleans indicate you may be performing unnecessary electric field self-embedding'
     !write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'warning: booleans indicate you may be performing unnecessary electric field self-embedding'
     !write(*,*)   'program exited without making changes to existing files'
     !stop
  endif

  if ((.NOT. findElecField).AND. (writeABCoefficientsPW.OR.writeABCoefficientsUF.OR.findJWDplanewave.OR.findJWDunitflux)) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: cannot write A,B or J,W without finding electric field ',&
                    findElecField,writeABCoefficientsUF,writeABCoefficientsPW,findJWDplanewave,findJWDunitflux
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: cannot write A,B or J,W without finding electric field ',&
                    findElecField,writeABCoefficientsUF,writeABCoefficientsPW,findJWDplanewave,findJWDunitflux
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: cannot write A,B or J,W without finding electric field ',&
                    findElecField,writeABCoefficientsUF,writeABCoefficientsPW,findJWDplanewave,findJWDunitflux
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

  if ((findJWDplanewave.OR.findJWDplanewave).AND.(.NOT.(writeSummedJW.OR.writeEachJW.OR.writeufJWDchan.OR.writepwJWDchan))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: booleans indicate you may be performing unnecessary electric field self-embedding'
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: booleans indicate you may be performing unnecessary electric field self-embedding'
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: booleans indicate you may be performing unnecessary electric field self-embedding'
     write(*,*)   'program exited without making changes to existing files'
     stop
  endif

  if ((.NOT. findElecField) .AND. writeABCoefficientsUF.OR.writeABCoefficientsPW) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'warning: cannot find AB,JW without finding electric field, so reset findElecField to .true.'
     close(133)
     write(125,*) 'warning: cannot find AB,JW without finding electric field, so reset findElecField to .true.'
     close(125) 
     write(*,*) 'warning: cannot find AB,JW without finding electric field, so reset findElecField to .true.'
     findElecField = .true.
  endif

  if (n_rlz.gt.99999999) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: n_rlz exceeds 8 digits, thus output files will not be recorded ',n_rlz
     write(133,*) 'program exited without making changes to existing files'
     close(133)
     write(125,*) 'WARNING: n_rlz exceeds 8 digits, thus output files will not be recorded ',n_rlz
     write(125,*) 'program exited without making changes to existing files'
     close(125)     
     write(*,*)   'WARNING: n_rlz exceeds 8 digits, thus output files will not be recorded ',n_rlz
     write(*,*)   'program exited without making changes to existing files'
     !== note: in practice, when n_rlz>99999999, the string becomes ":000", ":001", ... and 
     !== then one just renames the output file. Alternatively, we could increase number
     !== of converted digits to 6 if needed
     stop
  endif


  !===== print to "out_screen.dat" the variables used in calculations =====
  write(125,*) '***** start of initialization numbers *****'
  write(125,*) "L (system length)                               = ",L
  write(125,*) "W (width of the waveguide)                      = ",W
  write(125,*) "scattering strength of point scatterers         = ",alphas
  write(125,*) "initial and final active value                  = ",alphaai,alphaaf,n_alphaa
  write(125,*) "initial and final frequencies, number of steps  = ",omgi,omgf,n_omg
 !write(125,*) "error                                           = ",error
  write(125,*) "M (total number of scatterers)                  = ",M
  write(125,*) "lambda (wavelength ratio (unit-less wavelength)) = ",lambda
  write(125,*) "se_step (number of scatterers in a se step)     = ",se_step 
  write(125,*) "number of open   channels                       = ",n_open     
  write(125,*) "number of closed channels                       = ",n_closed
  write(125,*) "number of realizations                          = ",n_rlz
  write(125,*) "seed number                                     = ",seed

  !===== Parameter definitions ========
  !== M_se = number of self-embedding steps
  M_se=(M/se_step) !=(total number of scatterers)/(number of scatterers in a self-embedding step)
  write(125,*) "M_se (number of self-embedding steps)           = ", M_se
  Lmin = 0.5 * sqrt(L*W/M)   !== Minimum allowable distance between any two scatterers
  write(125,*) 'Average distance between scatterers             = ',sqrt(L*W/M)
  write(125,*) 'Minimum distance between scatterers             = ',Lmin
  write(125,*) "n_open                                          = ", n_open
  write(125,*) "Total # of channels nmax                        = ", nmax
  write(125,*) '************ end of initialization numbers **************'
  write(125,*) ' '

end subroutine inputvariables

!*****************************************************************************
!** OUTPUT FILES =======
!== note: if you open a file and write nothing to it, it's empty. Thus we 
!== avoid empty files by using if statements
!== FREQUENCY DEPENDENT only
!== CALLED BY main body
!== DEPENDS ON int2String
subroutine outputfiles(rank,ch_n,rlz,n_rlz,writeRunningStatus,writeFrequency,&
                         writeReflectionMatrix,writeTransmissionMatrix,&
                         writeTmRmDet,writeAccuracy, writeScatFreeMatrix, &
                         writeMedium,writeDeterminats,writeElecField, &   
                         writeSelfEmbedTmRm, writeGain, writeFinish,&
                         findElecField,writeABCoefficientsPW,writeABCoefficientsUF,writeElecFieldError,&
                         findJWDplanewave,findjwdunitflux,writeSummedJW,writeEachJW,pausable,&
                         percntwriteJWD,writeTEuf,writegEpw,&
                         writeTEufdist,writecriticalgain,writedistcriticalgain,&
                         writeconductionufdist,writeconductionpwdist,writeaveTab,&
                         writeufJWDchan,writepwJWDchan,findcriticalgain,writegEuf,find_dzy,dzy_gain_resolved)

  implicit none !== all of the variables should be of an explicit type
!  include 'mpif.h' !== MPIF90
  !== external function
  character*8 :: int2String
  !== arguments, scalar
  integer*4 :: rank
  integer :: rlz,&                      !== realization loop index, for filename
             n_rlz
  real*8 :: percntwriteJWD
  !== arguments, logical
  logical :: writeRunningStatus,&
             writeFrequency,&
             writeReflectionMatrix,& 
             writeTransmissionMatrix,& 
             writeTmRmDet,& 
             writeAccuracy,& 
             writeScatFreeMatrix,& 
             writeMedium,& 
             writeDeterminats,& 
             writeSelfEmbedTmRm,& 
             writeGain,& 
             writeFinish,&
             writeElecFieldError,&
             writeElecField,&
             findElecField,&
             writeABCoefficientsUF,&
             writeABCoefficientsPW,&
             findJWDplanewave,findjwdunitflux,&
             writeSummedJW,&
             writeEachJW,&
             pausable,&
             findcriticalgain,&
             writegEpw,&
             writeTEuf,&
             writegEuf,&
             writeTEufdist,&
             writecriticalgain,&
             writedistcriticalgain,&
             writeconductionufdist,&
             writeconductionpwdist,&
             writeaveTab,&
             writeufJWDchan,&
             writepwJWDchan,&
             find_dzy,&
             dzy_gain_resolved
  !== boolean "varyStrength" should not be included
  character*8              :: ch_n  !== note: if number of digits is increased, then variable size must change accordingly
  !== local, scalar
  character(len=70)        :: file_name

  writegEuf = .true.
!== file name conventions:
!== -output files start with "out_*"
!== -if unit flux OR plane wave output is available, specify "uf" or "pw" in the name
!== -if the output is a distribution, filename should start with "out_distribution(s)_*"
!== -if file output is every "percntwriteJWD", then the filename should have how many rlz that is
!== Note: it is desirable to use "percntwriteJWD" for some files (such as out_pw_gE.dat or out_uf_TE_chan*),
!== but if that is enabled, then only the realization where the percntwriteJWD condition is satisifed gets written

  !== formerly written to the terminal, now it gets its own output file
  if (findcriticalgain) then
     open(111,file='out_critical_gain_screen.dat',POSITION='APPEND')
  endif

  if (pausable) then
    file_name = "out_num_ran2_calls.dat"
    open(130,file=file_name,POSITION='APPEND')
  endif

  if (writeElecFieldError) then
     file_name = "out_"//int2string(rlz)//"_elecField_error.dat"
     open(136,file=file_name,POSITION='APPEND')
  endif

  if (writeABCoefficientsUF.OR.writeABCoefficientsPW) then
     file_name = "out_"//int2string(rlz)//"_ABCoefficients_fwd.dat"
     open(140,file=file_name,POSITION='APPEND')
!     file_name = "out_"//int2string(rlz)//"_ABCoefficients_bck.dat"
!     open(142,file=file_name,POSITION='APPEND')
  endif

  if (writeElecField) then
     file_name = "out_"//int2string(rlz)//"_elecField_fwd.dat"
     open(129,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2string(rlz)//"_running_status.dat"
  if (writeRunningStatus) then
     open(127,file=file_name,POSITION='APPEND')
  endif 

  file_name = "out_frequencies.dat"
  if(writeFrequency) then
     open(135,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2string(rlz)//"_MatrixReflection.dat"
  if (writeReflectionMatrix) then
     open(121,file=file_name,POSITION='APPEND')
  endif
  
  file_name = "out_"//int2string(rlz)//"_MatrixTransmission.dat"
  if (writeTransmissionMatrix) then
     open(122,file=file_name,POSITION='APPEND')
  endif
  
  file_name = "out_"//int2String(rlz)//"_detTmRm.dat"
  if (writeTmRmDet) then
     open(124,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2String(rlz)//"_accuracy.dat"
  if (writeAccuracy) then
     open(126,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2String(rlz)//"_MatrixScattering.dat"
  if (writeScatFreeMatrix) then
     open(131,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2String(rlz)//"_MatrixFree.dat"
  if (writeScatFreeMatrix) then
     open(138,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2String(rlz)//"_MatrixChunk.dat"
  if (writeScatFreeMatrix) then
     open(139,file=file_name,POSITION='APPEND')
  endif

  file_name = "out_"//int2String(rlz)//"_Medium.dat"
  if (writeMedium) then
    open (100, file=file_name,POSITION='APPEND')
  endif

  if (writeDeterminats) then
  file_name = "out_"//int2String(rlz)//"_detFreeScat.dat"
     open(101, file=file_name,POSITION='APPEND')
  file_name = "out_"//int2String(rlz)//"_detChunk.dat"
     open(102, file=file_name,POSITION='APPEND')
  endif

  if (writeSelfEmbedTmRm) then
  file_name = "out_"//int2string(rlz)//"_Tm.dat"
     open(103, file=file_name,POSITION='APPEND')
  file_name = "out_"//int2string(rlz)//"_Rm.dat"
     open(104, file=file_name,POSITION='APPEND')
  endif

  file_name = "out_gains.dat"
  if(writeGain) then
     open(137,file=file_name,POSITION='APPEND')
  endif

!   file_name = "out_"//int2string(rlz)//"_uf_JWD.dat"
!   if (writeEachJW) then
!      open(143,file=file_name,POSITION='APPEND')
!   endif

  if (writeSummedJW.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_uf_JWD_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(144,file=file_name,POSITION='APPEND')
!!     file_name="out_uf_JWD_averaged_over_"//int2string(rlz)//"_rlz.dat"
!!     open(145,file=file_name,POSITION='APPEND')
  endif

  if (writeTEuf) then !.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_uf_TE_channel_resolved.dat"
     open(146,file=file_name,POSITION='APPEND')
  endif

  if (writegEuf) then !.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_uf_gE.dat"
     open(161,file=file_name,POSITION='APPEND')
  endif

  if (writegEpw) then !.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_pw_gE.dat"
     open(141,file=file_name,POSITION='APPEND')
  endif

  if (writeTEufdist.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_distributions_"//int2string(rlz)//"_rlz_uf_TE_input_chan.dat"
     open(147,file=file_name,POSITION='APPEND')
  endif

  if (writecriticalgain.AND.findcriticalgain) then
     file_name="out_critical_gains.dat"
     open(148,file=file_name,POSITION='APPEND')
  endif

  if (writedistcriticalgain.AND.findcriticalgain.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_distribution_"//int2string(rlz)//"_rlz_critical_gains.dat"
     open(149,file=file_name,POSITION='APPEND')
  endif

  if (writeconductionufdist.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_distributions_"//int2string(rlz)//"_rlz_uf_gE.dat"
     open(150,file=file_name,POSITION='APPEND')
  endif

  !== note: there is not going to be an input channel resolved pw distribution,
  !== since g, E are calculated only after summing over input channels.
  if (writeconductionpwdist.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_distributions_"//int2string(rlz)//"_rlz_pw_gE.dat"
     open(160,file=file_name,POSITION='APPEND')
  endif

  if (writeaveTab) then
     file_name="out_aveTab.dat"
     open(151,file=file_name,POSITION='APPEND')
  endif

  if (writeufJWDchan.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_uf_J_channels_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(152,file=file_name,POSITION='APPEND')
!!     file_name="out_uf_J_channels_averaged_over_"//int2string(rlz)//"_rlz.dat"
!!     open(153,file=file_name,POSITION='APPEND')
     file_name="out_uf_W_channels_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(154,file=file_name,POSITION='APPEND')
!!     file_name="out_uf_W_channels_averaged_over_"//int2string(rlz)//"_rlz.dat"
!!     open(155,file=file_name,POSITION='APPEND')
  endif

  if (writepwJWDchan.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_pw_J_channels_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(156,file=file_name,POSITION='APPEND')
!!     file_name="out_pw_J_channels_averaged_over_"//int2string(rlz)//"_rlz.dat"
!!     open(157,file=file_name,POSITION='APPEND')
     file_name="out_pw_W_channels_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(158,file=file_name,POSITION='APPEND')
!!     file_name="out_pw_W_channels_averaged_over_"//int2string(rlz)//"_rlz.dat"
!!     open(159,file=file_name,POSITION='APPEND')
     file_name="out_pw_energy_corr_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(162,file=file_name,POSITION='APPEND')
  endif

  if (find_dzy.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
     file_name="out_pw_WJyz_summed_over_"//int2string(rlz)//"_rlz.dat"
     open(164,file=file_name,POSITION='APPEND')
  endif

  if (.true.) then
     file_name="out_eigenvalues.dat"
     open(166,file=file_name,POSITION='APPEND')
  endif

end subroutine outputfiles

!*****************************************************************************
!** Close all open files (end of program only)
!== CALLED BY main body
!== DEPENDS ON NONE
subroutine closeAllFiles(a)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer :: a  !== for some reason, a value must be passed. Cannot do ()

  close(100) !== out_%rlz%_Medium.dat
  close(101) !== out_%rlz%_detFree.dat
  close(102) !== out_%rlz%_detChunk.dat
  close(103) !== out_%rlz%_Tm.dat
  close(104) !== out_%rlz%_Rm.dat
  close(121) !== out_%rlz%_Reflection_matrix.dat
  close(122) !== out_%rlz%_Transmission_matrix.dat
  close(123) !== seed.input
  close(124) !== out_%rlz%_TmRm_det.dat
  close(125) !== out_screen.dat
  close(126) !== out_%rlz%_accuracy.dat
  close(127) !== out_running_status.dat
  close(128) !== quasi1d_rect_files.input
  close(129) !== out_%rlz%_elecField_fwd.dat
  close(130) !== out_num_ran2_calls.dat
  close(131) !== out_%rlz%_MatrixScattering.dat
  close(132) !== quasi1d_rect_parameters.input
  close(133) !== out_ERROR.log
  close(134) !== out_finished.status
  close(135) !== out_%rlz%_frequencies.dat
  close(136) !== out_%rlz%_elecField_error.dat
  close(137) !== out_%rlz%_gains.dat
  close(138) !== out_%rlz%_MatrixFree.dat
  close(139) !== out_%rlz%_MatrixChunk.dat
  close(140) !== out_%rlz%_ABCoefficients_fwd.dat
  close(141) !== out_pw_gE.dat
  close(142) !== out_%rlz%_ABCoefficients_bck.dat
  close(143) !== out_uf_%rlz%_JWD.dat
  close(144) !== out_uf_JWD_summed.dat
!!  close(145) !== out_uf_JWD_averaged_over_nrlz.dat
  close(146) !== out_uf_TE_channel_resolved.dat
  close(147) !== out_distributions_%rlz%_rlz_uf_TE_input_chan.dat
  close(148) !== out_%rlz%_critical_gain.dat
  close(149) !== out_%rlz%_critical_gain_distributions.dat
  close(150) !== out_distributions_%rlz%_rlz_uf_gE.dat
  close(151) !== out_aveTab.dat
  close(152) !== out_uf_J_channels_summed.dat
!!  close(153) !== out_uf_J_channels_averaged_over_nrlz.dat
  close(154) !== out_uf_W_channels_summed.dat
!!  close(155) !== out_uf_W_channels_averaged_over_nrlz.dat
  close(156) !== out_pw_J_channels_summed.dat
!!  close(157) !== out_pw_J_channels_averaged_over_nrlz.dat
  close(158) !== out_pw_W_channels_summed.dat
!!  close(159) !== out_pw_W_channels_averaged_over_nrlz.dat
  close(160) !== out_distributions_%rlz%_rlz_pw_gE.dat
  close(161) !== uf_gE
  close(162) !== correlationE_sum
  close(163) !== critical gain reasoning
  close(164) !== wdist_zy, jflux_y, jflux_z
  close(165) !== how write jw
  close(166) !== out_eigenvalues
end subroutine closeAllFiles

!**********************************
! J(z), W(z) to find D(z) for unit flux input
!== CALLED BY main body
!== DEPENDS ON int2String, flush
subroutine JWDunitflux(rank,ierr,vector,numscatterersplustwo,ntot,nmax,kpara,inputchan,pi,k2,W,&
                        scatr_indx,expon,L,ABstep,n_omg,fli,gli,n_open,writeufJWDchan,&
                        kparaMat, n_alphaa,alphaa,rlz,n_rlz,percntwriteJWD,M,Medium,&
                        writeEachJW,writeSummedJW,&
                        Jflux_r_sum,Jflux_l_sum,Wdist_r_sum,Wdist_l_sum,&
                        Wdist_r_n_sum,Wdist_l_n_sum,Jflux_r_n_sum,Jflux_l_n_sum,&
                        tebinsuf,condbinsuf,binsize,writeTEuf,writeTEufdist,gainbin,findcriticalgain,&
                        conduf,total_energyuf,total_passive_energyuf,writeconductionufdist)
  implicit none !== all of the variables should be of an explicit type
!  include 'mpif.h' !== MPIF90
  !== external functions
  character*8  :: int2String
  !== arguments,scalar
  integer*4  :: rank,ierr !== for MPI
  integer :: numscatterersplustwo,ntot,nmax,inputchan,n_omg,n_open,&
             fli,gli,n_alphaa,rlz,n_rlz,M,binsize
  real*8  :: pi,k2,L,ABstep,alphaa,percntwriteJWD
  real*8 :: W,Emin,Emax,Tmin,Tmax,Eamin,Eamax,Tamin,Tamax
  real*8 :: conduf,total_energyuf, total_passive_energyuf
  !== arguments, array
  complex*16  :: vector(1:ntot,1:numscatterersplustwo)
  integer :: scatr_indx(1:(int(L*ABstep)+1)),condbinsuf(1:8,1:binsize,1:n_alphaa)
  real*8 :: Kpara(nmax),kparaMat(1:nmax,1),energyPassiveuf(1:n_open)
  complex*16  :: expon(1:nmax,1:(int(L*ABstep)+1))
  real        :: Medium(M,2)
  real*8 :: Jflux_r_sum(1:(int(L*ABstep)+1),1:n_alphaa),&
            Jflux_l_sum(1:(int(L*ABstep)+1),1:n_alphaa),&
!            Wdist_sum(1:(int(L*ABstep)+1),1:n_alphaa),&
            Wdist_r_sum(1:(int(L*ABstep)+1),1:n_alphaa),&
            Wdist_l_sum(1:(int(L*ABstep)+1),1:n_alphaa),&
            Jflux_l_n_sum(1:n_open,1:n_open,1:(int(L*ABstep)+1),1:n_alphaa),&
            Jflux_r_n_sum(1:n_open,1:n_open,1:(int(L*ABstep)+1),1:n_alphaa),&
            Wdist_l_n_sum(1:n_open,1:nmax  ,1:(int(L*ABstep)+1),1:n_alphaa),&
            Wdist_r_n_sum(1:n_open,1:nmax  ,1:(int(L*ABstep)+1),1:n_alphaa)
  integer :: tebinsuf(1:8,1:binsize,1:n_open,1:n_alphaa)
  integer :: gainbin(1:binsize)
  !== arguments, logical
  logical :: writeufJWDchan,writeEachJW,writeSummedJW,writeconductionufdist,writeTEufdist,&
              writeTEuf,findcriticalgain,writegEuf
  !== local
  real*8  :: ufchannorm,k
  real*8 :: transuf,energyuf
  real   :: deltax
  complex*16  :: A_raw(1:nmax,1:numscatterersplustwo),&
                 B_raw(1:nmax,1:numscatterersplustwo),&
                 A_inside(1:nmax,1:(int(L*ABstep)+1)),&
                 B_inside(1:nmax,1:(int(L*ABstep)+1))
  real*8      :: Auf2(1:nmax,1:(int(L*ABstep)+1)),&
                 Buf2(1:nmax,1:(int(L*ABstep)+1)),&
                 Aufscat2(1:nmax,1:numscatterersplustwo),&
                 Bufscat2(1:nmax,1:numscatterersplustwo)
  integer :: a,b,c,d,f,h
  k = sqrt(k2)
  writegEuf = .true.
  ufchannorm = (1.0d+0)/(k*kpara(inputchan))**0.5
  !write(*,*) 'inputchan=',inputchan,' k=',k,' pi=',pi,'ufchannorm = ',ufchannorm

  A_raw(1:nmax,1:numscatterersplustwo) = (0.5)*( vector(1:nmax,1:numscatterersplustwo)-&
                      (0.0d+0,1.0d+0)*vector(nmax+1:2*nmax,1:numscatterersplustwo) )*(ufchannorm*(1.0d+0,0.0d+0))
  B_raw(1:nmax,1:numscatterersplustwo) = (0.5)*( vector(1:nmax,1:numscatterersplustwo)+&
                      (0.0d+0,1.0d+0)*vector(nmax+1:2*nmax,1:numscatterersplustwo) )*(ufchannorm*(1.0d+0,0.0d+0))
  !== NOTE: including the (1,0) for pwchan norm does not affect A,B_raw. Both are still complex

  !== find distributions of T_a,E_a,g,E

  if (writeTEuf.OR.writegEuf.OR.writeTEufdist.OR.writeconductionufdist) then
     Aufscat2 = cdabs(A_raw)**2
     Bufscat2 = cdabs(B_raw)**2
  
     do c = 1,numscatterersplustwo !== for every scatterer
    
        if (c.eq.1) then
           energyuf=0   
           DeltaX = Medium(1,1)
!           write(*,*) 'energyuf,k2,W,Deltax',energyuf,k2,W,Deltax
           energyuf = energyuf + k2*sum(Aufscat2(1:nmax,c)+Bufscat2(1:nmax,c))*W*DeltaX
        elseif (c.lt.numscatterersplustwo-1) then
           DeltaX = Medium(c,1) - Medium(c-1,1)
           if (DeltaX.lt.0) then
              write(133,*) rlz,gli,fli,DeltaX,c
           endif
!           write(*,*) 'energyuf,k2,W,Deltax',energyuf,k2,W,Deltax
           energyuf = energyuf + k2*sum(Aufscat2(1:nmax,c)+Bufscat2(1:nmax,c))*W*DeltaX
        elseif (c.eq.numscatterersplustwo-1) then
           DeltaX = L-Medium(M,1)
           energyuf = energyuf + k2*sum(Aufscat2(1:nmax,c)+Bufscat2(1:nmax,c))*W*DeltaX
        else
           !== note: there is no energyDensity associated with c=numscatterersplustwo
           !== only need to find conductance once
           transuf = k*sum(Aufscat2(1:n_open,numscatterersplustwo)*kparaMat(1:n_open,1))

           if (sum(Bufscat2(1:n_open,numscatterersplustwo)).gt.(1E-10)) then
              open(133,file='out_ERROR.log',POSITION='APPEND')
              write(133,*) 'WARNING: incident flux on right side !?'
              close(133)
              write(125,*) 'WARNING: incident flux on right side !?'
              close(125)
              write(*,*) 'WARNING: incident flux on right side !?'
          !if (rlz.lt.(int(.5*n_rlz))) then
!*              call closeAllFiles(1)
!*              stop ! die
          !endif
           endif

           !== reality check
           if ((transuf.lt.0).OR.(energyuf.lt.0)) then !.AND.(alphaa.ge.0)) then 
              open(133,file='out_ERROR.log',POSITION='APPEND')
              write(133,*) 'WARNING: E or T <0',transuf,energyuf,rlz,fli,gli
              close(133)
              write(125,*) 'WARNING: E or T <0',transuf,energyuf,rlz,fli,gli
              close(125)
              write(*,*) 'WARNING: energy or transmission Less than 0',transuf,energyuf,fli,gli
              !if (rlz.lt.(int(.5*n_rlz))) then
!*              call closeAllFiles(1)
!*              stop !== die
              !endif
           endif
            !== normalize the energy per channel by empty waveguide
           energyuf = energyuf/(k*L*W/kparamat(inputchan,1))
            if (writeTEuf) then !.AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz)))) then
              write(146,*) transuf,energyuf,gli,fli  !== per channel
           endif
           if (writeTEufdist.OR.writegEuf.OR.writeconductionufdist) then
              if (alphaa.eq.0) then
                 energyPassiveuf(inputchan) = energyuf
              endif
              !== the following distribution is input channel resolved [unit flux only]
!              call distributionsTEuf(Energyuf,transuf,EnergyPassiveuf(inputchan),tebinsuf,binsize,inputchan,n_open,n_alphaa,gli)  
              if (inputchan.eq.1) then
                 conduf = transuf
                 total_energyuf = energyuf
              else                 
                 conduf = conduf + transuf
                 total_energyuf = total_energyuf + energyuf
              endif
              if (alphaa.eq.0) then
                 total_passive_energyuf = total_energyuf
              endif
              if ((inputchan.eq.n_open).AND.writegEuf) then
                 write(161,*) conduf,total_energyuf,gli,fli
              endif
!               if ((inputchan.eq.n_open).AND.writeconductionufdist) then
!                  !== the following distribution is summed over input channels [unit flux case]
!                  call distributionsgEuf(conduf,total_energyuf,total_passive_energyuf,condbinsuf,binsize,n_alphaa,gli)
!               endif
           endif
        endif  !== c
     enddo !== c
  endif !== writeTE

  A_inside(1:nmax,1:(int(L*ABstep)+1)) = &
       A_raw(1:nmax,scatr_indx(1:(int(L*ABstep)+1)))*expon(1:nmax,1:(int(L*ABstep)+1))
  B_inside(1:nmax,1:(int(L*ABstep)+1)) = &
       B_raw(1:nmax,scatr_indx(1:(int(L*ABstep)+1)))*(1/expon(1:nmax,1:(int(L*ABstep)+1)))  !== this should be element-wise inversion. Tested & works

  if(.true.) then
    !open(933,file="out_AB_freq_"//int2String(fli)//"_chan_"//int2String(inputchan)//".dat",POSITION='APPEND')
    open(933,file="out_AB_freq_"//int2String(fli)//".dat",POSITION='APPEND')
    do a=1,(int(L*ABstep)+1)
      write(933,"("//int2String(4*nmax)//"ES24.15)") dreal(A_inside(1:nmax,a)), dimag(A_inside(1:nmax,a)), &
                                                      dreal(B_inside(1:nmax,a)), dimag(B_inside(1:nmax,a))
    enddo
  endif

  Auf2(1:nmax,1:(int(L*ABstep)+1))=cdabs(A_inside)**2
  Buf2(1:nmax,1:(int(L*ABstep)+1))=cdabs(B_inside)**2 
!==    J_r = J_
!==    J_l = J_
!==    J = J_r - J_l
!==    J = k c \sum_{n=1}^{N_{open}} (|A_n|^2 - |B_n|^2) k_{\para_n}
!==    W = k^2 \sum_{n=1}^{N_{open}+N_{closed}} (|A_n|^2 + |B_n|^2)
  do a = 1,(int(L*ABstep)+1)

     if ((alphaa.eq.0.0).AND.(writeTEufdist.OR.writeconductionufdist).AND.& !== passive only, need transuf
        (k*(sum(Auf2(1:n_open,a)*kparaMat(1:n_open,1))-sum(Buf2(1:n_open,a)*kparaMat(1:n_open,1)))-transuf).gt.(1E-7)) then 
        open(133,file='out_ERROR.log',POSITION='APPEND')
        write(133,*) 'WARNING: J is not constant within the sample ',gli,fli,rlz,inputchan,a!,&
         !(sum(Auf2(1:n_open,a)*kparaMat(1:n_open,1))-sum(Buf2(1:n_open,a)*kparaMat(1:n_open,1)))-transuf,transuf
        close(133)
        write(125,*) 'WARNING: J is not constant within the sample ',gli,fli,rlz,inputchan,a!,&
         !(sum(Auf2(1:n_open,a)*kparaMat(1:n_open,1))-sum(Buf2(1:n_open,a)*kparaMat(1:n_open,1)))-transuf,transuf
        close(125)
        write(*,*) 'WARNING: J is not constant within the sample ',gli,fli,rlz,inputchan,a,&
         (sum(Auf2(1:n_open,a)*kparaMat(1:n_open,1))-sum(Buf2(1:n_open,a)*kparaMat(1:n_open,1)))-transuf,transuf
!*        call closeAllFiles(1)
!*        stop ! die
!     elseif ((alphaa.eq.0.0).AND.(writeTEufdist.OR.writeconductionufdist)) then
!        write(*,*) 'passed',a,inputchan
       !== by passing this check, we can tell that the error occurs either 
        !==in the sum over incident channels, 
        !==in the sum over realizations, 
        !==or in writing to file
     endif
     if (writeEachJW.OR.writeSummedJW) then
        Jflux_r_sum(a,gli) = Jflux_r_sum(a,gli) + k*sum(Auf2(1:n_open,a)*kparaMat(1:n_open,1))
        Jflux_l_sum(a,gli) = Jflux_l_sum(a,gli) + k*sum(Buf2(1:n_open,a)*kparaMat(1:n_open,1))
        !Wdist_sum(a,gli)   = Wdist_sum(a,gli)   + k2*sum(Auf2(1:nmax,a)+Buf2(1:nmax,a))
        Wdist_r_sum(a,gli) = Wdist_r_sum(a,gli)   + k2*sum(Auf2(1:nmax,a))
        Wdist_l_sum(a,gli) = Wdist_l_sum(a,gli)   + k2*sum(Buf2(1:nmax,a))
     endif
     if (writeufJWDchan) then
  Jflux_r_n_sum(inputchan,1:nmax,a,gli)=Jflux_r_n_sum(inputchan,1:nmax,a,gli)+k*Auf2(1:nmax,a)*kparaMat(1:nmax,1)
  Jflux_l_n_sum(inputchan,1:nmax,a,gli)=Jflux_l_n_sum(inputchan,1:nmax,a,gli)+k*Buf2(1:nmax,a)*kparaMat(1:nmax,1)
  Wdist_r_n_sum(inputchan,1:nmax,a,gli)=Wdist_r_n_sum(inputchan,1:nmax,a,gli)+k2*Auf2(1:nmax,a)
  Wdist_l_n_sum(inputchan,1:nmax,a,gli)=Wdist_l_n_sum(inputchan,1:nmax,a,gli)+k2*Buf2(1:nmax,a)
     endif

     if ((alphaa.eq.0.0).AND.(a.gt.4).AND.(writeTEufdist.OR.writeconductionufdist).AND.& !== passive only, need transuf
   (int(10000*(Jflux_r_sum(a,gli) - Jflux_l_sum(a,gli))).ne.int(10000*(Jflux_r_sum(4,gli) - Jflux_l_sum(4,gli))))) then
        open(133,file='out_ERROR.log',POSITION='APPEND')
        write(133,*) 'WARNING: J_sum is not constant within the sample ',gli,fli,rlz,inputchan,a,&
         (Jflux_r_sum(a,gli) - Jflux_l_sum(a,gli)),transuf
        close(133)
        write(125,*) 'WARNING: J_sum is not constant within the sample ',gli,fli,rlz,inputchan,a,&
         (Jflux_r_sum(a,gli) - Jflux_l_sum(a,gli)),transuf
       ! close(125)
        write(*,*) 'WARNING: J_sum is not constant within the sample ',gli,fli,rlz,inputchan,a,&
         (Jflux_r_sum(a,gli) - Jflux_l_sum(a,gli)),(Jflux_r_sum(4,gli) - Jflux_l_sum(4,gli))
       ! call closeAllFiles(1)
       ! stop ! die
       !== if no failure, then J was constant throughout the sample
     endif
  enddo  !== position a

  !== reality check: compares incident to exiting fluxes (T+R=1)
 ! if ((inputchan.eq.n_open).AND.(alphaa.eq.(0.0)).AND.&  !== only works in passive media
  if ((alphaa.eq.(0.0)).AND.&  !== only works in passive media
!     .AND.((Jflux_l(1,gli)+Jflux_r((int(L*ABstep)+1),gli)-k).gt.(1E-10))) then 
  ((sum(Buf2(1:n_open,1)*kparaMat(1:n_open,1)) + sum(Auf2(1:n_open,(int(L*ABstep)+1))*kparaMat(1:n_open,1)) -1).gt.(1E-10))) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: Jflux_l(x=0) + Jflux_r(x=L) != k ',gli,fli,rlz,inputchan
     close(133)
     write(125,*) 'WARNING: Jflux_l(x=0) + Jflux_r(x=L) != k ',gli,fli,rlz,inputchan
     close(125)
     write(*,*) 'WARNING: Jflux_l(x=0) + Jflux_r(x=L) != k ',gli,fli,rlz,inputchan
     !if (rlz.lt.(int(.5*n_rlz))) then
!*     call closeAllFiles(1)
!*       stop ! die
     !endif
  endif

  !write to file once per realization
!   if ((fli.eq.n_omg).AND.(gli.eq.n_alphaa).AND.writeEachJW) then
!      do c = 1,n_alphaa
!         do a = 1,(int(L*ABstep)+1)  !== the following statement is written to file n_open times
!            write(143,*) Wdist(a,c), Jflux_r(a,c), Jflux_l(a,c)
!         enddo  
!      enddo
!   endif  

  !== NOTE: percntwriteJWD*n_rlz must occasionally be an integer!
  if ((fli.eq.n_omg).AND.(gli.eq.n_alphaa).AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz))) &
        .AND.(inputchan.eq.n_open).AND.writeSummedJW) then  !== write J_+, J_-, W to file
     !d = n_rlz*n_open*n_omg !== input channel, realization
!      do c = 1,n_alphaa  
!            !Jflux(:)=Jflux_r(:,b,c)-Jflux_l(:,b,c)  !== do not write to file, save disk space and write time
!         do a = 1,(int(L*ABstep)+1)
!            write(144,'(4ES24.15)') Wdist_l_sum(a,c),Wdist_r_sum(a,c),Jflux_l_sum(a,c),Jflux_r_sum(a,c)   !== sum 
!            !== 20090901: Ben removed /(n_open*1.0)
! !!  write(145,*) Wdist_sum(a,c)/(d*1.0),     Jflux_r_sum(a,c)/(d*1.0),     Jflux_l_sum(a,c)/(d*1.0)  !== average
!         enddo  
!      enddo  
  endif

  if ((fli.eq.n_omg).AND.(gli.eq.n_alphaa).AND.(0.eq.mod(rlz,int(percntwriteJWD*n_rlz))) &
        .AND.(inputchan.eq.n_open).AND.writeufJWDchan) then !== write J_+, J_-, W_+,W_- for each channel to file
!      do c = 1,n_alphaa
!         do f = 1,n_open !== input
!            do h = 1,n_open  !== interior
!               do a = 1,(int(L*ABstep)+1)
!                  write(152,*) Jflux_l_n_sum(f,h,a,c),Jflux_r_n_sum(f,h,a,c)  !== sum
! !!  write(153,*) Jflux_r_n_sum(f,h,a,c)/(n_rlz*n_omg*1.0),Jflux_l_n_sum(f,h,a,c)/(n_rlz*n_omg*1.0)  !== average
!               enddo  
!            enddo  
!            do h = 1,nmax  !== interior
!               do a = 1,(int(L*ABstep)+1)
!                  write(154,*) Wdist_l_n_sum(f,h,a,c),Wdist_r_n_sum(f,h,a,c)  !== sum
! !!  write(155,*) Wdist_r_n_sum(f,h,a,c)/(n_rlz*n_omg*1.0),Wdist_l_n_sum(f,h,a,c)/(n_rlz*n_omg*1.0)  !== average
!               enddo  
!            enddo  
!         enddo
!      enddo
  endif
  call flush(143)
  call flush(144)
  call flush(145)
  call flush(152)
  call flush(153)
  call flush(154)
  call flush(155)
end subroutine JWDunitflux

!********************
! distribution of critical gain values, not frequency resolved
!== for creating histogram of critical gain
!== CALLED BY 
!== DEPENDS ON NONE
subroutine distributionscriticalgain(cg,gainbin,binsize)
  implicit none
  !== argument, scalar
  integer :: binsize
  real*8 :: cg
  !== argument, array
  integer :: gainbin(1:binsize)
  !== local, scalar
  integer :: binmax,binmin,indx,bininc
  real*8 :: cgmax,cgmin,&
            stepsize !,step
!  write(*,*) 'critical gain is ',cg
  cgmax = 0 
  cgmin = -1
  stepsize = ((cgmax-cgmin)*(1.0d+0))/(binsize*(1.0d+0))
!   do indx = 1,binsize
!      step = cgmin + indx*stepsize
!      if ((cg.lt.(step+stepsize)).AND.(cg.gt.step)) then
!         gainbin(indx) = gainbin(indx)+1
!      endif
!   enddo
  bininc = int((cg-cgmin)/stepsize)
  if (bininc.lt.1) then
    bininc=1  !== if below the limits of the distribution, increment the first bin
  elseif (bininc.gt.binsize) then
    bininc=binsize !== if beyond the limits of the distribution, increment the last bin
  endif
  gainbin(bininc) = gainbin(bininc)+1

end subroutine distributionscriticalgain

!== convert integer to string
!== CALLED BY main body, inputvariables
!== DEPENDS ON NONE   
function int2string(inputinteger)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer :: inputinteger
  !== return value, scalar
  character*8 :: int2string
  !== local scalar
  integer*4  :: chr0,chr1,chr2,chr3,chr4,chr5,chr6,chr7

  if (inputinteger.lt.0) then
     open(133,file='out_ERROR.log',POSITION='APPEND')
     write(133,*) 'WARNING: someone tried to pass a negative integer to int2string()',inputinteger
     close(133)
     write(125,*) 'WARNING: someone tried to pass a negative integer to int2string()',inputinteger
     close(125)
     write(*,*) 'WARNING: someone tried to pass a negative integer to int2string()',inputinteger
     !call mpi_finalize(ierr)  !== to enable these, need to pass ierr to function
     !call mpi_abort(MPI_COMM_WORLD,ierr)
!*     call closeAllFiles(1)
!*     stop ! die
  endif

  !== works, but does not scale well
  chr0=int(    int(1.0*inputinteger+0.1)/10000000)
  chr1=int(mod(int(1.0*inputinteger+0.1),10000000)/1000000)
  chr2=int(mod(int(1.0*inputinteger+0.1),1000000)/100000)
  chr3=int(mod(int(1.0*inputinteger+0.1),100000)/10000)
  chr4=int(mod(int(1.0*inputinteger+0.1),10000)/1000)
  chr5=int(mod(int(1.0*inputinteger+0.1),1000 )/100)
  chr6=int(mod(int(1.0*inputinteger+0.1),100  )/10)
  chr7=    mod(int(1.0*inputinteger+0.1),10   )

  !== append all four characters to a string
  int2string=char(chr0+48)//char(chr1+48)//char(chr2+48)//&
               char(chr3+48)//char(chr4+48)//char(chr5+48)//&
               char(chr6+48)//char(chr7+48)

  !== alternative
!  write(outputstring,'(i8.8)') inputinteger 
!  write(*,'(a)') trim(outputstring) ! should output "00000455" 
  return
end function int2String

!*****************************
!== given a separation between scatterers, create the "free matrix" that propagates
!== electric field and its derivative over dz. Matrix values are real
!== CALLED BY main body
!== DEPENDS ON NONE
subroutine freematrixcreation(FreeMatrix,ntot,nmax,Kpara,DeltaX,n_open)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer         :: ntot,&
                     nmax,&
                     n_open
  real :: deltax  !== distance between two scatterers
  !== argument, array
  real*8   :: FreeMatrix(1:ntot,1:ntot),&
                     Kpara(nmax)
  !== local, scalar
  integer ::    g,&      !== local: loop index
                     p,&
                     q,&
                     Remainder
  real*8          ::  T  !== argument of trig functions

  FreeMatrix(:,:) = 0.0d+0
  do g=1, ntot  !== fills out the "FreeMatrix"
     Remainder = mod(g,nmax)
     if (Remainder.eq.0) then
        Remainder = nmax
     endif
     T = Kpara(Remainder) * DeltaX       !
     p = g + nmax                        !makes code easier
     q = g - nmax                        !
     if (g.le.n_open) then
        FreeMatrix(g,g) = COS(T)
        FreeMatrix(g,p) = SIN(T)
     else if (g.le.nmax) then
        FreeMatrix(g,g) = COSH(T)
        FreeMatrix(g,p) = SINH(T)
     else if (g.le.nmax+n_open) then
        FreeMatrix(g,q) =-SIN(T)
        FreeMatrix(g,g) = COS(T)
     else
        FreeMatrix(g,q) = SINH(T)
        FreeMatrix(g,g) = COSH(T)
     endif
  enddo !== end of the "g" index
end subroutine freematrixcreation

!== for a specific scatterer coordinate, create a scattering matrix that propagates
!== electric field and scrambles derivative of electric field into all channels
!== CALLED BY main body
!== DEPENDS ON ran2
subroutine scatmatrixcreation(ScatteringMatrix,varyStrength,seed,&
                              alphasign,seed_counter,nmax,ntot,alphas,&
                              alphaa,omg,W,Kpara,Kperp,Medium,M,a,firstpropagation)
  implicit none !== all of the variables should be of an explicit type
  !== external function
  real          :: ran2
  !== argument, logical
  logical       :: varyStrength,firstpropagation
  !== argument, scalar
  integer       :: seed,&
                   seed_counter,&
                   nmax,&
                   ntot,&
                   M,&   !== number of scatterers
                   alphasign(M)
  real*8        :: alphas,& !== strength
                   alphaa,& !== active component
                   omg,&
                   W
  !== argument, array
  real ::           Medium(M,2)
  complex*16    :: ScatteringMatrix(1:ntot,1:ntot)
  real*8 ::         Kpara(nmax),&
                   Kperp(nmax)
  !== local, scalar
  integer ::     c,& !== loop index for scatterer elements
                   a,& !== scatterer number
                   d !== loop index for scatterer elements

  !== will scattering strength be constant and non-zero average (varyStrength == false)
  !== OR (+/-)alpha  so that average is zero (vacuum-like refractive index) [non-physical for light!]
  if (varyStrength.AND.firstpropagation) then   !== varying alpha sign is in effect --> electronic
     !== SIGN( A, B ) returns the absolute value of A times the sign of B.
     alphasign(a) = sign(1.0,2*ran2(seed)-1)
     seed_counter = seed_counter +1  !== note that the seed was called. Necessary to track for pausing
     !write(*,*) alphasign(a), a
  elseif (.NOT.varyStrength) then  
     !== varying strength is not in effect --> photonic
     alphasign = -1
  endif
  do c=nmax+1, ntot
     !do d=1, nmax
        ScatteringMatrix(c,1:nmax) = ((alphasign(a)*alphas*(1.0d+0))-(0.0d+0,1.0d+0)*alphaa)*&
                                (omg**2)*2/W/Kpara(c-nmax)*&
                                SIN(Kperp(c-nmax)*Medium(a,2))*&
                                SIN(Kperp(1:nmax)*Medium(a,2))
     !enddo
  enddo !== end "c" indexed loop 
end subroutine scatmatrixcreation

!*****************************************************************************
!** Calculation of determinant of real matrix by LU decomposition **
!== CALLED BY main body
!== DEPENDS ON ludcmp
function determinant(m,size_of_m)  ! returns a real scalar
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer :: size_of_m
  !== argument, array
  real*8  :: m(1:size_of_m,1:size_of_m)
  !== local, scalar
  integer :: j,ID,ICODE
  real*8  :: determinant
  !== local, array
  integer :: INDX(1:size_of_m)
  real*8  :: A(1:size_of_m,1:size_of_m)

  determinant = 0
  A=m
  call LUDCMP(A,size_of_m,INDX,ID,ICODE)  
  determinant=dfloat(ID)
  do j=1,size_of_m
     determinant=determinant*A(j,j)
  enddo
  ! The return value should be stored in a variable with the same name as the function.
  ! Functions are terminated by the return statement instead of stop. 
  return
end function determinant

!*****************************************************************************
!** Calculation of determinant of complex matrix by LU decomposition **
!********************************************************************
!== CALLED BY 
!== DEPENDS ON ludcmpc
function determinantc(m,size_of_m)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer     :: size_of_m
  !== argument, array
  complex*16  :: m(1:size_of_m,1:size_of_m)
  !== local, scalar
  integer :: j,ID,ICODE
  complex*16  :: determinantc
  !== local, array
  integer     :: INDX(1:size_of_m)
  complex*16  :: A(1:size_of_m,1:size_of_m)
  determinantc = 0
  A=m
  call LUDCMPC(A,size_of_m,size_of_m,INDX,ID,ICODE)  
  determinantc=dcmplx(dfloat(ID),(0.0d+0))
  do j=1,size_of_m
     determinantc=determinantc*A(j,j)
  enddo
  return
end function determinantc

!== called during self-embedding
!== CALLED BY main body
!== DEPENDS ON  lubksbc, ludcmpc
subroutine inverse_c_matrix(m,size_of_m,mi)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  integer    :: size_of_m
  !== argument, array
  complex*16 :: mi(1:size_of_m,1:size_of_m)
  complex*16 :: m (1:size_of_m,1:size_of_m)
  !== local, scalar
  integer :: j,ID,ICODE
  !== local, scalar
  integer    :: INDX(1:size_of_m)
  complex*16 :: A (1:size_of_m,1:size_of_m)
  mi=(0.0d+0)
  do j=1,size_of_m
     mi(j,j)=(1.0d+0)
  enddo
  A=m
  call LUDCMPC(A,size_of_m,size_of_m,INDX,ID,ICODE)  
  do j=1,size_of_m
     mi(j,j)=(1.0d+0)
     call LUBKSBC(A,size_of_m,size_of_m,INDX,mi(1,j))
  enddo
  !m = mi
end subroutine inverse_c_matrix


!****************************************************************************
!** LU decomposition of a real matrix **
!*************************************
!== the reason the code is formatted using CAPS is probably for F77 compatibility
!== CALLED BY 
!== DEPENDS ON NONE
Subroutine LUDCMP(A,N,INDX,D,CODE)
  implicit none !== all of the variables should be of an explicit type
  INTEGER  :: N,NMAX
  REAL*8   :: TINY
  PARAMETER(NMAX=1000,TINY=1.5D-16)
  REAL*8   :: AMAX,DUM, SUM, A(N,N),VV(NMAX)
  INTEGER  :: CODE, D, INDX(N)
  INTEGER  :: I,IMAX,J,K
  D=1; CODE=0
  DO I=1,N
     AMAX=0.d0
     DO J=1,N
        IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
     enddo ! j loop
     IF(AMAX.LT.TINY) THEN
        CODE = 1
        RETURN
     endif
     VV(I) = 1.d0 / AMAX
  enddo ! i loop
  DO J=1,N
     DO I=1,J-1
        SUM = A(I,J)
        DO K=1,I-1
           SUM = SUM - A(I,K)*A(K,J) 
        enddo ! k loop
        A(I,J) = SUM
     enddo ! i loop
     AMAX = 0.d0
     DO I=J,N
        SUM = A(I,J)
        DO K=1,J-1
           SUM = SUM - A(I,K)*A(K,J) 
        enddo ! k loop
        A(I,J) = SUM
        DUM = VV(I)*DABS(SUM)
        IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
        endif
     enddo ! i loop  
     IF(J.NE.IMAX) THEN
        DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
        enddo ! k loop
        D = -D
        VV(IMAX) = VV(J)
     endif
     INDX(J) = IMAX
     IF(DABS(A(J,J)) < TINY) THEN 
        A(J,J) = TINY 
     ENDIF
     IF(J.NE.N) THEN
        DUM = 1.d0 / A(J,J)
        DO I=J+1,N
           A(I,J) = A(I,J)*DUM
        enddo ! i loop
     endif
  enddo ! j loop
  RETURN
END subroutine LUDCMP

!*****************************************************************************
!** some LU matrix trickery (?) **
!** called by inverse_c_matrix()
!*********************************
!== the reason the code is formatted using CAPS is probably for F77 compatibility
!== CALLED BY 
!== DEPENDS ON NONE
SUBROUTINE LUBKSB(A,N,NP,INDX,B)
  implicit none !== all of the variables should be of an explicit type
  integer :: N,NP,INDX(N)
  real*8  :: A(NP,NP),B(N)
  real*8  :: SUM
  integer :: I,II,J,LL
  II=0
  DO I=1,N
     LL=INDX(I)
     SUM=B(LL)
     B(LL)=B(I)
     IF (II.NE.0)THEN
        DO J=II,I-1
           SUM=SUM-A(I,J)*B(J)
        enddo
     ELSE IF (SUM.NE.0.) THEN
        II=I
     ENDIF
     B(I)=SUM
  enddo
  DO  I=N,1,-1
     SUM=B(I)
!     IF(I.LT.N)THEN
        DO J=I+1,N
           SUM=SUM-A(I,J)*B(J)
        enddo
!     ENDIF
     B(I)=SUM/A(I,I)
  enddo
  RETURN
END SUBROUTINE LUBKSB

!*****************************************************************************
! LUDCMPC(A,N,NP,INDX,D) **
!*************************************
!     Given a general complex matrix A, this routine replaces it by its 
!     LU decomposition of a row-wise permutation of itself.
!     This routine is used in combination with LUBKSBC(), a complex
!     extension of the routine LUBKSB() (DOUBLE COMPLEX).
!     For further details, refer to textbook (see below).
!     ========================================================
!     Source: Own adaption/extension to complex matrix of the
!             Subroutine LUDCMP() taken from
!             Press et al, "Numerical Recipes in Fortran"
!             The adaption follows the statements given in section 2.3
!             of the textbook "N.R. in C", following Eq.(2.3.16): 
!             - definition of variables, vector and matrix elements
!               as complex variables (use of complex arithmetic does
!               not necessitate any adaption in fortran).
!             - complex modulus instead of absolute values in the
!               construction of the vector vv and in the search for the
!               largest pivot elements.
!
!== see http://books.google.com/books?id=lZJD2zYGbNwC
!        http://books.google.com/books?id=1aAOdzK3FegC&pg=PA52&lpg=PA52&dq=%22lu+decomposition+of+a+rowwise+permutation+of+itself%22&source=web&ots=3gVoKeGrid&sig=K89Yq8qoHiaIgUzEEHEiFwHsLg4
!== page 52
!
!  Version: 28.08.2000   
!====================================================================
!== the reason the code is formatted using CAPS is probably for F77 compatibility
!== CALLED BY 
!== DEPENDS ON NONE
SUBROUTINE LUDCMPC(A,N,NP,INDX,D,CODE)
  implicit none !== all of the variables should be of an explicit type
  INTEGER    :: N,NP,NMAX
  REAL*8     :: TINY
  PARAMETER(NMAX=1000,TINY=1.5D-16)
  COMPLEX*16 :: A(NP,NP),SUM,DUMC,TINYC
  INTEGER    :: INDX(N),D,CODE
  REAL*8     :: VV(NMAX),DUM,AAMAX
  INTEGER    :: I,J,K,IMAX
  D=1; CODE=0
  TINYC=DCMPLX(TINY,0.D0)
  DO I=1,N
     AAMAX=0.
     DO J=1,N
        IF (CDABS(A(I,J)).GT.AAMAX) AAMAX=CDABS(A(I,J))
     enddo
     IF (AAMAX < TINY) then
        CODE = 1
        RETURN
     endif
     VV(I)=1./AAMAX
  enddo
  DO J=1,N
     DO I=1,J-1
        SUM=A(I,J)
        IF (I.GT.1)THEN
           DO K=1,I-1
              SUM=SUM-A(I,K)*A(K,J)
           enddo
           A(I,J)=SUM
        ENDIF
     enddo
     AAMAX = 0.d0
     DO I=J,N
        SUM=A(I,J)
        DO K=1,J-1
           SUM=SUM-A(I,K)*A(K,J)
        enddo
        A(I,J)=SUM
        DUM=VV(I)*CDABS(SUM)
        IF (DUM.GE.AAMAX) THEN
           IMAX=I
           AAMAX=DUM
        ENDIF
     enddo
     IF (J.NE.IMAX)THEN
        DO K=1,N
           DUMC=A(IMAX,K)
           A(IMAX,K)=A(J,K)
           A(J,K)=DUMC
        enddo
        D = -D
        VV(IMAX)=VV(J)
     ENDIF
     INDX(J)=IMAX
     IF(J.NE.N)THEN
        IF(CDABS(A(J,J)) < TINY) THEN 
           A(J,J)=TINYC 
        ENDIF
        DUMC=1./A(J,J)
        DO I=J+1,N
           A(I,J)=A(I,J)*DUMC
        enddo
     ENDIF
  enddo
  IF(CDABS(A(N,N))<TINY)THEN
     A(N,N)=TINYC
  endif
  RETURN
END SUBROUTINE LUDCMPC

!*****************************************************************************
!     Solution of the set of linear equations A' * x = b where 
!     A is input not as the original matrix, but as a LU decomposition 
!     of some original matrix A' as determined by the subroutine
!     LUDCMPC() (matrix and vectors being of type DOUBLE COMPLEX). 
!     INDX() is input as the permutation vector returned by LUDCMPC().
!     B() is input as the right hand side vector b of the Eqn. to solve
!     and returns with the solution vector x. 
!     A, N and INDX are not modified by this routine and can be left in
!     place for successive calls with different right-hand-sides b.
!     For further details, refer to textbook (see below).
!================================================================
!     Source: Own adaption/extension to complex matrix of the
!             Subroutine LUBKSB() taken from
!             Press et al, "Numerical Recipes in Fortran"
!             The adaption follows the statements given in section 2.3
!             of the textbook "N.R. in C", following Eq.(2.3.16).
!
!== see http://books.google.com/books?id=lZJD2zYGbNwC&pg=PA1017
!
!Version: 28.08.2000
!===============================================================
!== the reason the code is formatted using CAPS is probably for F77 compatibility
!== CALLED BY 
!== DEPENDS ON NONE
SUBROUTINE LUBKSBC(A,N,NP,INDX,B)
  implicit none !== all of the variables should be of an explicit type
  !== argument, scalar
  INTEGER    :: N,NP
  !== argument, array
  integer :: INDX(N)
  COMPLEX*16 :: A(NP,NP),B(N)
  COMPLEX*16 :: SUM,ZEROC
  !== local, scalar
  INTEGER    :: II,LL,I,J
  II=0
  ZEROC=DCMPLX(0.D0,0.D0)
  DO I=1,N
     LL=INDX(I)
     SUM=B(LL)
     B(LL)=B(I)
     IF (II.NE.0)THEN
        DO J=II,I-1
           SUM=SUM-A(I,J)*B(J)
        enddo
     ELSE IF (SUM.NE.ZEROC) THEN
        II=I
     ENDIF
     B(I)=SUM
  enddo
  DO I=N,1,-1
     SUM=B(I)
!     IF(I.LT.N)THEN
        DO J=I+1,N
           SUM=SUM-A(I,J)*B(J)
        enddo
!     ENDIF
     B(I)=SUM/A(I,I)
  enddo
  RETURN
END SUBROUTINE LUBKSBC

!== end of subroutines section

!== end of Fortran90 code
!== start of embedded comments

!=============================================================================
!== Quasi 1-d random media simulator
!== Written by:
!== Fortran90
!== Authors: Alexey Yamilov, Tom Mahler, Ben Payne
!== Last updated 20091204
!== Purpose: model quasi 1-D rectangular media for light incident from left 
!== side, propagating to the right using transfer matrix method and 
!== self-embedding technique with delta-function scatterers, randomly placed
!== 
!== addition: originally the delta-function scatterer was passive (Real alpha)
!== now alpha can complex, allowing for gain and absorption ("active" media) 
!== Caveat: Kramer-Kronig will come into play with frequency
!== 
!== addition: instead of fixed frequency of unity, the ability to vary 
!== frequency over a specified range
!==
!== Note: there is are associated input files needed:
!== "quasi1d_rect_files.input" to know which files to write and 
!== "quasi1d_rect_parameters.input" for the numerical values
!== 
!== conventions: 
!== -use "!== " to start descriptive comments
!== -use "!" to comment out a f90 line not in use
!== -try not to exceed 79 character columns, as that's the printer width
!==  (over 79 character breaks to next line or does not get printed)
!== -program is indented 2 spaces, loops and if statements are 3 indents
!== -the reason "argument, local; argument, scalar" are explicitly labeled is
!== (1) decrease confusion about where the variable is from
!== (2) prevent declaring array before specifying array size from argument
!== -for code analysis, use ftnchek -f90
!==
!=============================================================================

!== eof

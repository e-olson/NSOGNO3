MODULE p4zsms
   !!======================================================================
   !!                         ***  MODULE p4zsms  ***
   !! TOP :   SMELT Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             2014 (E. Olson) adapted to SMELT
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4zsms         :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE trcdta
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zbio          !  Biological model
   USE p4zsbc          !  External source of nutrients
   USE p4zriv          ! river input of nutrients
   USE p4zint          !  time interpolation
   USE iom             !  I/O manager
   USE trd_oce         !  Ocean trends variables
   USE trdtrc          !  TOP trends variables
   USE prtctl_trc      !  print control for debugging
   USE trcbc, only : trc_bc_read
   USE lib_mpp, ONLY: ctl_warn

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sms_init   ! called in p4zsms.F90
   PUBLIC   p4z_sms        ! called in p4zsms.F90


   REAL(wp) :: xfact1, xfact2, xfact3

   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: ???


   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsms.F90 3320 2012-03-05 16:37:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sms  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   ji, jj, jk, jnt, jn, jl
      REAL(wp) ::  ztra
      CHARACTER (len=25) :: charout
      CHARACTER (len=20) :: strsms
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sms')
      !
      IF( kt == nittrc000 ) THEN
        !
        ALLOCATE( xnegtr(jpi,jpj,jpk) )
        !
      ENDIF

      CALL trc_bc_read  ( kt )       ! tracers: surface and lateral Boundary Conditions [EO: this might be a problem if use time splitting]
      !   set time step size (Euler/Leapfrog)
      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN   ;    rfact = rdttrc(1)     !  at nittrc000
      ELSEIF( kt <= nittrc000 + nn_dttrc )                          THEN   ;    rfact = 2. * rdttrc(1)   ! at nittrc000 or nittrc000+nn_dttrc (Leapfrog)
      ENDIF
      !
      
      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN
         DO jn = jp_pcs0, jp_pcs1              !   SMS on tracer without Asselin time-filter
            trb(:,:,:,jn) = trn(:,:,:,jn)
         END DO
      ENDIF

      
      CALL p4z_int( kt ) ! computation of Q10s EVERY TIME STEP, BUT NOT SPLIT STEPS
      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio( kt, jnt )   ! Biology
      END DO

      !
      !
      !
      !!!EO: moved to put surface boundary effects in with physics to make river inputs consistent (as for T,S)
      !!! river inputs are not included in time splitting
      IF( ll_sbc )     CALL p4z_sbc( kt )   ! external sources of nutrients 
      IF (ln_river)    CALL p4z_riv( kt )   ! River input; in p4zsed.F90
      
      !
      IF( lrst_trc )  CALL p4z_rst( kt, 'WRITE' )  !* Write SMELT information in restart file 
      !

      IF ( lwm .AND. kt == nittrc000 ) CALL FLUSH    ( numonp )     ! flush output namelist bio
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sms')
      !
      !
   END SUBROUTINE p4z_sms

   SUBROUTINE p4z_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  p4z_sms_init  ***  
      !!
      !! ** Purpose :   read namelist
      !!----------------------------------------------------------------------
      NAMELIST/nampisbio/ nrdttrc 
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : SMELT variables
      READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : SMELT variables
      READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisbio )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' Namelist : nampisbio'
         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
         WRITE(numout,*) ' nrdttrc not equal 1 not yet supported in SMELT'
      ENDIF

   END SUBROUTINE p4z_sms_init


   SUBROUTINE p4z_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_rst  ***
      !!
      !!  ** Purpose : Read or write variables in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!---------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
   END SUBROUTINE p4z_rst

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_sms: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_sms
#endif 

   !!======================================================================
END MODULE p4zsms 

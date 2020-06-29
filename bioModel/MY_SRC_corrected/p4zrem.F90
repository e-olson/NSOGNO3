MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   SMELT Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   2014 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   USE p4zopt          !  optical model
   !USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE trdtrc
   USE trd_oce

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   zz_remin_NH      ! 1/s 2 months at 10 degrees (Alain 0.1 d-1 at 10degrees)
   REAL(wp), PUBLIC ::   zz_remin_D_DON   ! 1/s DON detritus remineralization rate
   REAL(wp), PUBLIC ::   zz_remin_D_PON   ! 1/s PON detritus remineralization rate
   REAL(wp), PUBLIC ::   zz_remin_D_bSi   ! 1/s ammonium remineralization rate
   
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zz_NH, zz_D_DON, zz_D_PON, zz_D_bSi, zz_NH_oxid
      REAL(wp) ::   zz_remin_NH_DON, zz_remin_NH_PON, zz_Si_remin
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d0      ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
       
      CALL wrk_alloc( jpi , jpj, jpk , z3d0 )
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
      
      z3d0(:,:,:)=0.0_wp
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zz_NH = trn(ji,jj,jk,jpnh4)
                zz_D_DON = trn(ji,jj,jk,jpdon)
                zz_D_PON = trn(ji,jj,jk,jppon)

                ! Remineralization:
                !
                ! Bacterial oxidation of NH3 to NO pool proportional to NH^2
                zz_NH_oxid = zz_remin_NH * zz_NH**2 * tgfunc(ji,jj,jk)
                ! Dissolved organic nitrogen
                zz_remin_NH_DON = zz_remin_D_DON * zz_D_DON * tgfunc(ji,jj,jk)
                ! Particulate organic nitrogen
                zz_remin_NH_PON = zz_remin_D_PON * zz_D_PON * tgfunc(ji,jj,jk)

                IF (zz_remin_NH < 0.) THEN
                   PRINT "(A)","zz_remin_NH < 0. in derivs.f90"
                   PRINT *, zz_remin_NH
                   CALL EXIT(1)
                END IF

                z3d0(ji,jj,jk)=zz_NH_oxid*cvol(ji,jj,jk)
                tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zz_NH_oxid 
                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zz_remin_NH_DON + zz_remin_NH_PON - zz_NH_oxid) 
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zz_remin_NH_DON 
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) - zz_remin_NH_PON 
            END DO
         END DO
      END DO
      
      IF ( iom_use("NITR" )) THEN 
         CALL lbc_lnk( z3d0, 'T', 1. )
         CALL iom_put( "NITR", z3d0 )
      ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                ! Biogenic silicon
                zz_D_bSi = trn(ji,jj,jk,jpdsi)
                zz_Si_remin = zz_remin_D_bSi * zz_D_bSi * tgfunc(ji,jj,jk)
                tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zz_Si_remin
                tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zz_Si_remin
            END DO
         END DO
      END DO

      IF( l_trdtrc ) THEN
         strsms="REM_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF

      !
      CALL wrk_dealloc( jpi , jpj, jpk , z3d0 )
      IF( nn_timing == 1 )  CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ zz_remin_NH, zz_remin_D_DON, zz_remin_D_PON, zz_remin_D_bSi
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )      ! Namelist nampisrem in reference namelist : remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )

      REWIND( numnatp_cfg )      ! Namelist nampisrem in configuration namelist : remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, nampisrem'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    remineralisation rate detritus->NH    zz_remin_NH    =', zz_remin_NH
         WRITE(numout,*) '    remineralization rate detritus->DON   zz_remin_D_DON =', zz_remin_D_DON
         WRITE(numout,*) '    remineralization rate detritus->PON   zz_remin_D_PON =', zz_remin_D_PON
         WRITE(numout,*) '    remineralization rate detritus->Si    zz_remin_D_bSi =', zz_remin_D_bSi
      ENDIF

   END SUBROUTINE p4z_rem_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_rem                    ! Empty routine
   END SUBROUTINE p4z_rem
#endif 

   !!======================================================================
END MODULE p4zrem

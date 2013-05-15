PROGRAM LESSSIMPLE

  !set up modules
  USE sps_vars; USE nrtype; USE sps_utils
  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90
  INTEGER :: i, j, k
  REAL :: zave
  !define variable for SSP spectrum
  REAL, DIMENSION(ntfull,nspec)  :: spec_ssp
  !define variables for Mass and Lbol info
  REAL, DIMENSION(ntfull)    :: mass_ssp,lbol_ssp
  CHARACTER(100) :: file1='', file2=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp

  !-----------------------------------------------------------!
  
  ! Now we're going to show you how to use full  
  ! metallicity-dependent info                   
  
  imf_type = 0    ! Salpeter IMF
  dust_type = 0
  pset%sfh = 4    ! compute delayed SFH
  pset%dust_index =-0.7

  !here we have to read in all the librarries
  CALL SPS_SETUP(-1)

  !compute all SSPs (i.e. at all Zs)
  !nz and the various *ssp_zz arrays are stored 
  !in the common block set up in sps_vars.f90
  DO i=1,nz
     pset%zmet = i
     CALL SSP_GEN(pset,mass_ssp,&
          lbol_ssp_zz(i,:),spec_ssp_zz(i,:,:))

     DO j=1,ntau
        pset%tau_sf=my_tau[j]
        DO k=0,1
           pset%dust1=0.7*k
           pset%dust2=0.3*k

           write(file1) 'fspsSalp_delay'
           CALL COMPSP(3,1,file2,mass_ssp_zz(i,:),lbol_ssp_zz(i:),&
                spec_ssp_zz(i,:),pset,ocompsp) 



  ENDDO


END PROGRAM LESSSIMPLE

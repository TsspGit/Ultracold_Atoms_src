!***********************************************************************
  PROGRAM plottingci
!***********************************************************************

! plots of the wave functions and trap

! Program plottingci comprises:

!  SUBROUTINES:     initialization_tad2hm
!  -----------      read_data_tagtd2hm
!                   calcul_knot_seq
!                   finish_programm

!=======================================================================

!   ===
    USE b_spline
    USE quadrature
    USE global_data_math
    USE global_data_physics
    USE global_data_otad2h ! fabio 
    USE global_data_citad2h  ! fabio 
    USE global_data_plotting
    USE global_data_optical_lattice
!   ===

    IMPLICIT NONE
!   -------------

    CHARACTER( 256 )        :: str_operation              ! auxiliary string that keeps temporary string values

!----------------------------------------------------------------------------------------

  CALL global_initialize
  CALL read_data_tagtd2hm ! fabio                             ! reading data from input files
!  CALL harmonic_values
!write(6,*) 'w_ho : ', w_ho
!pause
  CALL getenv('COORDPLOT', str_operation)
  str_operation = ADJUSTL(str_operation)

  IF ( TRIM(str_operation) .EQ. 'trap' ) THEN

     CALL trap_plot

  ELSE IF (str_operation(1:2) .EQ. 'E_') THEN 

     CALL read_orbitals                           ! in this procedure size of the CI matrices is calculated
     CALL create_ci_vector
     CALL energy_values

  ELSE

     CALL calcul_knot_seq                         ! calculation of knot sequences
     CALL read_orbitals                           ! in this procedure size of the CI matrices is calculated
     CALL create_ci_vector

     IF (      str_operation(1:2)  .EQ. 'ci' ) THEN
        CALL plot_ci
     ELSE IF ( str_operation(1:2)  .EQ. 'co' ) THEN
        CALL plot_config
     ELSE IF ( str_operation(1:2)  .EQ. 'cu' ) THEN
        CALL plot_ci_cut

     ELSE IF (TRIM(ADJUSTL(str_operation)) .EQ. 'raddens' ) THEN

        CALL plot_raddens

     END IF

  END IF

  CALL finish_programm

  CONTAINS
! --------

!**************************************************************************
 SUBROUTINE error_sym(i,j)
!**************************************************************************

! % -----------------------------------------------%
! |   if this subroutine is called then xxx.dci    |
! |   file is constructed incorrectly              |
! %------------------------------------------------%

   
   IMPLICIT NONE
!  -------------

   INTEGER, INTENT(IN)   :: i, j

!--------------------------------------------------------------------------

   WRITE(*,*)
   WRITE(*,*) " %--------------------%"
   WRITE(*,*) " |ERROR in input file |"
   WRITE(*,*) " %--------------------%"
   WRITE(*,*)
   WRITE(*,*) sym_ci, " CI vector is configured incorrectly."
   WRITE(*,*) "----------------------------------------"
   WRITE(*,*) "The multiplication table rule violation!"
   WRITE(*,*) "In input file the pair multiplier eigther"
   WRITE(*,*) "absent or the wrong one is given:"
   WRITE(*,*) "            combination "
   WRITE(*,*) "      ",sym_ci, " = ", symmetry_set(i),"(rm) x ", symmetry_set(j), "(CM)"
   WRITE(*,*) "           does not exist or"
   WRITE(*,*) "only one of ", symmetry_set(i), "(rm) or ", symmetry_set(j), "(CM) is zero"
   WRITE(*,*)
   STOP

!**************************************************************************
 END SUBROUTINE error_sym
!**************************************************************************


!**************************************************************************
 SUBROUTINE read_orbitals
!**************************************************************************

! % -----------------------------------------------%
! |  This subroutine gets information from         |
! |  configuration file and creates configurations |
! %------------------------------------------------%


   IMPLICIT NONE
!  -------------

   CHARACTER(LEN=*), PARAMETER      :: group_name6   = 'ACTIVE ORBITALS:' 
   CHARACTER(LEN=*), PARAMETER      :: group_name_rm = 'RELATIVE COORDINATE:'
   CHARACTER(LEN=*), PARAMETER      :: group_name_CM = 'CENTER-OF-MASS COORDINATE:'

   CHARACTER( 512 )                 :: str
   CHARACTER( 512 )                 :: file_config, name_config

   INTEGER                          :: i, j, k, n_i, n_f                    ! auxiliary
   INTEGER                          :: ini_lev, fin_lev                     ! indices of orbitals
   INTEGER                          :: tmp_ini_n, tmp_fin_n                 ! auxiliary
   INTEGER                          :: count_create
   INTEGER                          :: start_position, final_position       ! config. fill in array gradually in every cicle starting in fin. pos.
   INTEGER                          :: config_alpha_dim, config_beta_dim    ! maximum value of n "quantum number" used in input config. file

!---------------------------------------------------------------------------


   ALLOCATE ( sym_ActOrbit_rm( 0:8 ), STAT = istatus )
   IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_ActOrbit_rm' )

   ALLOCATE ( sym_ActOrbit_CM( 0:8 ), STAT = istatus )
   IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_ActOrbit_CM' )

   CALL getenv( 'CONFIGFILE', file_config)
   OPEN( UNIT=iconfig, FILE=file_config, STATUS="OLD", &
         & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

   sym_ActOrbit_rm = 0
   sym_ActOrbit_CM = 0


   DO count_create = 1, 2

      DO

         READ( iconfig, '(A)' ) str

         IF ( TRIM( ADJUSTL(str) ) .EQ. group_name6 ) THEN
            EXIT
         END IF

      END DO

      n_AO_rm = 0
      start_position = 1

!          %--------------------%
      DO ! |RELATIVE COORDINATE:|
!          %--------------------%

         READ( iconfig, '(A)' ) str
         str = ADJUSTL(str)

         IF ( TRIM(str) .EQ. group_name_CM ) THEN
            EXIT
         END IF

         IF ( (str(1:1) .EQ. 'a') .OR. (str(1:1) .EQ. 'b') ) THEN

            IF ( str(1:1) .EQ. 'a' ) THEN
               sym_name = str(1:2)
            ELSE IF ( str(1:1) .EQ. 'b' ) THEN
               sym_name = str(1:3)
            END IF

            SELECT CASE(TRIM(ADJUSTL(sym_name)))

                   CASE("ag")

                      i_sym    =  1

                      a_D      = aD(1)
                      b_D      = bD(1)
                      D_01     = 0
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

!                     We compute the maximum value of m 
!                     because -m_max <= m <= m_max <= l_max                        ! fabio
!
!                     We have to account for m_start in order to calculate well the
!                     orbitals dimensions                                          ! fabio
                      m_start  = 0                                                 ! fabio
                      CALL mmax( max_lrm, max_mrm, max_m ) ! only even values of m ! fabio 
                      max_mrm_D = max_m - MOD( max_m, 2)                           ! fabio 
     
                      IF( MOD( max_lrm_D + max_m, 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                    ! fabio 
                      ELSE                                                         ! fabio 
                          max_lrm_D_ref = max_m - 1                                ! fabio 
                      END IF                                                       ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b1g") 

                      i_sym    =  2

                      a_D  = aD(2)
                      b_D  = bD(2)
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start   = 2                                                 ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m ) ! only even values of m  ! fabio 
                      max_mrm_D = max_m - MOD( max_m, 2)                            ! fabio 
     
                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b2g")

                      i_sym    =  3

                      a_D      = aD(3)
                      b_D      = bD(3)
                      D_01     =  1
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start   = 1                                                ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m ) ! only odd values of m  ! fabio 
                      max_mrm_D = max_m - MOD( max_m + 1, 2)                       ! fabio 

                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                   ! fabio 
                          max_lrm_D_ref = max_m                                    ! fabio 
                      ELSE                                                         ! fabio 
                          max_lrm_D_ref = max_m - 1                                ! fabio 
                      END IF                                                       ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b3g")

                      i_sym    =  4

                      a_D      = aD(4)
                      b_D      = bD(4)
                      D_01     =  1
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start  = 1                                                  ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m )  ! only odd values of m  ! fabio 
                      max_mrm_D = max_m - MOD( max_m + 1, 2)                        ! fabio 

                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("au")

                      i_sym    =  5

                      a_D  = aD(5)
                      b_D  = bD(5)
                      D_01 =  1
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start  = 2                                                  ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m )  ! only even values of m ! fabio 
                      max_mrm_D = max_m - MOD( max_m, 2)                            ! fabio 
     
                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b1u")

                      i_sym    =  6

                      a_D  =  aD(6)
                      b_D  =  bD(6)
                      D_01 =  1
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start  = 0                                                  ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m ) ! only even values of m  ! fabio 
                      max_mrm_D = max_m - MOD( max_m, 2 )                           ! fabio 

                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b2u")

                      i_sym    =  7

                      a_D  =  aD(7)
                      b_D  =  bD(7)
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start  = 1                                                 ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m )  ! only odd values of m  ! fabio 
                      max_mrm_D = max_m - MOD( max_m + 1, 2)                        ! fabio 

                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

                   CASE("b3u")

                      i_sym    =  8

                      a_D  =  aD(8)
                      b_D  =  bD(8)
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start  = 1                                                  ! fabio 
                      CALL mmax( max_lrm, max_mrm, max_m )   ! only odd values of m ! fabio 
                      max_mrm_D = max_m - MOD( max_m + 1, 2)                        ! fabio 

                      IF( MOD( max_lrm_D + max_m , 2 ) == 0)THEN                    ! fabio 
                          max_lrm_D_ref = max_m                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_lrm_D_ref = max_m - 1                                 ! fabio 
                      END IF                                                        ! fabio 

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_rm(i_sym) = dim_rm
                      maxlrmD(i_sym)    = max_lrm_D

            END SELECT

            CALL get_data(str,str)

            config_alpha_dim = -1

            tmp_ini_n = 0
            tmp_fin_n = 0

            DO n_i = 1, 256
            
               IF ( str(n_i:n_i) .EQ. "[" ) THEN
                  DO n_f = n_i, 256  
                     IF ( str(n_f:n_f) .EQ. "-" ) THEN
                        READ( str(n_i+1:n_f-1), * ) ini_lev
                     
                        IF ( ini_lev .LE. tmp_fin_n ) THEN
                           WRITE(*,*) 'THE CALL OF THE CONFIGURATION IS ORGANIZED NON-CORRECTLY'
                           WRITE(*,*) 'in [a->b][c->d] for (r) either'
                           WRITE(*,*) 'a<1 or b<=c'
                           STOP
                        END IF
                     
                        EXIT
                     
                     END IF
                  END DO
               END IF
               
               IF ( str(n_i:n_i) .EQ. ">" ) THEN
                  DO n_f = n_i, 256  
                     IF ( str(n_f:n_f) .EQ. "]" ) THEN
                        READ( str(n_i+1:n_f-1), * ) fin_lev
                        
                        IF ( fin_lev .LT. ini_lev  ) THEN
                           WRITE(*,*) 'THE CALL OF THE CONFIGURATION IS ORGANIZED NON-CORRECTLY'
                           WRITE(*,*) 'in [a->b] for (r) MUST NOT BE b<a'
                           STOP
                        END IF

                        tmp_ini_n = ini_lev
                        tmp_fin_n = fin_lev

                        IF ( fin_lev .GT. config_alpha_dim ) THEN
                           config_alpha_dim = fin_lev
                        END IF

                        IF ( count_create .EQ. 1 ) THEN
                           n_AO_rm = n_AO_rm + (fin_lev - ini_lev) + 1
                        ELSE
                           n_AO_rm = n_AO_rm + (fin_lev - ini_lev) + 1
                           CALL create_ref_vector( ini_lev, fin_lev, start_position, final_position, "rm" )
                           start_position = final_position+1
                        END IF

                        EXIT

                     END IF
                  END DO

               END IF ! >

            END DO ! n_i

            IF ( config_alpha_dim .GT. dim_rm ) THEN

               WRITE(*,*) 'ONE-COORDINATE INPUT BASIS FILE'
               WRITE(*,*) 'DOES NOT CONTAIN rm ORBITALS WHICH ARE CALLED'
               WRITE(*,*) 'IN CONFIGURATION FILE'
               STOP
            END IF

         END IF ! sym_name

      END DO ! RELATIVE COORDINATE


      IF ( count_create .EQ. 1 ) THEN            
         ALLOCATE ( ref_ActOrbit_rm( 1:n_AO_rm ), STAT = istatus )
         IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'ref_ActOrbit_rm' )
      END IF


      n_AO_CM = 0
      start_position = 1

!          %-------------------------%
      DO ! |CENTER-OF-MASS COORDINATE|
!          %-------------------------%

         READ(iconfig,'(A)') str
         str = ADJUSTL(str)

         IF ( TRIM(ADJUSTL(str)) .EQ. "END" ) THEN
            EXIT
         END IF

         IF ( (str(1:1) .EQ. 'A') .OR. (str(1:1) .EQ. 'B') ) THEN

            IF ( str(1:1) .EQ. 'A' ) THEN
               sym_name = str(1:2)
            ELSE IF (str(1:1) .EQ. 'B') THEN
               sym_name = str(1:3)
            END IF

            SELECT CASE(TRIM(ADJUSTL(sym_name)))

                   CASE("Ag")
                  
                      i_sym    =  1

                      a_D      = 2
                      b_D      = 0
                      D_01     = 0
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start   = 0                                                ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M ) ! only even values of M ! fabio 
                      max_MCM_D = max_M - MOD( max_M, 2)                           ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                    ! fabio 
                           max_LCM_D_ref = max_M                                   ! fabio 
                      ELSE                                                         ! fabio 
                           max_LCM_D_ref = max_M - 1                               ! fabio 
                      END IF                                                       ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B1g") 

                      i_sym    =  2
                  
                      a_D  = -2
                      b_D  = -8
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start  = 2                                                  ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only even values of M ! fabio 
                      max_MCM_D = max_M - MOD( max_M, 2)                            ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio 

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B2g") 

                      i_sym    =  3

                      a_D      = -2
                      b_D      = -4
                      D_01     =  1
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start   = 1                                                 ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only odd values of M  ! fabio 
                      max_MCM_D = max_M - MOD( max_M + 1, 2)                        ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B3g") 

                      i_sym    =  4

                      a_D      = -2
                      b_D      = -4
                      D_01     =  1
                      max_lrm_D = max_lrm - MOD(max_lrm,2) ! only even values
                      max_LCM_D = max_LCM - MOD(max_LCM,2)

                      m_start   = 1                                                 ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only odd values of M  ! fabio 
                      max_MCM_D = max_M - MOD( max_M + 1, 2)                        ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                           max_LCM_D_ref = max_M                                    ! fabio 
                      ELSE                                                          ! fabio 
                           max_LCM_D_ref = max_M - 1                                ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("Au") 

                      i_sym    =  5

                      a_D  = -4
                      b_D  = -5
                      D_01 =  1
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start  = 2                                                  ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only even values of M ! fabio 
                      max_MCM_D = max_M - MOD( max_M, 2)                            ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B1u")

                      i_sym    =  6

                      a_D  =  0
                      b_D  = -1
                      D_01 =  1
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start   = 0                                                 ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only even values of M ! fabio 
                      max_MCM_D = max_M - MOD( max_M, 2 )                           ! fabio 
     
                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B2u") 
                      
                      i_sym    =  7

                      a_D  =  0
                      b_D  = -5
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start   = 1                                                 ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )  ! only odd values of M  ! fabio 
                      max_MCM_D = max_M - MOD( max_M + 1, 2)                        ! fabio 

                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 
                      
                      sym_dim_CM(i_sym) = dim_CM
                      maxLCMD(i_sym)    = max_LCM_D

                   CASE("B3u") 

                      i_sym    =  8
                  
                      a_D  =  0
                      b_D  = -5
                      D_01 =  0
                      max_lrm_D = max_lrm - MOD(max_lrm+1,2) ! only ODD values
                      max_LCM_D = max_LCM - MOD(max_LCM+1,2)

                      m_start   = 1                                                 ! fabio 
                      CALL mmax( max_LCM, max_MCM, max_M )   ! only odd values of M ! fabio 
                      max_MCM_D = max_M - MOD( max_M + 1, 2)                        ! fabio 
     
                      IF( MOD( max_LCM_D + max_M, 2 ) == 0)THEN                     ! fabio 
                          max_LCM_D_ref = max_M                                     ! fabio 
                      ELSE                                                          ! fabio 
                          max_LCM_D_ref = max_M - 1                                 ! fabio 
                      END IF                                                        ! fabio

                      CALL initialization_tad2hm ! fabio 
!                     CALL initialization        ! fabio 

                      sym_dim_CM(i_sym) = dim_CM 
                      maxLCMD(i_sym)    = max_LCM_D

                   END SELECT

                END IF


                CALL get_data(str,str)

                config_beta_dim = -1

                tmp_ini_n = 0
                tmp_fin_n = 0

                DO n_i = 1, 256

                   IF ( str(n_i:n_i) .EQ. "[" ) THEN
                      DO n_f = n_i, 256
                         IF ( str(n_f:n_f) .EQ. "-" ) THEN
                            READ( str(n_i+1:n_f-1), * ) ini_lev

                            IF ( ini_lev .LE. tmp_fin_n ) THEN
                               WRITE(*,*) 'THE CALL OF THE CONFIGURATION IS ORGANIZED NON-CORRECTLY'
                               WRITE(*,*) 'in [a->b][c->d] for (r) either'
                               WRITE(*,*) 'a<1 or b<=c'
                               STOP
                            END IF

                            EXIT

                         END IF
                      END DO
                   END IF

                   IF ( str(n_i:n_i) .EQ. ">" ) THEN
                      DO n_f = n_i, 256
                         IF ( str(n_f:n_f) .EQ. "]" ) THEN
                            READ( str(n_i+1:n_f-1), * ) fin_lev

                            IF ( fin_lev .LT. ini_lev  ) THEN
                               WRITE(*,*) 'THE CALL OF THE CONFIGURATION IS ORGANIZED NON-CORRECTLY'
                               WRITE(*,*) 'in [a->b] for (r) MUST NOT BE b<a'
                               STOP
                            END IF

                            tmp_ini_n = ini_lev
                            tmp_fin_n = fin_lev

                            IF ( fin_lev .GT. config_beta_dim ) THEN
                               config_beta_dim = fin_lev
                            END IF

                            IF ( count_create .EQ. 1 ) THEN
                               n_AO_CM = n_AO_CM + (fin_lev - ini_lev) + 1
                            ELSE
                               n_AO_CM = n_AO_CM + (fin_lev - ini_lev) + 1
                               CALL create_ref_vector( ini_lev, fin_lev, start_position, final_position, "CM" )
                               start_position = final_position+1
                            END IF

                            EXIT

                         END IF
                      END DO

                   END IF ! >

                END DO ! n_i

                IF ( config_beta_dim .GT. dim_CM ) THEN
                   WRITE(*,*) 'ONE-COORDINATE INPUT BASIS FILE'
                   WRITE(*,*) 'DOES NOT CONTAIN CM ORBITALS WHICH ARE CALLED'
                   WRITE(*,*) 'IN CONFIGURATION FILE'
                   STOP
                END IF

             END DO ! CENTER-OF-MASS COORDINATE


             IF ( count_create .EQ. 1 ) THEN
                ALLOCATE ( ref_ActOrbit_CM( 1:n_AO_CM ), STAT = istatus )
                IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'ref_ActOrbit_CM' )
             END IF

             REWIND( iconfig )

          END DO ! count_create


! %------------------------------%
! | check the correctness of the |
! |   CI input file content      |
! %------------------------------%

   CALL getenv( 'NAMECONFIG', name_config) 

   i = INDEX(name_config,'_')
   sym_CI = name_config(1:i-1)


   SELECT CASE( TRIM(ADJUSTL(sym_CI)) )

          CASE("Ag")

             DO i = 1, 8

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(i) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(i)) .NE. 0) ) THEN
                   CALL error_sym(i,i)
                END IF

             END DO
             
          CASE("B1g")

             DO i = 1, 8

                j = i+(-1)**(i+1)

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("B2g")
             
             DO i = 1, 8

                j  = i + 2*(-1)**((i-MOD(i-1,2))/2)

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("B3g")


             DO i = 1, 8

                j  = i + (3-2*MOD(i/2,2))*(-1)**(i+1+i/2)

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("Au")

             DO i = 1, 8

                j  = i + 4*SIGN(1,4-i)

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("B1u")

             DO i = 1, 8

                j  = 9 - ( i + 2*(-1)**((i-MOD(i-1,2))/2) )

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("B2u")

             DO i = 1, 8

                j  = 9 - ( i+(-1)**(i+1) )

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE("B3u")

             DO i = 1, 8

                j  = 9 - i

                IF ( (sym_ActOrbit_rm(i)*sym_ActOrbit_CM(j) .EQ. 0) .AND. ((sym_ActOrbit_rm(i)+sym_ActOrbit_CM(j)) .NE. 0) ) THEN
                   CALL error_sym(i,j)
                END IF

             END DO

          CASE DEFAULT  ! the CI file name is given non-correctly

             WRITE(*,*) "FORMAT OF THE INPUT FILE IS D_xxx.dci"
             WRITE(*,*) " D_ IS CASE SENSITIVE AND HERE IS"
             WRITE(*,*) "GIVEN NON-CORRECTLY"
             STOP

   END SELECT

   CLOSE(iconfig)

!**************************************************************************
 END SUBROUTINE read_orbitals
!**************************************************************************

!***************************************************************************************
 SUBROUTINE create_ref_vector( ini_lev, fin_lev, start_position, final_position, coord )
!***************************************************************************************

!  %-----------------------------------------------------------------%
!  | This subroutine creates vector which contains numbers of active |
!  | orbitals involved in configuration storing them in ascending    |
!  | order                                                           |
!  %-----------------------------------------------------------------%


   IMPLICIT NONE
!  -------------

   CHARACTER( * )                      :: coord                                           ! "rm" or "CM"

   INTEGER                             :: i, j                                            ! auxiliary
   INTEGER, INTENT( IN )               :: ini_lev, fin_lev, start_position                ! quantum numbers defyning active orbitals
   INTEGER, INTENT( OUT )              :: final_position                                  ! together with start_position keep sequence of storage

!---------------------------------------------------------------------------------------------------------


   IF ( coord .EQ. "rm" ) THEN

      j = 0
      DO i = start_position, start_position + (fin_lev-ini_lev)
         ref_ActOrbit_rm(i) = ini_lev + j
         j = j + 1
      END DO

      final_position = start_position + (fin_lev-ini_lev)

      sym_ActOrbit_rm(i_sym) = final_position
!     ----------------------

   ELSE IF ( coord .EQ. "CM" ) THEN

      j = 0
      DO i = start_position, start_position + (fin_lev-ini_lev)
         ref_ActOrbit_CM(i) = ini_lev + j
         j = j + 1
      END DO

      final_position = start_position + (fin_lev-ini_lev)

      sym_ActOrbit_CM(i_sym) = final_position
!     ----------------------

   END IF

!*********************************************************************************************************
 END SUBROUTINE create_ref_vector
!*********************************************************************************************************

!*********************************************************************************************************
  SUBROUTINE create_ci_vector
!*********************************************************************************************************

! % ---------------------------------%
! | This subroutine creates CI vector|
! | using ordered active orbitals    |
! %----------------------------------%

    IMPLICIT NONE
!   -------------

    CHARACTER( 512 )                          :: fileEVArm, fileEVACM  ! files with 
    CHARACTER( 512 )                          :: fileEVErm, fileEVECM  !  orbital solutions

    INTEGER                                   :: i, j, k, n
    INTEGER                                   :: C_ActOrbit_rm_dim     ! dimention of C_ActOrbit_rm 
    INTEGER                                   :: C_ActOrbit_CM_dim     ! dimention of C_ActOrbit_CM
    INTEGER                                   :: non0_IREP             ! number of non-zero irredicible presentations in rm&cm must be equal

    REAL( rk ), ALLOCATABLE, DIMENSION (:)    :: orbit_EVA             ! eigenvalues of one-coordinate solutions
    REAL( rk ), ALLOCATABLE, DIMENSION (:)    :: orbit_EVE             ! eigenvectors

!---------------------------------------------------------------------------------------------------------

! %--------------------------------------------------------------------------------%
! |  Some IREP can be not used therefore sym_ActOrbit_rm(:) and sym_ActOrbit_CM(:) |
! |  arrays are allocated again with only active IREP. Number of non-zero IREP     |
! |  will be stored in variable non0_IREP                                          |
! |                                                                                |
! | START NON-ZERO IREP:                                                           |
! |                                                                                |
! | --> sym_istart_C_AO_rm and sym_istart_C_AO_CM are artificially used for        |
! |     strorage sym_ActOrbit_rm(:) and sym_ActOrbit_CM(:) non-zero values         |
! |     but the functcion of these matrices after  NON-ZERO IREP block             |
! |     will be different as given in note in variable declaration part            |
! %--------------------------------------------------------------------------------%


    ALLOCATE ( sym_istart_C_AO_rm( 0:8 ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_istart_C_AO_rm' )

    ALLOCATE ( sym_istart_C_AO_CM( 0:8 ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_istart_C_AO_CM' )


    i_sym_rm  = 0
    i_sym_CM  = 0

    non0_IREP = 0

    DO i = 1, 8

       IF ( sym_ActOrbit_rm(i) .NE. 0 ) THEN

          non0_IREP                     = non0_IREP + 1
          i_sym_rm(non0_IREP)           = i
          sym_istart_C_AO_rm(non0_IREP) = sym_ActOrbit_rm(i)
       END IF

    END DO


    non0_IREP = 0

    DO i = 1, 8

       IF ( sym_ActOrbit_CM(i) .NE. 0 ) THEN

          non0_IREP                     = non0_IREP + 1
          i_sym_CM(non0_IREP)           = i
          sym_istart_C_AO_CM(non0_IREP) = sym_ActOrbit_CM(i)

       END IF

    END DO


    DEALLOCATE ( sym_ActOrbit_rm, STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'sym_ActOrbit_rm')

    DEALLOCATE ( sym_ActOrbit_CM, STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'sym_ActOrbit_CM')

    ALLOCATE ( sym_ActOrbit_rm( 0:non0_IREP ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_ActOrbit_rm' )

    ALLOCATE ( sym_ActOrbit_CM( 0:non0_IREP ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_ActOrbit_CM' )


    DO i = 1, non0_IREP

       sym_ActOrbit_rm(i) = sym_istart_C_AO_rm(i)
       sym_ActOrbit_CM(i) = sym_istart_C_AO_CM(i)

    END DO

    sym_ActOrbit_rm(0) = 0
    sym_ActOrbit_CM(0) = 0


    DEALLOCATE ( sym_istart_C_AO_rm, STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'sym_istart_C_AO_rm')

    DEALLOCATE ( sym_istart_C_AO_CM, STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'sym_istart_C_AO_CM')

    ALLOCATE ( sym_istart_C_AO_rm( 0:non0_IREP ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_istart_C_AO_rm' )

    ALLOCATE ( sym_istart_C_AO_CM( 0:non0_IREP ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'sym_istart_C_AO_CM' )

! END NON_ZERO IREP <<


    C_ActOrbit_rm_dim     = 0
    C_ActOrbit_CM_dim     = 0
    sym_istart_C_AO_rm(0) = 0
    sym_istart_C_AO_CM(0) = 0

! write(*,*) sym_dim_rm(1),sym_dim_rm(2),sym_dim_rm(3),sym_dim_rm(4)

    DO i = 1, non0_IREP

       C_ActOrbit_rm_dim     = C_ActOrbit_rm_dim + (sym_ActOrbit_rm(i)-sym_ActOrbit_rm(i-1))* sym_dim_rm( i_sym_rm(i) )
       sym_istart_C_AO_rm(i) = C_ActOrbit_rm_dim

       C_ActOrbit_CM_dim     = C_ActOrbit_CM_dim + (sym_ActOrbit_CM(i)-sym_ActOrbit_CM(i-1))* sym_dim_CM( i_sym_CM(i) )
       sym_istart_C_AO_CM(i) = C_ActOrbit_CM_dim

    END DO

! orbit EVA/EVE:
    ALLOCATE ( E_ActOrbit_rm( 1:UBOUND(ref_ActOrbit_rm,1) ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'E_ActOrbit_rm' )

    ALLOCATE ( E_ActOrbit_CM( 1:UBOUND(ref_ActOrbit_CM,1) ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'E_ActOrbit_CM' )

    ALLOCATE ( C_ActOrbit_rm( 1:C_ActOrbit_rm_dim ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'C_ActOrbit_rm' )

    ALLOCATE ( C_ActOrbit_CM( 1:C_ActOrbit_CM_dim ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'C_ActOrbit_CM' )

! pointer:
    ALLOCATE ( IREP_a( 1:UBOUND(ref_ActOrbit_rm,1) ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'aa_p_IREP' )

    ALLOCATE ( IREP_b( 1:UBOUND(ref_ActOrbit_CM,1) ), STAT = istatus )
    IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'bb_p_IREP' )


    DO n = 1, non0_IREP

!  %---------------------------------%
!  | relative motion active orbitals:|
!  %---------------------------------%

       i_sym = i_sym_rm(n)

       sym_name = ADJUSTL( symmetry_set( i_sym ) )

       CALL getenv( 'fileEVArm'//TRIM(sym_name), fileEVArm )
       CALL getenv( 'fileEVErm'//TRIM(sym_name), fileEVErm )


       OPEN( UNIT=1000, FILE = fileEVArm, STATUS='OLD', FORM = 'UNFORMATTED' )
       OPEN( UNIT=1001, FILE = fileEVErm, STATUS='OLD', FORM = 'UNFORMATTED' )

       ALLOCATE ( orbit_EVA( 1:sym_dim_rm(i_sym) ), STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_EVA' )

       ALLOCATE ( orbit_EVE( 1:sym_dim_rm(i_sym) ), STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_EVE' )


       READ( 1000 ) orbit_EVA

       DO i = sym_ActOrbit_rm(n-1)+1, sym_ActOrbit_rm(n)

          E_ActOrbit_rm( i ) = orbit_EVA( ref_ActOrbit_rm(i) )

          IREP_a(i) = n
!         -------------
       END DO

! write(*,*) 'sym_ActOrbit_rm(n)',sym_ActOrbit_rm(n),n
! pause
       DO i = sym_ActOrbit_rm(n-1)+1, sym_ActOrbit_rm(n)

          DO j = 1, ref_ActOrbit_rm( i )
             READ( 1001 )  orbit_EVE
          END DO

          REWIND( 1001 )

          k = i - sym_ActOrbit_rm(n-1) - 1

          DO j = 1, sym_dim_rm(i_sym)
!write(*,*) 'ini_AO ??',sym_istart_C_AO_rm(n-1) + sym_dim_rm(i_sym) * k + j

             C_ActOrbit_rm( sym_istart_C_AO_rm(n-1) + sym_dim_rm(i_sym) * k + j ) = orbit_EVE( j )
          END DO

       END DO

       DEALLOCATE ( orbit_EVA , STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'orbit_EVA')

       DEALLOCATE ( orbit_EVE , STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'orbit_EVE')

       CLOSE( 1000 )
       CLOSE( 1001 )


!  %---------------------------------%
!  | center-of-mass active orbitals: |
!  %---------------------------------%

       i_sym = i_sym_CM(n)

       sym_name = ADJUSTL( symmetry_set( i_sym ) )

       CALL getenv( 'fileEVACM'//TRIM(sym_name), fileEVACM )
       CALL getenv( 'fileEVECM'//TRIM(sym_name), fileEVECM )

       OPEN( UNIT=1000, FILE = fileEVACM, STATUS='OLD', FORM = 'UNFORMATTED' )
       OPEN( UNIT=1001, FILE = fileEVECM, STATUS='OLD', FORM = 'UNFORMATTED' )


       ALLOCATE ( orbit_EVA( 1:sym_dim_CM(i_sym) ), STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_EVA' )

       ALLOCATE ( orbit_EVE( 1:sym_dim_CM(i_sym) ), STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_EVE' )

       READ( 1000 ) orbit_EVA

       DO i = sym_ActOrbit_CM(n-1)+1, sym_ActOrbit_CM(n)

          E_ActOrbit_CM( i ) = orbit_EVA( ref_ActOrbit_CM(i) )

          IREP_b(i) = n
!         -------------
       END DO

       DO i = sym_ActOrbit_CM(n-1)+1, sym_ActOrbit_CM(n)

          DO j = 1, ref_ActOrbit_CM( i )
             READ( 1001 )  orbit_EVE
          END DO

          REWIND( 1001 )

          k = i - sym_ActOrbit_CM(n-1) - 1

          DO j = 1, sym_dim_CM(i_sym)
             C_ActOrbit_CM( sym_istart_C_AO_CM(n-1) + sym_dim_CM(i_sym) * k + j ) = orbit_EVE( j )
          END DO

       END DO

       DEALLOCATE ( orbit_EVA , STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'orbit_EVA')

       DEALLOCATE ( orbit_EVE , STAT = istatus )
       IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'orbit_EVE')

       CLOSE( 1000 )
       CLOSE( 1001 )

    END DO
! pause
! write(*,*) "dim" ,sym_dim_rm(1),sym_dim_rm(2),sym_dim_rm(3),sym_dim_rm(4)
! pause
! %--------------------------------------------------------------------------------%
! |  D2h product table:                                                            |
! |  ------------------                                                            |
! |                                                                                |
! | <n>-indext of IREP together with construct "sym_ActOrbit_rm(n-1)+1 to          |
! | sym_ActOrbit_rm(n)" point to active orbitals of <n>-number non-zero IREP or    |
! | rm orbitals. To form IREP of CI this <n>-IREP directly product to <n'> from    |
! | COM according to the multiplication tabel. <n'> which correspond to <n> are    |
! | stored in D2h_PT(n) = n' <n'>-indext of IREP together with construct           |
! | "sym_ActOrbit_CM(n'-1)+1 to sym_ActOrbit_CM(n')" point to active orbitals of   |
! | <n'>-number non-zero IREP or CM orbitals.                                      |
! %--------------------------------------------------------------------------------%

    ALLOCATE ( D2h_PT( 1:non0_IREP ), STAT = istatus )
    IF ( istatus /= 0 ) CALL alloc_error ( ilog, 0, istatus, 'D2h_PT' )

    dim_ci = 0
!   ------
    DO n = 1, non0_IREP

       i = i_sym_rm(n)

       SELECT CASE( TRIM(ADJUSTL(sym_CI)) )

              CASE("Ag")
                 j  = i
              CASE("B1g")
                 j  = i+(-1)**(i+1)
              CASE("B2g")
                 j  = i + 2*(-1)**((i-MOD(i-1,2))/2)
              CASE("B3g")
                 j  = i + (3-2*MOD(i/2,2))*(-1)**(i+1+i/2)
              CASE("Au")
                 j  = i + 4*SIGN(1,4-i)
              CASE("B1u")
                 j  = 9 - ( i + 2*(-1)**((i-MOD(i-1,2))/2) )
              CASE("B2u")
                 j  = 9 - ( i+(-1)**(i+1) )
              CASE("B3u")
                 j  = 9 - i

       END SELECT

       DO k = 1, non0_IREP

          IF ( i_sym_cm(k) .EQ. j ) THEN
             D2h_PT(n) = k
             EXIT
          END IF

       END DO

       dim_ci =  dim_ci + ( sym_ActOrbit_rm(n) - sym_ActOrbit_rm(n-1) ) * ( sym_ActOrbit_CM(k) - sym_ActOrbit_CM(k-1) )

    END DO

    WRITE(*,*) "CI VECTOR LENGTH: ", dim_ci
    WRITE(*,*) "NUMBER OF rm ACTIVE ORBITALS: ", n_AO_rm
    WRITE(*,*) "NUMBER OF CM ACTIVE ORBITALS: ", n_AO_CM

!*********************************************************************************************************
  END SUBROUTINE create_ci_vector
!*********************************************************************************************************
 
!*****************************************
  SUBROUTINE read_coord_plot
!*****************************************

!  this subroutine reads points
!  where cut pseudo 3D plot will
!  be produced

 IMPLICIT NONE
!-------------

 CHARACTER( 512 )                :: plot_input      ! points in format (c*,A,B,C,D,E,F) or (c*,A,B,C,D,E,F,cut coordinate)

 INTEGER                         :: i, j, k         ! aux

 REAL(rk)                        :: point

!-------------------------------------------------------------------------

 CALL getenv( 'COORDPLOT', plot_input )

 k = 0 ! if not Ci-N then this is active coordinate

 DO i = 4, 256

    IF ( plot_input(i:i+1) .EQ. 'x1' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(1) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 1
       END IF

    ELSE IF ( plot_input(i:i+1) .EQ. 'x2' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(2) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 2
       END IF

    ELSE IF ( plot_input(i:i+1) .EQ. 'y1' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(3) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 3
       END IF

    ELSE IF ( plot_input(i:i+1) .EQ. 'y2' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(4) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 4
       END IF

    ELSE IF ( plot_input(i:i+1) .EQ. 'z1' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(5) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 5
       END IF

    ELSE IF ( plot_input(i:i+1) .EQ. 'z2' ) THEN

       IF ( plot_input(i+2:i+2) .EQ. '-' ) THEN
          DO j = i+3, 256
             IF ( (plot_input(j:j) .EQ. ',') .OR. (plot_input(j:j) .EQ. ' ') ) THEN
                READ( plot_input(i+3:j-1), * ) point
                coordinates(6) = point
                EXIT
             END IF
          END DO
       ELSE IF ( (plot_input(i+2:i+2) .EQ. ',') .OR. (plot_input(i+2:i+2) .EQ. ' ') ) THEN
          k = k + 1
          active_coord(k) = 6
       END IF

       EXIT

    END IF

  END DO

  IF ( plot_input(1:2) .EQ. 'cu' ) THEN
     cut_coord = plot_input(j+1:j+2)
  END IF

!*****************************************
  END SUBROUTINE read_coord_plot
!*****************************************

!************************************************************************************
  SUBROUTINE abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)
!************************************************************************************

! this subroutine transfers absoulute coordines into COM/REL motion 
! coordinates in spherical frame

 IMPLICIT NONE
!-------------

 REAL( rk ), INTENT(OUT)      :: rr_rm, theta_rm, phi_rm
 REAL( rk ), INTENT(OUT)      :: RR_CM, Theta_CM, Phi_CM ! COM/REL in spherical frame

 REAL( rk )                   :: REL_x, REL_y, REL_z     ! relative cartesian frame
 REAL( rk )                   :: COM_x, COM_y, COM_z     ! center-of-mass cartesian frame

!-------------------------------------------------------------------------

 REL_x = coordinates(1) - coordinates(2)
 REL_y = coordinates(3) - coordinates(4)
 REL_z = coordinates(5) - coordinates(6)

 COM_x = mu_s(1)*coordinates(1)+mu_s(2)*coordinates(2)
 COM_y = mu_s(1)*coordinates(3)+mu_s(2)*coordinates(4)
 COM_z = mu_s(1)*coordinates(5)+mu_s(2)*coordinates(6)

! %-----%
! | REL |
! %-----%

 IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_y .GT. 0.0_rk)) .AND. (REL_z .NE. 0.0_rk) ) THEN

    rr_rm    = Sqrt(REL_y**2.0_rk+REL_z**2.0_rk)
    phi_rm   = pi/2.0_rk
    theta_rm = aTan( ABS(REL_y)/REL_z )

 ELSE IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_y .LT. 0.0_rk)) .AND. (REL_z .NE. 0.0_rk) ) THEN

    rr_rm    = Sqrt(REL_y**2.0_rk+REL_z**2.0_rk)
    phi_rm   = 1.5_rk*pi
    theta_rm = aTan( ABS(REL_y)/REL_z )

 ELSE IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_y .EQ. 0.0_rk)) .AND. (REL_z .GT. 0.0_rk) ) THEN

    rr_rm  = ABS(REL_z)
    phi_rm = 0.0_rk     ! can be any value
    theta_rm = 0.0_rk

 ELSE IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_y .EQ. 0.0_rk)) .AND. (REL_z .LT. 0.0_rk) ) THEN

    rr_rm  = ABS(REL_z)
    phi_rm = 0.0_rk     ! can be any value
    theta_rm = pi

 ELSE IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_z .EQ. 0.0_rk)) .AND. (REL_y .GT. 0.0_rk) ) THEN

    rr_rm = ABS(REL_y)
    theta_rm = pi/2.0_rk
    phi_rm = pi/2.0_rk

 ELSE IF ( ((REL_x .EQ. 0.0_rk) .AND. (REL_z .EQ. 0.0_rk)) .AND. (REL_y .LT. 0.0_rk) ) THEN

    rr_rm = ABS(REL_y)
    theta_rm = pi/2.0_rk
    phi_rm = 1.5_rk*pi

 ELSE IF ( ((REL_x .GT. 0.0_rk) .AND. (REL_y .EQ. 0.0_rk)) .AND. (REL_z .NE. 0.0_rk) ) THEN

    rr_rm  = Sqrt(REL_x**2.0_rk+REL_z**2.0_rk)
    phi_rm = 0.0_rk
    theta_rm = aTan(ABS(REL_x)/REL_z)

 ELSE IF ( ((REL_x .LT. 0.0_rk) .AND. (REL_y .EQ. 0.0_rk)) .AND. (REL_z .NE. 0.0_rk) ) THEN

    rr_rm  = Sqrt(REL_x**2.0_rk+REL_z**2.0_rk)
    phi_rm = pi
    theta_rm = aTan(ABS(REL_x)/REL_z)

 ELSE IF (  ((REL_y .EQ. 0.0_rk) .AND. (REL_z .EQ. 0.0_rk)) .AND. (REL_x .GT. 0.0_rk) ) THEN

    rr_rm    = ABS(REL_x)
    theta_rm = pi/2.0_rk
    phi_rm   = 0.0_rk

 ELSE IF (  ((REL_y .EQ. 0.0_rk) .AND. (REL_z .EQ. 0.0_rk)) .AND. (REL_x .LT. 0.0_rk) ) THEN

    rr_rm    = ABS(REL_x)
    theta_rm = pi/2.0_rk
    phi_rm   = pi

 ELSE IF ( ((REL_x .NE. 0.0_rk) .AND. (REL_y .NE. 0.0_rk)) .AND. (REL_z .EQ. 0.0_rk) ) THEN

    rr_rm    = Sqrt(REL_x**2.0_rk+REL_y**2.0_rk)
    phi_rm   = aTan(REL_y/REL_x)
    theta_rm = pi/2.0_rk

 ELSE

    rr_rm    = Sqrt(REL_x**2.0_rk+REL_y**2.0_rk+REL_z**2.0_rk)
    phi_rm   = aTan(REL_y/REL_x)
    theta_rm = aTan(Sqrt(REL_x**2.0_rk+REL_y**2.0_rk)/REL_z)

 END IF


! %-----%
! | COM |
! %-----%

 IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_y .GT. 0.0_rk)) .AND. (COM_z .NE. 0.0_rk) ) THEN

    RR_CM    = Sqrt(COM_y**2.0_rk+COM_z**2.0_rk)
    Phi_CM   = pi/2.0_rk
    Theta_CM = aTan( ABS(COM_y)/COM_z )

 ELSE IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_y .LT. 0.0_rk)) .AND. (COM_z .NE. 0.0_rk) ) THEN

    RR_CM    = Sqrt(COM_y**2.0_rk+COM_z**2.0_rk)
    Phi_CM   = 1.5_rk*pi
    Theta_CM = aTan( ABS(COM_y)/COM_z )

 ELSE IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_y .EQ. 0.0_rk)) .AND. (COM_z .GT. 0.0_rk) ) THEN

    RR_CM  = ABS(COM_z)
    Phi_CM = 0.0_rk     ! can be any value
    Theta_CM = 0.0_rk

 ELSE IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_y .EQ. 0.0_rk)) .AND. (COM_z .LT. 0.0_rk) ) THEN

    RR_CM  = ABS(COM_z)
    Phi_CM = 0.0_rk     ! can be any value
    Theta_CM = pi

 ELSE IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_z .EQ. 0.0_rk)) .AND. (COM_y .GT. 0.0_rk) ) THEN

    RR_CM = ABS(COM_y)
    Theta_CM = pi/2.0_rk
    Phi_CM = pi/2.0_rk

 ELSE IF ( ((COM_x .EQ. 0.0_rk) .AND. (COM_z .EQ. 0.0_rk)) .AND. (COM_y .LT. 0.0_rk) ) THEN

    RR_CM = ABS(COM_y)
    Theta_CM = pi/2.0_rk
    Phi_CM = 1.5_rk*pi

 ELSE IF ( ((COM_x .GT. 0.0_rk) .AND. (COM_y .EQ. 0.0_rk)) .AND. (COM_z .NE. 0.0_rk) ) THEN

    RR_CM  = Sqrt(COM_x**2.0_rk+COM_z**2.0_rk)
    Phi_CM = 0.0_rk
    Theta_CM = aTan(ABS(COM_x)/COM_z)

 ELSE IF ( ((COM_x .LT. 0.0_rk) .AND. (COM_y .EQ. 0.0_rk)) .AND. (COM_z .NE. 0.0_rk) ) THEN

    RR_CM  = Sqrt(COM_x**2.0_rk+COM_z**2.0_rk)
    Phi_CM = pi
    Theta_CM = aTan(ABS(COM_x)/COM_z)

 ELSE IF (  ((COM_y .EQ. 0.0_rk) .AND. (COM_z .EQ. 0.0_rk)) .AND. (COM_x .GT. 0.0_rk) ) THEN

    RR_CM    = ABS(COM_x)
    Theta_CM = pi/2.0_rk
    Phi_CM   = 0.0_rk

 ELSE IF (  ((COM_y .EQ. 0.0_rk) .AND. (COM_z .EQ. 0.0_rk)) .AND. (COM_x .LT. 0.0_rk) ) THEN

    RR_CM    = ABS(COM_x)
    Theta_CM = pi/2.0_rk
    Phi_CM   = pi

 ELSE IF ( ((COM_x .NE. 0.0_rk) .AND. (COM_y .NE. 0.0_rk)) .AND. (COM_z .EQ. 0.0_rk) ) THEN

    RR_CM    = Sqrt(COM_x**2.0_rk+COM_y**2.0_rk)
    Phi_CM   = aTan(COM_y/COM_x)
    Theta_CM = pi/2.0_rk

 ELSE

    RR_CM    = Sqrt(COM_x**2.0_rk+COM_y**2.0_rk+COM_z**2.0_rk)
    Phi_CM   = aTan(COM_y/COM_x)
    Theta_CM = aTan(Sqrt(COM_x**2.0_rk+COM_y**2.0_rk)/COM_z)

 END IF
    

! %-----------------------------------------------%
! | correction due to aTan is defined in (-pi,pi] |
! %-----------------------------------------------%

 IF ( theta_rm .LT. 0.0_rk ) THEN
    theta_rm = pi + theta_rm
 END IF

 IF ( REL_x .LT. 0.0_rk ) THEN                                     ! second/third quarter
    phi_rm = pi + phi_rm
 ELSE IF ( (REL_x .GT. 0.0_rk) .AND. (REL_y .LT. 0.0_rk) ) THEN    ! fourth quarter
    phi_rm = twopi + phi_rm
 END IF


 IF ( Theta_CM .LT. 0.0_rk ) THEN
    Theta_CM = pi + Theta_CM
 END IF

 IF ( COM_x .LT. 0.0_rk ) THEN                                     ! second/third quarter
    Phi_CM = pi + Phi_CM
 ELSE IF ( (COM_x .GT. 0.0_rk) .AND. (COM_y .LT. 0.0_rk) ) THEN    ! fourth quarter
    Phi_CM = twopi + Phi_CM
 END IF

!***************************************************************************
  END SUBROUTINE abscart_to_RELCOMspher
!***************************************************************************

!***************************************************************************
 SUBROUTINE ActiveCI
!***************************************************************************

! this subroutine determines CI vector where orbital of rm and CM
! form a configuration which is dominant in this vector

 IMPLICIT NONE
!-------------

 CHARACTER( 512 )          :: str                 ! aux
 CHARACTER( 3 )            :: sym_rm, sym_CM      ! symetry names

 INTEGER                   :: i_EVE               ! sequence number of CI solution with a given active orbitals
 INTEGER                   :: n_EVE_active        ! total number of CI solutions with a given AOs ( if > 1 then it is not clear what to plot)
 INTEGER                   :: i, j, k, i_ci       ! aux
 INTEGER                   :: a_AO, b_AO          ! active orbitals numbers
 INTEGER                   :: non0_IREP           ! number of non-zero irredicible presentations in rm&cm must be equal
 INTEGER                   :: leadO_rm, leadO_CM  ! number of orbitals with leading contribution
 INTEGER                   :: staterm, stateCM    ! active orbitals which must be leading in CI vector
 INTEGER                   :: stateCI		  ! CI state to plot if setted
 INTEGER                   :: nstartrm, nstartCM  ! staterm and stateCM must be recalculatend according to the XXX_x.dci content
 INTEGER                   :: nn_rm, nn_CM

 REAL( 8 )                 :: C_max               ! maximal coefficient
 REAL( 8 )                 :: suma, shannon_entropy
!-------------------------------------------------------------------------

! %------------------------------%
! | number of active orbitals in |
! | accordance with numberring   |
! | in CI input file             |
! %------------------------------%

  CALL getenv('STATErm', str)
  DO i = 1, 256

     IF ( (str(i:i) .EQ. 'a') .OR. (str(i:i) .EQ. 'b') ) THEN
        READ(str(1:i-1),*) staterm
        READ(str(i:i+2),*) sym_rm
        EXIT
     END IF

  END DO

  CALL getenv('STATEcm', str)
  DO i = 1, 256

     IF ( (str(i:i) .EQ. 'A') .OR. (str(i:i) .EQ. 'B') ) THEN
        READ(str(1:i-1),*) stateCM
        READ(str(i:i+2),*) sym_CM
        EXIT
     END IF

  END DO

  CALL getenv('STATEci', str)
  READ(str,*) stateCI

  CALL getenv( 'ciEVE', str )
  OPEN( UNIT=1001, FILE = str, STATUS='OLD', FORM = 'UNFORMATTED' )

  CALL AOstart( 'rm', sym_rm, nstartrm )
  CALL AOstart( 'CM', sym_CM, nstartcm )

  IF ( nstartrm .GT. staterm ) THEN
     WRITE(*,*) 'starting AO of rm in input file is higher than a given one'
     STOP
  ELSE IF ( nstartCM .GT. stateCM ) THEN
     WRITE(*,*) 'starting AO of rm in input file is higher than a given one'
     STOP
  END IF

  staterm = staterm - nstartrm + 1
  stateCM = stateCM - nstartCM + 1

! %-----------------------%
! | determination of the  |
! |  active CI vector     |
! %-----------------------%

  ALLOCATE ( EVA( 1:dim_ci ), STAT=istatus )                      ! eigen vectors EVE(i,level) but for one "level" only
  IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'EVA' )

  non0_IREP = UBOUND(sym_ActOrbit_rm,1)

  n_EVE_active = 0

!*********************************
! FABIO  
147 continue
  CLOSE(1001)
  OPEN( UNIT=1001, FILE = str, STATUS='OLD', FORM = 'UNFORMATTED' )
!   write(6,*) str
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*) 'i_EVE?'
  read(5,*) i_EVE
  
    DO k = 1, i_EVE
     READ(1001) EVA
!     do i=1,dim_ci
!       write(6,*) k, i, EVA(i)
!     end do
!     pause     
  END DO

  suma = 0d0
  C_max = 0.0d0
  shannon_entropy = 0.d0
  do k = 1, dim_ci
     suma = suma + EVA( k )**2   
     IF( EVA( k )**2 .gt. 1d-14 ) THEN
       shannon_entropy = shannon_entropy + EVA( k )**2 * dlog( EVA( k )**2 )     
!      ELSE
!        shannon_entropy = shannon_entropy + dlog( EVA( k )**(2d0*EVA( k )**2))
     END IF
     if( EVA( k )**2 .gt. C_max ) then
       C_max = EVA( k )**2
     end if
  end do
  shannon_entropy = - shannon_entropy / dlog( dfloat( dim_ci ) )
  write(6,*) 'C_max            : ', C_max  
  write(6,*) 'suma             : ', suma 
  write(6,*) 'shannon_entropy  : ', shannon_entropy 
  write(6,*) 
  do k = 1, dim_ci
     if( EVA( k )**2 * 100d0 / suma .gt. 15d0 ) then
!        if( k .le. n_AO_CM ) then
!          nn_rm = 1 + int(dfloat(k)/dfloat(n_AO_CM))
!          nn_CM = k - ( int(dfloat(k)/dfloat(n_AO_CM)) ) * n_AO_CM     
!        else
!          if( int(dfloat(k)/dfloat(n_AO_CM)) .eq. k/n_AO_CM ) then
!            nn_rm = int(dfloat(k)/dfloat(n_AO_CM))
!            nn_CM = k + 1 - ( int(dfloat(k)/dfloat(n_AO_CM)) ) * n_AO_CM     
!          else
!            nn_rm = 1 + int(dfloat(k)/dfloat(n_AO_CM))
!            nn_CM = k - ( int(dfloat(k)/dfloat(n_AO_CM)) ) * n_AO_CM               
!          end  if       
!        end if
         nn_rm = 1 + int(dfloat(k)/dfloat(n_AO_CM))
         nn_CM = k - ( int(dfloat(k)/dfloat(n_AO_CM)) ) * n_AO_CM     
       write(6,*) k, 1 + int(dfloat(k)/dfloat(n_AO_CM)), k - ( int(dfloat(k)/dfloat(n_AO_CM)) ) * n_AO_CM
       write(6,*) 'rm : ', nn_rm, E_ActOrbit_rm( nn_rm )
       write(6,*) 'CM : ', nn_CM, E_ActOrbit_CM( nn_CM )
       write(6,*) 'ci : ', k, EVA( k )**2 * 100d0 / suma, EVA( k )**2 / C_max
       write(6,*)
       write(6,*)
       write(6,*)       
     end if 
  end do
  goto 147
  stop
!*********************************
  CLOSE(1001)

  WRITE(*,*) "Number of CI vector ", i_EVE

!************************************************************************
 END SUBROUTINE ActiveCI
!************************************************************************

!************************************************************************
 SUBROUTINE AOstart( CMrm, strsym, nstart )
!************************************************************************

! the counting of active orbitals in CI code does not start from the 
! very first one  This subroutine returns the starting number of the 
! active orbital

 IMPLICIT NONE
!-------------

 CHARACTER(512)             :: file_config
 CHARACTER(2)               :: CMrm         ! CM or rm
 CHARACTER(3)               :: strsym       ! symetry name
 CHARACTER(512)             :: str          ! aux

 INTEGER                    :: n_i, n_f     ! aux

 INTEGER, INTENT(OUT)       :: nstart       ! starting active orbital

!-------------------------------------------------------------------------

  CALL getenv( 'CONFIGFILE', file_config )
  OPEN( UNIT=iconfig, FILE=file_config, STATUS="OLD", &
      & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  IF ( CMrm .EQ. 'rm' ) THEN

!         %--------------------%
     DO ! |RELATIVE COORDINATE:|
!         %--------------------%

        READ( iconfig, '(A)' ) str
        str = ADJUSTL(str)

        IF ( (str(1:1) .EQ. 'a') .OR. (str(1:1) .EQ. 'b') ) THEN

           IF ( str(1:1) .EQ. 'a' ) THEN
              sym_name = str(1:2)
           ELSE IF ( str(1:1) .EQ. 'b' ) THEN
              sym_name = str(1:3)
           END IF


           IF ( sym_name .EQ. TRIM(ADJUSTL(strsym)) ) THEN

              DO n_i = 1, 256

                 IF ( str(n_i:n_i) .EQ. "[" ) THEN
                    DO n_f = n_i, 256  
                       IF ( str(n_f:n_f) .EQ. "-" ) THEN
                          READ( str(n_i+1:n_f-1), * ) nstart
                          EXIT
                       END IF
                    END DO
                    EXIT
                 END IF
                 
              END DO

              EXIT

           END IF

        END IF

     END DO

  ELSE IF ( CMrm .EQ. 'CM' ) THEN

!         %-------------------------%
     DO ! |CENTER-OF-MASS COORDINATE|
!         %-------------------------%

        READ( iconfig, '(A)' ) str
        str = ADJUSTL(str)
        
        IF ( (str(1:1) .EQ. 'A') .OR. (str(1:1) .EQ. 'B') ) THEN

           IF ( str(1:1) .EQ. 'A' ) THEN
              sym_name = str(1:2)
           ELSE IF ( str(1:1) .EQ. 'B' ) THEN
              sym_name = str(1:3)
           END IF


           IF ( sym_name .EQ. TRIM(ADJUSTL(strsym)) ) THEN

              DO n_i = 1, 256

                 IF ( str(n_i:n_i) .EQ. "[" ) THEN
                    DO n_f = n_i, 256  
                       IF ( str(n_f:n_f) .EQ. "-" ) THEN
                          READ( str(n_i+1:n_f-1), * ) nstart
                          EXIT
                       END IF
                    END DO
                    EXIT
                 END IF

              END DO

              EXIT

           END IF

        END IF

     END DO

  END IF

  CLOSE(iconfig)

!***********************************************************************
 END SUBROUTINE AOstart
!***********************************************************************

!************************************************************************
 SUBROUTINE plot_ci
!************************************************************************

! this subroutine produces points in absolute coordinate and
! store conditional wave function density into a file

 IMPLICIT NONE
!-------------

 CHARACTER(LEN=512)        :: str                 ! auxiliary string that keeps temporary string values
 CHARACTER(LEN=512)        :: twoDinfo            ! information for xmatrix program e.g. # array2D(1,999,1,999)
 CHARACTER(LEN=512)        :: file_plot           ! file where points for plotting will be stored

 INTEGER                   :: i, j
 INTEGER                   :: nstep               ! number of points in X direction
 INTEGER                   :: act_1, act_2        ! number of active coordinate in "coordinates" array

 REAL(rk)                  :: box                 ! active zone to be plotted
 REAL(rk)                  :: step                ! distanse between two points
 REAL(rk)                  :: act_co1, act_co2    ! value of active coordinate

 REAL(rk)                  :: new_box_rm          ! due to pecularities of the spherical coordinates
 REAL(rk)                  :: new_box_CM          !  for some points solution exists beyond XYZ setteled in input
 REAL( rk )                :: REL_x, REL_y, REL_z ! relative cartesian frame
 REAL( rk )                :: COM_x, COM_y, COM_z ! center-of-mass cartesian frame

 REAL(rk)                  :: rr_rm, RR_CM        ! spherical frame
 REAL(rk)                  :: phi_rm, Phi_CM      !                  COM /
 REAL(rk)                  :: theta_rm, Theta_CM  !                        REL

 REAL(rk)                  :: WFval               ! wave function density

! %-----------------------------------------------%
! | orbital code calculates scaled functions      |
! | therefore the total ci function is devided to |
! | radial parts To avoid numerical devision to   |
! | zero zeroes are excluded                      |
! %-----------------------------------------------%

!-------------------------------------------------------------------------

 CALL ActiveCI
 CALL read_coord_plot

 norm_pref(0) = 0.5_rk
 norm_pref(1) = Sqrt(0.5_rk)

 ALLOCATE ( orbit_aAO( 1:UBOUND(ref_ActOrbit_rm,1) ), STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_aAO' )

 ALLOCATE ( orbit_bAO( 1:UBOUND(ref_ActOrbit_CM,1) ), STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_bAO' )

 DEALLOCATE ( ref_ActOrbit_rm, STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'ref_ActOrbit_rm' )

 DEALLOCATE ( ref_ActOrbit_CM, STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'ref_ActOrbit_CM' )


! %-------------------%
! | {y1,y2,z1,z2} = 0 |
! %-------------------%

  CALL getenv('plotdir', str)
  CALL getenv('UNIQUENAME',file_plot)

  file_plot = TRIM(ADJUSTL(str))//'/'//TRIM(ADJUSTL(file_plot))

  CALL getenv('STATErm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'_'//TRIM(ADJUSTL(str))

  CALL getenv('STATEcm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'CM'//TRIM(ADJUSTL(str))

  CALL getenv('STATEci', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'_LEV'//TRIM(ADJUSTL(str))//'.dat'

  write(*,*) 'file: ', TRIM( file_plot ) ! fabio 
  
  OPEN(UNIT = 666, FILE=file_plot, STATUS="REPLACE", &
     & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  CALL getenv('ACTIVEZONE', str)
  READ(str,*) box
  box = box - 1.0_rk

  REL_x = coordinates(1) - coordinates(2)
  REL_y = coordinates(3) - coordinates(4)
  REL_z = coordinates(5) - coordinates(6)

  COM_x = mu_s(1)*coordinates(1)+mu_s(2)*coordinates(2)
  COM_y = mu_s(1)*coordinates(3)+mu_s(2)*coordinates(4)
  COM_z = mu_s(1)*coordinates(5)+mu_s(2)*coordinates(6)

  new_box_rm = Sqrt( (2.0*box)**2.0 + REL_x**2.0 + REL_y**2.0 + REL_z**2.0 )
  new_box_CM = Sqrt( box**2.0 + COM_x**2.0 + COM_y**2.0 + COM_z**2.0 )

  IF ( new_box_rm .GT. max_rm ) THEN
     WRITE(*,*) 'zone', new_box_rm,'is beyond', max_rm, 'is not included in REL motion calculations'
     WRITE(*,*) 'ask to plot with smaller box than', box
     STOP
  END IF

  IF ( new_box_CM .GT. max_CM ) THEN
     WRITE(*,*) 'zone', new_box_CM,'is beyond', max_CM, 'is not included in COM motion calculations'
     WRITE(*,*) 'ask to plot with smaller box than', box
     STOP
  END IF

  CALL getenv('NPOINTS', str)
  READ(str,*) nstep

  WRITE(str,*) nstep
  twoDinfo = '# array2D(1,'//TRIM(ADJUSTL(str))//',1,'//TRIM(ADJUSTL(str))//')'

!  WRITE(666,*) TRIM(ADJUSTL(twoDinfo))


! %--------------%
! | plotting ... |
! %--------------%

  step = 2.0_rk*box/(nstep-1)

  act_1 = active_coord(1)
  act_2 = active_coord(2)

  act_co1 = - box - step

  DO i = 1, nstep

     act_co1 = act_co1 + step

     act_co2 = -box - step

     DO j = 1, nstep

        act_co2 = act_co2 + step

        coordinates(act_1) = act_co1
        coordinates(act_2) = act_co2

        CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)

        IF ( (rr_rm .EQ. 0.0_rk) .OR. (RR_CM .EQ. 0.0_rk) ) THEN

           act_co2 = act_co2 + 0.001*step
           coordinates(act_2) = act_co2
           CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)
           act_co2 = act_co2 - 0.001*step

        END IF
        CALL ciWFval(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM, WFval)
        WRITE(666,*) act_co1, act_co2, REAL(WFval,4)*REAL(WFval,4)

     END DO
        WRITE(666,*)
  END DO

  CLOSE(666)

!************************************************************************
 END SUBROUTINE plot_ci
!************************************************************************

!************************************************************************
 SUBROUTINE plot_ci_cut
!************************************************************************

! this subroutine produces points in absolute coordinate and
! store conditional wave function density into a file

 IMPLICIT NONE
!-------------

 CHARACTER(LEN=512)        :: str                ! auxiliary string that keeps temporary string values
 CHARACTER(LEN=512)        :: twoDinfo           ! information for xmatrix program e.g. # array2D(1,999,1,999)
 CHARACTER(LEN=512)        :: file_plot          ! file where points for plotting will be stored
 CHARACTER(LEN=512)        :: name_fixed         ! name_fixed = fixed&fixed_val

 INTEGER                   :: i, j
 INTEGER                   :: nstep              ! number of points in X direction
 INTEGER                   :: act_1, act_2       ! number of active coordinate in "coordinates" array

 REAL(rk)                  :: act_co2            ! value of active coordinate
 REAL(rk)                  :: rr_rm, RR_CM       ! spherical frame
 REAL(rk)                  :: phi_rm, Phi_CM     !                  COM/
 REAL(rk)                  :: theta_rm, Theta_CM !                      REL

 REAL(rk)                  :: box                ! active zone to be plotted
 REAL(rk)                  :: step               ! distanse between two points
 REAL(rk)                  :: fixed_val          ! value of the fixed coordinate

 REAL(rk)                  :: WFval              ! wave function density

! %-----------------------------------------------%
! | orbital code calculates scaled functions      |
! | therefore the total ci function is devided to |
! | radial parts To avoid numerical devision to   |
! | zero zeroes are excluded                      |
! %-----------------------------------------------%

!-------------------------------------------------------------------------

 CALL ActiveCI
 CALL read_coord_plot

 norm_pref(0) = 0.5_rk
 norm_pref(1) = Sqrt(0.5_rk)

 ALLOCATE ( orbit_aAO( 1:UBOUND(ref_ActOrbit_rm,1) ), STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_aAO' )

 ALLOCATE ( orbit_bAO( 1:UBOUND(ref_ActOrbit_CM,1) ), STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'orbit_bAO' )

 DEALLOCATE ( ref_ActOrbit_rm, STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'ref_ActOrbit_rm' )

 DEALLOCATE ( ref_ActOrbit_CM, STAT = istatus )
 IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'ref_ActOrbit_CM' )


! %---------------%
! | fix the fixed |
! %---------------%

  SELECT CASE( cut_coord )

         CASE('x1')
            act_1 = 1
         CASE('x2')
            act_1 = 2
         CASE('y1')
            act_1 = 3
         CASE('y2')
            act_1 = 4
         CASE('z1')
            act_1 = 5
         CASE('z2')
            act_1 = 6

  END SELECT

  fixed_val = coordinates(act_1)

  CALL getenv('plotdir', str)
  CALL getenv('UNIQUENAME',file_plot)

  WRITE(name_fixed,*) fixed_val
  name_fixed = cut_coord//TRIM(ADJUSTL(name_fixed))

  file_plot = TRIM(ADJUSTL(str))//'/'//TRIM(ADJUSTL(name_fixed))//TRIM(ADJUSTL(file_plot))

  CALL getenv('STATErm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'_'//TRIM(ADJUSTL(str))

  CALL getenv('STATEcm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'CM'//TRIM(ADJUSTL(str))//'.dat'

  OPEN(UNIT = 666, FILE=file_plot, STATUS="REPLACE", &
     & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  CALL getenv('ACTIVEZONE', str)
  READ(str,*) box

  IF ( box .GT. max_rm/2.0_rk ) THEN
     WRITE(*,*) 'zone beyound', max_rm, 'is not included in calculations'
     WRITE(*,*) 'ask to plot with smaller box than', box
     STOP
  END IF

  CALL getenv('NPOINTS', str)
  READ(str,*) nstep

  WRITE(str,*) nstep
  twoDinfo = '# array2D(1,'//TRIM(ADJUSTL(str))//',1,'//TRIM(ADJUSTL(str))//')'

  WRITE(666,*) TRIM(ADJUSTL(twoDinfo))

! %--------------%
! | plotting ... |
! %--------------%

  step = 2.0_rk*box/(nstep-1)

  act_2 = active_coord(1)

  act_co2 = -box - step

  DO j = 1, nstep

     act_co2 = act_co2 + step

     coordinates(act_2) = act_co2

     CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)

     IF ( (rr_rm .EQ. 0.0_rk) .OR. (RR_CM .EQ. 0.0_rk) ) THEN

        act_co2 = act_co2 + 0.001*step
        coordinates(act_2) = act_co2
        CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)
        act_co2 = act_co2 - 0.001*step

     END IF

     CALL ciWFval(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM, WFval)
     WRITE(666,*) act_co2, WFval

  END DO

  CLOSE(666)

!************************************************************************
 END SUBROUTINE plot_ci_cut
!************************************************************************

!************************************************************************
 SUBROUTINE ciWFval(rm, theta_rm, phi_rm, CM, Theta_CM, Phi_CM, WFval)
!************************************************************************

! this subroutine calculates ci wave function in a given point

 IMPLICIT NONE
!-------------

 INTEGER                     ::  l, m                          ! quantum numbers
 INTEGER                     ::  n, t, nb                      ! aux
 INTEGER                     ::  a_AO, b_AO                    ! active orbitals number
 INTEGER                     ::  step_back                     ! this variable helps to determine correct span number
 INTEGER                     ::  iAO, i_CI, ini_b_AO, fin_b_AO ! aux
 INTEGER                     ::  ispan_rm, ispan_CM         ! number of the span
 INTEGER                     ::  st_rm, st_CM, fn_rm, fn_CM ! to satisfy boundary conditions 1&n B-S set to zerro
 INTEGER                     ::  alpha, beta      ! aux for index precomputation
 INTEGER                     ::  alpha_part, beta_part      ! aux for index precomputation
 INTEGER                     ::  part_i1, part_i, ini_AO    ! PC of the orbital index

 REAL( rk )                  ::  orbit                      ! aux sum of orbital
 REAL( rk )                  ::  m_rk                       ! REAL(m,rk)
 REAL( rk )                  ::  pl_mi                      ! + or -
 REAL( rk )                  ::  costheta_rm, costheta_CM   ! Cos(theta)

 REAL( rk ), INTENT(IN)      ::  rm, theta_rm, phi_rm       ! point
 REAL( rk ), INTENT(IN)      ::  CM, Theta_CM, Phi_CM       !       in space

 REAL( rk ), EXTERNAL        ::  LegendreP                  ! associated Legendre function of the first kind
 REAL( rk ), EXTERNAL        ::  FactorialRatio             ! n!/m!

 REAL( rk ), INTENT(OUT)     ::  WFval
 REAL( rk ), DIMENSION(k_rm) ::  b_s_rm                     ! B-Splines
 REAL( rk ), DIMENSION(k_CM) ::  b_s_CM                     ! B-Splines

!-------------------------------------------------------------------------

! %---------------------%
! | quest of the span   |
! | for the radial part |
! %---------------------%

   ! initialization needed
   step_back = 1

   DO ispan_rm = k_rm, n_rm

      IF ( knots_rm(ispan_rm) .EQ. rm ) THEN
         step_back = 0
         EXIT
      END IF

      IF ( knots_rm(ispan_rm) .GT. rm ) THEN
         step_back = 1
         EXIT
      END IF

   END DO


   ispan_rm = ispan_rm - step_back

   IF ( ispan_rm .EQ. k_rm ) THEN  ! due to boundary conditions C_1,b = 0
      st_rm = lower_rm
   ELSE
      st_rm = 1
   END IF

   IF ( ispan_rm .EQ. n_rm ) THEN  ! due to boundary condition basis has n_rm-1 B-S and not n_rm B(alpha=n_rm) = 0
      fn_rm = k_rm-1
   ELSE
      fn_rm = k_rm
   END IF


   DO ispan_CM = k_CM, n_CM

      IF ( knots_CM(ispan_CM) .EQ. CM ) THEN
         step_back = 0
         EXIT
      END IF

      IF ( knots_CM(ispan_CM) .GT. CM ) THEN
         step_back = 1
         EXIT
      END IF

   END DO

   ispan_CM = ispan_CM - step_back

   IF ( ispan_CM .EQ. k_CM ) THEN  ! due to boundary conditions C_1,b = 0
      st_CM = lower_CM
   ELSE
      st_CM = 1
   END IF

   IF ( ispan_CM .EQ. n_CM ) THEN  ! due to boundary condition basis has n_CM-1 B-S and not n_CM B(alpha=n_CM) = 0
      fn_CM = k_CM-1
   ELSE
      fn_CM = k_CM
   END IF

   CALL b_spl ( k_rm, n_rm, ispan_rm, rm, knots_rm, b_s_rm )
   CALL b_spl ( k_CM, n_CM, ispan_CM, CM, knots_CM, b_s_CM )


! %-----------------%
! | relative motion |
! %-----------------%

   CosTheta_rm = Cos(theta_rm)
   alpha_part = ispan_rm - k_rm - shift_rm

! write(*,*) sym_istart_C_AO_rm(0),sym_istart_C_AO_rm(1),sym_istart_C_AO_rm(2),sym_istart_C_AO_rm(3),sym_istart_C_AO_rm(4)
! write(*,*) sym_ActOrbit_rm(0),sym_ActOrbit_rm(1),sym_ActOrbit_rm(2),sym_ActOrbit_rm(3),sym_ActOrbit_rm(4)

   orbit_aAO = 0.0_rk
! write(*,*) 'n_AO_rm',n_AO_rm
   DO a_AO = 1, n_AO_rm

      n      = IREP_a(a_AO)
      i_sym  = i_sym_rm(n)
      ini_AO = sym_istart_C_AO_rm(n-1) + sym_dim_rm(i_sym) * ( a_AO - sym_ActOrbit_rm(n-1) - 1 )

!if (a_AO < 30) 
!write(*,*) "ini_AO,n,a_AO",ini_AO,n,a_AO
!if (a_AO == 30) pause

! write(*,*)  'sym_dim_rm(i_sym)' ,sym_dim_rm(i_sym)
! write(*,*)  'sym_ActOrbit_rm(n-1)',sym_ActOrbit_rm(n-1)

      a_D       = aD(i_sym)
      b_D       = bD(i_sym)
      l_start   = lstart(i_sym)
      m_start   = mstart(i_sym)
      max_lrm_D = maxlrmD(i_sym)

      CALL mmax(max_lrm_D, max_mrm, max_m)                                                 ! fabio 
      IF(MOD(max_lrm_D+max_m, 2) == 0)THEN                                                 ! fabio 
         max_lrm_D_ref   = max_m                                                           ! fabio 
      ELSE                                                                                 ! fabio 
         max_lrm_D_ref   = max_m - 1                                                       ! fabio 
      END IF !MOD(max_lrm_D_ref+max_m, 2) = 0                                              ! fabio 
    
! if (i_sym == 1 ) write(*,*) 'max_mrm',max_mrm 
! if (i_sym == 1 ) write(*,*) 'm_start',m_start 
      IF(MOD(max_mrm+m_start, 2) == 0)THEN                                               ! fabio 
          max_mrm_D   = max_mrm                                                           ! fabio 
      ELSE                                                                                 ! fabio 
          max_mrm_D   = max_mrm - 1                                                        ! fabio 
      END IF !MOD(max_mrm_D_ref+max_lrm_D, 2) .NE. 0                                       ! fabio

      pl_mi     = ReIm_sym(i_sym)

      DO l = l_start, max_lrm_D, 2

!         part_i1 = l**2 + a_D*l + b_D
         CALL mmax( l, max_mrm, max_m )                        !fabio 
!if(i_sym==1) write(*,*) 'max_m',max_m 
         DO m = m_start, max_m, 2                    !fabio
!         DO m = m_start, l, 2

            IF( l <= max_mrm_D ) THEN							  	  ! bruno
              part_i = ( ( l**2 + a_D*l + b_D + 4*m   )/8 ) * alpha_dim!+ alpha_part			  	  ! bruno
            ELSE									 	  	  ! bruno
              part_i = ((max_lrm_D_ref**2 + a_D*max_lrm_D_ref + b_D + 4*( max_lrm_D_ref))/8 + &	  ! bruno
                  & (((l-max_lrm_D_ref)/2.0_rk)) * (1+((max_mrm_D-m_start)/2.0_rk)) + &	  ! bruno
                  & (ABS(m)-max_mrm_D)/2) * alpha_dim!+ alpha_part  					  	  ! bruno
            ENDIF									  	  	  ! bruno

!if (i_sym == 1 ) write(*,*) 'lm_dim_rm', part_i/alpha_dim,max_mrm_D,max_lrm_D_ref,max_mrm
!            part_i = ( ( part_i1 + 4*m )/8 ) * alpha_dim + alpha_part
            m_rk   = REAL(m,rk)

            orbit = 0.0_rk

            DO t = st_rm, fn_rm

               alpha = ispan_rm - k_rm + t

! write(*,*) "part_i,alpha,shift_rm", part_i,alpha,shift_rm,t,ispan_rm
! pause
               iAO = part_i + ( alpha - shift_rm )
!               iAO = part_i + t

!if ( iAO > sym_dim_rm(i_sym)) write(*,*) 'iAO, sym_dim_rm(i_sym)',iAO,sym_dim_rm(i_sym),i_sym
!   write(*,*)  'ini_AO', ini_AO
!   write(*,*)  "iAO   " ,iAO
!   write(*,*)  "" 
               orbit = orbit + C_ActOrbit_rm( ini_AO + iAO ) * b_s_rm(t)

            END DO

            orbit_aAO(a_AO) = orbit_aAO(a_AO)                               + &
                            & orbit * Sqrt( (REAL(2*l+1,rk)/(4.0_rk * pi) ) * &
                            & FactorialRatio(l-m,l+m) )                     * &
                            & LegendreP( l, m, CosTheta_rm )                * &
                            & (Cos(m_rk*phi_rm)*(1+pl_mi*(-1.0_rk)**m_rk)   + &
                            & Sin(m_rk*phi_rm)*(1-pl_mi*(-1.0_rk)**m_rk) )  * &
                            & norm_pref( m/max(1,m) )

         END DO

      END DO

! write(*,*) 'iAO',iAO


   END DO


! %----------------%
! | center-of-mass |
! %----------------%

   CosTheta_CM = Cos(Theta_CM)
   beta_part = ispan_CM - k_CM - shift_CM

   orbit_bAO = 0.0_rk

   DO b_AO = 1, n_AO_CM

      n      = IREP_b(b_AO)
      i_sym  = i_sym_CM(n)
      ini_AO = sym_istart_C_AO_CM(n-1) + sym_dim_CM(i_sym) * ( b_AO - sym_ActOrbit_CM(n-1) - 1 )
      
      a_D       = aD(i_sym)
      b_D       = bD(i_sym)
      L_start   = Lstart(i_sym)
      M_start   = Mstart(i_sym)
      max_LCM_D = maxLCMD(i_sym)

      CALL mmax(max_LCM_D, max_MCM, max_M)                                                 ! fabio 
      IF(MOD(max_LCM_D+max_M, 2) == 0)THEN                                                 ! fabio 
         max_LCM_D_ref   = max_M                                                           ! fabio 
      ELSE                                                                                 ! fabio 
         max_LCM_D_ref   = max_M - 1                                                       ! fabio 
      END IF !MOD(max_LCM_D_ref+max_M, 2) = 0                                              ! fabio 
     
      IF(MOD(max_MCM+M_start, 2) == 0)THEN                                               ! fabio 
          max_MCM_D   = max_MCM                                                            ! fabio 
      ELSE                                                                                 ! fabio 
          max_MCM_D   = max_MCM - 1                                                        ! fabio 
      END IF !MOD(max_MCM_D_ref+max_LCM_D, 2) .NE. 0                                       ! fabio

      pl_mi     = ReIm_sym(i_sym)

      DO L = L_start, max_LCM_D, 2


!         part_i1 = L**2 + a_D*L + b_D
         CALL mmax( L, max_MCM, max_M )                        !fabio
         DO M = M_start, max_M, 2                              !fabio 
!         DO M = M_start, L, 2

         IF( L <= max_MCM_D ) THEN
           part_I = ( ( L**2 + A_D*L + B_D + 4*ABS(M)   )/8 ) * beta_dim + beta_part	           ! bruno
         ELSE									   ! bruno
           part_I = ((max_LCM_D_ref**2 + A_D*max_LCM_D_ref + B_D + 4*ABS(max_LCM_D_ref))/8 + &	   ! bruno
                & (((L-max_LCM_D_ref)/2.0_rk)) * (1+((max_MCM_D-M_start)/2.0_rk)) + &	   ! bruno
                & (ABS(M)-max_MCM_D)/2) * beta_dim + beta_part 						   ! bruno
         ENDIF     								   ! bruno

!            part_i = ( ( part_i1 + 4*M )/8 ) * beta_dim + beta_part
            M_rk   = REAL(M,rk)

            orbit = 0.0_rk

            DO t = st_CM, fn_CM

               iAO = part_I + t

               orbit = orbit + C_ActOrbit_CM( ini_AO + iAO ) * b_s_CM(t)

            END DO

            orbit_bAO(b_AO) = orbit_bAO(b_AO)                               + &
                            & orbit * Sqrt( (REAL(2*L+1,rk)/(4.0_rk * pi) ) * &
                            & FactorialRatio(L-M,L+M) )                     * &
                            & LegendreP( L, M, CosTheta_CM )                * &
                            & (Cos(M_rk*phi_CM)*(1+pl_mi*(-1.0_rk)**M_rk)   + &
                            & Sin(M_rk*phi_CM)*(1-pl_mi*(-1.0_rk)**M_rk) )  * &
                            & norm_pref(M/max(1,M))

         END DO

      END DO

   END DO


! %-----------------------%
! | value of the function |
! %-----------------------%

   I_ci  = 0
   WFval = 0.0_rk

   DO a_AO = 1, n_AO_rm

      nb       = D2h_PT(IREP_a(a_AO))

      ini_b_AO = sym_ActOrbit_CM(nb-1)+1
      fin_b_AO = sym_ActOrbit_CM(nb)

      DO b_AO = ini_b_AO, fin_b_AO

         I_ci = I_ci + 1

         WFval = WFval + EVA(I_ci) * orbit_aAO(a_AO) * orbit_bAO(b_AO)
!write(*,*) EVA(I_ci),orbit_aAO(a_AO),orbit_bAO(b_AO)
      END DO

   END DO

   WFval = WFval / ( rm * CM )

!***************************************************************************
 END SUBROUTINE ciWFval
!***************************************************************************

!************************************************************************
 SUBROUTINE plot_config
!************************************************************************

! this subroutine plots configuration in absolute coordinates


 IMPLICIT NONE
!-------------

 CHARACTER( 3 )            :: sym_rm, sym_CM      ! symetry names

 CHARACTER(LEN=256)        :: str                 ! auxiliary string that keeps temporary string values
 CHARACTER(LEN=256)        :: twoDinfo            ! information for xmatrix program e.g. # array2D(1,999,1,999)
 CHARACTER(LEN=256)        :: file_plot           ! file where points for plotting will be stored

 INTEGER                   :: i, j
 INTEGER                   :: nstep               ! number of points in X direction

 INTEGER                   :: act_rm, act_cm      ! number of active orbitals
 INTEGER                   :: act_1, act_2        ! number of active coordinate in "coordinates" array

 REAL(rk)                  :: box                 ! active zone to be plotted
 REAL(rk)                  :: step                ! distanse between two points
 REAL(rk)                  :: act_co1, act_co2    ! value of active coordinate

 REAL(rk)                  :: rr_rm, RR_CM        ! spherical frame
 REAL(rk)                  :: phi_rm, Phi_CM      !                  COM/
 REAL(rk)                  :: theta_rm, Theta_CM  !                      REL

 REAL(rk)                  :: config_result       ! configuration in a given point

!-------------------------------------------------------------------------

! %-----------------------------%
! | configuration specification |
! %-----------------------------%

  norm_pref(0) = 0.5_rk
  norm_pref(1) = Sqrt(0.5_rk)

  CALL getenv('STATErm', str)
  DO i = 1, 256

     IF ( (str(i:i) .EQ. 'a') .OR. (str(i:i) .EQ. 'b') ) THEN
        READ(str(1:i-1),*) act_rm
        READ(str(i:i+2),*) sym_rm
        IF ( sym_rm(1:1) .EQ. 'a' ) THEN
           sym_rm(1:1) = 'A'
        ELSE
           sym_rm(1:1) = 'B'
        END IF

        EXIT
     END IF

  END DO

  CALL getenv('STATEcm', str)
  DO i = 1, 256

     IF ( (str(i:i) .EQ. 'A') .OR. (str(i:i) .EQ. 'B') ) THEN
        READ(str(1:i-1),*) act_CM
        READ(str(i:i+2),*) sym_CM
        EXIT
     END IF

  END DO

  DO i = 1, 8

     IF ( TRIM(ADJUSTL(symmetry_set(i))) .EQ. TRIM(ADJUSTL(sym_rm)) ) THEN
        i_sym_rm_ordered = i
     END IF

     IF ( TRIM(ADJUSTL(symmetry_set(i))) .EQ. TRIM(ADJUSTL(sym_CM)) ) THEN
        i_sym_CM_ordered = i
     END IF

  END DO

! %-------------------%
! | {y1,y2,z1,z2} = 0 |
! %-------------------%

  CALL getenv('plotdir', str)
  CALL getenv('UNIQUENAME',file_plot)

  file_plot = TRIM(ADJUSTL(str))//'/'//TRIM(ADJUSTL(file_plot))

  CALL getenv('STATErm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'_'//TRIM(ADJUSTL(str))

  CALL getenv('STATEcm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'CM'//TRIM(ADJUSTL(str))//'.dat'

  OPEN(UNIT = 666, FILE=file_plot, STATUS="REPLACE", &
     & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  CALL getenv('ACTIVEZONE', str)
  READ(str,*) box

  IF ( box .GT. max_rm/2.0_rk ) THEN
     WRITE(*,*) 'zone beyound', max_rm, 'is not included in calculations'
     WRITE(*,*) 'ask to plot with smaller box than', box
     STOP
  END IF

! %-----------------------------------------------%
! | orbital code calculates scaled functions      |
! | therefore the total ci function is devided to |
! | radial parts To avoid numerical devision to   |
! | zero zeroes are excluded                      |
! %-----------------------------------------------%

  CALL getenv('NPOINTS', str)
  READ(str,*) nstep

  WRITE(str,*) nstep
  twoDinfo = '# array2D(1,'//TRIM(ADJUSTL(str))//',1,'//TRIM(ADJUSTL(str))//')'

  WRITE(666,*) TRIM(ADJUSTL(twoDinfo))


! %--------------%
! | plotting ... |
! %--------------%

  CALL read_coord_plot

  step = 2.0_rk*box/(nstep-1)

  act_1 = active_coord(1)
  act_2 = active_coord(2)

  act_co1 = - box - step

  DO i = 1, nstep

     act_co1 = act_co1 + step

     act_co2 = -box - step

     DO j = 1, nstep

        act_co2 = act_co2 + step

        coordinates(act_1) = act_co1
        coordinates(act_2) = act_co2

        CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)

        IF ( (rr_rm .EQ. 0.0_rk) .OR. (RR_CM .EQ. 0.0_rk) ) THEN

           act_co2 = act_co2 + 0.001*step
           coordinates(act_2) = act_co2
           CALL abscart_to_RELCOMspher(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM)
           act_co2 = act_co2 - 0.001*step

        END IF

        CALL config_val(rr_rm, theta_rm, phi_rm, RR_CM, Theta_CM, Phi_CM, act_rm, act_CM, config_result)
        WRITE(666,*) REAL(config_result,4)

     END DO

  END DO

  CLOSE(666)

!************************************************************************
 END SUBROUTINE plot_config
!************************************************************************

!************************************************************************************************
 SUBROUTINE config_val(rm, theta_rm, phi_rm, CM, Theta_CM, Phi_CM, act_rm, act_CM, config_result)
!************************************************************************************************

!  this subroutine calculates configuration in a given point

 IMPLICIT NONE
!-------------

 INTEGER                     ::  l, m                          ! quantum numbers
 INTEGER                     ::  n, t, nb                      ! aux
 INTEGER                     ::  step_back                     ! this variable helps to determine correct span number
 INTEGER                     ::  ispan_rm, ispan_CM            ! number of the span
 INTEGER                     ::  st_rm, st_CM, fn_rm, fn_CM    ! to satisfy boundary conditions 1&n B-S set to zerro
 INTEGER                     ::  alpha_part, beta_part         ! aux for index precomputation
 INTEGER                     ::  part_i1, part_i               ! PC of the orbital index

 INTEGER                     ::  iAO
 INTEGER                     ::  ini_AO                        ! initial active orbital
 INTEGER                     ::  non0_IREP                     ! number of non-zero irredicible presentations in rm&cm must be equal
 INTEGER, INTENT(IN)         ::  act_rm, act_CM
 

 REAL( rk )                  ::  orbit                         ! aux sum of orbital
 REAL( rk )                  ::  m_rk                          ! REAL(m,rk)
 REAL( rk )                  ::  pl_mi                         ! + or -
 REAL( rk )                  ::  costheta_rm, costheta_CM      ! Cos(theta)
 REAL( rk )                  ::  wf_rm, wf_CM
 REAL( rk )                  ::  radial_rm, radial_CM

 REAL( rk ), INTENT(IN)      ::  rm, theta_rm, phi_rm       ! point
 REAL( rk ), INTENT(IN)      ::  CM, Theta_CM, Phi_CM       !       in space

 REAL( rk ), EXTERNAL        ::  LegendreP                  ! associated Legendre function of the first kind
 REAL( rk ), EXTERNAL        ::  FactorialRatio             ! n!/m!

 REAL( rk ), INTENT(OUT)     ::  config_result
 REAL( rk ), DIMENSION(k_rm) ::  b_s_rm                     ! B-Splines
 REAL( rk ), DIMENSION(k_CM) ::  b_s_CM                     ! B-Splines

!-------------------------------------------------------------------------

! %---------------------%
! | quest of the span   |
! | for the radial part |
! %---------------------%

   ! initialization needed
   step_back = 1

   DO ispan_rm = k_rm, n_rm

      IF ( knots_rm(ispan_rm) .EQ. rm ) THEN
         step_back = 0
         EXIT
      END IF

      IF ( knots_rm(ispan_rm) .GT. rm ) THEN
         step_back = 1
         EXIT
      END IF

   END DO


   ispan_rm = ispan_rm - step_back

   IF ( ispan_rm .EQ. k_rm ) THEN  ! due to boundary conditions C_1,b = 0
      st_rm = lower_rm
   ELSE
      st_rm = 1
   END IF

   IF ( ispan_rm .EQ. n_rm ) THEN  ! due to boundary condition basis has n_rm-1 B-S and not n_rm B(alpha=n_rm) = 0
      fn_rm = k_rm-1
   ELSE
      fn_rm = k_rm
   END IF


   DO ispan_CM = k_CM, n_CM

      IF ( knots_CM(ispan_CM) .EQ. CM ) THEN
         step_back = 0
         EXIT
      END IF

      IF ( knots_CM(ispan_CM) .GT. CM ) THEN
         step_back = 1
         EXIT
      END IF

   END DO

   ispan_CM = ispan_CM - step_back

   IF ( ispan_CM .EQ. k_CM ) THEN  ! due to boundary conditions C_1,b = 0
      st_CM = lower_CM
   ELSE
      st_CM = 1
   END IF

   IF ( ispan_CM .EQ. n_CM ) THEN  ! due to boundary condition basis has n_CM-1 B-S and not n_CM B(alpha=n_CM) = 0
      fn_CM = k_CM-1
   ELSE
      fn_CM = k_CM
   END IF

   CALL b_spl ( k_rm, n_rm, ispan_rm, rm, knots_rm, b_s_rm )
   CALL b_spl ( k_CM, n_CM, ispan_CM, CM, knots_CM, b_s_CM )


! %-----------------%
! | relative motion |
! %-----------------%

   non0_IREP = UBOUND(sym_ActOrbit_rm,1)

   CosTheta_rm = Cos(theta_rm)
   alpha_part = ispan_rm - k_rm - shift_rm


   DO n = 1, non0_IREP

      i_sym = i_sym_rm(n)
      
      IF ( i_sym .EQ. i_sym_rm_ordered ) THEN
         EXIT
      END IF

   END DO

   ini_AO = sym_istart_C_AO_rm(n-1) + sym_dim_rm(i_sym) * ( act_rm - 1 )

   a_D       = aD(i_sym)
   b_D       = bD(i_sym)
   l_start   = lstart(i_sym)
   m_start   = mstart(i_sym)
   max_lrm_D = maxlrmD(i_sym)

   pl_mi     = ReIm_sym(i_sym)


   wf_rm = 0.0_rk

   DO l = l_start, max_lrm_D, 2

      part_i1 = l**2 + a_D*l + b_D

      DO m = m_start, l, 2

         part_i = ( ( part_i1 + 4*m )/8 ) * alpha_dim + alpha_part
         m_rk   = REAL(m,rk)


         DO t = st_rm, fn_rm

            iAO = part_i + t

            wf_rm = wf_rm + C_ActOrbit_rm( ini_AO + iAO ) * b_s_rm(t)

         END DO

         wf_rm = wf_rm                                         * &
               & Sqrt( (REAL(2*l+1,rk)/(4.0_rk * pi) )         * &
               & FactorialRatio(l-m,l+m) )                     * &
               & LegendreP( l, m, CosTheta_rm )                * &
               & (Cos(m_rk*phi_rm)*(1+pl_mi*(-1.0_rk)**m_rk)   + &
               & Sin(m_rk*phi_rm)*(1-pl_mi*(-1.0_rk)**m_rk) )  * &
               & norm_pref( m/max(1,m) )

      END DO

   END DO


! %----------------%
! | center-of-mass |
! %----------------%

   CosTheta_CM = Cos(Theta_CM)
   beta_part = ispan_CM - k_CM - shift_CM

   non0_IREP = UBOUND(sym_ActOrbit_CM,1)

   DO n = 1, non0_IREP

      i_sym = i_sym_CM(n)

      IF ( i_sym .EQ. i_sym_CM_ordered ) THEN
         EXIT
      END IF

   END DO

   ini_AO = sym_istart_C_AO_CM(n-1) + sym_dim_CM(i_sym) * ( act_cm - 1 )

   a_D       = aD(i_sym)
   b_D       = bD(i_sym)
   L_start   = Lstart(i_sym)
   M_start   = Mstart(i_sym)
   max_LCM_D = maxLCMD(i_sym)

   pl_mi     = ReIm_sym(i_sym)


   wf_CM = 0.0_rk

   DO L = L_start, max_LCM_D, 2

      part_i1 = L**2 + a_D*L + b_D

      DO M = M_start, L, 2

         part_i = ( ( part_i1 + 4*M )/8 ) * beta_dim + beta_part
         M_rk   = REAL(M,rk)

         radial_CM = 0.0_rk
         DO t = st_CM, fn_CM

            iAO = part_i + t

            radial_CM = radial_CM + C_ActOrbit_CM( ini_AO + iAO ) * b_s_CM(t)

         END DO

         radial_CM = radial_CM                                 * &
               & Sqrt( (REAL(2*L+1,rk)/(4.0_rk * pi) )         * &
               & FactorialRatio(L-M,L+M) )                     * &
               & LegendreP( L, M, CosTheta_CM )                * &
               & (Cos(M_rk*phi_CM)*(1+pl_mi*(-1.0_rk)**M_rk)   + &
               & Sin(M_rk*phi_CM)*(1-pl_mi*(-1.0_rk)**M_rk) )  * &
               & norm_pref(M/max(1,M))

         wf_CM = wf_CM + radial_CM

      END DO

   END DO


! %-----------------------%
! | value of the function |
! %-----------------------%

   config_result = wf_rm * wf_CM / ( rm * CM )

!**********************************************************************************************
 END SUBROUTINE config_val
!**********************************************************************************************

!**********************************************************************************************
 SUBROUTINE plot_raddens
!**********************************************************************************************

! this sobroutine plots true radial dencities
! of the full solution (COM is integrated out)

 IMPLICIT NONE
!-------------

 CHARACTER( 256 )                           ::  str                              ! aux
 CHARACTER( 256 )                           ::  file_plot
 LOGICAL                                    ::  OK_flag

 INTEGER                                    ::  t, tp
 INTEGER                                    ::  gqp_rm                           ! numbers of the Gaussiab quadrature point
 INTEGER                                    ::  ispan_rm                         ! number of the span variables
 INTEGER                                    ::  alpha_part
 INTEGER                                    ::  alpha, alpha_p                   ! b-spline indices
 INTEGER                                    ::  l, m, l_p, m_p
 INTEGER                                    ::  st_rm, fn_rm                     ! to satisfy boundary conditions 1&n B-S set to zerro

 INTEGER                                    ::  IAO, JAO
 INTEGER                                    ::  I_ci, J_ci
 INTEGER                                    ::  na, nap, nb
 INTEGER                                    ::  ini_AO_a, ini_AO_ap 
 INTEGER                                    ::  l_start_a, m_start_a
 INTEGER                                    ::  l_start_ap, m_start_ap
 INTEGER                                    ::  a_AO, ap_AO, b_AO, bp_AO
 INTEGER                                    ::  max_lrm_D_a, max_lrm_D_ap
 INTEGER                                    ::  ini_b_AO, fin_b_AO
 INTEGER                                    ::  ini_bp_AO, fin_bp_AO
 INTEGER                                    ::  a_D_a, b_D_a, a_D_ap, b_D_ap
 INTEGER                                    ::  part_i, part_j, part_i1, part_j1
 INTEGER, ALLOCATABLE, DIMENSION(:,:)       ::  ci_index

 REAL( rk )                                 ::  rm
 REAL( rk )                                 ::  rmDEN
 REAL( rk )                                 ::  l_span_rm                         ! length of the span
 REAL( rk )                                 ::  orbit_part 

 REAL( rk ), DIMENSION(k_rm)                ::  b_s_rm
 REAL( rk ), DIMENSION(n_gqp_rm)            ::  root_rm, weight_rm                ! roots and weights for Gaussian point

!------------------------------------------------------------------------

! %------------------------% 
! | configurations indices |
! %------------------------%

  ALLOCATE ( ci_index( n_AO_rm, n_AO_CM ), STAT=istatus )                      ! eigen vectors EVE(i,level) but for one "level" only
  IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'ci_index' )

  I_ci = 0
  DO a_AO = 1, n_AO_rm

     na       = IREP_a(a_AO)
     nb       = D2h_PT(na)

     ini_b_AO = sym_ActOrbit_CM(nb-1)+1
     fin_b_AO = sym_ActOrbit_CM(nb)

     DO b_AO = ini_b_AO, fin_b_AO

        I_ci = I_ci + 1

        ci_index( a_AO, b_AO ) = I_ci

     END DO

  END DO

! %------------------%
! | storage location |
! %------------------%

  CALL getenv('plotdir', str)
  CALL getenv('UNIQUENAME',file_plot)

  file_plot = 'rmD_'//TRIM(ADJUSTL(file_plot))
  file_plot = TRIM(ADJUSTL(str))//'/'//TRIM(ADJUSTL(file_plot))

  CALL getenv('STATErm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'_'//TRIM(ADJUSTL(str))

  CALL getenv('STATEcm', str)
  file_plot = TRIM(ADJUSTL(file_plot))//'CM'//TRIM(ADJUSTL(str))//'.dat'

  OPEN(UNIT = 666, FILE=file_plot, STATUS="REPLACE", &
       & ACCESS="SEQUENTIAL", FORM="FORMATTED" )


! %-----------------------%
! | densities calculation |
! %-----------------------%

 CALL ActiveCI

 CALL gauss( n_gqp_rm, root_rm, weight_rm, ilog )

 DO ispan_rm = k_rm, n_rm           ! run over sections of "rm"

    l_span_rm = knots_rm( ispan_rm+1 ) - knots_rm( ispan_rm ) ! calculation of the length of the span

    IF ( ispan_rm .EQ. k_rm ) THEN  ! due to boundary conditions C_1,b = 0
       st_rm = lower_rm
    ELSE
       st_rm = 1
    END IF

    IF ( ispan_rm .EQ. n_rm ) THEN  ! due to boundary condition basis has n_rm-1 B-S and not n_rm B(alpha=n_rm) = 0
       fn_rm = k_rm-1
    ELSE
       fn_rm = k_rm
    END IF

    alpha_part = ispan_rm - k_rm - shift_rm

    DO gqp_rm = 1, n_gqp_rm         ! run over Gaussian points

       rm = knots_rm( ispan_rm ) + root_rm( gqp_rm ) * l_span_rm     ! this can be found in description of the code manual

       CALL b_spl ( k_rm, n_rm, ispan_rm, rm, knots_rm, b_s_rm )     ! calculates all non-zero B-s in point "r" of the span ispan_rm


       rmDEN = 0.0_rk

       DO a_AO = 1, n_AO_rm

          na       = IREP_a(a_AO)
          i_sym    = i_sym_rm(na)
          ini_AO_a = sym_istart_C_AO_rm(na-1) + sym_dim_rm(i_sym) * ( a_AO - sym_ActOrbit_rm(na-1) - 1 )

          a_D_a  = aD(i_sym)
          b_D_a  = bD(i_sym)

          l_start_a   = lstart(i_sym)
          m_start_a   = mstart(i_sym)
          max_lrm_D_a = maxlrmD(i_sym)

          nb       = D2h_PT(na)
          ini_b_AO = sym_ActOrbit_CM(nb-1)+1
          fin_b_AO = sym_ActOrbit_CM(nb)

          DO ap_AO = 1, n_AO_rm

             nap     = IREP_a(ap_AO)
             i_sym  = i_sym_rm(nap)
             ini_AO_ap = sym_istart_C_AO_rm(nap-1) + sym_dim_rm(i_sym) * ( ap_AO - sym_ActOrbit_rm(nap-1) - 1 )

             a_D_ap  = aD(i_sym)
             b_D_ap  = bD(i_sym)

             l_start_ap   = lstart(i_sym)
             m_start_ap   = mstart(i_sym)
             max_lrm_D_ap = maxlrmD(i_sym)

             nb        = D2h_PT(nap)
             ini_bp_AO = sym_ActOrbit_CM(nb-1)+1
             fin_bp_AO = sym_ActOrbit_CM(nb)


             DO b_AO = ini_b_AO, fin_b_AO

                OK_flag = .FALSE.
                do  bp_AO = ini_bp_AO, fin_bp_AO

                   if ( bp_AO .EQ. b_AO ) then
                      OK_flag = .TRUE.
                      EXIT
                   end if

                end do

                IF ( OK_flag .EQ. .FALSE. ) THEN
                   EXIT
                END IF


                orbit_part = 0.0_rk

                DO l = l_start_a, max_lrm_D_a, 2

                   OK_flag = .FALSE.
                   do l_p = l_start_ap, max_lrm_D_ap, 2

                      if ( l_p .EQ. l ) then
                         OK_flag = .TRUE.
                         EXIT
                      end if

                   end do

                   IF ( OK_flag .EQ. .FALSE. ) THEN
                      EXIT
                   END IF

                   part_i1 = l**2 + a_D_a*l  + b_D_a
                   part_j1 = l**2 + a_D_ap*l + b_D_ap

                   DO m = m_start_a, l, 2

                      OK_flag = .FALSE.
                      DO m_p = m_start_ap, l, 2

                         IF ( m_p .EQ. m ) THEN
                            OK_flag = .TRUE.
                            EXIT
                         END IF

                      END DO

                      IF ( OK_flag .EQ. .FALSE. ) THEN
                         EXIT
                      END IF

                      part_i = ( ( part_i1 + 4*m )/8 ) * alpha_dim + alpha_part
                      part_j = ( ( part_j1 + 4*m )/8 ) * alpha_dim + alpha_part


                      DO t = st_rm, fn_rm

                         iAO = part_i + t

                         DO tp = st_rm, fn_rm

                            jAO = part_j + tp

                            orbit_part = orbit_part                                                         + &
                                       & C_ActOrbit_rm( ini_AO_a + iAO ) * C_ActOrbit_rm( ini_AO_ap + jAO ) * &
                                       & b_s_rm(t) * b_s_rm(tp)

                         END DO ! alpha_p

                      END DO ! alpha

                   END DO ! m

                END DO ! l
                
                I_ci = ci_index( a_AO,  b_AO )
                J_ci = ci_index( ap_AO, b_AO )

                rmDEN = rmDEN + EVA(I_ci)*EVA(J_ci)*orbit_part

             END DO ! b_AO

          END DO ! ap_AO

       END DO ! a_AO

       WRITE(666,*) rm, rmDEN

    END DO ! gqp_rm

 END DO ! ispan_rm

!**********************************************************************************************
 END SUBROUTINE plot_raddens
!**********************************************************************************************

!**********************************************************************************************
  SUBROUTINE energy_values
!**********************************************************************************************

! %---------------------------------------%
! | this routine reads energies of the CI |
! | and output them on a screen           |
! %---------------------------------------%

 IMPLICIT NONE
!-------------

 CHARACTER( 256 )                :: str             ! aux
 CHARACTER( 256 )                :: nE_range        ! enviromantal variable with information about range of the energy level

 INTEGER                         :: i, j            ! aux

 INTEGER                         :: nEini, nEfin    ! initial level energy number, final level energy number

!------------------------------------------------------------------------

  CALL getenv('COORDPLOT', nE_range)
  str_operation = ADJUSTL(nE_range) 
  

  DO i = 1, 256

     IF ( nE_range(i:i) .EQ. 'I') THEN

        DO j = i, 256

           IF ( nE_range(j:j) .EQ. 'F') THEN

              READ(nE_range(i+1:j-2),*) nEini
              READ(nE_range(j+1:256),*) nEfin

              EXIT

           END IF

        END DO

        EXIT

     END IF

  END DO


  CALL getenv( 'ciEVA', str )
  OPEN( UNIT=1111, FILE = str, STATUS='OLD', FORM = 'UNFORMATTED' )

  ALLOCATE ( EVA( 1:dim_ci ), STAT=istatus )                      ! eigen vectors EVE(i,level) but for one "level" only
  IF (istatus /= 0) CALL alloc_error ( ilog, 0, istatus, 'EVA' )

  READ(1111) EVA

  DO i = nEini, nEfin
     WRITE(*,*) i, EVA(i), EVA(i)/w_1kHz
  END DO
  
  DEALLOCATE ( EVA, STAT = istatus )
  IF (istatus /= 0) CALL alloc_error ( ilog, 1, istatus, 'EVA')
  
  CLOSE(1111)

!************************************************************************
  END SUBROUTINE energy_values
!************************************************************************

!************************************************************************
 SUBROUTINE trap_plot
!************************************************************************

! this subroutine make conditional plot of the trap
! in absolute coordinates

 IMPLICIT NONE
!-------------

 CHARACTER(LEN=256)        :: str                    ! auxiliary string that keeps temporary string values
 CHARACTER(LEN=256)        :: file_plot, file_plot2  ! file where points for plotting will be stored
 CHARACTER(LEN=256)        :: twoDinfo               ! information for xmatrix program e.g. # array2D(1,999,1,999)

 INTEGER                   :: i, j
 INTEGER                   :: nstep_x                ! number of points in X direction

 REAL(rk)                  :: Vx1, Vx2               ! lattice depths
 REAL(rk)                  :: box_x                  ! active zone to be plotted
 REAL(rk)                  :: x1, x2, y1             ! point
 REAL(rk)                  :: y2, z1, z2             !       in space
 REAL(rk)                  :: rm_x, CM_x             ! cartesian coordinates
 REAL(rk)                  :: rr_rm, RR_CM           ! point in space
 REAL(rk)                  :: trap_expand            ! value of the trap in expanded situation
 REAL(rk)                  :: step_x                 ! distanse between two points

!------------------------------------------------------------------------

! %-------------------%
! | {y1,y2,z1,z2} = 0 |
! %-------------------%

  WRITE(*,*) 'trap plotting ...'

  CALL getenv('plotdir', str)
  CALL getenv('UNIQUENAME',file_plot)

  file_plot2 = file_plot

  file_plot  = TRIM(ADJUSTL(str))//'/trapReal_'//TRIM(ADJUSTL(file_plot))
  file_plot2 = TRIM(ADJUSTL(str))//'/trapExpand_'//TRIM(ADJUSTL(file_plot2))

  CALL getenv('STATErm', str)
  file_plot  = TRIM(ADJUSTL(file_plot))//'_'//TRIM(ADJUSTL(str))
  file_plot2 = TRIM(ADJUSTL(file_plot2))//'_'//TRIM(ADJUSTL(str))

  CALL getenv('STATEcm', str)
  file_plot  = TRIM(ADJUSTL(file_plot))//'CM'//TRIM(ADJUSTL(str))//'.dat'
  file_plot2 = TRIM(ADJUSTL(file_plot2))//'CM'//TRIM(ADJUSTL(str))//'.dat'

  OPEN(UNIT = 666, FILE=file_plot, STATUS="REPLACE", &
     & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  OPEN(UNIT = 777, FILE=file_plot2, STATUS="REPLACE", &
     & ACCESS="SEQUENTIAL", FORM="FORMATTED" )

  CALL getenv('NPOINTS', str)

  READ(str,*) nstep_x
  str = ADJUSTL(str)

  twoDinfo = '# array2D(1,'//TRIM(str)//',1,'//TRIM(str)//')'
  WRITE(666,*) TRIM(ADJUSTL(twoDinfo))
  WRITE(777,*) TRIM(ADJUSTL(twoDinfo))

  CALL getenv('ACTIVEZONE', str)
  READ(str,*) box_x

  IF ( box_x .GT. max_rm/2.0_rk ) THEN
     WRITE(*,*) 'zone beyound', max_rm/2.0, 'is not included in calculations'
     WRITE(*,*) 'ask to plot with smaller box than', max_rm/2.0
     STOP
  END IF

  Vx1 = Intensity(1) * Polarizability(1)
  Vx2 = Intensity(1) * Polarizability(2)


  step_x = 2.0_rk*box_x/(nstep_x-1)
  x1 = - box_x - step_x

  DO i = 1, nstep_x

     x1 = x1 + step_x

     x2 = -box_x - step_x

     DO j = 1, nstep_x

        x2 = x2 + step_x

        WRITE(666,*) REAL( Vx1*Cos(kc(1)*x1)**2.0_rk + Vx2*Cos(kc(1)*x2)**2.0_rk, 4 )

        rr_rm    = x1-x2
        RR_CM    = mu_s(1)*x1 + mu_s(2)*x2


        CALL trap_expansion( rr_rm, RR_CM, 'x1x2', trap_expand )

        WRITE(777,*) REAL( trap_expand, 4 )

     END DO

  END DO

  CLOSE(666)
  CLOSE(777)


!************************************************************************
 END SUBROUTINE trap_plot
!************************************************************************

!************************************************************************
 SUBROUTINE trap_expansion(rr_rm, RR_CM, ci_cj, trap_expand)
!************************************************************************

! gives value of the trap in a given point in space

 IMPLICIT NONE
!-------------

 CHARACTER(4), INTENT(IN)   :: ci_cj               ! what to plot x1x2, x1y1, x2z1, ...

 INTEGER                    :: nx, eta_s
 INTEGER                    :: i, j, t, k, s

 REAL( rk )                 :: constV              ! look in
 REAL( rk )                 :: sin_term            !         formulas for 
 REAL( rk )                 :: cos_term            !                     COM/RM in absolute frame
 
 REAL( rk ), INTENT(IN)     :: rr_rm, RR_CM        ! center of mass and relative motion point
 REAL( rk ), INTENT(OUT)    :: trap_expand         ! value of the trap in a given points

 REAL( rk ), EXTERNAL       ::  Factorial

!-------------------------------------------------------------------------


 SELECT CASE(ci_cj)

        CASE('x1x2')

           nx = Tnx/2

           sin_term = 0.0_rk
           cos_term = 0.0_rk
              
           DO s = 1, 2

              eta_s = s+(-1)**(s-1)

              DO j = 0, nx-1

                 DO i = 0, nx-1-j

                    sin_term = sin_term                                                        + &
                             & Intensity(1) * Polarizability(s)                                * &
                             & SinCosSign(1) * (-1)**eta_s                                        * &
                             & (-1.0_rk)**Real(i+j,rk) / (Factorial(2*i+1) * Factorial(2*j+1)) * &
                             & (2.0_rk * kc(1) * RR_CM )**REAL(2*i+1,rk)                       * &
                             & (2.0_rk * kc(1)*mu_s(eta_s) * rr_rm)**REAL(2*j+1,rk)

                 END DO ! i

              END DO ! j

              
              DO t = 0, nx

                 DO k = 0, nx-t

                    cos_term = cos_term                                                   + &
                             & Intensity(1) * Polarizability(s)                           * &
                             & (-1.0_rk)**Real(k+t,rk)/ (Factorial(2*k) * Factorial(2*t)) * &
                             & (2.0_rk * kc(1) * RR_CM )**REAL(2*k,rk)                    * &
                             & (2.0_rk * kc(1)*mu_s(eta_s) * rr_rm)**REAL(2*t,rk)

                 END DO

              END DO


           END DO ! s

           constV = Intensity(1) * Polarizability(1) + Intensity(1) * Polarizability(2)
           trap_expand =  0.5 * (constV + sin_term - SinCosSign(1) * cos_term )

 END SELECT

!***************************************************************************************
 END SUBROUTINE trap_expansion
!***************************************************************************************

!***********************************************************************
 SUBROUTINE  finish_programm
!***********************************************************************

! This subroutine writes final lines in log file and closes it
! This is THE END ...

!=======================================================================

  WRITE(ilog,1050) get_time(0)

  1050 FORMAT(6X,"Final time : ",A)

  CLOSE(ilog)

!=======================================================================
 END SUBROUTINE  finish_programm
!=======================================================================


 END PROGRAM plottingci

!***********************************************************************
!***********************************************************************
!                                         0ooo
!                                 ooo0    (   )
!                                (   )     ) /
!                                 \ (     (_/
!                                  \_)

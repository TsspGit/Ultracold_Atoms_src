#!/bin/tcsh
#  C shell program ci2at3d.csh:
#  ============================
#
#  Date:
#  Author:
#
#  EXAMPLE of the program run: 
#                              ./eigenvector_citaold2hm_upm.csh Li7Li7_b_x15000_y15000_z15000_152rm8g2l34m14_102CM8g1L34M14 LiLi_a3Sup_0.70970 Li7a200_Li7a200_kx1064_ky1064_kz1064_ix2273_iy2273_iz2273 sin sin sin 6 6 6  Ag_vsLiLi_610-810_75 610ag 1Ag ci,x1,x2,y1-0,y2-0,z1-0,z2-0 100 8500 D m64 X clea 5
#							   ./eigenvector_citagsd2hm_upm.csh Li7Li7_x18000_y18000_z18000_140rm8g2l60m8_80CM8g1L60M8 LiLi_a3Sup_0.71330 Li7Li7_nx6ny6nz6_nx6ny6nz6_ix4993iy4993iz4993_ix4993iy4993iz50 Ag_vsLiLi_3d1d_rm250_CM50_dd 1007ag 1Ag ci,x1-0,x2-0,y1-0,y2-0,z1,z2 100 8500 D m64 X clea 1002
#                               NOTE:    ci,x1,x2,y1-0,y2-0,z1-0,z2-0       full solution
#                                        co,x1,x2,y1-0,y2-0,z1-0,z2-0       product of orbitals (configuration)
#                                        cu,x1,x2-0,y1-0,y2-0,z1-0,z2-0,x2  cut of full solution for x2 is fixed
#
#======================================================================================================
#
   date
   set echo
#
   set MAINCODE     = eigenvector_citagsd2hm
   set GLOBDATAOP   = global_data_otad2h
   set GLOBDATACI   = global_data_citad2h
   set MAINDIR      = $PWD
   set GENERAL_LIB  = $MAINDIR/../../../../GENERAL_LIB
   set OUTPUTDIR    = ~/TwoAtInOL/d2h/scratch
#
   set CODE        = $MAINDIR/code
   set INPUT       = $MAINDIR/../../config/input
   set COM_LIB     = $MAINDIR/../../common_lib
   set INPUTORBIT  = $MAINDIR/../../orbit/input
   set ORBITOUTPUT = $OUTPUTDIR/otagsd2h
   set OUTPUT      = $OUTPUTDIR/citagsd2h
   set PLOTOUTPUT  = $MAINDIR/out/ci
   set EXEC        = $MAINDIR/exec
   set LOG         = $MAINDIR/log
   set SCR         = $MAINDIR/scratch/
   set path        = ( $path $CODE $EXEC )
   set trapdir     = $COM_LIB/trap
#
#
#      Store arguments in variables:
#      -----------------------------
#
   set BASATDES     = $1  # name of the input file that contains information about BASIS (characteristics of
                          #  number and order of the b-splines, box size) and PARAMETERS OF THE ATOMS
   set MOLECPOTFILE = $2  # file contains numerically or analytically specifyed molecular potential
   set FILECOEF     = $3  # input file (without extension) that contains characteristics of the external confinement
   set CONFIG       = $4 # configuration file
   set AORM         = $5 # relative motion active orbital
   set AOCM         = $6 # center-of-mass active orbital 
   set COORDPLOT    = $7 # ci/co,x1,x2,y1,y2,z1,z2->ci,config
                          # cu,x130,x2,y1,y2,z1,z2->ci-cut (must not contain blanks)
                          # raddens->true radial dens
                          # E_Ixxx_Fyyy just energy values of the CI xxx, yyy are number of levels
   set NPOINTS      = $8 # number of points to be plotted NxN (contour plot)
   set ACTIVEZONE   = $9 # box can not exceed the actual used in calculations
   set PRECISION    = $10 # precision
                          #  Q - quadrupole
                          #  D - double
   set BITMODE      = $11 # ifort32 or ifort64
   set COMPOPTIO    = $12 # compilation options
                          #  E - execute with agressive optimization -O3
                          #  C - execute with Fortran diagnostic
                          #  X - run executable file (previosli compiled with E option)
   set COPYCLEAN    = $13 # delete outputs of all previous calculations or copy current outputs into all claster machines
   set LEV	    = $14 # set CI level to plot directly; Set to -1 if to be neglected
#
   set SINCOSx = sin #fabio
   set SINCOSy = sin #fabio
   set SINCOSz = sin #fabio
   set TAYLORNX = 0   #fabio
   set TAYLORNY = 0   #fabio
   set TAYLORNZ = 0   #fabio
#
   set PLOTOUTPUT = output/ci/                # fabio upm
   set PLOTOUTPUT = output/ci/                # fabio upm
#
   if ( $BITMODE == m32 ) then
      set GENERAL_LIB = $GENERAL_LIB/m32
      set COMPILERBIN = $INTEL_FORTRAN_COMPILER_32
      set COMP_OPT    = "$INTEL_FORTRAN_WORK_OPT_32"
      set ADD_TO_NAME = '_32'
   else if ( $BITMODE == m64 ) then
      set GENERAL_LIB = $GENERAL_LIB/m64
      set COMPILERBIN = '/sw/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/bin/intel64/ifort'       # fabio upm
      set COMP_OPT    = ""                           # fabio upm
#     set COMPILERBIN = $INTEL_FORTRAN_COMPILER_64   # fabio upm
#     set COMP_OPT    = "$INTEL_FORTRAN_WORK_OPT_64" # fabio upm
      set ADD_TO_NAME = '_64'
   endif
#
   unset echo
 if ($COPYCLEAN == 'clean') then
#
    rm -rf $PLOTOUTPUT/*
    exit
#
 endif
# 
   setenv COORDPLOT     $COORDPLOT
#
   setenv STATErm       $AORM
   setenv STATEcm       $AOCM
   setenv STATEci       $LEV            #test
   setenv NPOINTS       $NPOINTS
   setenv ACTIVEZONE    $ACTIVEZONE
   setenv plotdir       $PLOTOUTPUT
#
   setenv SINCOSx       $SINCOSx
   setenv SINCOSy       $SINCOSy
   setenv SINCOSz       $SINCOSz
   setenv TAYLORNX      $TAYLORNX
   setenv TAYLORNY      $TAYLORNY
   setenv TAYLORNZ      $TAYLORNZ
   setenv SCRDIR        $SCR
   setenv NAMEBASATDES  $BASATDES.dba
   setenv NAMECONFIG    $CONFIG.dci
   setenv nameCOEF      $INPUTORBIT/trap/generic_trap_coeff/$FILECOEF.coeff # Tomas
   setenv NAMEMOLECPOT  $MOLECPOTFILE.pot
   setenv BASATDESFILE  $INPUTORBIT/basis/$NAMEBASATDES                     # .dba - data basis atom
   setenv MOLECPOT      $INPUTORBIT/interatompot/$NAMEMOLECPOT              # .pot - potential
   setenv MOLECPOTFILE  $INPUTORBIT/interatompot/$NAMEMOLECPOT              # .pot - potential   
   setenv LOGFILE       $LOG/$CONFIG$BASATDES$MOLECPOTFILE.log  
 #setenv EXTCONFFILE   $INPUTORBIT/trap/$NAMEEXTCONF                       # .dtr - data trap
   setenv CONFIGFILE    $INPUT/$NAMECONFIG                                  # .dci - data ci
#
   setenv RESULTS       $MAINDIR/../analysis/results/$MOLECPOTFILE
#
# S Y M M E T R Y:
#
   set ORBITPART    = $BASATDES'_'$MOLECPOTFILE'_'$FILECOEF
   set ORBITPARTCM  = $BASATDES'_'$FILECOEF
   setenv  UNIQUENAME $SCR/$ORBITPART$CONFIG
#
   setenv NONSYMORBITrm  $ORBITOUTPUT'/eva/rm/'$ORBITPART'.eva'
   setenv NONSYMORBITCM  $ORBITOUTPUT'/eva/CM/'$ORBITPART'.eva'
#
   set symmetry = (Ag B1g B2g B3g Au B1u B2u B3u)
#
    while ($#symmetry > 0 )
#
	   set SYM = $symmetry[1]
#
           set nameEVArm    = $ORBITOUTPUT'/eva/rm/'$SYM'_'$ORBITPART'.eva'
           set nameEVACM    = $ORBITOUTPUT'/eva/CM/'$SYM'_'$ORBITPARTCM'.eva'
           set nameEVErm    = $ORBITOUTPUT'/eve/rm/'$SYM'_'$ORBITPART'.eve'
           set nameEVECM    = $ORBITOUTPUT'/eve/CM/'$SYM'_'$ORBITPARTCM'.eve'
#
           setenv fileEVArm$SYM   $nameEVArm
           setenv fileEVACM$SYM   $nameEVACM
           setenv fileEVErm$SYM   $nameEVErm
           setenv fileEVECM$SYM   $nameEVECM
#
	  shift symmetry # after this comand array is destroied
#
    end # while
#
   set EVAdir = $OUTPUT'/eva/'$ORBITPART
   set EVEdir = $OUTPUT'/eve/'$ORBITPART
#
   setenv ciEVA $EVAdir/$CONFIG.eva
   setenv ciEVE $EVEdir/$CONFIG.eve
   setenv EVErm $nameEVErm
   setenv EVECM $nameEVECM
#
# modules for Quadrupole and Double precision are in different folders
#
  if ( $PRECISION == 'Q') then
   set MODULES      =   $GENERAL_LIB/modules/quadrupole
  endif
#
  if ( $PRECISION == 'D') then
   set MODULES      =   $GENERAL_LIB/modules/double
  endif
#
   set CODELIB      = ( $MAINDIR/../../config/code/$GLOBDATACI.f90 $MAINDIR/../../orbit/code/$GLOBDATAOP.f90 $MAINDIR/code/global_data_plotting.f90)
   set TRAPLIB      = ( $trapdir/global_data_optical_lattice.f90 $trapdir/expansion.f90 )
   set SHAREDLIB1   = ( $GENERAL_LIB/$PRECISION'_vital_math_physics.so' )
   set SHAREDLIB2   = ( $GENERAL_LIB/blas.so $GENERAL_LIB/lapack.so $GENERAL_LIB/machine_dependent.so)
   set SHAREDLIB    = ( $SHAREDLIB1 $SHAREDLIB2 )
   set AUXFILES     = ( $SHAREDLIB $CODELIB $TRAPLIB $COM_LIB/*.f90 )
#
   cd $SCR
#
 if ($COMPOPTIO == 'C') then
    cd $SCR
    rm -f $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
    $COMPILERBIN -C -O0 $COMP_OPT -o $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x $AUXFILES $CODE/$MAINCODE.f90 -I$MODULES -L/usr/X11R6/lib #-lX11
    unset echo
    rm -f $SCR/*
 endif
#
 if ($COMPOPTIO == 'E') then
    cd $SCR
    rm -f $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
    $COMPILERBIN -w $COMP_OPT -o $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x $AUXFILES $CODE/$MAINCODE.f90 -I$MODULES -L/usr/X11R6/lib #-lX11
    echo $COMPILERBIN -w $COMP_OPT -o $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x $AUXFILES $CODE/$MAINCODE.f90 -I$MODULES -L/usr/X11R6/lib
 
    unset echo
    rm -f $SCR/*
 endif
#
 if ($COMPOPTIO == 'X') then
    unset echo
    echo $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
    $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
 endif
#
#
  if ($COMPOPTIO == 'S') then
    unset echo
    cd $MAINDIR
    set QSUBFILE = $BASATDES$MOLECPOTFILE$FILECOEF$SINCOS$TAYLORNX$TAYLORNY$TAYLORNZ$CONFIG$AORM$AOCM$PRECISION$COMPILER
    set OUTPATH     = $OUTDIR/$QSUBFILE.out
    set FULLOUTPATH = frontend.cluster.physik.hu-berlin.de:$OUTPATH
# - - - -
    echo \# > $QSUBFILE.csh
    echo cd \$PBS_O_WORKDIR >> $QSUBFILE.csh
    echo rm -f $OUTPATH >> $QSUBFILE.csh
    echo ci2at3d.csh $BASATDES $MOLECPOTFILE $FILECOEF $SINCOS $TAYLORNX $TAYLORNY $TAYLORNZ $CONFIG $AORM $AOCM\
                   \ $PRECISION $COMPILER X $COPYCLEAN >> $QSUBFILE.csh
    echo rm -f $QSUBFILE.csh >> $QSUBFILE.csh
# - - - -
    chmod u+x $QSUBFILE.csh
    qsub -o $FULLOUTPATH.o -e $FULLOUTPATH.e $QSUBFILE.csh
#
  endif
#
#
#============== END ===================

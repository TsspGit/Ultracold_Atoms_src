#!/bin/csh
#  C shell program citagsd2hm.csh:
#  ============================
#
#  Date:
#  Author:
#
#  Date:   May 2015
#  Author: Fabio Revuelta (based on Serguei's original ci2at3d.f90 program)
#
#  EXAMPLE of the program run: 
#
#  citagsd2hm.csh Li7Li7_b_x10000_y10000_z10000_100rm8g1l6m4_68CM8g1L6M4 LiLi_noint_no_dipole Li7Li7sinsinsin2 Ag_vsLiLi_1-10 D V m64 X clea 
#
# Fabio Revuelta
# fabio.revuelta@upm.es, fabiorevuelta@hotmail.com
#======================================================================================================
#
   date
   set echo 
#
   export  LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2019/linux/lib/intel64_lin/:$LD_LIBRARY_PATH
#
   set MAINCODE    = citagsd2hm         #fabio
   set GLOBDATAOP  = global_data_otad2h
   set GLOBDATACI  = global_data_citad2h
   set MAINDIR     = $PWD
   set GENERAL_LIB  = $MAINDIR/../../../../GENERAL_LIB
   set OUTPUTDIR    = ~/TwoAtInOL/d2h/scratch  # fabio upm
#  set OUTPUTDIR    = $OUTPUTDIR/TwoAtInOL/d2h # fabio upm
#
   set CODE        = $MAINDIR/code
   set INPUT       = $MAINDIR/input
   set COM_LIB     = $MAINDIR/../common_lib 
   set INPUTORBIT  = $MAINDIR/../orbit/input
   set ORBITOUTPUT = $OUTPUTDIR/otagsd2h
   set OUTPUT      = $OUTPUTDIR/citagsd2h
   set EXEC        = $MAINDIR/exec
   set LOG         = $MAINDIR/log
   set SCR      = ~/TwoAtInOL/d2h/scratch/scr/ # fabio upm
#  set SCR      = $SCRATCH                     # fabio upm
   set path        = ( $path $CODE $EXEC ) #fabio
#   set path        = \( $path $CODE $EXEC \) #fabio
   set trapdir     = $COM_LIB/trap
#
#
#      Store arguments in variables:
#      -----------------------------
#
   set BASATDES     = $1   # name of the input file that contains information about BASIS (characteristics of
                           #  number and order of the b-splines, box size) and PARAMETERS OF THE ATOMS
   set MOLECPOTFILE = $2   # file contains numerically or analytically specifyed molecular potential
#   set EXTCONFFILE  = $3   # input file (without extension) that contains characteristics of the external confinement # fabio 
#   set SINCOSx      = $4  # type of the trap "sin" or "cos"
#   set SINCOSy      = $5  # type of the trap "sin" or "cos"
#   set SINCOSz      = $6  # type of the trap "sin" or "cos"
    set FILECOEF     = $3  # fabio: name of the file containing the expansion coefficients of the potential
#    set ATOMION      = $5  #fabio: atom-ion interaction constant
#   set TAYLORNX     = $5  # maximum order of Taylor expansion of external confinement in X direction #fabio $7->$5
#   set TAYLORNY     = $6  # maximum order of Taylor expansion of external confinement in Y direction #fabio $8->$6
#   set TAYLORNZ     = $7  # maximum order of Taylor expansion of external confinement in Z direction #fabio $9->$7
   set CONFIG       = $4   #fabio $10->$8
   set PRECISION    = $5   # precision
                             #  Q - quadrupole
                             #  D - double  #fabio $11->$9
   set JOBZ         = $6   #  N - calculate eigenvalues only
                           #  V - calculate both eigenvalues and eigenvectors  #fabio $12->$10
   set BITMODE      = $7   # ifort32 or ifort64   #fabio $13->$11
   set COMPOPTIO    = $8   # compilation options  #fabio $14->$12
                           #  E - execute with agressive optimization -O3
                           #  C - execute with Fortran diagnostic
                           #  X - run executable file (previosli compiled with E option)
   set COPYCLEAN    = $9  # delete outputs of all previous calculations or copy current outputs into all claster machines  #fabio $15->$13
#
set SINCOSx = sin #fabio
set SINCOSy = sin #fabio
set SINCOSz = sin #fabio
set TAYLORNX = 0   #fabio
set TAYLORNY = 0   #fabio
set TAYLORNZ = 0   #fabio

#setenv ATOM_ION_COEFF  $ATOMION   #fabio

   if ( $BITMODE == m32 ) then
      set GENERAL_LIB = $GENERAL_LIB/m32
      set COMPILERBIN = $INTEL_FORTRAN_COMPILER_32
      set COMP_OPT    = "/sw/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/bin/intel64/ifort"
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
 if ($COPYCLEAN == 'clean') then
#
    rm -rf $OUTPUT/eve/CM/*
    rm -rf $OUTPUT/eva/CM/*
    rm -rf $OUTPUT/eve/rm/*
    rm -rf $OUTPUT/eva/rm/*
    exit
#
 endif
#
 unset echo
#
   setenv SCRDIR        $SCR
   setenv NAMEBASATDES  $BASATDES.dba
   setenv NAMECONFIG    $CONFIG.dci
#    setenv NAMEEXTCONF   $EXTCONFILE.dtr # fabio 
   setenv NAMEMOLECPOT  $MOLECPOTFILE.pot
#   setenv LOGFILE       $LOG/$CONFIG$BASATDES$MOLECPOTFILE$EXTCONFFILE'_Tnx'$TAYLORNX'_Tny'$TAYLORNY'_Tnz'$TAYLORNZ.log #fabio
   setenv LOGFILE       $LOG/$CONFIG$BASATDES$MOLECPOTFILE.log              #fabio
   setenv BASATDESFILE  $INPUTORBIT/basis/$NAMEBASATDES                     # .dba - data basis atom
   setenv MOLECPOT      $INPUTORBIT/interatompot/$NAMEMOLECPOT              # .pot - potential
   setenv MOLECPOTFILE  $INPUTORBIT/interatompot/$NAMEMOLECPOT              # .pot - potential   
#   setenv EXTCONFFILE   $INPUTORBIT/trap/$NAMEEXTCONF                      # .dtr - data trap  # fabio 
   setenv CONFIGFILE    $INPUT/$NAMECONFIG                                  # .dci - data ci
#
   setenv SINCOSx       $SINCOSx
   setenv SINCOSy       $SINCOSy
   setenv SINCOSz       $SINCOSz
   setenv TAYLORNX      $TAYLORNX
   setenv TAYLORNY      $TAYLORNY
   setenv TAYLORNZ      $TAYLORNZ
#
   setenv JOBZ          $JOBZ
   setenv RESULTS       $MAINDIR/../analysis/results/$MOLECPOTFILE
#
# S Y M M E T R Y:
#
#   set ORBITPART   = $BASATDES'_'$EXTCONFFILE'_'$MOLECPOTFILE'_'$SINCOSx'Tnx'$TAYLORNX'_'$SINCOSy'Tny'$TAYLORNY'_'$SINCOSz'Tnz'$TAYLORNZ
#   set ORBITPARTCM = $BASATDES'_'$EXTCONFFILE'_'$SINCOSx'Tnx'$TAYLORNX'_'$SINCOSy'Tny'$TAYLORNY'_'$SINCOSz'Tnz'$TAYLORNZ
#   set ORBITPART   = $BASATDES'_'$EXTCONFFILE'_'$MOLECPOTFILE'_'$FILECOEF'_Tnx'$TAYLORNX'_Tny'$TAYLORNY'_Tnz'$TAYLORNZ #fabio
#   set ORBITPARTCM = $BASATDES'_'$EXTCONFFILE'_'$FILECOEF'_Tnx'$TAYLORNX'_Tny'$TAYLORNY'_Tnz'$TAYLORNZ #fabio
#
   set ORBITPART   = $BASATDES'_'$MOLECPOTFILE'_'$FILECOEF #fabio
   set ORBITPARTCM = $BASATDES'_'$FILECOEF                 #fabio

#   set ORBITPART   = $BASATDES'_'$EXTCONFFILE'_'$MOLECPOTFILE'_'$FILECOEF'_'$ATOMION #fabio
#   set ORBITPARTCM = $BASATDES'_'$EXTCONFFILE'_'$FILECOEF'_'$ATOMION                 #fabio

   setenv  UNIQUENAME    $SCR/$ORBITPART$CONFIG
#
   setenv NONSYMORBITrm  $ORBITOUTPUT'/eva/rm/'$ORBITPART'.eva'
   setenv NONSYMORBITCM  $ORBITOUTPUT'/eva/CM/'$ORBITPART'.eva'
#
#  set nameCOEF   = $INPUT'/../../orbit/input/coef/'$FILECOEF'.dat' # fabio
   set nameCOEF   = $INPUT'/../../orbit/input/trap/generic_trap_coeff/'$FILECOEF'.coeff' # fabio
   setenv nameCOEF  $nameCOEF

   set symmetry = (Ag B1g B2g B3g Au B1u B2u B3u)    #fabio
#   set symmetry = \( Ag B1g B2g B3g Au B1u B2u B3u \) #fabio
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
#   rm -rf $EVAdir
#   rm -rf $EVEdir
   mkdir -p $EVAdir
   mkdir -p $EVEdir
   setenv ciEVA $EVAdir/$CONFIG.eva
   setenv ciEVE $EVEdir/$CONFIG.eve
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
#----------------------------------
 #fabio
   set CODELIB      = ( $CODE/$GLOBDATACI.f90 $MAINDIR/../orbit/code/$GLOBDATAOP.f90 )
   set TRAPLIB      = ( $trapdir/global_data_optical_lattice.f90 $trapdir/expansion.f90 )
   set SHAREDLIB1   = ( $GENERAL_LIB/$PRECISION'_vital_math_physics.so' )
   set SHAREDLIB2   = ( $GENERAL_LIB/blas.so $GENERAL_LIB/lapack.so $GENERAL_LIB/machine_dependent.so)
   set SHAREDLIB    = ( $SHAREDLIB1 $SHAREDLIB2 )
   set AUXFILES     = ( $SHAREDLIB $CODELIB $TRAPLIB $COM_LIB/*.f90 )
#   set CODELIB      = \( $CODE/$GLOBDATACI.f90 $MAINDIR/../orbit/code/$GLOBDATAOP.f90 \)
#   set TRAPLIB      = \( $trapdir/global_data_optical_lattice.f90 $trapdir/expansion.f90 \)
#   set SHAREDLIB1   = \( $GENERAL_LIB/$PRECISION'_vital_math_physics.so' \)
#   set SHAREDLIB2   = \( $GENERAL_LIB/blas.so $GENERAL_LIB/lapack.so $GENERAL_LIB/machine_dependent.so \)
#   set SHAREDLIB    = \( $SHAREDLIB1 $SHAREDLIB2 \)
#   set AUXFILES     = \( $SHAREDLIB $CODELIB $TRAPLIB $COM_LIB/*.f90 \)
#----------------------------------   
#
   cd $SCR
#
 if ($COMPOPTIO == 'C') then
    cd $SCR
    rm -f $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
    $COMPILERBIN -w $COMP_OPT -C -o $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x $AUXFILES $CODE/$MAINCODE.f90 -I$MODULES -L/usr/X11R6/lib #-lX11
    unset echo
    rm -f $SCR/*
 endif
#
 if ($COMPOPTIO == 'E') then
    cd $SCR
    rm -f $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
    $COMPILERBIN -w $COMP_OPT -o $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x $AUXFILES $CODE/$MAINCODE.f90 -I$MODULES -L/usr/X11R6/lib #-lX11
    unset echo
    rm -f $SCR/*
 endif
#
 if ($COMPOPTIO == 'X') then
    unset echo
    $EXEC/$PRECISION'_'$MAINCODE$ADD_TO_NAME.x
 endif
#
#
  if ($COMPOPTIO == 'S') then
    unset echo
    cd $MAINDIR
#    set QSUBFILE = $BASATDES$MOLECPOTFILE$EXTCONFFILE$SINCOSx$SINCOSy$SINCOSz$TAYLORNX$TAYLORNY$TAYLORNZ$CONFIG$PRECISION$JOBZ$COMPILER # fabio 
#
    set QSUBFILE = $BASATDES$MOLECPOTFILE$FILECOEF$CONFIG                        # fabio 
#    set QSUBFILE = $BASATDES$MOLECPOTFILE$EXTCONFFILE$FILECOEF$ATOMION$CONFIG   # fabio 
    set OUTPATH     = $OUTDIR/$QSUBFILE.out
    set FULLOUTPATH = frontend.cluster.physik.hu-berlin.de:$OUTPATH
# - - - -
    echo \# > $QSUBFILE.csh
    echo cd \$PBS_O_WORKDIR >> $QSUBFILE.csh
    echo rm -f $OUTPATH >> $QSUBFILE.csh
#    echo ci2at3d_gen.csh $BASATDES $MOLECPOTFILE $EXTCONFFILE $SINCOSx $SINCOSy $SINCOSz $TAYLORNX $TAYLORNY $TAYLORNZ $CONFIG \ # fabio 
#                   \ $PRECISION $JOBZ $COMPILER X $COPYCLEAN >> $QSUBFILE.csh # fabio 
    echo citagsd2hm.csh $BASATDES $MOLECPOTFILE $FILECOEF $CONFIG \
                   \ $PRECISION $JOBZ $COMPILER X $COPYCLEAN >> $QSUBFILE.csh # fabio 
    echo rm -f $QSUBFILE.csh >> $QSUBFILE.csh
# - - - -
    chmod u+x $QSUBFILE.csh
    qsub -o $FULLOUTPATH.o -e $FULLOUTPATH.e $QSUBFILE.csh
#
  endif
#
#
#============== END ===================


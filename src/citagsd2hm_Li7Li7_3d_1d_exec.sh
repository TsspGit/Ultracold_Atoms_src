#!/bin/csh
##----------------------- Start job description -----------------------
#SBATCH --partition=standard
#SBATCH --job-name=config4493_3d_1d
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=12000
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tssanchezpastor@gmail.com
##------------------------ End job description ------------------------
# 18.000 a.u
srun ./citagsd2hm_upm.csh Li7Li7_x18000_y18000_z18000_140rm8g2l60m8_80CM8g1L60M8 LiLi_a3Sup_0."$1" Li7Li7_nx6ny6nz6_nx6ny6nz6_ix4993iy4993iz4993_ix4993iy4993iz50 Ag_vsLiLi_3d1d_rm250_CM50_dd D V m64 X clea
#srun ./citagsd2hm_upm.csh Li7Li7_x18000_y18000_z18000_140rm8g2l60m8_80CM8g1L60M8 LiLi_noint_no_dipole Li7Li7_nx6ny6nz6_nx6ny6nz6_ix4993iy4993iz4993_ix4993iy4993iz50 Ag_vsLiLi_1-10_dd D N m64 X clea

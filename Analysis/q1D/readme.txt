
# Correr los programas con Li7Li7  para un potencial de sin^2 hasta sexto orden en una trampa anisótropa:

# Sin interaccion
#################
# Orbitales:
nohup ./otaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_100rm8g2l50m14_68CM8g1L50M14 LiLi_noint_no_dipole Li7a200_Li7a200_kx1000_ky1000_kz1000_ix6604_iy4993_iz50 sin sin sin 6 6 6 D V m64 X clea Ag &

# Configuraciones:
./citaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_100rm8g2l50m14_68CM8g1L50M14 LiLi_noint_no_dipole Li7a200_Li7a200_kx1000_ky1000_kz1000_ix6604_iy4993_iz50 sin sin sin 6 6 6 Ag_vsLiLi_1-10_dd D V m64 X clea Ag

# Con interaccion:
##################
# Orbitales: 20.60
nohup ./otaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_120rm8g2l60m14_68CM8g1L60M14 LiLi_a3Sup071 Li7a200_Li7a200_kx1000_ky1000_kz1000_ix6604_iy4993_iz50 sin sin sin 6 6 6 D N m64 X clea Ag &

# Configuraciones:
nohup ./citaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_120rm8g2l60m14_68CM8g1L60M14 LiLi_a3Sup071 Li7a200_Li7a200_kx1000_ky1000_kz1000_ix6604_iy4993_iz50 sin sin sin 6 6 6 Ag_vsLiLi_int_rm200_CM50_dd D V m64 X clea Ag &

#############################################################################################
potpoints_many:                                                                             #
#############################################################################################
./potpoints_many_upm.csh LiLi LiLi_a3Sup 0.716 0.730 0.002 D m64 X

#############################################################################################
scatlength_many:                                                                            #
#############################################################################################
./scatlength_many_upm.csh ~/TwoAtInOL/d2h/orbit/input/interatompot/ LiLi_a3Sup 7.0160040 7.0160040 3.0 0.70777 0.70831 0.00002 rep D m64 X

# Una vez hecho podemos ejecutar a lo bestia unos cuantos orbitales así:
for i in `seq 70700 5 70750`; do echo "LiLi_a3Sup_0."$i""; nohup ./otaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_120rm8g2l60m14_68CM8g1L60M14 LiLi_a3Sup_0."$i" Li7a200_Li7a200_kx1000_ky1000_kz1000_ix6604_iy4993_iz50 sin sin sin 6 6 6 D V m64 X clea Ag rm & sleep 3; done;

# En el servidor del cevisma:
for i in `seq 70700 5 71200`; do echo "LiLi_a3Sup_0."$i""; sbatch otaold2hm_Li7Li7_120rm8g2l60m14_68CM8g1L60M14_LiLi_a3Sup_ix6604_iy4993_iz50_exec.sh $i; done;

for i in `seq 70700 5 71200`; do echo "LiLi_a3Sup_0."$i""; sbatch citaold2hm_Li7Li7_120rm8g2l60m14_68CM8g1L60M14_LiLi_a3Sup_ix6604_iy4993_iz50_exec.sh $i; done;

########################
plotiing_citaold2hm_upm#
########################
./plotting_citaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_120rm8g2l60m14_68CM8g1L60M14 LiLi_a3Sup_0.70780 Li7a200_Li7a200_kx1000_ky1000_kz1000_ix7190_iy4993_iz50 sin sin sin 6 6 6 Ag_vsLiLi_int_rm200_CM50_dd 1007ag 1Ag ci,x1,x2,y1-0,y2-0,z1-0,z2-0.1 100 8500 D m64 X clea 5

###########################
eigenvector_citaold2hm_upm#
###########################
./eigenvector_citaold2hm_upm.csh Li7Li7_x10000_y10000_z10000_120rm8g2l60m14_68CM8g1L60M14 LiLi_a3Sup_0.71115 Li7a200_Li7a200_kx1000_ky1000_kz1000_ix7190_iy4993_iz50 sin sin sin 6 6 6 Ag_vsLiLi_int_rm200_CM50_dd 1007ag 1Ag ci,x1-0,x2-0,y1-0,y2-0,z1,z2 100 8500 D m64 X clea 1002

################
Acerca de 7802:#
################
Hay que ejecutar primero los potenciales de 5 en 5 desde 70700 hasta 70120. Después se ejecutan de dos en dos desde 70777 hasta 70831.

##################################################################################################
Ayuda para ver las autofunciones!!!                                                              #
Usar LiLi_a3Sup_0.0790, provee un asc de -1E-2, que es muy pequeña y se verá mucho mejor, x=-25e3#
################################################################################################## 

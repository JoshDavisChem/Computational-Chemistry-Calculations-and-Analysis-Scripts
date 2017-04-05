#!/bin/bash
############################################################################
#
# Kpoint test for T-LiFeO2, two lithum per calculation
# ECUT = 130
# RHOCUT = 1040
# 
# edited for spelling errors
###########################################################################
for BASE_K in 1 2 3 4 5 6 7
do
X_DIR=$((1 + $BASE_K))
Y_DIR=$BASE_K
Z_DIR=$((2 + $BASE_K))

# make directory for each removed lithium
if [ ! -d Li1-2xFeO2_KP$BASE_K ]; then
mkdir Li1-2xFeO2_KP$BASE_K
fi

cd Li1-2xFeO2_KP$BASE_K

# make the initial skeleton file for lithium to be removed
cat >  Li1-2xFeO2_KP$BASE_K.in << EOF
&CONTROL
   title = "LiFeO2_Li$BASE_K",
   restart_mode = 'from_scratch',
!  restart_mode = 'restart',  
   calculation = 'scf',
   pseudo_dir = '/mnt/home/davis101/pot',
   outdir = './scratch',
   prefix = "lifeo2_$BASE_K",
   etot_conv_thr = 1.0D-5,
   forc_conv_thr = 1.0D-4
   verbosity = 'high',
   wf_collect = .true.,
   max_seconds = 13500
 /

 &SYSTEM
   ibrav = 0,
   nat = 14,
   ntyp = 6,
   starting_magnetization(1) = 0,
   starting_magnetization(2) = 1,
   starting_magnetization(3) = -1,
   starting_magnetization(4) = -1,
   starting_magnetization(5) = 1,
   starting_magnetization(6) = 0,
   ecutwfc = 130,
   ecutrho = 1040,
   nspin = 2,
   occupations = 'smearing',
   smearing = 'gaussian',
   degauss = 2.0D-3
 /

 &ELECTRONS
   mixing_beta = 0.7,
   electron_maxstep = 200,
   conv_thr = 1.D-8
 /

! &IONS
!   ion_dynamics = 'bfgs'
! /

! &CELL
!   cell_dynamics = 'bfgs',
!   cell_dofree = 'xyz'
! /

ATOMIC_SPECIES
 Li  6.94 Li.pbe-s-kjpaw_psl.0.2.1.UPF
 Fe1 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe2 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe3 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe4 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 O   15.999 O.pbe-n-kjpaw_psl.0.1.UPF

CELL_PARAMETERS (angstrom)
   5.331350825   0.000000000   0.000000000
   0.000000000   6.695312159   0.000000000
   0.000000000   0.000000000   5.162466421

ATOMIC_POSITIONS (crystal)
Li       0.593649091   0.863314807   0.000477450
Li       0.906341919   0.363321766   0.500454101
Fe1      0.081662996   0.117212822   0.010290051
Fe2      0.918288162   0.862888670   0.502984310
Fe3      0.581723941   0.362874221   0.002996094
Fe4      0.418338163   0.617121777   0.510323028
O        0.010977569   0.095230577   0.359247630
O        0.939704953   0.904211108   0.840607741
O        0.560279460   0.404133307   0.340639488
O        0.489010125   0.595223534   0.859285161
O        0.422181667   0.133619209   0.949175725
O        0.609986978   0.825268456   0.381284119
O        0.890027180   0.325339969   0.881229803
O        0.077827781   0.633639747   0.449173146

K_POINTS (automatic)
$X_DIR $Y_DIR $Z_DIR 0 0 0

EOF
VOLK=$(($X_DIR * $Y_DIR * $Z_DIR))
numnodes=$(($VOLK / 2))
gbyte="$((4 * $numnodes))gb"
npvalue=$((4 * $numnodes))
# designate working directory
fold=$(pwd)

cat > Li1-2xFeO2_KP$BASE_K.qsub << EOF
#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=$numnodes:ppn=4,mem=$gbyte
#PBS -N Li1-2xFeO2_KP$BASE_K
#PBS -m abe

# Source Comands
source /mnt/home/davis101/.bashrc
module swap GNU Intel
module load FFTW
module list

# change to the working directory where your code is located
cd $fold

# set the number of OpenMP threads
export OMP_NUM_THREADS=4

mpirun -np $npvalue pw.x -nk $numnodes -nt 4 <Li1-2xFeO2_KP$BASE_K.in>>Li1-2xFeO2_KP$BASE_K.out

EOF

#submit file

qsub Li1-2xFeO2_KP$BASE_K.qsub

cd ..

done

#!/bin/bash
############################################################################
#
# adding lithium to T-LiFeO2, two lithum per calculation
# ECUT = 130
# RHOCUT = 1040
# 210 KPT
# edited for spelling errors
############################################################################
fun1() {
	# make directory for each removed lithium
if [ ! -d Li15FeO2_AddLi"$atom1""_""$atom2"_CalcFold ]; then
mkdir Li15FeO2_AddLi"$atom1""_""$atom2"_CalcFold
fi
	# make the initial skeleton file for lithium to be removed
cat > Li15FeO2_Li_"$atom1""_""$atom2"_PreEdit.in << EOF
&CONTROL
   title = "Li15FeO2_AddLi$atom1$atom2",
   restart_mode = 'from_scratch',
!  restart_mode = 'restart',  
   calculation = 'vc-relax',
   pseudo_dir = '/mnt/home/davis101/pot',
   outdir = './scratch',
   prefix = "li15feo2_$atom1$atom2",
   etot_conv_thr = 1.0D-5,
   forc_conv_thr = 1.0D-4
   verbosity = 'high',
   wf_collect = .true.,
   max_seconds = 12600
 /

 &SYSTEM
   ibrav = 0,
   nat = 18,
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

 &IONS
   ion_dynamics = 'bfgs'
 /

 &CELL
   cell_dynamics = 'bfgs',
   cell_dofree = 'xyz'
 /

ATOMIC_SPECIES
 Li  6.94 Li.pbe-s-kjpaw_psl.0.2.1.UPF
 Fe1 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe2 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe3 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 Fe4 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 O   15.999 O.pbe-n-kjpaw_psl.0.1.UPF

CELL_PARAMETERS {angstrom}
5.51600 0.00000 0.00000
0.00000 6.41390 0.00000
0.00000 0.00000 5.07890

ATOMIC_POSITIONS (crystal)
Li     0.251930000         0.881220000         0.124450000
Li     0.248070000         0.381210000         0.624450000
Li     0.751930000         0.618790000         0.124450000
Li     0.748070000         0.118780000         0.624450000
Li     0.413300008         0.198300004         0.546299994
Li     0.586699963         0.801699996         0.046299994
Li     0.913300037         0.301699996         0.546299994
Li     0.086699992         0.698300004         0.046299994
Fe1    0.077100001         0.116949998         0.000000000
Fe2    0.922900021         0.883050025         0.500000000
Fe3    0.577099979         0.383049995         0.000000000
Fe4    0.422899991         0.616949975         0.500000000
O      0.067199998         0.133811995         0.358011991
O      0.932799995         0.866187990         0.858011961
O      0.567200005         0.366187990         0.358011991
O      0.432799995         0.633812010         0.858011961
O      0.421200007         0.153510004         0.890879989
O      0.578799963         0.846490026         0.390879989
O      0.921200037         0.346489996         0.890879989
O      0.078799993         0.653509974         0.390879989


K_POINTS (automatic)
6 5 7 0 0 0

EOF

# change to correct folder

cd ./Li15FeO2_AddLi"$atom1""_""$atom2"_CalcFold/

#edit the file to remove a specific lithium

sed "$((62 + $atom1))d;$((62 + $atom2))d" ../Li15FeO2_Li_"$atom1""_""$atom2"_PreEdit.in >  Li15FeO2_AddLi_"$atom1"-"$atom2".in

# designate working directory
fold=$(pwd)

cat > Li15FeO2_AddLi_"$atom1"-"$atom2".qsub << EOF
#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=105:ppn=4,mem=420gb
#PBS -N Li15FeO2_Li_$atom1-$atom2
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

mpirun -np 420 pw.x -nk 105 -nt 4 <Li15FeO2_AddLi_$atom1-$atom2.in>>Li15FeO2_AddLi_$atom1-$atom2.out

EOF

#submit it
qsub Li15FeO2_AddLi_"$atom1"-"$atom2".qsub
cd ..
rm Li15FeO2_Li_"$atom1""_""$atom2"_PreEdit.in
}


for atom1 in 1 2 3 4; do
  if [ "$atom1" -eq 1 ]; then
		for atom2 in 2 3 4; do
       fun1
		   done
 	fi	
  if [ "$atom1" -eq 2 ]; then
		for atom2 in 3 4; do
       fun1
	     done
  fi
  if [ "$atom1" -eq 3 ]; then
		for atom2 in 4; do
       fun1
	     done
  fi


done

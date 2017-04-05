#!/bin/bash
############################################################################
#
#removing lithium from T-LiFeO2, only one lithum per calculation
# ECUT = 130
# RHOCUT = 1040
# 125 KPT
# edited for spelling errors
############################################################################
fun1() {
# make directory for each removed lithium
if [ ! -d LiFeO2_Li"$atom1""_""$atom2"_KPUpdate ]; then
mkdir LiFeO2_Li"$atom1""_""$atom2"_KPUpdate
fi

# change to correct folder
cd ./LiFeO2_Li"$atom1""_""$atom2"_KPUpdate/

# make the initial skeleton file for lithium to be removed
cat > LiFeO2_Li_"$atom1"-"$atom2".in<< EOF
&CONTROL
   title = "LiFeO2_Li$atom1$atom2",
   restart_mode = 'from_scratch',
!  restart_mode = 'restart',  
   calculation = 'vc-relax',
   pseudo_dir = '/mnt/home/davis101/pot',
   outdir = './scratch',
   prefix = "lifeo2_$atom1$atom2",
   etot_conv_thr = 1.0D-5,
   forc_conv_thr = 1.0D-4
   verbosity = 'high',
   wf_collect = .true.,
   max_seconds = 84600
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

K_POINTS (automatic)
6 5 7 0 0 0




EOF

#edit the file to remove a specific lithium

#sed "$((62 + $atom1))d;$((62 + $atom2))d" LiFeO2_Li_"$atom1""_""$atom2"_PreEdit.in >  LiFeO2_Li_"$atom1"-"$atom2".in

# designate working directory
fold=$(pwd)

cat > LiFeO2_Li_"$atom1"-"$atom2".qsub << EOF
#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=105:ppn=4,mem=420gb
#PBS -N LiFeO2_Li_$atom1-$atom2

# Source Comands
source /mnt/home/davis101/.bashrc
module swap GNU Intel
module load FFTW
module list

# change to the working directory where your code is located
cd $fold

# set the number of OpenMP threads
export OMP_NUM_THREADS=4

mpirun -np 420 pw.x -nk 105 -nt 4 <LiFeO2_Li_$atom1-$atom2.in>>LiFeO2_Li_$atom1-$atom2.out

EOF

#submit it
#qsub LiFeO2_Li_"$atom1"-"$atom2".qsub

#rm LiFeO2_Li_"$atom1""_""$atom2"_PreEdit.in
cd ..
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

#!/bin/bash
##################################################
# This is a k-points convergence 
# of LiFeO2 Ferromagnetic
##################################################

for BASE_K in 1 2 3 4 5 6 7 8 
do
rm -rf ./scratch/*
X_DIR=$((1 + $BASE_K))
Y_DIR=$BASE_K
Z_DIR=$((2 + $BASE_K))
cat > LiFeO2_kp$BASE_K.in << EOF
 &CONTROL
   title = 'LiFeO2ECUTtest',
   calculation = 'vc-relax',
   pseudo_dir = '../pot',
   outdir = './scratch',
   prefix = 'lifeo2ecut',
   etot_conv_thr = 1.0D-5,
   forc_conv_thr = 1.0D-4
 /

 &SYSTEM
   ibrav = 8,
   a = 5.51600,
   b = 6.41390,
   c = 5.07890,
   nat = 16,
   ntyp = 3,
   starting_magnetization(1) = 0,
   starting_magnetization(2) = 6,
   starting_magnetization(3) = 0,
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
 Li 6.94 Li.pbe-s-kjpaw_psl.0.2.1.UPF
 Fe 55.845 Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
 O 15.999 O.pbe-n-kjpaw_psl.0.1.UPF

ATOMIC_POSITIONS (crystal)
Li     0.413300008         0.198300004         0.546299994
Li     0.586699963         0.801699996         0.046299994
Li     0.913300037         0.301699996         0.546299994
Li     0.086699992         0.698300004         0.046299994
Fe     0.077100001         0.116949998         0.000000000
Fe     0.922900021         0.883050025         0.500000000
Fe     0.577099979         0.383049995         0.000000000
Fe     0.422899991         0.616949975         0.500000000
O      0.067199998         0.133811995         0.358011991
O      0.932799995         0.866187990         0.858011961
O      0.567200005         0.366187990         0.358011991
O      0.432799995         0.633812010         0.858011961
O      0.421200007         0.153510004         0.890879989
O      0.578799963         0.846490026         0.390879989
O      0.921200037         0.346489996         0.890879989
O      0.078799993         0.653509974         0.390879989


K_POINTS (automatic)
$X_DIR $Y_DIR $Z_DIR 0 0 0

EOF

echo $BASE_K Started...

mpirun -np 8 pw.x -nt 8 -nd 8 < LiFeO2_kp$BASE_K.in > LiFeO2_kp$BASE_K.out

echo $BASE_K Done

done

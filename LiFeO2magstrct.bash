#!/bin/bash
##################################################
# This is a magnetic structure test for
# of LiFeO2
##################################################
for FE1 in 1 -1
do
for FE2 in 1 -1
do
for FE3 in 1 -1 
do
for FE4 in 1 -1
do
rm -rf ./scratch/*
cat > LiFeO2magstrct_"$FE1"_"$FE2"_"$FE3"_"$FE4"_.in << EOF
 &CONTROL
   title = 'LiFeO2magstruct',
   calculation = 'vc-relax',
   pseudo_dir = '../pot',
   outdir = './scratch',
   prefix = 'lifeo2mag',
   etot_conv_thr = 1.0D-5,
   forc_conv_thr = 1.0D-4
 /

 &SYSTEM
   ibrav = 0,
   nat = 16,
   ntyp = 6,
   starting_magnetization(1) = 0,
   starting_magnetization(2) = $FE1,
   starting_magnetization(3) = $FE2,
   starting_magnetization(4) = $FE3,
   starting_magnetization(5) = $FE4,
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
4 3 5 0 0 0

EOF
echo "$FE1"_"$FE2"_"$FE3"_"$FE4" Started...

mpirun -np 24 pw.x -nk 6 -nt 4 -nd 16 < LiFeO2magstrct_"$FE1"_"$FE2"_"$FE3"_"$FE4"_.in > LiFeO2magstrct_"$FE1"_"$FE2"_"$FE3"_"$FE4"_.out

echo "$FE1"_"$FE2"_"$FE3"_"$FE4" Finished

done
done
done
done

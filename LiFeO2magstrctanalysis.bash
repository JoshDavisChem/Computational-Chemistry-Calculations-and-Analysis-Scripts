#!/bin/bash
##################################################
# This is a magnetic structure test for
# of LiFeO2 (Analysis
##################################################
for FE1 in 1 -1
do
for FE2 in 1 -1
do
for FE3 in 1 -1
do
for FE4 in 1 -1
do

echo "$FE1"_"$FE2"_"$FE3"_"$FE4"; tac LiFeO2magstrct_"$FE1"_"$FE2"_"$FE3"_"$FE4"_.out | grep -m1 !

done
done
done
done

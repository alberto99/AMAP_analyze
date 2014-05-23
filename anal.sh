#echo -n ALA
#/cavern/alberto/Ffield/OPENMM/Scripts/LR.py ALA/phi.txt ALA/psi.txt  |grep "List" |tail -1
for a in  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  SER  THR  TRP  TYR  VAL PRO
do
echo -n $a
/cavern/alberto/Ffield/OPENMM/Scripts/LR_mix.py $a/phi.txt $a/psi.txt  |grep "List" |tail -1
done

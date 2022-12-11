#!/bin/bash
### Date: 2022-09
### Author: Junpeng Yuan
### Address: Zhengzhou China
### Email: junpeng_yuan@163.com
### Function: The local atomic descriptor is obtained
### Usage: sh getlocalAtom.sh

cellLengthARange=1
cellLengthBRange=1
cellLengthCRange=1

###
location=$(pwd)

adsorbentName=PorousMOF.cif
adsorbateName=C6H6.pdb


# STEP 1
############################################################
cifName=PorousMOF
moleculeName=C6H6

###  UnitCell
HeadCif(){

headCifNumber=$(cat -n $adsorbentName |grep "_atom_site_charge" |awk '{print $1}')

head -n $headCifNumber $adsorbentName > $location/tempcifxyz1

cifLengthA=$(grep "_cell_length_a" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')
cifLengthB=$(grep "_cell_length_b" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')
cifLengthC=$(grep "_cell_length_c" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')

cifAngleAlpha=$(grep "_cell_angle_alpha" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')
cifAngleBeta=$(grep "_cell_angle_beta" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')
cifAngleGamma=$(grep "_cell_angle_gamma" "$cifName".cif |awk '{print $2}' |tr -d "()" |tr -d '\r')

}

###  UnitCell
TailCif(){

CifNumber=$(cat $adsorbentName |wc -l)

TailCifNumber=$(echo $CifNumber $headCifNumber |awk '{print $1 - $2}')

tail -n $TailCifNumber $adsorbentName |grep "0." > $location/tempcifxyz2

}

###  UnitCell
InfoSupercell(){

supercellLengthA=$(echo $cifLengthA $cellLengthARange |awk '{print $1 * $2}')
supercellLengthB=$(echo $cifLengthB $cellLengthBRange |awk '{print $1 * $2}')
supercellLengthC=$(echo $cifLengthC $cellLengthCRange |awk '{print $1 * $2}')

cat > "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-1.cif <<!
data_$cifName

_audit_creation_method RASPA-1.0
_audit_creation_date 2022-9-1
_audit_author_name ''

_cell_length_a    $supercellLengthA
_cell_length_b    $supercellLengthB
_cell_length_c    $supercellLengthC
_cell_angle_alpha $cifAngleAlpha
_cell_angle_beta  $cifAngleBeta
_cell_angle_gamma $cifAngleGamma

_symmetry_cell_setting          triclinic
_symmetry_space_group_name_Hall 'P 1'
_symmetry_space_group_name_H-M  'P 1'
_symmetry_Int_Tables_number     1

_symmetry_equiv_pos_as_xyz 'x,y,z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_charge
!

}

###  UnitCell
Unitcell2Supercell(){

xyzNumber=$(cat "$location"/$1 |wc -l)


for ((n=1; n<=$xyzNumber; n++))
do

for ((i=1; i<=$cellLengthARange; i++))
do

for ((j=1; j<=$cellLengthBRange; j++))
do

for ((k=1; k<=$cellLengthCRange; k++))
do

label=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $1}')
symbol=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $2}')
cifXfraction=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $3}')
cifYfraction=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $4}')
cifZfraction=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $5}')
charge=$(head -n $n "$location"/$1 |tail -n 1 |awk '{print $6}')

supercellXfraction=$(echo $i 1 $cifXfraction $cellLengthARange |awk '{print (($1 - $2) / $4) + ($3 / $4)}')
supercellYfraction=$(echo $j 1 $cifYfraction $cellLengthBRange |awk '{print (($1 - $2) / $4) + ($3 / $4)}')
supercellZfraction=$(echo $k 1 $cifZfraction $cellLengthCRange |awk '{print (($1 - $2) / $4) + ($3 / $4)}')

echo "$label $symbol $supercellXfraction $supercellYfraction $supercellZfraction $charge" |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-2.cif

done

done

done

echo "$cifName: $n"

done
}


###
Cif2Xyz(){

xyzNumber=$(cat "$location"/$1 |wc -l)

for ((j=1; j<=$xyzNumber; j++))
do
label=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $1}')
symbol=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $2}')
xfraction=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $3}')
yfraction=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $4}')
zfraction=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $5}')
charge=$(head -n $j "$location"/$1 |tail -n 1 |awk '{print $6}')

cat > "$location"/cifxyz2xyz.py <<!
from math import *

aLength=$cifLengthA
bLength=$cifLengthB
cLength=$cifLengthC
alphaAngle=$cifAngleAlpha / 180 * pi
betaAngle=$cifAngleBeta / 180 * pi
gammaAngle=$cifAngleGamma / 180 * pi

ax = aLength
ay = 0.0000000000
az = 0.0000000000
bx = bLength*(cos(gammaAngle))
by = bLength*(sin(gammaAngle))
bz = 0.0000000000
cx = cLength*(cos(betaAngle))
cy = (cLength*(cos(alphaAngle)) - (cos(gammaAngle))*cLength*(cos(betaAngle)))/(sin(gammaAngle))
cz = sqrt(cLength*cLength - cLength*(cos(betaAngle))*cLength*(cos(betaAngle)) - \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle))))* \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle)))))

x = $xfraction * ax + $yfraction * bx + $zfraction * cx
y = $xfraction * ay + $yfraction * by + $zfraction * cy
z = $xfraction * az + $yfraction * bz + $zfraction * cz

print ('%.9f' % x, ' ', '%.9f' % y, '%.9f' % z)
!

done

}

### Molecule
HeadPdb(){

pdbNumber=$(cat $adsorbateName |wc -l)

headPdbNumber=$(cat $adsorbateName -n |grep "ENDMDL" |tail -n 2 |head -n 1 |awk '{print $1}')

tailPdbNumber=$(echo $pdbNumber $headPdbNumber |awk '{print $1 - $2}')

tail -n $tailPdbNumber $adsorbateName > $location/"$moleculeName"1.pdb

grep "CRYST1" $location/"$moleculeName"1.pdb > $location/temppdbxyz1

cifLengthA=$(cat $location/temppdbxyz1 |awk '{print $2}')
cifLengthB=$(cat $location/temppdbxyz1 |awk '{print $3}')
cifLengthC=$(cat $location/temppdbxyz1 |awk '{print $4}')

cifAngleAlpha=$(cat $location/temppdbxyz1 |awk '{print $5}')
cifAngleBeta=$(cat $location/temppdbxyz1 |awk '{print $6}')
cifAngleGamma=$(cat $location/temppdbxyz1 |awk '{print $7}')

}

### Molecule
TailPdb(){

grep "ATOM" $location/"$moleculeName"1.pdb |awk '{print $3, $3, $5, $6, $7, $9}' |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' > $location/temppdbxyz2

}

### Molecule
Pdb2cif(){

xyzNumber=$(cat "$location"/$1 |wc -l)

for ((k=1; k<=$xyzNumber; k++))
do
label=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $1}')
symbol=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $2}')
xcartesian=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $3}')
ycartesian=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $4}')
zcartesian=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $5}')
charge=$(head -n $k "$location"/$1 |tail -n 1 |awk '{print $6}')

cat > "$location"/pdb2cif.py <<!
from math import *

import numpy as np

aLength=$cifLengthA
bLength=$cifLengthB
cLength=$cifLengthC
alphaAngle=$cifAngleAlpha / 180 * pi
betaAngle=$cifAngleBeta / 180 * pi
gammaAngle=$cifAngleGamma / 180 * pi

ax = aLength
ay = 0.0000000000
az = 0.0000000000
bx = bLength*(cos(gammaAngle))
by = bLength*(sin(gammaAngle))
bz = 0.0000000000
cx = cLength*(cos(betaAngle))
cy = (cLength*(cos(alphaAngle)) - (cos(gammaAngle))*cLength*(cos(betaAngle)))/(sin(gammaAngle))
cz = sqrt(cLength*cLength - cLength*(cos(betaAngle))*cLength*(cos(betaAngle)) - \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle))))* \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle)))))

list1 = [[ax, bx, cx], [ay, by, cy], [az, bz, cz]]

#matrix1 = np.matrix(list1)

imatrix1 = np.linalg.inv(list1)

xyz0 = [$xcartesian, $ycartesian, $zcartesian]

xyz1 = np.dot(imatrix1, xyz0)

print('%.9f' % xyz1[0], '%.9f' % xyz1[1], '%.9f' % xyz1[2])
!

python "$location"/pdb2cif.py >> tempcifxyz3.2

done

cat temppdbxyz2 |awk '{print $1, $2}' > tempcifxyz3.1
cat temppdbxyz2 |awk '{print $6}' > tempcifxyz3.3
paste tempcifxyz3.1 tempcifxyz3.2 tempcifxyz3.3 |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' > $location/tempcifxyz3
rm tempcifxyz3.1 tempcifxyz3.2 tempcifxyz3.3

}

### Molecule
Molecule2Supercell(){

cat > molecule2supercell.py <<!
import numpy as np
from math import *

###
f = open('temppdbxyz1', mode='r')
line = f.readline()
list1 = []
while line:
    LengthAngle = line.split()
    LengthA = float(LengthAngle[1])
    LengthB = float(LengthAngle[2])
    LengthC = float(LengthAngle[3])
    AngleAlpha = float(LengthAngle[4])
    AngleBeta = float(LengthAngle[5])
    AngleGamma = float(LengthAngle[6])
    line = f.readline()
f.close()

###
f = open('temppdbxyz2', mode='r')
line = f.readline()
list2 = []
x=0
y=0
z=0
while line:
    a = line.split()
    if (a[0] == 'C'):
        x = float(a[2]) + x
        y = float(a[3]) + y
        z = float(a[4]) + z
    line = f.readline()
mx = x/6
my = y/6
mz = z/6
f.close()

###
aLength=LengthA
bLength=LengthB
cLength=LengthC
alphaAngle=AngleAlpha / 180 * pi
betaAngle=AngleBeta / 180 * pi
gammaAngle=AngleGamma / 180 * pi

ax = aLength
ay = 0.0000000000
az = 0.0000000000
bx = bLength*(cos(gammaAngle))
by = bLength*(sin(gammaAngle))
bz = 0.0000000000
cx = cLength*(cos(betaAngle))
cy = (cLength*(cos(alphaAngle)) - (cos(gammaAngle))*cLength*(cos(betaAngle)))/(sin(gammaAngle))
cz = sqrt(cLength*cLength - cLength*(cos(betaAngle))*cLength*(cos(betaAngle)) - \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle))))* \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle)))))

list3 = [[ax, bx, cx], [ay, by, cy], [az, bz, cz]]

ucA = [ax, bx, cx]
ucB = [ay, by, cy]
ucC = [az, bz, cz]
mxyz = [mx, my, mz]

benzenedistance=7

## AB
axial = np.cross(ucA,ucB).tolist()
c2 = sqrt(np.sum([i ** 2 for i in axial]))
eaxial = [x/c2 for x in axial]
distance = np.dot(eaxial, ucC)
distanceC = np.dot(eaxial, mxyz)

if (distanceC >= benzenedistance and distanceC <= (distance - benzenedistance)):
    cellLengthCRange=1
    k=1
    print("cellLengthCRange: 1")
elif (distanceC < 7):
    cellLengthCRange=2
    k=2
    print("cellLengthCRange: 2")
else:
    cellLengthCRange=2
    k=2
    print("cellLengthCRange: 2")

## BC
axial = np.cross(ucB,ucC).tolist()
c2 = sqrt(np.sum([i ** 2 for i in axial]))
eaxial = [x/c2 for x in axial]
distance = np.dot(eaxial, ucA)
distanceA = np.dot(eaxial, mxyz)

if (distanceA >= benzenedistance and distanceA <= (distance - benzenedistance)):
    cellLengthARange=1
    i=1
    print("cellLengthARange: 1")
elif (distanceA < 7):
    cellLengthARange=2
    i=2
    print("cellLengthARange: 2")
else:
    cellLengthARange=2
    i=2
    print("cellLengthARange: 2")

## CA
axial = np.cross(ucC,ucA).tolist()
c2 = sqrt(np.sum([i ** 2 for i in axial]))
eaxial = [x/c2 for x in axial]
distance = np.dot(eaxial, ucB)
distanceB = np.dot(eaxial, mxyz)

if (distanceB >= benzenedistance and distanceB <= (distance - benzenedistance)):
    cellLengthBRange=1
    j=1
    print("cellLengthBRange: 1")
elif (distanceB < 7):
    cellLengthBRange=2
    j=2
    print("cellLengthBRange: 2")
else:
    cellLengthBRange=2
    j=2
    print("cellLengthBRange: 2") 

###
fileNameCif=f"$moleculeName{i}{j}{k}.cif"
fileNamePdb=f"$moleculeName{i}{j}{k}.pdb"
file1 = open('tempcifxyz3', mode='r')
line = file1.readline()
while line:
    mclist = line.split()
    cifmcx = ((i - 1) / cellLengthARange) + (float(mclist[2]) / cellLengthARange)
    cifmcy = ((j - 1) / cellLengthBRange) + (float(mclist[3]) / cellLengthBRange)
    cifmcz = ((k - 1) / cellLengthCRange) + (float(mclist[4]) / cellLengthCRange)
    list4=[mclist[0],mclist[1],str(cifmcx),str(cifmcy),str(cifmcz),mclist[5]]
    with open(fileNameCif,'a') as file2:
        file2.write("     ".join(list4))
        file2.write("\n")
    pdbmcx = (cifmcx * ax) + (cifmcy * bx) + (cifmcz * cx)
    pdbmcy = (cifmcx * ay) + (cifmcy * by) + (cifmcz * cy)
    pdbmcz = (cifmcx * az) + (cifmcy * bz) + (cifmcz * cz)
    list5=[mclist[0],mclist[1],str(pdbmcx),str(pdbmcy),str(pdbmcz),mclist[5]]
    with open(fileNamePdb,'a') as file3:
        file3.write("     ".join(list5))
        file3.write("\n")

    line = file1.readline()
file1.close()
!

python molecule2supercell.py > cellSupercellInfo

}

#for(())
#do
step1Unitcell(){
echo "UnitCell"
HeadCif
TailCif

cellLengthARange=$(grep "cellLengthARange" cellSupercellInfo |awk '{print $2}')
cellLengthBRange=$(grep "cellLengthBRange" cellSupercellInfo |awk '{print $2}')
cellLengthCRange=$(grep "cellLengthCRange" cellSupercellInfo |awk '{print $2}')

InfoSupercell
Unitcell2Supercell tempcifxyz2
}
#done

step1Molecule(){
echo "Molecule"
HeadPdb
TailPdb
Pdb2cif temppdbxyz2
Molecule2Supercell
}

step1MoleculeUnitcell(){
echo "step1: yes"
cat "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-1.cif > step1"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
cat "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-2.cif >> step1"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step1"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
}

step1Molecule
step1Unitcell
step1MoleculeUnitcell


# STEP 2
############################################################
step2InteractLocalAtom(){
cat > localAtom.py <<!

import numpy as np
from math import *

f = open('$cifName$cellLengthARange$cellLengthBRange$cellLengthCRange-1.cif', mode='r')
line = f.readline()
counts = 1

while line:
    LengthAngle = line.split()
    if counts == 7:
        LengthA = float(LengthAngle[1])
    if counts == 8:
        LengthB = float(LengthAngle[1])
    if counts == 9:
        LengthC = float(LengthAngle[1])
    if counts == 10:
        AngleAlpha = float(LengthAngle[1])
    if counts == 11:
        AngleBeta = float(LengthAngle[1])
    if counts == 12:
        AngleGamma = float(LengthAngle[1])
    if counts >= 13:
        break
    
    line = f.readline()
    counts = counts + 1

f.close()

aLength=LengthA
bLength=LengthB
cLength=LengthC
alphaAngle=AngleAlpha / 180 * pi
betaAngle=AngleBeta / 180 * pi
gammaAngle=AngleGamma / 180 * pi

ax = aLength
ay = 0.0000000000
az = 0.0000000000
bx = bLength*(cos(gammaAngle))
by = bLength*(sin(gammaAngle))
bz = 0.0000000000
cx = cLength*(cos(betaAngle))
cy = (cLength*(cos(alphaAngle)) - (cos(gammaAngle))*cLength*(cos(betaAngle)))/(sin(gammaAngle))
cz = sqrt(cLength*cLength - cLength*(cos(betaAngle))*cLength*(cos(betaAngle)) - \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle))))* \
    ((cLength*bLength*(cos(alphaAngle)) - bLength*(cos(gammaAngle))*cLength*(cos(betaAngle)))/(bLength*(sin(gammaAngle)))))

###
f = open('$moleculeName$cellLengthARange$cellLengthBRange$cellLengthCRange.cif', mode='r')
line = f.readline()
x=0
y=0
z=0
while line:
    mclist = line.split()
    if (mclist[0] == 'C'):
        x = float(mclist[2]) + x
        y = float(mclist[3]) + y
        z = float(mclist[4]) + z
    line = f.readline()
cifmcx = x/6
cifmcy = y/6
cifmcz = z/6
f.close()

pdbmcx = (cifmcx * ax) + (cifmcy * bx) + (cifmcz * cx)
pdbmcy = (cifmcx * ay) + (cifmcy * by) + (cifmcz * cy)
pdbmcz = (cifmcx * az) + (cifmcy * bz) + (cifmcz * cz)

###
benzeneDistance = 7
fileNameCif='$cifName-LocalAtom.cif'
fileNamePdb='$cifName-LocalAtom.pdb'
file1 = open('$cifName$cellLengthARange$cellLengthBRange$cellLengthCRange-2.cif', mode='r')
line = file1.readline()
while line:
    uclist = line.split()
    cifucx = float(uclist[2])
    cifucy = float(uclist[3])
    cifucz = float(uclist[4])
    pdbucx = (cifucx * ax) + (cifucy * bx) + (cifucz * cx)
    pdbucy = (cifucx * ay) + (cifucy * by) + (cifucz * cy)
    pdbucz = (cifucx * az) + (cifucy * bz) + (cifucz * cz)
    
    r2 = (pdbucx - pdbmcx)**2 + (pdbucy - pdbmcy)**2 + (pdbucz - pdbmcz)**2 
    
    if (r2 <= (benzeneDistance**2)):
        list5=[uclist[0],uclist[1],str(cifucx),str(cifucy),str(cifucz),uclist[5]]
        with open(fileNameCif,'a') as file3:
            file3.write("     ".join(list5))
            file3.write("\n")
        list6=[uclist[0],uclist[1],str(pdbucx),str(pdbucy),str(pdbucz),uclist[5]]
        with open(fileNamePdb,'a') as file4:
            file4.write("     ".join(list6))
            file4.write("\n")
    line = file1.readline()
file1.close()
!

python localAtom.py
}

step2UnitcelltoMolecule(){

cellLengthARange=$(grep "cellLengthARange" cellSupercellInfo |awk '{print $2}')
cellLengthBRange=$(grep "cellLengthBRange" cellSupercellInfo |awk '{print $2}')
cellLengthCRange=$(grep "cellLengthCRange" cellSupercellInfo |awk '{print $2}')

step2InteractLocalAtom

}

step2MoleculeUnitcell(){
echo "step2: yes"
cat "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-1.cif > step2"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
cat $cifName-LocalAtom.cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step2"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step2"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName".cif
}

step2OutputData(){
cat > AtomRadius.set <<!
H 0.38  
Li 0.86 
Be 0.53 
B 1.01  
C 0.88  
N 0.86  
O 0.89  
F 0.82  
Na 1.15 
Mg 1.28 
Al 1.53 
Si 1.38 
P 1.28  
S 1.20  
Cl 1.17 
K 1.44  
Ca 1.17 
Sc 1.62 
Ti 1.65 
V 1.51  
Cr 1.53 
Mn 1.53 
Fe 1.43 
Co 1.31 
Ni 1.33 
Cu 1.31 
Zn 1.41 
Ga 1.40 
Ge 1.35 
As 1.39 
Se 1.40 
Br 1.39 
Rb 1.65 
Sr 1.30 
Y 1.84  
Zr 1.73 
Nb 1.66 
Mo 1.57 
Ru 1.58 
Rh 1.63 
Pd 1.68 
Ag 1.56 
Cd 1.56 
In 1.53 
Sn 1.64 
Sb 1.64 
Te 1.65 
I  1.58 
Cs 1.85 
Ba 1.52 
La 1.91 
Ce 1.98 
Pr 1.75 
Nd 1.92 
Sm 1.89 
Eu 1.83 
Gd 1.79 
Tb 1.82
Dy 1.79
Ho 1.63
Er 1.80
Tm 1.84
Yb 1.80
Lu 1.86
Hf 1.73
W 1.33
Re 1.29
Ir 1.50
Pt 1.66
Au 1.68
Hg 1.88
Pb 1.72
Bi 1.72
Th 1.97
U 1.76
Np 1.73
Pu 1.71
!

cat > AtomMass.set <<!
H,1.008
He,4.0026
Li,6.94
Be,9.0122
B,10.81
C,12.011
N,14.007
O,15.999
F,18.998
Ne,20.180
Na,22.990
Mg,24.305
Al,26.982
Si,28.085
P,30.974
S,32.06
Cl,35.45
Ar,39.948
K,39.098
Ca,40.078
Sc,44.956
Ti,47.867
V,50.942
Cr,51.996
Mn,54.938
Fe,55.845
Co,58.933
Ni,58.693
Cu,63.546
Zn,65.38
Ga,69.723
Ge,72.630
As,74.922
Se,78.971
Br,79.904
Kr,83.798
Rb,85.468
Sr,87.62
Y,88.906
Zr,91.224
Nb,92.906
Mo,95.95
Tc,98
Ru,101.07
Rh,102.91
Pd,106.42
Ag,107.87
Cd,112.41
In,114.82
Sn,118.71
Sb,121.76
Te,127.60
I,126.90
Xe,131.29
Cs,132.91
Ba,137.33
La,138.91
Ce,140.12
Pr,140.91
Nd,144.24
Pm,145
Sm,150.36
Eu,151.96
Gd,157.25
Tb,158.93
Dy,162.50
Ho,164.93
Er,167.26
Tm,168.93
Yb,173.05
Lu,174.97
Hf,178.49
Ta,180.95
W,183.84
Re,186.21
Os,190.23
Ir,192.22
Pt,195.08
Au,196.97
Hg,200.59
Tl,204.38
Pb,207.2
Bi,208.98
Po,209
At,210
Rn,222
Fr,223
Ra,226
Ac,227
Th,232.04
Pa,231.04
U,238.03
Np,237
Pu,244
Am,243
Cm,247
Bk,247
Cf,251
Es,252
Fm,257
Md,258
No,259
Lr,266
!

echo "descriptor" > LocalAtomDescriptor.dat
echo "geometric" >> LocalAtomDescriptor.dat
TotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |wc -l)
TotalAtomTypeNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |sort -u |wc -l)
NTotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "N" |wc -l)
OTotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "O" |wc -l)
FTotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "F" |wc -l)
STotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "S" |wc -l)
ClTotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "Cl" |wc -l)
BrTotalAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |grep -w "Br" |wc -l)
ElectroNegativityTotalAtomNumber=$(echo $NTotalAtomNumber $OTotalAtomNumber $FTotalAtomNumber $STotalAtomNumber $ClTotalAtomNumber $BrTotalAtomNumber |awk '{print $1 + $2 + $3 + $4 + $5 + $6}')
ElectroNegativityRatio=$(echo $ElectroNegativityTotalAtomNumber $TotalAtomNumber |awk '{print $1 / $2}')
echo "TotalAtomNumber: $TotalAtomNumber" >> LocalAtomDescriptor.dat
echo "TotalAtomTypeNumber: $TotalAtomTypeNumber" >> LocalAtomDescriptor.dat
echo "NTotalAtomNumber: $NTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "OTotalAtomNumber: $OTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "FTotalAtomNumber: $FTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "STotalAtomNumber: $STotalAtomNumber" >> LocalAtomDescriptor.dat
echo "ClTotalAtomNumber: $ClTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "BrTotalAtomNumber: $BrTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "ElectroNegativityTotalAtomNumber: $ElectroNegativityTotalAtomNumber" >> LocalAtomDescriptor.dat
echo "ElectroNegativityRatio: $ElectroNegativityRatio" >> LocalAtomDescriptor.dat

C12kg=1.66e-27
pai=3.141592654
radiim=7
ai2m=1e-10

LocalAtomMass=0
LocalAtomVolume=0

for((i=1; i<="$TotalAtomTypeNumber"; i++))
do

OneAtomName=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |sort |uniq -c |head -n $i |tail -n 1 |awk '{print $2}')
OneAtomNumber=$(cat "$cifName"-LocalAtom.pdb |awk '{print $1}' |sort |uniq -c |head -n $i |tail -n 1 |awk '{print $1}')

OneAtomMass=$(grep -w "$OneAtomName" AtomMass.set |awk -F "," '{print $2}')
OneAtomTypeMass=$(echo $OneAtomMass $C12kg $OneAtomNumber |awk '{print $1 * $2 * $3}')

OneAtomRadius=$(grep -w "$OneAtomName" AtomRadius.set |awk '{print $2}')
OneAtomTypeVolume=$(echo 4 3 $pai $OneAtomRadius $OneAtomNumber |awk '{print (($1 / $2) * $3) * ($4 * $4 * $4) * $5}')

cat > MassVolume.py <<!
LocalAtomMass=$LocalAtomMass
OneAtomTypeMass=$OneAtomTypeMass
LocalAtomVolume=$LocalAtomVolume
OneAtomTypeVolume=$OneAtomTypeVolume
LocalAtomMass=LocalAtomMass+OneAtomTypeMass
LocalAtomVolume=LocalAtomVolume+OneAtomTypeVolume
print("LocalAtomMass: %e" %LocalAtomMass)
print("LocalAtomVolume: %e" %LocalAtomVolume)
!

python MassVolume.py > MassVolume.dat

LocalAtomMass=$(grep "LocalAtomMass" MassVolume.dat |awk '{print $2}')
LocalAtomVolume=$(grep "LocalAtomVolume" MassVolume.dat |awk '{print $2}')

done

LocalAtomNumberDensityA_3=$(echo $TotalAtomNumber 4 3 $pai $radiim |awk '{print $1 / ((($2 / $3) * $4) * ($5 * $5 * $5))}')
LocalAtomDensitykgm3=$(echo $LocalAtomMass 4 3 $pai $ai2m $radiim |awk '{print $1 / ((($2 / $3) * $4) * ($5 * $6 * $5 * $6 * $5 * $6))}')
LocalAtomDensitygcm3=$(echo $LocalAtomDensitykgm3 1000 |awk '{print $1 / $2}')
LocalFreeVolumeA3=$(echo 4 3 $pai $radiim $LocalAtomVolume |awk '{print (($1 / $2) * $3) * ($4 * $4 * $4) - $5}')

echo "LocalAtomNumberDensityA_3: $LocalAtomNumberDensityA_3" >> LocalAtomDescriptor.dat
echo "LocalAtomDensitygcm3: $LocalAtomDensitygcm3" >> LocalAtomDescriptor.dat
echo "LocalFreeVolumeA3: $LocalFreeVolumeA3" >> LocalAtomDescriptor.dat

}

step2UnitcelltoMolecule
step2MoleculeUnitcell
step2OutputData


# STEP 3
############################################################
LocalAtomBCNOSiPS(){
grep -w 'B' $cifName-LocalAtom.cif > LocalAtom.cif
grep -w 'C' $cifName-LocalAtom.cif >> LocalAtom.cif
grep -w 'N' $cifName-LocalAtom.cif >> LocalAtom.cif
grep -w 'O' $cifName-LocalAtom.cif >> LocalAtom.cif
grep -w 'Si' $cifName-LocalAtom.cif >> LocalAtom.cif
grep -w 'P' $cifName-LocalAtom.cif >> LocalAtom.cif
grep -w 'S' $cifName-LocalAtom.cif >> LocalAtom.cif

grep -w 'B' $cifName-LocalAtom.pdb > LocalAtom.pdb
grep -w 'C' $cifName-LocalAtom.pdb >> LocalAtom.pdb
grep -w 'N' $cifName-LocalAtom.pdb >> LocalAtom.pdb
grep -w 'O' $cifName-LocalAtom.pdb >> LocalAtom.pdb
grep -w 'Si' $cifName-LocalAtom.pdb >> LocalAtom.pdb
grep -w 'P' $cifName-LocalAtom.pdb >> LocalAtom.pdb
grep -w 'S' $cifName-LocalAtom.pdb >> LocalAtom.pdb
}

clusterSplit(){
LocalAtomNumber=$(cat LocalAtom.pdb |wc -l)
moleculeN=1

while (($LocalAtomNumber > 0))
do

head -n 1 LocalAtom.pdb > tempNode

tempNodeNumber=$(cat tempNode |wc -l)

while (($tempNodeNumber > 0))
do

atomName=$(head -n 1 tempNode |awk '{print $1}')
atomX0=$(head -n 1 tempNode |awk '{print $3}' |tr -d '\r')
atomY0=$(head -n 1 tempNode |awk '{print $4}' |tr -d '\r')
atomZ0=$(head -n 1 tempNode |awk '{print $5}' |tr -d '\r')

for((i=1; i<="$LocalAtomNumber"; i++))
do

atomName=$(head -n $i LocalAtom.pdb |tail -n 1 |awk '{print $1}')
atomX=$(head -n $i LocalAtom.pdb |tail -n 1 |awk '{print $3}' |tr -d '\r')
atomY=$(head -n $i LocalAtom.pdb |tail -n 1 |awk '{print $4}' |tr -d '\r')
atomZ=$(head -n $i LocalAtom.pdb |tail -n 1 |awk '{print $5}' |tr -d '\r')

cat > R2XYZatom.py <<!
R2 = ($atomX - $atomX0)**2 + ($atomY - $atomY0)**2 + ($atomZ - $atomZ0)**2
if (R2 <= (1.6**2)):
    print("yes")
else:
    print("no")
!

python R2XYZatom.py > atomYes2No

atomState=$(cat atomYes2No)

if [[ $atomState = yes ]]
then
head -n $i LocalAtom.pdb |tail -n 1 >> tempNode
head -n $i LocalAtom.pdb |tail -n 1 >> molecule"$moleculeN".pdb
head -n $i LocalAtom.cif |tail -n 1 >> molecule"$moleculeN".cif
sed -i "$i d" LocalAtom.pdb
sed -i "$i d" LocalAtom.cif
((i--))
((LocalAtomNumber--))
fi
done

sed -i '1d' tempNode

tempNodeNumber=$(cat tempNode |wc -l)

done

echo "molecule$moleculeN: yes"

LocalAtomNumber=$(cat LocalAtom.pdb |wc -l)
((moleculeN++))

done
}

step3AtomSplit(){
LocalAtomBCNOSiPS
clusterSplit
}

step3MoleculeUnitcell(){
echo "step3: yes"

cellLengthARange=$(grep "cellLengthARange" cellSupercellInfo |awk '{print $2}')
cellLengthBRange=$(grep "cellLengthBRange" cellSupercellInfo |awk '{print $2}')
cellLengthCRange=$(grep "cellLengthCRange" cellSupercellInfo |awk '{print $2}')

moleculeNumber=$(ls molecule*.cif |wc -l)

for((i=1; i<="$moleculeNumber"; i++))
do
cat "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-1.cif > step3"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName"molecule"$i".cif
cat molecule"$i".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step3"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName"molecule"$i".cif
cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step3"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName"molecule"$i".cif
done

}

step3AtomSplit
step3MoleculeUnitcell


# STEP 4
############################################################
step4Get6MR(){
moleculeNumber=$(ls molecule*.pdb  |wc -l)

for((mName=1; mName<="$moleculeNumber"; mName++))
do

rm6=1

LocalAtomNumber=$(cat molecule"$mName".pdb |wc -l)

if [[ "$LocalAtomNumber" -lt 6 ]]
then

echo "No 6MR"

else

while (($LocalAtomNumber > 0))
do

atomName0=$(head -n 1 molecule"$mName".pdb |awk '{print $1}')
atomX0=$(head -n 1 molecule"$mName".pdb |awk '{print $3}' |tr -d '\r')
atomY0=$(head -n 1 molecule"$mName".pdb |awk '{print $4}' |tr -d '\r')
atomZ0=$(head -n 1 molecule"$mName".pdb |awk '{print $5}' |tr -d '\r')

for((i=1; i<="$LocalAtomNumber"; i++))
do

atomName=$(head -n $i molecule"$mName".pdb |tail -n 1 |awk '{print $1}')
atomX=$(head -n $i molecule"$mName".pdb |tail -n 1 |awk '{print $3}' |tr -d '\r')
atomY=$(head -n $i molecule"$mName".pdb |tail -n 1 |awk '{print $4}' |tr -d '\r')
atomZ=$(head -n $i molecule"$mName".pdb |tail -n 1 |awk '{print $5}' |tr -d '\r')

cat > R2XYZatom.py <<!
R2 = ($atomX - $atomX0)**2 + ($atomY - $atomY0)**2 + ($atomZ - $atomZ0)**2
if ((2.6**2) <= R2 <= (2.9**2)):
    print("yes")
else:
    print("no")
!

python R2XYZatom.py > atomYes2No

atomState=$(cat atomYes2No)

if [[ $atomState = yes ]]
then

for((j=1; j<="$LocalAtomNumber"; j++))
do

atomName_2=$(head -n $j molecule"$mName".pdb |tail -n 1 |awk '{print $1}')
atomX_2=$(head -n $j molecule"$mName".pdb |tail -n 1 |awk '{print $3}' |tr -d '\r')
atomY_2=$(head -n $j molecule"$mName".pdb |tail -n 1 |awk '{print $4}' |tr -d '\r')
atomZ_2=$(head -n $j molecule"$mName".pdb |tail -n 1 |awk '{print $5}' |tr -d '\r')

cat > R2XYZatom_2.py <<!
Xm=($atomX+$atomX0)/2
Ym=($atomY+$atomY0)/2
Zm=($atomZ+$atomZ0)/2

R2 = ($atomX_2 - Xm)**2 + ($atomY_2 - Ym)**2 + ($atomZ_2 - Zm)**2

if ((1.25**2) <= R2 <= (1.55**2)):
    print("yes")
else:
    print("no")
!

python R2XYZatom_2.py > atomYes2No_2

atomState_2=$(cat atomYes2No_2)

if [[ $atomState_2 = yes ]]
then

head -n $j molecule"$mName".pdb |tail -n 1 >> molecule"$mName"-6MR"$rm6".pdb
head -n $j molecule"$mName".cif |tail -n 1 >> molecule"$mName"-6MR"$rm6".cif

fi

done

mr6Number=$(cat molecule"$mName"-6MR"$rm6".pdb |wc -l)

if [[ $mr6Number -ne 6 ]]
then
rm molecule"$mName"-6MR"$rm6".pdb
rm molecule"$mName"-6MR"$rm6".cif
echo "molecule$mName-6MR$rm6: No"
else
echo "molecule$mName-6MR$rm6: Yes"
((rm6++))
fi

fi

done

sed -i "1 d" molecule"$mName".pdb
sed -i "1 d" molecule"$mName".cif

LocalAtomNumber=$(cat molecule"$mName".pdb |wc -l)
done

fi

done
}


step4MoleculeUnitcell(){
echo "step4: yes"

cellLengthARange=$(grep "cellLengthARange" cellSupercellInfo |awk '{print $2}')
cellLengthBRange=$(grep "cellLengthBRange" cellSupercellInfo |awk '{print $2}')
cellLengthCRange=$(grep "cellLengthCRange" cellSupercellInfo |awk '{print $2}')

#if [[ -f molecule[0-9]-6MR[0-9].cif ]]
#then
#SixMemberedRingNumber=$(ls molecule[0-9]-6MR[0-9].cif |wc -l)
#else
#SixMemberedRingNumber=0
#fi

SixMemberedRingNumber=$(ls molecule*-6MR[0-9].cif 2>/dev/null |wc -l)

for((i=1; i<="$SixMemberedRingNumber"; i++))
do

SixMemberedRingName=$(ls molecule*-6MR[0-9].cif |head -n $i |tail -n 1 |awk -F ".cif" '{print $1}')

cat "$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange"-1.cif > step4"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName""$SixMemberedRingName".cif
cat "$SixMemberedRingName".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step4"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName""$SixMemberedRingName".cif
cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".cif |awk '{printf "%-8s %-3s %17.12f %18.12f %18.12f %12.6f\n", $1,$2,$3,$4,$5,$6}' >> step4"$cifName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange""$moleculeName""$SixMemberedRingName".cif

done
}

step4OutputData(){
echo "paiStacking" >> LocalAtomDescriptor.dat
SandwichTypeNumber=0
TShapeTypeNumber=0
ParallelDisplacedTypeNumber=0

bc1x=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $3}')
bc1y=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $4}')
bc1z=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $5}')

bc2x=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $3}')
bc2y=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $4}')
bc2z=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $5}')

bc3x=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $3}')
bc3y=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $4}')
bc3z=$(cat "$moleculeName""$cellLengthARange""$cellLengthBRange""$cellLengthCRange".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $5}')

SixMemberedRingNumber=$(ls molecule*-6MR[0-9].pdb 2>/dev/null |wc -l)

if [[ "$SixMemberedRingNumber" -ge 1 ]]
then
for((i=1; i<="$SixMemberedRingNumber"; i++))
do

SixMemberedRingName=$(ls molecule*-6MR[0-9].pdb |head -n $i |tail -n 1 |awk -F ".pdb" '{print $1}')

rc1x=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $3}')
rc1y=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $4}')
rc1z=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 1 |tail -n 1 |awk '{print $5}')

rc2x=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $3}')
rc2y=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $4}')
rc2z=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 2 |tail -n 1 |awk '{print $5}')

rc3x=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $3}')
rc3y=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $4}')
rc3z=$(cat "$SixMemberedRingName".pdb |grep -w "C" |head -n 3 |tail -n 1 |awk '{print $5}')

cat > "$SixMemberedRingName".py <<!
import numpy as np
from math import *

f = open('$moleculeName$cellLengthARange$cellLengthBRange$cellLengthCRange.pdb', mode='r')
line = f.readline()
list1 = []
x=0
y=0
z=0
while line:
    a = line.split()
    if (a[0] == 'C'):
        x = float(a[2]) + x
        y = float(a[3]) + y
        z = float(a[4]) + z
    line = f.readline()
f.close()

#
bmx = x/6
bmy = y/6
bmz = z/6
#
bA = [$bc2x-$bc1x,$bc2y-$bc1y,$bc2z-$bc1z]
bB = [$bc3x-$bc1x,$bc3y-$bc1y,$bc3z-$bc1z]
bC = np.cross(bA,bB)
bC2 = sqrt(np.sum([i ** 2 for i in bC]))
be = bC / bC2

f = open('$SixMemberedRingName.pdb', mode='r')
line = f.readline()
list2 = []
x=0
y=0
z=0
while line:
    a = line.split()
    x = float(a[2]) + x
    y = float(a[3]) + y
    z = float(a[4]) + z
    line = f.readline()
f.close()

#
rmx = x/6
rmy = y/6
rmz = z/6
#
rA=[$rc2x-$rc1x,$rc2y-$rc1y,$rc2z-$rc1z]
rB=[$rc3x-$rc1x,$rc3y-$rc1y,$rc3z-$rc1z]
rC = np.cross(rA,rB)
rC2 = sqrt(np.sum([i ** 2 for i in rC]))
re = rC / rC2

#
rbrCoordinate = [rmx-bmx,rmy-bmy,rmz-bmz]
rbr2 = np.sum([i ** 2 for i in rbrCoordinate])
rbh2 = np.dot(be, rbrCoordinate) ** 2
rbdistance = sqrt(rbr2 - rbh2)
rbanglecos = np.dot(be, re) / (sqrt(np.sum([i ** 2 for i in be])) * sqrt(np.sum([i ** 2 for i in re])))
rbanglearccos = np.arccos(rbanglecos)
rbangle = np.degrees(rbanglearccos)

# c-c: 1.392
# h-c-c-h: 2.785
# c-c-c:  

if (0 <= rbdistance <= 2.785 and (0 <= rbangle <= 45 or 135 <= rbangle <= 180)):
    print ("Sandwich")
elif (rbdistance > 2.785 and (0 <= rbangle <= 45 or 135 <= rbangle <= 180)):    
    print("ParallelDisplaced")
else:
    print("TShape")
!

python "$SixMemberedRingName".py >> PaiStackingType

done

fi

SandwichTypeNumber=$(grep "Sandwich" PaiStackingType |wc -l)
ParallelDisplacedTypeNumber=$(grep "ParallelDisplaced" PaiStackingType |wc -l)
TShapeTypeNumber=$(grep "TShape" PaiStackingType |wc -l)
PaiStackingTypeNumber=$(echo $SandwichTypeNumber $ParallelDisplacedTypeNumber $TShapeTypeNumber |awk '{print $1 + $2 + $3}')
if [ $PaiStackingTypeNumber -ne 0 ]
then
SandwichTypeRatio=$(echo $SandwichTypeNumber $PaiStackingTypeNumber |awk '{print $1 / $2}')
ParallelDisplacedTypeRatio=$(echo $ParallelDisplacedTypeNumber $PaiStackingTypeNumber |awk '{print $1 / $2}')
TShapeTypeRatio=$(echo $TShapeTypeNumber $PaiStackingTypeNumber |awk '{print $1 / $2}')
else
SandwichTypeRatio=0
ParallelDisplacedTypeRatio=0
TShapeTypeRatio=0
fi

echo "PaiStackingTypeNumber: $PaiStackingTypeNumber" >> LocalAtomDescriptor.dat
echo "SandwichTypeNumber: $SandwichTypeNumber" >> LocalAtomDescriptor.dat
echo "ParallelDisplacedTypeNumber: $ParallelDisplacedTypeNumber" >> LocalAtomDescriptor.dat
echo "TShapeTypeNumber: $TShapeTypeNumber" >> LocalAtomDescriptor.dat
echo "SandwichTypeRatio: $SandwichTypeRatio" >> LocalAtomDescriptor.dat
echo "ParallelDisplacedTypeRatio: $ParallelDisplacedTypeRatio" >> LocalAtomDescriptor.dat
echo "TShapeTypeRatio: $TShapeTypeRatio" >> LocalAtomDescriptor.dat
}

step4Get6MR
step4MoleculeUnitcell
step4OutputData


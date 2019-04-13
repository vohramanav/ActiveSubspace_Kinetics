#!/bin/bash

path_kinetics=$PWD
path_h2=$path_kinetics/TChem/example/ign-c/run/h2

for k in `seq 1 10000`;
do
#----------------------------UPDATE chem.inp---------------------------------

  y=$(grep H+O2 chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a1
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $1}') # new pce4D_19D_60ue of a1
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a1
  y=$(grep H+O2 chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a1
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $20}') # new pce4D_19D_60ue of a1
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a1

  y=$(grep O+H2 chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a2
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $2}') # new pce4D_19D_60ue of a2
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a2
  y=$(grep O+H2 chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a2
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $21}') # new pce4D_19D_60ue of a2
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a2
 
  y1=$(grep H2+OH chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a3
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $3}') # new pce4D_19D_60ue of a3
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a3
  y1=$(grep H2+OH chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a3
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $22}') # new pce4D_19D_60ue of a3
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a3
 
  y=$(grep OH+OH chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a4
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $4}') # new pce4D_19D_60ue of a4
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a4
  y=$(grep OH+OH chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a4
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $23}') # new pce4D_19D_60ue of a4
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a4
 
  y1=$(grep H2+M chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a5
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $5}') # new pce4D_19D_60ue of a5
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a5
  y1=$(grep H2+M chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a5
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $24}') # new pce4D_19D_60ue of a5
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a5
 
  y1=$(grep O+O+M chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a6
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $6}') # new pce4D_19D_60ue of a6
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a6
  y1=$(grep O+O+M chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a6
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $25}') # new pce4D_19D_60ue of a6
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a6
 
  y1=$(grep O+H+M chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a7
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $7}') # new pce4D_19D_60ue of a7
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a7
  y1=$(grep O+H+M chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a7
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $26}') # new pce4D_19D_60ue of a7
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a7
 
  y=$(grep H+OH+M chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a8
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $8}') # new pce4D_19D_60ue of a8
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a8
  y=$(grep H+OH+M chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a8
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $27}') # new pce4D_19D_60ue of a8
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a8

  y1=$(grep H+O2+M chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a9
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $9}') # new pce4D_19D_60ue of a9
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a9
  y1=$(grep H+O2+M chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a9
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $28}') # new pce4D_19D_60ue of a9
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a9
 
  y=$(grep HO2+H chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a10
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $10}') # new pce4D_19D_60ue of a10
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a10
  y=$(grep HO2+H chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a10
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $29}') # new pce4D_19D_60ue of a10
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a10
 
  y=$(grep HO2+H chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a11
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $11}') # new pce4D_19D_60ue of a11
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a11
  y=$(grep HO2+H chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a11
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $30}') # new pce4D_19D_60ue of a11
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a11
 
  y=$(grep HO2+O chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a12
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $12}') # new pce4D_19D_60ue of a12
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a12
  y=$(grep HO2+O chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a12
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $31}') # new pce4D_19D_60ue of a12
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a12
 
  y1=$(grep HO2+OH chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a13
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $13}') # new pce4D_19D_60ue of a13
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a13
  y1=$(grep HO2+OH chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a13
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $32}') # new pce4D_19D_60ue of a13
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a13
 
  y1=$(grep HO2+HO2 chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a14
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $14}') # new pce4D_19D_60ue of a14
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a14
  y1=$(grep HO2+HO2 chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a14
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $33}') # new pce4D_19D_60ue of a14
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a14
 
  y1=$(grep H2O2+M chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a15
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $15}') # new pce4D_19D_60ue of a15
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a15
  y1=$(grep H2O2+M chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a15
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $34}') # new pce4D_19D_60ue of a15
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a15
 
  y=$(grep H2O2+H chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a16
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $16}') # new pce4D_19D_60ue of a16
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a16
  y=$(grep H2O2+H chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $1}')  # old pce4D_19D_60ue of a16
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $35}') # new pce4D_19D_60ue of a16
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a16
 
  y=$(grep H2O2+H chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a17
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $17}') # new pce4D_19D_60ue of a17
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a17
  y=$(grep H2O2+H chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a17
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $36}') # new pce4D_19D_60ue of a17
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a17
 
  y=$(grep H2O2+O chem.inp | awk '{print $2}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a18
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $18}') # new pce4D_19D_60ue of a18
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a18
  y=$(grep H2O2+O chem.inp | awk '{print $4}')
  y1=$(echo $y | awk '{print $2}')  # old pce4D_19D_60ue of a18
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $37}') # new pce4D_19D_60ue of a18
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a18
 
  y1=$(grep H2O2+OH chem.inp | awk '{print $2}') # old pce4D_19D_60ue of a19
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $19}') # new pce4D_19D_60ue of a19
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a19
  y1=$(grep H2O2+OH chem.inp | awk '{print $4}') # old pce4D_19D_60ue of a19
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $38}') # new pce4D_19D_60ue of a19
  sed -i "s/$y1/$y2/g" "chem.inp" # update pce4D_19D_60ue of a19
  
  dos2unix chem.inp # to remove arriage return symbol (^M)

#------------------------------Update input.dat------------------------------------

  y1=$(grep pfac input.dat | awk '{print $2}') # current P
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $39}') # replacement value for P
  sed -i "s/$y1/$y2/g" "input.dat" # update P in input.dat

  # update T
  y1=$(grep Tini input.dat | awk '{print $2}')
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $40}')
  sed -i "s/$y1/$y2/g" "input.dat" # update T in input.dat

  # update H2 molar ratio
  y=$(grep spec input.dat | awk '{print $3}')
  y1=$(echo $y | awk '{print $1}')
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $41}')
  sed -i "s/$y1/$y2/g" "input.dat"
  
  # update O2 molar ratio
  y=$(grep spec input.dat | awk '{print $3}')
  y1=$(echo $y | awk '{print $2}')
  y2=$(head -$k pts_N1e4.txt | tail -1 | awk '{print $42}')
  sed -i "s/$y1/$y2/g" "input.dat"

  echo $k

  cp chem.inp $path_h2/. # copy chem.inp to target dir
  cp input.dat $path_h2/. # copy updated input.dat to target dir
  cd $path_h2 # cd to targer dir
  ./ign > outfile.dat # run the case and store output in outfile.dat
  id=$(tail -1 outfile.dat) # extract ignition delay from outfile.dat
  rm outfile.dat
  echo $id >> record_id.txt # record ignition delay in a text file record_id.txt
  cd $path_kinetics  # return to the starting point
done

mv $path_h2/record_id.txt $path_kinetics/record_id_N1e4.txt # copy record_id.txt to starting point

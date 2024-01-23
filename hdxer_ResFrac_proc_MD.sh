#!/bin/bash

rm prep_*.pdb *.scrpt temp *.log 

        #HDXer is particular on the PDB input- you can't have headers, and anything that doesn't match the topology also doesn't work

read -e -p "Input Trajectory: " trajin
read -e -p "Input prmtop: " prmtop
read -p "Frame interval (n): " frameint

cframe=1	#start frame count at 1
let nframe=cframe+frameint
#frameint=11
eframe=`cpptraj -p ${prmtop} -y ${trajin} -tl | awk '{print $2}'`
#eframe=5

	conda env list > envs  #pull out the HDXER_ENV variable
	hdxenv=`grep HDXER_ENV envs | tr -d " \t\n\r"`

while  [[ ${cframe} -le ${eframe} ]]; do 
	#HDXer portion

cat>HDXer_${cframe}.sh<<EOF
	conda init bash
	conda activate ${hdxenv}
	
	python $HDXER_PATH/HDXer/calc_hdx.py -t ${trajin} -p ${prmtop} -s ${cframe} -e ${nframe} -m BestVendruscolo -log test.log -out init -seg res_seg.txt -mopt "{ 'save_detailed' : True }" -dt 0.5 1 3 10 30 60 120
EOF

	chmod 755 HDXer_${cframe}.sh
	bash -i HDXer_${cframe}.sh
	wait
		mv initSegment_average_fractions.dat Segment_average_${cframe}
	#awk '{print $1}' initProtection_factors.dat | awk '{for(i=p+1; i<$1; i++) print i} {p=$1}'
		wait
		rm *.tmp *.scrpt temp *.log prep*.pdb init*
			let cframe=cframe+frameint
			wait
			let nframe=cframe+frameint

done

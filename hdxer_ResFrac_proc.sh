#!/bin/bash

rm prep_*.pdb *.scrpt temp *.log 

        #HDXer is particular on the PDB input- you can't have headers, and anything that doesn't match the topology also doesn't work
sed -n -E -i '/ATOM/,$ p' *.pdb

inpdb=`ls *.pdb | head -n 1`

	#convert the files by adding H and making topology 
cat>tleap.scrpt<<EOF
source leaprc.protein.ff14SB
m1 = loadpdb ${inpdb}
saveamberparm m1 model.prmtop model.inpcrd
quit
EOF

	tleap -f tleap.scrpt 

for i in *.pdb; do
cat>convert.scrpt<<EOF
	source leaprc.protein.ff14SB
	m1 = loadpdb ${i}
	savepdb m1 prep_${i}
	quit
EOF
	tleap -f convert.scrpt
done

	#HDXer portion
	
	echo "conda info --envs"  #pull out the HDXER_ENV variable
	hdxenv=`grep HDXER_ENV envs | tr -d " \t\n\r"`
	 
cat>HDXer.sh<<EOF
	conda init bash
	conda activate ${hdxenv}
	
	python $HDXER_PATH/HDXer/calc_hdx.py -t prep_*.pdb -p *.prmtop -m BestVendruscolo -log test.log -out init -seg res_seg.txt -mopt "{ 'save_detailed' : True }" -dt 0.5 1 3 10 30 60 120
EOF

chmod 755 HDXer.sh
bash -i HDXer.sh


        #awk '{print $1}' initProtection_factors.dat | awk '{for(i=p+1; i<$1; i++) print i} {p=$1}'

fragment=hdx_res_seg.txt

echo "" > residue_perc_time

peplist=`wc -l ${fragment} | awk '{print $1}'`  #total count of peptide fragments 
spep=1  #start for the counter

while  [[ $spep -le $peplist ]]; do
        res1=`awk -v a=${spep} 'NR==a {print $1-1}' ${fragment}`
                        #this is done as (i-1, i) because any deuterium labeling of the first residue in each HDX-MS peptide segment is completely back-exchanged to hydrogen after the quenching of labeling and digestion of the intact protein, as the hydrolyzed amide is converted to an extremely labile, N-terminal amine
        res2=`awk -v a=${spep} 'NR==a {print $2}' ${fragment}`  #pull the last residue of the peptide fragment 
        echo $res1 $res2
        wait
                awk 'NR>1 {print $0}' *Residue_fractions.dat | awk -v x=${res1} -v y=${res2} 'NR>=x&&NR<=y' | awk '{$1=""}1'| awk '{for(i=1; i<=NF; i++) {a[i]+=$i; if($i!="") b[i]++}}; END {for(i=1; i<=NF; i++) printf "%s%s", (a[i]/b[i])*100, (i==NF?ORS:OFS)}' >> residue_perc_time
                wait
                        let spep=spep+1
done

        tr -s ' '  '\n' < residue_perc_time > formatted | awk '($1 > 0.01)' formatted > temp

rm *.tmp *.scrpt temp *.log prep*.pdb
#!/bin/bash
export AMBERHOME=/home/sisig/bin/amber16
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AMBERHOME/lib
source /home/sisig/bin/amber16/amber.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-8.0/lib64
export CUDA_VISIBLE_DEVICES=1

rm 05_prod*.in 05_prod*.out *.nc 05_prod*.mdinfo moves R_MC.R

#counter for # ns
final=2600000
begin=1
init=0

#counter for position to track temp
start=1
#counter for move up or down or stationary
#move=0

#know when the list ends
endtemp=`wc -l temp_list | awk '{print $1}'`

while [ ${begin} -le ${final} ]; do
	
	#initial temperature; assume starts at 1
	tempi=`awk -v b=${start} 'NR==b {print $1}' temp_list`
	#determine next position
	#sum=`python -c "print (${start}+${move})"`
	#temp to go or stay at
	temp0=`awk -v b=${start} 'NR==b {print $1}' temp_list`

	#this runs for 0.5 ns steps
cat >05_prod_${begin}.in<<EOF
Production 0.5 ns MD 
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntxo= 2,
  barostat = 2,
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 9.0, 
  ntc = 2, 
  ntf = 2,
  tempi = ${tempi}, 
  temp0 = ${temp0},
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  iwrap=1,
  nstlim = 2500, dt = 0.002,
  ntpr = 1000, ntwx = 1000, ntwr = 25000,
  ntr = 0,
  ntwprt = 2289
 /
EOF

	wait

	#start from some initial point, ie 0 for lowest temp rst
	${AMBERHOME}/bin/pmemd.cuda -O -i 05_prod_${begin}.in -o 05_prod_${begin}.out -p *.prmtop -c 05_prod_${init}.rst -r 05_prod_${begin}.rst -x 05_prod_${begin}.nc -inf 05_prod_${begin}.mdinfo -ref 05_prod_${init}.rst

	#mpirun -np 2 ${AMBERHOME}/bin/sander.MPI -O -i 05_prod_${begin}.in -o 05_prod_${begin}.out -p *.prmtop -c 05_prod_${init}.rst -r 05_prod_${begin}.rst -x 05_prod_${begin}.nc -inf 05_prod_${begin}.mdinfo

	wait

	#echo $start $begin $temp0

	#pull pot energies and append to longer list
	grep "EPtot" 05_prod_${begin}.out | head -n-2 | awk '{print $9}' | awk -v a=${temp0} '{print $1,$2=a}' >> pote_list

	wait

	#the size+lenght of the pot energies accumulated seems to bottleneck it a bit
	#so instead have external loop pull avg potE, last known potE, and last known weight

	#it isn't really great so might come back to later
	
	listgencount=`wc -l temp_list | awk '{print $1}'`
	listgenstart=1

	#remove the old lookuplist
	rm lookuplist

	while [ ${listgenstart} -le ${listgencount} ]; do
		#set lookup temp
		lookupt=`awk -v t=${listgenstart} 'NR==t {print $1}' temp_list | awk '{print $1}'`
		
		#pull last known pote and avg pote
		potavg=`awk -v a=${lookupt} '$2 == a {print  $1}' pote_list | awk '{sum += $1} END {if (NR>0) print sum/NR}'`
		lastpot=`awk -v a=${lookupt} '$2 == a {print  $1}' pote_list | tail -n 1`

		#pull last known weight- is this actually faster here? 
		lastweight=`awk -v a=${lookupt} '$2 == a {print  $1}' initw | tail -n 1`

		#structure of lookup is 
		# avgPotE lastPotE lastWei Temp
		#it does weird stuff otherwise with the formatting 
		echo "${potavg} ${lastpot} ${lastweight} ${lookupt}" >> lookuplist

		let listgenstart=listgenstart+1

	done

	wait

	rm R_MC.R

#run the long R script to determine move/direction/etc
cat >R_MC.R<<EOF
reweight <- function (a,b,c,d,e,f) {
	kb=0.008314
	x1 <<- (1/(kb*c))
	x2 <<- (1/(kb*a))
	x <<- (1/(kb*c))-(1/(kb*a))
	#use avg energies
	y <<- (b+d)/2
	#weight of previous; assume 0 at T1
	w1 <- e
	##check <<- x
	##check2 <<- y
	##check3 <<- e
	w2 <<- w1+x*y
	#print(weight)
	z <<- -1*((x)*f-(w2-w1))
	pmc <<- min(1,exp(z))
	#print(pmc)
}

isEmpty <- function(x) {
    return(length(x)==0)
}

#this got replaced by the above awk loop for a single lookup list
#pelist=NULL
#pelist <- as.matrix(read.table("pote_list"))
#weight=NULL
#weight <- as.matrix(read.table("initw"))
#tlist=NULL
#tlist <- as.matrix(read.table("temp_list"))

lookup <- as.matrix(read.table("lookuplist"))

start=${start}
tlist=${endtemp}
moves=NULL
calcw=NULL
record=NULL

#print("errortrack1")

	if (start == 1) {

		coin <- rbinom(n=1, size=1, prob=0.5)

			if (coin == 1) {
				#attempt move up
				#order does matter

				#grab input- temp is in col 4 
				t1 <- lookup[start,4]
				#make sure there's something there- mean potE is in col 1
					if (isEmpty(lookup[lookup[,4]==t1,1])) {
						print("Missing mean energy value1")
						quit()
					} else {
						#avg potential energies for temp
						m1 <- lookup[lookup[,4]==t1,1]
						#pull last known energy- col 2
						e1 <- lookup[lookup[,4]==t1,2]
					}
				#grab next set from list
				t2 <- lookup[(start+1),4]
				#make sure there's something there
					if (isEmpty(lookup[lookup[,4]==t2,1])) {
						print("Missing mean energy value2")
						quit()
					} else {
						#avg potential energies for temp- col 1
						e2 <- lookup[lookup[,4]==t2,1]
					}			
				#grab last known weight for T1- col 3
				w1 <- lookup[lookup[,4]==t1,3]

				reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
				calcw <- rbind(calcw,cbind(w2,t2))
				write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

					#check if pass possible
					if (pmc == 1) { 

						rdnum=0
						#record stuff
						record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,+1))
						#move for external counter
						moves <- rbind(moves,+1)

						#write out
						write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")
						#print out move for updating
						write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

					} else {

						#make random roll
						rdnum <- runif(1, min=0, max=1)
				
							if (pmc >= rdnum) {

								record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,+1))
								#write out
								write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

								#move for external counter
								moves <- rbind(moves,+1)
								#print out move for updating
								write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

							} else {

								record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,0))
								#write out
								write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")			

								#move for external counter
								moves <- rbind(moves,0)
								#print out move for updating
								write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")						

							}

						}
			} else {

				#grab input
				t1 <- lookup[start,4]
				#make sure there's something there
					if (isEmpty(lookup[lookup[,4]==t1,1])) {
						print("Missing mean energy value3")
						quit()
					} else {
						#avg potential energies for temp
						m1 <- lookup[lookup[,4]==t1,1]
						#pull last known energy- col 2
						e1 <- lookup[lookup[,4]==t1,2]
					}
				#grab next set from list
				t2 <- 0
				e2 <- 0
				#grab last known weight for previous
				w1 <- lookup[lookup[,4]==t1,3]
				w2 <- 0			

				#important- this will not update weights 

				rdnum <- 1
				record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,0,rdnum,start,0))

				#write out
				write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

				#move for external counter
				moves <- rbind(moves,0)
				#print out move for updating
				write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

			}

	#at end can only go down
	} else if (start == tlist) {

		coin <- rbinom(n=1, size=1, prob=0.5)

			if (coin == 1) {

			#order does matter

			#grab input
			t1 <- lookup[start,4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t1,1])) {
					print("Missing mean energy value5")
					quit()
				} else {
					#avg potential energies for temp
					m1 <- lookup[lookup[,4]==t1,1]
					#pull last known energy- col 2
					e1 <- lookup[lookup[,4]==t1,2]
				}
			#grab next set from list
			t2 <- lookup[(start-1),4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t2,1])) {
					print("Missing mean energy value6")
					quit()
				} else {
					#avg potential energies for temp- col 1
					e2 <- lookup[lookup[,4]==t2,1]
				}			
			#grab last known weight for T1- col 3
			w1 <- lookup[lookup[,4]==t1,3]

			reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
			calcw <- rbind(calcw,cbind(w2,t2))
			write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				#check if pass possible
				if (pmc == 1) { 

					rdnum=0
					#record stuff
					record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,+start,-1))
					#move down for external counter
					moves <- rbind(moves,-1)

					#write out
					write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")
					#print out move for updating
					write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

				} else {

					#make random roll
					rdnum <- runif(1, min=0, max=1)
					
						if (pmc >= rdnum) {

							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,-1))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

							#move down for external counter
							moves <- rbind(moves,-1)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

						} else {

							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,0))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")			

							#move for external counter
							moves <- rbind(moves,0)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")						

						}

					}	

			} else {

				#grab input- temp is in col 4 
				t1 <- lookup[start,4]
				#make sure there's something there
					if (isEmpty(lookup[lookup[,4]==t1,1])) {
						print("Missing mean energy value7")
						quit()
					} else {
						#avg potential energies for temp
						m1 <- lookup[lookup[,4]==t1,1]
						#pull last known energy- col 2
						e1 <- lookup[lookup[,4]==t1,2]
					}
				#grab next set from list
				t2 <- 0
				e2 <- 0			
				#grab last known weight for T1- col 3
				w1 <- lookup[lookup[,4]==t1,3]
				w2 <- 0			

				#important- this will not update weights 

				rdnum <- 1
				record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,0,rdnum,start,0))
				#write out
				write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

				#move for external counter
				moves <- rbind(moves,0)
				#print out move for updating
				write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

			}

	} else {

	#flip coin to make down or up
	coin <- rbinom(n=1, size=1, prob=0.5)

		if (coin == 1) {

			#attempt move up

			#grab input- temp is in col 4 
			t1 <- lookup[start,4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t1,1])) {
					print("Missing mean energy value9")
					quit()
				} else {
					#avg potential energies for temp
					m1 <- lookup[lookup[,4]==t1,1]
					#pull last known energy- col 2
					e1 <- lookup[lookup[,4]==t1,2]
				}
			#grab next set from list
			t2 <- lookup[(start+1),4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t2,1])) {
					print("Missing energy value10")
					quit()
				} else {
					#avg potential energies for temp- col 1
					e2 <- lookup[lookup[,4]==t2,1]
				}			
			#grab last known weight for previous
			w1 <- lookup[lookup[,4]==t1,3]

			reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
			calcw <- rbind(calcw,cbind(w2,t2))
			write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				#check if pass possible
				if (pmc == 1) { 

					rdnum=0
					#move for external counter
					moves <- rbind(moves,+1)

					record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,+1))
					#write out
					write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")
					
					#print out move for updating
					write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

				} else {

					#make random roll
					rdnum <- runif(1, min=0, max=1)
					
						if (pmc >= rdnum) {

							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,+1))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

							#move for external counter
							moves <- rbind(moves,+1)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

						} else {
			
							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,0))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

							#move for external counter
							moves <- rbind(moves,0)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")						

						}

					}
			
		} else {

			#attempt move down

			#order does matter

			#grab input- temp is in col 4 
			t1 <- lookup[start,4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t1,1])) {
					print("Missing mean energy value11")
					quit()
				} else {
					#avg potential energies for temp
					m1 <- lookup[lookup[,4]==t1,1]
					#pull last known energy- col 2
					e1 <- lookup[lookup[,4]==t1,2]
				}
			#grab next set from list
			t2 <- lookup[(start-1),4]
			#make sure there's something there
				if (isEmpty(lookup[lookup[,4]==t2,1])) {
					print("Missing energy value12")
					quit()
				} else {
					#avg potential energies for temp- col 1
					e2 <- lookup[lookup[,4]==t2,1]
				}			
			#grab last known weight for T1- col 3
			w1 <- lookup[lookup[,4]==t1,3]

			reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
			calcw <- rbind(calcw,cbind(w2,t2))
			write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				#check if pass possible
				if (pmc == 1) { 

					rdnum=0
					#move down for external counter
					moves <- rbind(moves,-1)

					record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,-1))
					#write out
					write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

					#print out move for updating
					write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

				} else {

					#make random roll
					rdnum <- runif(1, min=0, max=1)
					
						if (pmc >= rdnum) {

							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,-1))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")

							#move down for external counter
							moves <- rbind(moves,-1)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

						} else {

							record <- rbind(record,cbind(t1,m1*4.184,t2,e2*4.184,e1*4.184,w1,w2,pmc,rdnum,start,0))
							#write out
							write.table(record,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="record")
			
							#move for external counter
							moves <- rbind(moves,0)
							#print out move for updating
							write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")						

						}

					}		

		}

	}
EOF

	#attempt to update weight part out of R for speed
	#if (( ${begin} % 10 == 0 )); then

		#combine stored weight with initial just in case temp isn't reached yet
		cat initw swei | tac | awk '!seen[$2]++' | sort -k2 > temp
		cp initw temp_w
		cat temp_w temp > initw

		#echo "errortrack2"

	#fi

	wait

	#echo "errortrack3"
	
	Rscript R_MC.R

	wait

	move=`awk '{print $1}' moves`

	if [[ ${move} -ne 0 ]]; then

		new=`python -c "print (${start}+${move})"`
		#first temp
		temp0=`awk -v b=${start} 'NR==b {print $1}' temp_list`
		#new temp
		tempi=`awk -v c=${new} 'NR==c {print $1}' temp_list`
		#sqrt(Tnew/Told)
		scale=`python -c "print ((float(${tempi})/${temp0})**(0.5))"`

cat >rescale.ptraj<<EOF
parm *.prmtop
trajin 05_prod_${begin}.rst
setvelocity scale factor ${scale}
trajout test.rst
EOF

	/home/sisig/bin/cpptraj-master/bin/cpptraj rescale.ptraj > ptraj_log

	wait

		rm 05_prod_${begin}.rst
		mv test.rst 05_prod_${begin}.rst

	fi

	start=`python -c "print (${start}+${move})"`

	#echo $move $temp0 $tempi $start >> check_log

	if ! (( ${begin} % 25 )) ; then

		#extra cleanup for now every N steps
		rm *.nc *.in *.mdinfo 

	fi

	let begin=begin+1
	let init=init+1

done

#avgpote=`grep "EPtot" 05_prod${iter}.out | awk '{print $9}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`

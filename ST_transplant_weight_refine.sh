#!/bin/bash
export AMBERHOME=/home/sisig/bin/amber16
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AMBERHOME/lib
source /home/sisig/bin/amber16/amber.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-8.0/lib64
export CUDA_VISIBLE_DEVICES=0

rm 05_prod*.in 05_prod*.out *.nc 05_prod*.mdinfo moves R_MC.R

#counter for # ns
final=15000
begin=1
init=0

#counter for position to track temp
start=1
#counter for move up or down or stationary
#move=0

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
  ntx = 5,
  ntb = 1,
  ntxo = 2,
  cut = 9.0,
  ntc = 2, 
  ntf = 2,
  tempi = ${tempi}, 
  temp0 = ${temp0},
  ntt = 3,
  ioutfm = 1,
  gamma_ln = 1.0,
  ig = -1,
  iwrap = 1,
  nstlim = 2500, dt = 0.002,
  ntpr = 500, ntwx = 500, ntwr = 25000,
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

pelist=NULL
pelist <- as.matrix(read.table("pote_list"))
weight=NULL
weight <- as.matrix(read.table("initw"))
tlist=NULL
tlist <- as.matrix(read.table("temp_list"))
start=${start}
moves=NULL
calcw=NULL
record=NULL

	if (start == 1) {

		coin <- 1

			if (coin == 1) {
				#attempt move up
				#order does matter

				#grab input
				t1 <- tlist[start,1]
				#make sure there's something there
					if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
						print("Missing energy value1")
						quit()
					} else {
					#avg potential energies for temp
						m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
						e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
					}
				#grab next set from list
				t2 <- tlist[(start+1),1]
				#make sure there's something there
					if (isEmpty(pelist[pelist[,2]==tlist[start+1,1],1])) {
						print("Missing energy value2")
						quit()
					} else {
					#avg potential energies for temp
						e2 <- mean(pelist[pelist[,2]==tlist[start+1,1],1])
					}			
				#grab last known weight for previous
				w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)

				#just update weight
				reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
				calcw <- rbind(calcw,cbind(w2,t2))
				write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

					if (${begin}%%30 == 0) {

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

						#this does not do anything except update weights
						#this will not generate a record entry and should not
						moves <- rbind(moves,0)
						#print out move for updating
						write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

					}

			} else {

				#grab input
				t1 <- tlist[start,1]
				#make sure there's something there
					if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
						print("Missing energy value3")
						quit()
					} else {
					#avg potential energies for temp
						m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
						e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
					}
				#grab next set from list
				t2 <- 0
				e2 <- 0
				#grab last known weight for previous
				w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)
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
	} else if (start == nrow(tlist)) {

		coin <- 1

			if (coin == 1) {

			#order does matter

			#grab input
			t1 <- tlist[start,1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
					print("Missing energy value5")
					quit()
				} else {
				#avg potential energies for temp
					m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
					e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
				}
			#grab next set from list
			t2 <- tlist[(start-1),1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start-1,1],1])) {
					print("Missing energy value6")
					quit()
				} else {
				#avg potential energies for temp
					e2 <- mean(pelist[pelist[,2]==tlist[start-1,1],1])
				}			
			#grab last known weight for previous
			w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)

			#just update weight
			reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
			calcw <- rbind(calcw,cbind(w2,t2))
			write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				if (${begin}%%30 == 0) {

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

						#move for external counter
						moves <- rbind(moves,0)
						#print out move for updating
						write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")
					}

			} else {

				#grab input
				t1 <- tlist[start,1]
				#make sure there's something there
					if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
						print("Missing energy value7")
						quit()
					} else {
					#avg potential energies for temp
						m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
						e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
					}
				#grab next set from list
				t2 <- 0
				e2 <- 0			
				#grab last known weight for previous
				w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)
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

			#grab input
			t1 <- tlist[start,1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
					print("Missing energy value9")
					quit()
				} else {
				#avg potential energies for temp
					m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
					e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
				}
			#grab next set from list
			t2 <- tlist[(start+1),1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start+1,1],1])) {
					print("Missing energy value10")
					quit()
				} else {
				#avg potential energies for temp
					e2 <- mean(pelist[pelist[,2]==tlist[start+1,1],1])
				}			
			#grab last known weight for previous
			w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)

				reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
				calcw <- rbind(calcw,cbind(w2,t2))
				write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				if (${begin}%%30 == 0) {

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

						#move for external counter
						moves <- rbind(moves,0)
						#print out move for updating
						write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")	

					}
			
		} else {

			#attempt move down

			#order does matter

			#grab input
			t1 <- tlist[start,1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start,1],1])) {
					print("Missing energy value11")
					quit()
				} else {
				#avg potential energies for temp
					m1 <- mean(pelist[pelist[,2]==tlist[start,1],1])
					e1 <- tail(pelist[pelist[,2]==tlist[start,1],1],n=1)
				}
			#grab next set from list
			t2 <- tlist[(start-1),1]
			#make sure there's something there
				if (isEmpty(pelist[pelist[,2]==tlist[start-1,1],1])) {
					print("Missing energy value12")
					quit()
				} else {
				#avg potential energies for temp
					e2 <- mean(pelist[pelist[,2]==tlist[start-1,1],1])
				}			
			#grab last known weight for previous
			w1 <- tail(weight[weight[,2]==tlist[start,1],1],n=1)

			reweight(t1,m1*4.184,t2,e2*4.184,w1,e1*4.184)
			calcw <- rbind(calcw,cbind(w2,t2))
			write.table(calcw,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,file="swei")

				if (${begin}%%30 == 0) {

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

				} else {

					#move for external counter
					moves <- rbind(moves,0)
					#print out move for updating
					write.table(moves,row.names=FALSE,col.names=FALSE,quote=FALSE,file="moves")

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

	#fi

	wait
	
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

	let begin=begin+1
	let init=init+1

done

#avgpote=`grep "EPtot" 05_prod${iter}.out | awk '{print $9}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`

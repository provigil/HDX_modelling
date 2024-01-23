#!/bin/bash

    #for determining if the difference between two sets of HDX data is more than the noise/uncertainty
    #based on Reliable identification of significant differences in differential HX-MS measurements using a hybrid significance testing approach [doi/10.1021/acs.analchem.9b01325]

    #the input HDX data should be a CSV with the structure of 
    #TimePoint  PF(HDX) PF(+/-StdDev)   RepCount    PF(HDX) PF(+/-StdDev)   RepCount



read -e -p "Avaliable nput data: 
    1. HDX
    2. Cross linking
    3. AlphaFold (log.txt)
    4. 
" var1   #prompt for input file
read -p "Confidence Interval (ie, 0.05): " pval #the threshold to use in the T test

cat >hdx_sigcalc.R<<EOF    #make R script to then do the math

hdx_confint <- function(hdx_data, i) {

    if (!is.na(as.numeric(hdx_data[i, 1]))) {   #an ugly check, but it works
            #pull out data
        sd1 <- as.numeric(hdx_data[i, 3])
        n1 <- as.numeric(hdx_data[i, 4])
        sd2 <- as.numeric(hdx_data[i, 6])
        n2 <- as.numeric(hdx_data[i, 7])

        hdx1 <- as.numeric(hdx_data[i, 2])
        hdx2 <- as.numeric(hdx_data[i, 5])

        poolsd <- sqrt((((n1-1)*sd1)+(sd2*(n2-1)))/(n1+n2-2))   #get pooled SD
        sem <- sqrt(((poolsd**2)/n1)+((poolsd**2)/n2))  #get std mean error

        tval <- abs(qt(p=${pval}, df=(n1+n2-2), lower.tail=TRUE))   #get T value; assumes one sided

        hdx_ci <- (tval*sem)    #find conf interval

        dhdx <- abs(as.numeric(abs(hdx1)-abs(hdx2)))  #absolute to keep things simple

        if (dhdx >= hdx_ci) {
            confmat <<- rbind(confmat, cbind(1, hdx_ci)) #when above CI, print 1
        } else if (dhdx < hdx_ci) {
            #when not sig, print 0
            confmat <<- rbind(confmat, cbind(0, hdx_ci))
        }
    } else {
        confmat <<- rbind(confmat, 'NA') #fills in for text spaces
    }

}

inhdx <- read.csv("${var1}")

confmat=NULL    #init blank matrix for significance recording
confmat <- matrix(, nrow = 1, ncol = 2)

for(i in 1:nrow(inhdx)) {
    hdx_confint(inhdx, i)
}

write.table(confmat, file ="sigtest", quote=FALSE, append=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

EOF

wait

Rscript hdx_sigcalc.R
wait
rm hdx_sigcalc.R
firstrun=18
lastrun=43
totalres=269
steps=25
sysname="run2"
analysis="hbond"
dir="data/"
ncores = detectCores()
firstres=21
lastres=40

#####   bin size for histogram
# binsizer <- 1
# binsizez <- 10

#####   output file
outname=paste("corrected",sysname,"_",analysis, "_",firstrun, "to", lastrun,"_",firstres,"to",lastres,".dat",sep="")


#####   hbond and angle and length cut off
hbl<-4
hba<-50 #in degree
hba<- hba*3.14/180 #in radial

#####   donor and acceptor list for each residue
reslist=c("CSH", "CSSC")
minCSSC <- 1
maxCSSC <- 145
acceptorlistCSSC<-c("N1","N2","N3","N4")
acceptorlistCSH<-c("N1","N2")
donorlistCSSC<-c("O1","O2","O3","O4")
donorlistCSH<-c("O1","O2")

#####   reference axis for angle calculation
vref<-c(0,0,1)
vref<-as.vector(vref)
vrefl<- sqrt(sum(vref^2))

library(plyr)
library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
registerDoParallel(cores = ncores)

#####   read pdb and dcd files
pdbname=paste(dir,sysname,".pdb",sep="")
mypdb<-read.pdb(pdbname)
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
print(ncores)
#####	for testing...
#nframes<-1
#totalres=40

#####   main function for getting hbond
gethbond <- function(firstresn,lastresn,j) {
  #####   get coordinate for each atom of interest
  cooracc<-as.vector(mydcd[j,selacc$xyz])
  coorh<-as.vector(mydcd[j,selh$xyz])
  coordor<-as.vector(mydcd[j,seldor$xyz])
  #####   get distance between acceptor and donor
  vhb <- coordor-cooracc
  #####   get norm
  vhbl <- sqrt(sum(vhb^2))
  #####   hbond length cutoff
  if (vhbl <= hbl) {
    #####   get vector between hydrogen and donor
    va <- coorh - cooracc
    #####   get norm
    val <- sqrt(sum(va^2))
    #####   get angle
	vhbacos<-sum(va * vhb)/(val * vhbl)
	vhba <- acos(vhbacos)
    #####   hbond angle cutoff
    if (vhba <= hba) {
      #####   get center coord of hbond
      coorhb <- (coorh+coordor)/2
      #####   get angle of hbond with reference axis
      vhbcos <- (sum(vhb * vref)/(vhbl * vrefl))
      #####   get moiety indicator
      if (acceptor=="N1" || acceptor=="N4") {
        if (donor=="O1" || donor=="O4") {
          indicator<-"AM1-AM1" 
        } else {
          indicator<-"AM1-AM2"
        }
      } else {
        if (donor=="O1" || donor=="O4") {
          indicator<-"AM1-AM2"
        } else {
          indicator<-"AM2-AM2"
        }
      }
      #####   get residue indicator
      if (k>=minCSSC && k<=maxCSSC) {
        if (i>=minCSSC && i<=maxCSSC) {
          mindicator<-"CSSC"
        } else {
          mindicator<-"CSSC-CSH"
        }
      } else {
        if (i>=minCSSC && i<=maxCSSC) {
          mindicator<-"CSSC-CSH"
        } else {
          mindicator<-"CSH"
        }
      }
      #####   add to output matrix
      c(coorhb,vhbacos,vhbcos,indicator,mindicator,vhbl)
   }
  }
}

    

#####   generate empty output data frame
out<-c()
for (i in firstres:lastres) {
  #####   check which residue we are dealing with
  if (i>=minCSSC && i<=maxCSSC) {
    acceptorlist<-acceptorlistCSSC
  } else {
    acceptorlist<-acceptorlistCSH
  }
  #####   loop through h-bond acceptors 
  for (acceptor in acceptorlist) {
    ######   get associated hydrogen list
    if (acceptor=="N1") {
      hydrogenlist<-c("H8","H9")
    } else if (acceptor=="N2") {
      hydrogenlist<-c("H7")
    } else if (acceptor=="N3") {
      hydrogenlist<-c("H16")
    } else {
      hydrogenlist<-c("H21","H22")
    }
    #####   select acceptor
    selacc<-atom.select(mypdb,elety = acceptor,resno =i)
    #####   loop through h-bond hydrogen
    for (hydrogen in hydrogenlist) {
      #####   select hydrogen
      selh<-atom.select(mypdb,elety = hydrogen,resno =i)
      for (k in 1:totalres) {
        #####   avoid same residue
        if (i==k) {
          next
        }
        #####   check donor residue
        if (k>=minCSSC && k<=maxCSSC) {
          donorlist<-donorlistCSSC
        } else {
          donorlist<-donorlistCSH
        }
        #####   loop through donor
        for (donor in donorlist) {
          #####   select donor
          seldor<-atom.select(mypdb,elety = donor,resno =k)
          #####   loop through frames
          out<-rbind(out,foreach (j=1:nframes,.combine=rbind) %dopar% {
                       gethbond(1,totalres,j)
                     })
 
        }
      }
    }
  }
  write.table(out, append=TRUE,file=outname,col.names = F,row.names=F,quote = F)
  out<-c()
  old<-Sys.time()
  new<-Sys.time()-old
  print(i)
  print(new)
}

#####   write table

#####   get r coordinate
#out <- out %>% mutate(r = sqrt(X1^2+X2^2))
#binnumr <- (0-min(particle$r)+max(particle$r))/binsizer - 1
#rdist <- as.data.frame(table(cut(particle$r,breaks =seq(min(particle$r),max(particle$r),binsizer))))
#out <- out[,-c(1,2)]




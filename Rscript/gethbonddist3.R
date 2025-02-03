# if two atoms are further than avoidcutoff, then avoid further computation of thse two residues 
outfreq<-40
avoidcutoff<-40
maxcssc<-411
totalres<-763
firstres=1
lastres=totalres-1
step<-10
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
#dcdname<-"test.dcd"
pdbname<-"data/mix_27A.pdb"
test<-FALSE
getdist<-FALSE
totpair<-4*4*(firstres+lastres)*lastres/2
#####   output file
#####   bin size for histogram
# binsizer <- 1
# binsizez <- 10
#outname=paste(sysname,"_",analysis, "_",firstrun, "to", lastrun,"_",firstres,"to",lastres,".dat",sep="")
outname<-"mix_27A_hbond_2D.rda"
#####   hbond and angle and length cut off
hbl<-3.2
hba<-30 #in degree
hba<- hba*3.14/180 #in radial
#####   donor and acceptor list for each residue
donorlistCSSC<-c("N1","N2","N3","N4")
donorlistCSH<-c("N1","N2")
acceptorlistCSSC<-c("O1","O2","O3","O4")
acceptorlistCSH<-c("O1","O2")

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
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
registerDoParallel(cores = detectCores())
print(detectCores())
#####   read pdb and dcd files
mypdb<-read.pdb(pdbname)
mydcd<-read.dcd(dcdname,big=TRUE)
mycell<-read.dcd(dcdname,cell=T)
nframes<-nrow(mydcd)

#####	for testing...
if(test==TRUE) {
  nframes<-2
  totalres=40
  lastres=totalres
}
# calculate the Euclidean distance between two points
euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))  # Calculate Euclidean distance
}
# implement periodic boundary conditions
pbc <- function(comi, cur, currbox) {
  # Loop through x, y, and z coordinates
  for (l in 1:3) {
    diff <- cur[l] - comi[l]  # Calculate the difference between the current coordinate and the center of mass coordinate
    absdiff <- abs(diff)  # Calculate the absolute value of the difference
    # Update coordinate only when necessary to account for periodic boundary conditions
    if (absdiff > currbox[l]/2) {
      cur[l] <- cur[l] - (diff/abs(diff)) * currbox[l] * ceiling((abs(diff) - currbox[l]/2) / currbox[l])
    } else {
      cur[l] <- cur[l]
    }
  }
  
  # Return updated current coordinates
  cur
}
#####   main function for getting hbond
gethbond <- function(j) {
  rawdistm<-big.matrix(totpair,3 ) #set up a matrix for cuthis function 
  counter<-1 #counter for which line in matrix to update
  avoidlist <- list()
  for (i in firstres:totalres) {
    avoidlist[[i]]<-c(0)
  }
  for (i in firstres:lastres) {
    #####   check which residue we are dealing with
    if ( i<=maxcssc) {
      donorlist<-donorlistCSSC
    } else {
      donorlist<-donorlistCSH
    }
    for (k in (i+1):totalres) {
      #####   avoid same residue
      if (i==k || k%in%avoidlist[[i]]==TRUE || i%in%avoidlist[[k]]==TRUE) {
        next
      }
      #####   check acceptor residue
      if ( k<=maxcssc) {
        acceptorlist<-acceptorlistCSSC
      } else {
        acceptorlist<-acceptorlistCSH
      }
      #####   loop through h-bond donors 
      for (donor in donorlist) {
        #####   avoid same residue
        if (k%in%avoidlist[[i]]==TRUE || i%in%avoidlist[[k]]==TRUE) {
          break
        }
        ######   get associated hydrogen list
        if (donor=="N1") {
          hydrogenlist<-c("H8","H9")
        } else if (donor=="N2") {
          hydrogenlist<-c("H7")
        } else if (donor=="N3") {
          hydrogenlist<-c("H16")
        } else {
          hydrogenlist<-c("H21","H22")
        }
        #####   select donor
        seldor<-atom.select(mypdb,elety = donor,resno =i)
        coordor<-as.vector(mydcd[j,seldor$xyz])
        for (acceptor in acceptorlist) {
          #####   select acceptor
          selacc<-atom.select(mypdb,elety = acceptor,resno =k)
          cooracc<-as.vector(mydcd[j,selacc$xyz])
          cooracc<-pbc(coordor,cooracc,mycell[j,])
          #####   get distance between acceptor and donor
          vhb <- cooracc-coordor
          #####   get norm, acceptor donor distance
          vhbl <- sqrt(sum(vhb^2))
          if (vhbl>avoidcutoff) {
            avoidlist[[i]]<-append(avoidlist[[i]],k)
            avoidlist[[k]]<-append(avoidlist[[k]],i)
            break
          }
          # Set an impossible OH distance
          OHdis<-mycell[j,1]*2  
          #####   loop through h-bond hydrogen
          for (hydrogen in hydrogenlist) {
            #####   select hydrogen
            selh<-atom.select(mypdb,elety = hydrogen,resno =i)
            coorh<-as.vector(mydcd[j,selh$xyz])
            tepOHdis<-euclidean(cooracc,coorh) # Get OH distance 
            if (tepOHdis<OHdis) {
              OHdis<-tepOHdis # Compare to set minimal OH distance
              coorhclose<-coorh # Save h coord for later
            }
          }
          #####   get vector between hydrogen and donor
          va <- coorhclose - coordor
          #####   get norm
          val <- sqrt(sum(va^2))
          #####   get angle
          vhbacos<-sum(va * vhb)/(val * vhbl)
          vhba <- acos(vhbacos)
          #####   hbond angle cutoff
          if (getdist==TRUE) {
            #####   get center coord of hbond
            coorhb <- (coorh+cooracc)/2
            #####   get angle of hbond with reference axis
            vhbcos <- (sum(vhb * vref)/(vhbl * vrefl))
            #####   get moiety indicator
            if (donor=="N1" || donor=="N4") {
              if (acceptor=="O1" || acceptor=="O4" || acceptor=="N1" || acceptor=="N4") {
                indicator<-"AM1-AM1" 
              } else {
                indicator<-"AM1-AM2"
              }
            } else {
              if (acceptor=="O1" || acceptor=="O4" || acceptor=="N1" || acceptor=="N4") {
                indicator<-"AM2-AM1"
              } else {
                indicator<-"AM2-AM2"
              }
            }
            #####   get residue indicator
            if ( k<=maxcssc) {
              if ( i<=maxcssc) {
                mindicator<-"CSSC-CSSC"
              } else {
                mindicator<-"CSSC-CSH"
              }
            } else {
              if ( i<=maxcssc) {
                mindicator<-"CSH-CSSC"
              } else {
                mindicator<-"CSH-CSH"
              }
            }
            #####   add to output matrix
            outj<-rbind(outj,c(coorhb,vhbacos,vhbcos,indicator,mindicator,vhbl))
          } else {
            rawdistm[counter,]<-c(vhbl, vhbacos, vhba) # Update matrix (distance, angle)
            counter<-counter+1 
          }
        }
      }
    }
  }
  rawdistm<-as.matrix(rawdistm[rowSums(!is.na(rawdistm[,])) > 0, ])
  cuthis(rawdistm,dishisbreak,anglehisbreak)
}
# Set the distance bin size max and min for distance
disbinsize <- 0.1
disbinmax <- 40# sqrt(2 * 100^2) / 2
disbinmin <- 0
# Set the distance bin size max and min for distance
anglebinsize <- 0.05
anglebinmax <- 1
anglebinmin <- -1
# Create an array of midpoint values for each distance 
dishistx <- seq(disbinmin + disbinsize / 2, disbinmax - disbinsize / 2, by = disbinsize)
# Create an array of distance bin 
dishisbreak <- c(dishistx - disbinsize / 2, disbinmax)
# Create an array of midpoint values 
anglehistx <- seq(anglebinmin + anglebinsize / 2, anglebinmax - anglebinsize / 2, by = anglebinsize)
# Create an array of angle bin 
anglehisbreak<-c(anglehistx-anglebinsize/2,anglebinmax)
# create 2D histogram table
cuthis <- function(coord, breakseq1, breakseq2) {
  # Use the cut() function to bin each variable
  var1_binned <- cut(coord[,1], breaks=breakseq1, include.lowest=TRUE)
  var2_binned <- cut(coord[,2], breaks=breakseq2, include.lowest=TRUE)
  
  # Use the table() function to create the 2D histogram table
  hist_table <- table(var1_binned, var2_binned)
  
  # Print the table
  hist_table
}


out<-c()
timelist<-c()
#####   generate empty output data frame
#####   loop through frames
l<-1
while (TRUE) {
  old<-Sys.time()
  m<-l+outfreq-1
  print(old)
  if (l<=nframes & m<=nframes) {
    m<-m
  } else if (l<=nframes) {
    m<-nframes
  } else {
    break
  }
  outtemp<-foreach (j=l:m) %dopar% {
    gethbond(j)
  }
  out<-append(out,outtemp)
  new<-Sys.time()-old
  l<-l+outfreq
  print(m/nframes)
  print(new)
  timelist<-c(timelist,as.numeric(gsub("([0.-9.]+).*$", "\\1", new)))
  print(paste("Time remain", mean(timelist)*(nframes-m)/outfreq, "progress", m/nframes,"frame", m))
  save(out,file=outname)
}



#####   write table

#####   get r coordinate
#out <- out %>% mutate(r = sqrt(X1^2+X2^2))
#binnumr <- (0-min(particle$r)+max(particle$r))/binsizer - 1
#rdist <- as.data.frame(table(cut(particle$r,breaks =seq(min(particle$r),max(particle$r),binsizer))))
#out <- out[,-c(1,2)]




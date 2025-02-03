maxcssc=411
dihetype<-"S"
library(bio3d)
library(sna)
library(doParallel)
library(dplyr)
corenum=24
corenum=detectCores()
registerDoParallel(cores = corenum)
nottip=FALSE
### for pdb & traj file names
sysname="mix_27A"
analysis="density"
dir="data/"
steps=1
firstrun=1
lastrun=9
### first frame to be considered
startframe<-2500
### for pdb & traj file names
pdbname=paste0("data/CSH_CSSC_mix_27A_cyl_water_ions_redo.pdb")
#pdbname=paste(dir,sysname,"_resonly.pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_all_every_", steps, ".dcd",sep="")
### how many residues in total
namelistphe1<-c("C6","C5","C4","N2")
namelistphe2<-c("C19","C16","C17","N3")
namelistS1<-c("S1","C3","C2","N2")
namelistS2<-c("S2","C11","C15","N3")
binsize<-1
binmax<-180
binmin<--180
nbins <- ceiling((binmax - binmin) / binsize)
histx<-seq(binmin + binsize / 2, binmax - binsize / 2, by = binsize)
getdihedS <- function(j,i) {
  myresids<-resids[[j]][[i]]
  myresids<-myresids[which(myresids<maxcssc & myresids%%2!=0)]
  if (length(myresids)>0) {
    xyzlist<-c()
    for (k in 1:length(myresids)) {
      tempres<-myresids[k]
      tempreslist<-c(tempres,tempres,tempres+1,tempres+1)
      for (l in 1:4) {
        sel<-atom.select(mypdb,resno = tempreslist[l],elety=namelist[l],operator="AND")
        xyzlist<-c(xyzlist,mydcd[j,sel$xyz])
      }
    }
    dihedtemp<-torsion.xyz(xyzlist, atm.inc = 4)
    dihedtemp<-as.data.frame(table(cut(dihedtemp,breaks =seq(-180,180,dihedbinsize))))
    dihedtemp<-dihedtemp %>% mutate(dihed=seq(-180,180-dihedbinsize,dihedbinsize))
    dihedtemp[,2]
  } else {
    c()
  }
}
getdihedCSSC <- function(j,namelist1,namelist2) {
  myresids<-resids[[j+totalframe]]
  if (nottip==TRUE) {
    myresids<-setdiff(seq(1, 763, 1), resids[[j+totalframe]])
  }
  myresids<-myresids[which(myresids<=maxcssc)]
  if (length(myresids)>0) {
    xyzlist<-c()
    for (k in 1:length(myresids)) {
      tempres<-myresids[k]
      for (l in 1:4) {
        sel<-atom.select(mypdb,resno = tempres,elety=namelist1[l],operator="AND")
        xyzlist<-c(xyzlist,mydcd[j,sel$xyz])
      }
      for (l in 1:4) {
        sel<-atom.select(mypdb,resno = tempres,elety=namelist2[l],operator="AND")
        xyzlist<-c(xyzlist,mydcd[j,sel$xyz])
      }
    }
    dihedtemp<-torsion.xyz(xyzlist, atm.inc = 4)
    dihedtemp<-as.data.frame(table(cut(dihedtemp,breaks =histx)))
    dihedtemp[,2]
    } else {
    c()
  }
}
getdihedCSH <- function(j,namelist1) {
  myresids<-resids[[j+totalframe]]
  if (nottip==TRUE) {
    myresids<-setdiff(seq(1, 763, 1), resids[[j+totalframe]])
  }
  myresids<-myresids[which(myresids>maxcssc)]
  if (length(myresids)>0) {
    xyzlist<-c()
    for (k in 1:length(myresids)) {
      tempres<-myresids[k]
      for (l in 1:4) {
        sel<-atom.select(mypdb,resno = tempres,elety=namelist1[l],operator="AND")
        xyzlist<-c(xyzlist,mydcd[j,sel$xyz])
        
      }
    }
    dihedtemp<-torsion.xyz(xyzlist, atm.inc = 4)
    dihedtemp<-as.data.frame(table(cut(dihedtemp,breaks =histx)))
    dihedtemp[,2]
    } else {
    c()
  }
}
update_histogram <- function(histogram, dis) {
  if (dis < binmin) {
    histogram[1, 2] <- histogram[1, 2] + 1
  } else if (dis > binmax) {
    histogram[nrow(histogram), 2] <- histogram[nrow(histogram), 2] + 1
  } else {
    bin <- floor((dis - binmin) / binsize) + 1
    histogram[bin, 2] <- histogram[bin, 2] + 1
  }
  histogram
}
# read the file
# read the file
lines <- readLines("tips_id.dat")
# create a list
resids <- lapply(lines, function(line) {
  # split each line into a vector of numbers
  as.numeric(unlist(strsplit(line, " ")))
})
mypdb<-read.pdb(pdbname,hex=TRUE)
#mydcd<-read.dcd(dcdname)
#nframes<-nrow(mydcd)
outcsscphe<-0
outcshphe<-0
outcsscS<-0
outcshS<-0
totalframe<-0

for (runnum in firstrun:lastrun) {
  if (startframe > totalframe) {
    totalframe<-totalframe+1000
    next
  } else {
    totalframe<-totalframe+1000
    restartframe<-startframe %% 1000
    mydcd<-read.dcd(paste0("../nvt0",runnum,".dcd"))
    endframe<-nrow(mydcd)
    outcsscphe1<-foreach (framenum=restartframe:endframe,.combine="+") %dopar% {
      getdihedCSSC(framenum,namelistphe1,namelistphe2) 
    }  
    outcsscphe<-outcsscphe+outcsscphe1
    # outcshphe1<-foreach (framenum=restartframe:endframe,.combine="+") %dopar% {
    #   getdihedCSH(framenum,namelistphe1) 
    # }  
    # outcshphe<-outcshphe+outcshphe1
    # outcsscS1<-foreach (framenum=restartframe:endframe,.combine="+") %dopar% {
    #   getdihedCSSC(framenum,namelistS1,namelistS2) 
    # }  
    # outcsscS<-outcsscS+outcsscS1
    # outcshS1<-foreach (framenum=restartframe:endframe,.combine="+") %dopar% {
    #   getdihedCSH(framenum,namelistS1) 
    # } 
    # outcshS1<-outcshS+outcshS11
  }
}
outcsscphe<-cbind(outcsscphe,histx)
#outcshphe<-cbind(outcshphe,histx)
#outcshphe <- cbind(outcshphe, "CSH")
#outcsscphe <- cbind(outcsscphe, "CSSC")
# merge the two data frames by histx
mergedphe <- merge(outcshphe, outcssc, by="histx")
# write merged data frame to file
if (nottip==TRUE) {
write.table(mergedphe, "nottip_dihed_phe.dat", sep="\t", row.names=FALSE)
} else {
  write.table(mergedphe, "tip_dihed_phe.dat", sep="\t", row.names=FALSE)
}
 
#outcsscS<-cbind(outcsscS,histx)
#outcshS<-cbind(outcshS,histx)
#outcshS <- cbind(outcshS, "CSH")
#outcsscS <- cbind(outcsscS, "CSSC")
# merge the two data frames by histx
mergedS <- merge(outcshS, outcsscS, by="histx")
# write merged data frame to file
if (nottip==TRUE) {
  write.table(mergedS, "nottip_dihed_CYS.dat", sep="\t", row.names=FALSE)
} else {
  write.table(mergedS, "tip_dihed_CYS.dat", sep="\t", row.names=FALSE)
}
if (0) {
# Read tip_dihed_phe.dat into a data frame
tip_dihed_phe <- read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/tip_dihed_phe.dat", header = TRUE)

# Create separate data frames for CSH and CSSC data
csh_data <- subset(tip_dihed_phe, Identifier == "CSH")
cssc_data <- subset(tip_dihed_phe, Identifier == "CSSC")

# Plot the data using ggplot
ggplot(tip_dihed_phe, aes(x = histx, y = outcsh, color = "CSH")) +
  geom_line() +
  geom_line(aes(y = outcssc, color = "CSSC")) +
  scale_color_manual(values = c("blue", "red"), name = "") +
  xlab("Dihedral Angle (degrees)") +
  ylab("Count") +
  ggtitle("Dihedral Angle Distribution Phe sidechain in tips") +
  theme_bw()
}





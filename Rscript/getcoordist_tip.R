
getrcoor=FALSE
centermass=FALSE
sysname="mix_27A"
analysis="density"
dir="data/"
steps=25
firstrun=6
lastrun=10
### for pdb & traj file names
pdbname=paste0("data/CSH_CSSC_mix_27A_cyl_water_ions_redo.pdb")
#pdbname=paste(dir,sysname,"_resonly.pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_all_every_", steps, ".dcd",sep="")
### first frame to be considered
startframe<-1
moietyname <- c("N1", "C1", "O1", "C2", "C3", "S", "N2", "C4", "O2", "C5", "C6", "C7", "C8", "C9", "C10","S1")
moietyname2 <- c("S2", "C11", "C12", "C13", "C14", "N3", "C15", "C16", "C17", "C18", "C19", "C20", "O3", "O4", "N4")
phemoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
am1moietyname <- c("O2", "C4", "N2")
am2moietyname <- c("O1", "C1", "N1")
cysmoietyname <- c("C3", "S1", "C2","S")
phemoietyname2 <- c("C12", "C13", "C14", "C16", "C18", "C19")
am1moietyname2 <- c("O3", "C17", "N3")
am2moietyname2 <- c("N4", "C20", "O4")
cysmoietyname2 <- c("C11", "S2", "C15")
totnohmoiety<-15*763+length(moietyname2)*411
totphemoiety<-length(phemoietyname)*763+length(phemoietyname2)*411
### how many residues in total
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
tipminz<--106
tipmaxz<-109
### analysis name for output name
analysis="coordist"
### binsize and max for normalization
if (getrcoor==TRUE) {
  binsize<-0.1
  binmax<-70
  binmin<-0
} else {
  binsize<-5
  binmax<-200
  binmin<--200
}
nbins <- ceiling((binmax - binmin) / binsize)
histx<-seq(binmin + binsize / 2, binmax - binsize / 2, by = binsize)
#####  Volume of sphere shell for rdf normalization
nvshell<-4*3.14/3
box<-97.7713223315*97.7713223315*321.607782715
### residue name for selection
resnamelist<-c("CSSC","CSH")
### output data file name
outname=paste0("data/",sysname, analysis, ".dat")
outname2=paste0("data/",sysname, analysis, "_normalized.dat")
### load all the library
library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)
### register the core
registerDoParallel(cores = detectCores())
### function for updating histogram
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
### function for cutting histogram
cuthis<-function(coord) {
  as.data.frame(table(cut(coord,breaks = c(histx-0.5,binmax) )))[,2]
}
### write tip index to file
write_tips <- function(tips, filename) {
  # open the file for writing
  fileConn <- file(filename, "w")
  for (i in seq_along(tips)) {
    line <- paste(tips[[i]], collapse = " ")
    # loop through each element of the tips list
    cat(line, "\n", file = fileConn)
  }
  # close the file
  close(fileConn)
}
### function for getting zcoordhist for different moiety 
getmoietyzcoorhist<-function(inputmoietyname,resname,framenum,getr=FALSE){
  # Create an empty histogram with zero counts
  temphistogram <- matrix(0, nrow = nbins, ncol = 2)
  # Set the bin centers
  temphistogram[, 1] <- histx
  if (resname=="CSSC") {
    firsti<-mincssc
    lasti<-maxcssc
  } else {
    firsti<-mincsh
    lasti<-maxcsh
  }
  #loop through resids
  if (centermass) {
    for (i in firsti:lasti) {
      #select moiety and resname
      #####   select moiety
      moietysele<-atom.select(mypdb,resno=i,resid=resname,elety=inputmoietyname,operator = "AND")
      #####   get coordinate
      moietyxyz<-mydcd[framenum,moietysele$xyz]
      #####   calculate COM of res i
      com<-as.vector(com.xyz(moietyxyz))-com_t[framenum,]
      #update histogram
      if (getr==FALSE) {
        temphistogram <- update_histogram(temphistogram, com[3])
      } else {
        temphistogram <- update_histogram(temphistogram, sqrt(com[1]^2+com[2]^2))
      }
    }
    #return the count column 
    temphistogram[,2]
  } else {
    moietyxyz<-atom.select(mypdb,resid=resname,elety=inputmoietyname,operator = "AND")
    if (getr==FALSE) {
      moietyz<-moietyxyz[seq(3, length(moietyxyz), 3)]-com_t[framenum,3]
      cuthis(moietyz)[,2]
    } else {
      moietyx<-moietyxyz[seq(1, length(moietyxyz)-2, 3)]-com_t[framenum,1]
      moietyy<-moietyxyz[seq(2, length(moietyxyz)-1, 3)]-com_t[framenum,2]
      
    }
  }
}
#####   read dcd files
mydcd<-read.dcd(dcdname)
mypdb<-read.pdb(pdbname,hex=TRUE)
#mydcd<-read.dcd("/Users/song/Documents/Research/MRSEC/CSH-CSSC/sumtraj/mix_27A_sum1to8_all_every_20.dcd")
#mypdb<-read.pdb("/Users/song/Documents/Research/MRSEC/CSH-CSSC/sumtraj/CSH_CSSC_mix_27A_cyl_water_ions.pdb",hex=TRUE)
#mypdb<-read.pdb("/Users/song/Documents/Research/MRSEC/CSH-CSSC/cyl_hpc/CSH_CSSC_3to7_cyl_water_ions.pdb")
nframes<-nrow(mydcd)
#loop through each res
if (1) {
  resnoh<-atom.select(mypdb,elety=c(moietyname,moietyname2))
  for (framenum in startframe:nframes) {
    if (framenum==startframe) {
      com_t<-matrix(nrow=nframes,ncol = 3)
      counter<-1
      coorlab<-c("x","y","z")
    }
    #####   get coordinate
    resnohxyz<-mydcd[framenum,resnoh$xyz]
    #####   calculate COM of res i
    comi<-as.vector(com.xyz(resnohxyz))
    ##### for tracking com shift
    for (i in 1:3) {
      com_t[framenum,]<-comi
      counter<-counter+1
    }
  }
  outm <- c()
  for (resname in resnamelist) {
    # Phe moiety
    pheoutm <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
      getmoietyzcoorhist(phemoietyname, resname, framenum,getrcoor)
    }
    # Reconstruct output
    pheoutm <- cbind(histx, pheoutm, resname, "PHE")
    # AM1 moiety
    am1outm <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
      getmoietyzcoorhist(am1moietyname, resname, framenum,getrcoor)
    }
    # Reconstruct output
    am1outm <- cbind(histx, am1outm, resname, "AM1")
    # AM2 moiety
    am2outm <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
      getmoietyzcoorhist(am2moietyname, resname, framenum,getrcoor)
    }
    # Reconstruct output
    am2outm <- cbind(histx, am2outm, resname, "AM2")
    # Cys moiety
    cysoutm <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
      getmoietyzcoorhist(cysmoietyname, resname, framenum,getrcoor)
    }
    # Reconstruct output
    cysoutm <- cbind(histx, cysoutm, resname, "CYS")
    #each moiety has a second list for CSSC
    if (resname == "CSSC") {
      # Phe moiety
      pheoutm1 <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
        getmoietyzcoorhist(phemoietyname2, resname, framenum,getrcoor)
      }
      # Reconstruct output
      pheoutm1 <- cbind(histx, as.numeric(pheoutm[, 2]) + pheoutm1, resname, "PHE")
      
      # AM1 moiety
      am1outm1 <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
        getmoietyzcoorhist(am1moietyname2, resname, framenum,getrcoor)
      }
      # Reconstruct output
      am1outm1 <- cbind(histx, as.numeric(am1outm[, 2]) + am1outm1, resname, "AM1")
      
      # AM2 moiety
      am2outm1 <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
        getmoietyzcoorhist(am2moietyname2, resname, framenum,getrcoor)
      }
      # Reconstruct output
      am2outm1 <- cbind(histx, as.numeric(am2outm[, 2]) + am2outm1, resname, "AM2")
      
      # Cys moiety
      cysoutm1 <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
        getmoietyzcoorhist(cysmoietyname2, resname, framenum,getrcoor)
      }
      # Reconstruct output
      cysoutm1 <- cbind(histx, as.numeric(cysoutm[, 2]) + cysoutm1, resname, "CYS")
      
      # Collect each output
      outm <- rbind(outm, pheoutm1, am1outm1, am2outm1, cysoutm1)
      
    } else {
      #collect each output
      # Combine all outputs
      outm <- rbind(outm,pheoutm, am1outm, am2outm, cysoutm)    }
  }
  if (getrcoor==TRUE) {
    colnames(outm)<-c("r","Freq","Resname","Type")
    write.table(outm, file=paste("moiety_rcoordist",firstrun, "to", lastrun,"every", steps,"dat",sep = "_"), row.names=FALSE,quote = F)
  } else {
    colnames(outm)<-c("z","Freq","Resname","Type")
    write.table(outm, file=paste("moiety_zcoordist",firstrun, "to", lastrun,"every", steps,"dat",sep = "_"), row.names=FALSE,quote = F)
  }
}
if (0) {
  resnoh<-atom.select(mypdb,elety=c(moietyname,moietyname2))
  cshnoh<-atom.select(mypdb,resid="CSH",elety=moietyname)
  csscnoh<-atom.select(mypdb,resid="CSSC",elety=c(moietyname,moietyname2))
  watnoh<-atom.select(mypdb,elety="OH2")
  #what are the index for resnoh, use resnoh$atom
  atomind<-resnoh$atom
  for (framenum in startframe:nframes) {
    if (framenum==startframe) {
      com_t<-matrix(nrow=nframes*3,ncol = 3)
      counter<-1
      coorlab<-c("x","y","z")
      if (getrcoor==TRUE) {
        cshrcoorhist<-0
      csscrcoorhist<-0
      watrcoorhist<-0
      } else {
        cshzcoorhist<-0
        cssczcoorhist<-0
        watzcoorhist<-0
      }
      tipind<-c()
      nottipind<-c()
    }
    #####   get coordinate
    resnohxyz<-mydcd[framenum,resnoh$xyz]
    #####   calculate COM of res i
    comi<-as.vector(com.xyz(resnohxyz))
    ##### for tracking com shift
    for (i in 1:3) {
      com_t[counter,]<-c(framenum,comi[i],coorlab[i])
      counter<-counter+1
    }
    # #####   get coordinate
    cshnohxyz<-mydcd[framenum,cshnoh$xyz]
    csscnohxyz<-mydcd[framenum,csscnoh$xyz]
    watresnohxyz<-mydcd[framenum,watnoh$xyz]
    ####   get zcoordinate and correct for COM
    #cshzcoord<-cshnohxyz[seq(3, length(cshnohxyz), 3)]-comi[3]
    #cssczcoord<-csscnohxyz[seq(3, length(csscnohxyz), 3)]-comi[3]
    #watzcoord <- watresnohxyz[seq(3, length(watresnohxyz), 3)]-comi[3]
    #watxcoord <- watresnohxyz[seq(1, length(watresnohxyz)-2, 3)]-comi[1]
    #watycoord <- watresnohxyz[seq(2, length(watresnohxyz)-1, 3)]-comi[2]
    #watrcoord<-sqrt(watxcoord^2+watycoord^2)
    if (getrcoor==TRUE) {
      # CSH  get x y z coordinate and correct accordingly and also r=sqrt(x^2+y^2)
      cshxcoord<-cshnohxyz[seq(1, length(cshnohxyz)-2, 3)]-comi[1]
      cshycoord<-cshnohxyz[seq(2, length(cshnohxyz)-1, 3)]-comi[2]
      cshrcoord<-sqrt(cshxcoord^2+cshycoord^2)[which(cshzcoord>tipminz & cshzcoord<tipmaxz)]
      # CSSC  get x y z coordinate and correct accordingly and also r=sqrt(x^2+y^2)
      csscxcoord<-csscnohxyz[seq(1, length(csscnohxyz)-2, 3)]-comi[1]
      csscycoord<-csscnohxyz[seq(2, length(csscnohxyz)-1, 3)]-comi[2]
      csscrcoord<-sqrt(csscxcoord^2+csscycoord^2)[which(cssczcoord>tipminz & cssczcoord<tipmaxz)]
      #####   get coordinate
      watrcoord<-watrcoord[which(watzcoord>tipminz & watzcoord<tipmaxz )]
      ##### update r histogram
      cshrcoorhist<-cshrcoorhist+cuthis(cshrcoord)
      csscrcoorhist<-csscrcoorhist+cuthis(csscrcoord)
      watrcoorhist<-watrcoorhist+cuthis(watrcoord)
    } else {
      # ##### update z histogram
      # watzcoorhisttemp<-cuthis(watzcoord[which(watrcoord<27.5)])
      # watzcoorhist<-watzcoorhist+watzcoorhisttemp
      # cshzcoorhist<-cshzcoorhist+cuthis(cshzcoord)
      # # ##### update histogram
      # cssczcoorhist<-cssczcoorhist+cuthis(cssczcoord)
    }
    if (0) {
      # #what z value to consider for water
    # watminz<-histx[max(which(resnohzcoorhist<avefreq/10 & histx<0))]
    # watmaxz<-histx[min(which(resnohzcoorhist<avefreq/10 & histx>0))]
    # #now to deal with water...
    # watrcoord<-sqrt(watxcoord^2+watycoord^2)
    #watinfibzcoord<-watzcoord#[which( watrcoord<27.5)]
    # watinfibrcoord<-watrcoord[which(watzcoord<watmaxz & watzcoord>watminz & watrcoord<27.5)]
    ##### update histogram
    #watzcoorhist<-watzcoorhist+cuthis(watinfibrcoord)
    # #### now lets figure out tips 
    #resnohzcoorhist<-cuthis(resnohzcoor)/sum(cuthis(resnohzcoor))
    #what's the average normalized freq resnohzcoorhist
    #avefreq<-(mean(watzcoorhisttemp[which(histx>-90 & histx<90)])+mean(watzcoorhisttemp[which(histx>120 & histx<150)]))/2
    # #wht's the first and last z value that matches this average normalized freq
    #tipminz<-histx[max(which(watzcoorhisttemp<avefreq & histx< -100 & histx > -120))]
    #tipmaxz<-histx[min(which(watzcoorhisttemp<avefreq & histx< 120 & histx > 100))]
    }
    resnohzcoor<-resnohxyz[seq(3, length(resnohxyz), 3)]-comi[3]
    #which resnohzcoor are outside the previous defined first and last z value, and find their corrsponding index in renoh$atom
    tipind[[framenum]]<-atomind[which(resnohzcoor< -85 | resnohzcoor > 85)]
    nottipind[[framenum]]<-atomind[which(resnohzcoor < 20 & resnohzcoor > -20)]
    #either store in a list, or write it down, or both
  }
  write_tips(tipind, paste("tips_index",firstrun, "to", lastrun,"every", steps,".dat",sep = "_"))
  write_tips(nottipind, paste("nottips_index",firstrun, "to", lastrun,"every", steps,".dat",sep = "_"))
  #save(tipind,file=paste("tips_index",firstrun, "to", lastrun,"every", steps,"rda",sep = "_"))
  ### write center deviation to the table
  #colnames(com_t)<-c("Frame","Deviation_from_center","Coor")
  #com_t<-as.data.frame(com_t)
  #write.table(com_t, file = "COM_shift.dat", row.names=FALSE,quote = F)
  if (getrcoor==TRUE) {
    ## label all the histgram file
    csscrcoorhist<-cbind(histx,csscrcoorhist,"CSSC")
    cshrcoorhist<-cbind(histx,cshrcoorhist,"CSH")
    watrcoorhist<-cbind(histx,watrcoorhist,"Water")
    totrcoordhist<-as.data.frame(rbind(csscrcoorhist,cshrcoorhist,watrcoorhist))
    colnames(totrcoordhist)<-c("r","Freq","Type")
    write.table(totrcoordhist, file=paste("rcoordist",firstrun, "to", lastrun,"every", steps,"dat",sep = "_"), row.names=FALSE,quote = F)
  } else {
    # ## label all the histgram file
    # cssczcoorhist<-cbind(histx,cssczcoorhist,"CSSC")
    # cshzcoorhist<-cbind(histx,cshzcoorhist,"CSH")
    # watzcoorhist<-cbind(histx,watzcoorhist,"Water")
    # totzcoordhist<-as.data.frame(rbind(cssczcoorhist,cshzcoorhist,watzcoorhist))
    # colnames(totzcoordhist)<-c("z","Freq","Type")
    # write.table(totzcoordhist, file=paste("zcoordist",firstrun, "to", lastrun,"every", steps,"dat",sep = "_"), row.names=FALSE,quote = F)
  }
}
if (0) {
  ##read com_t and totzcoordhist for plotting
  if (rzcoor) {
    com_t<-read.table(file="/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/COM_shift.dat",header = TRUE)
    ggplot(com_t,aes(x=Frame,y=Deviation_from_center,col=Coor))+geom_line()
    totzcoordhist<-read.table(file="/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/zcoordist_6_to_10_every_1_dat",header = TRUE)
    totzcoordhistnorm <- totzcoordhist %>%
      group_by(Type) %>%
      mutate(Freq_norm = Freq/((5000)*(27.5^2)*3.1415)) %>%
      ungroup()
    abstotzcoordhistnorm <- totzcoordhistnorm %>% 
      mutate(abs_z = abs(z)) %>% 
      group_by(abs_z, Type) %>% 
      summarize(mean_freq_norm = mean(Freq_norm))
    cylzcoordistplt<-{
      ggplot(abstotzcoordhistnorm, aes(x = abs_z, y = mean_freq_norm, col = Type)) +
        geom_line() +
        ggtitle("All noh")+
        #xlim(0, 160) +
        ylim(0,0.04)+
        xlab("|z| (Å)") +
        #ylab(paste0("Number density (",expression(Å^3),")"))+
        ylab(expression(paste("Number density (", Å^3, ")")))+
        
        scale_color_manual(values = c("CSSC" = "blue", "CSH" = "red", "Water" = "black")) +
        scale_x_continuous(breaks = seq(0,160,20),limits = c(0,160))+
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = c(0.5,0.9),
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )+
        geom_vline(xintercept = c(20, 85), linetype = "dashed",col="dark green")+
        geom_rect(xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, fill = "antiquewhite", alpha = 0.02,color=NA) +
        geom_rect(xmin = 85, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "cyan", alpha = 0.002,color=NA) 
      
    }
    ggsave("~/Documents/Research/MRSEC/abscylzcoordist.png",
           plot=cylzcoordistplt, dpi=300, 
           width=5.5*1.5, height=5,units="in")
    
    
    totrcoordhist<-read.table(file="/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/rcoordist_6_to_10_every_1_dat",header = TRUE)
    totrcoordhist$r<-as.numeric(totrcoordhist$r)
    totrcoordhist$Freq<-as.numeric(totrcoordhist$Freq)
    totrcoordhistnorm <- totrcoordhist %>%
      group_by(Type) %>%
      mutate(Freq_norm = Freq/(((r+binsize)^2-(r)^2)*3.1415 * (5000) * 215)) %>%
      ungroup()
    cylrcoordist<-{
      ggplot(totrcoordhist,aes(x=r,y=Freq/(((r+0.1)^2-(r)^2)*3.1415 * (5000) * 215),col=Type))+
        geom_line()+
        xlim(0, 40) +
        ylim(0,0.04)+
        xlab("r (Å)") +
        ylab(expression(paste("Number density (", Å^3, ")")))+
        scale_color_manual(values = c("CSSC" = "blue", "CSH" = "red", "Water" = "black")) +
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = c(0.5,0.9),
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )
      
    }
    ggsave("~/Documents/Research/MRSEC/cylrcoordist.png",
           plot=cylrcoordist, dpi=600, 
           width=5.5*1.5, height=5,units="in")
  }
  if (moietycoor) {
    moietycoordhist<-read.table(file="/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/moiety_zcoordist_6_to_10_every_25_dat",header = TRUE)
    colnames(moietycoordhist)<-c("z","Freq","Resname","Type")
    moietycoordhist <- moietycoordhist %>%
      mutate(Freq_norm = Freq/(5*(5000/25)*(27.5^2)*3.1415)) %>%
      ungroup()
    watzcoordhist<-abstotzcoordhistnorm[which(abstotzcoordhistnorm$Type=="Water"),]
    colnames(watzcoordhist)<-c("z","Type","Freq_norm")
    for (moiety in unique(moietycoordhist$Type)) {
      tempzcoordhist <- moietycoordhist[which(moietycoordhist$Type == moiety), ]
      print(sum(tempzcoordhist$Freq))
      abstempzcoordhist <- tempzcoordhist %>% 
        mutate(abs_z = abs(z)) %>% 
        group_by(abs_z, Resname) %>% 
        summarize(mean_freq_norm = mean(Freq_norm))
      moietyplt<- {
        ggplot() +
          geom_line(data = abstempzcoordhist, aes(x = abs_z, y = mean_freq_norm/6, col = Resname)) +
          #geom_line(data = watzcoordhist, aes(x = z, y = Freq_norm/5 )) +
          xlab("|z| (Å)") +
          ylab(expression(paste("Number density (", Å^3, ")")))+
          ggtitle(moiety) +
          ylim(0,0.0020)+
          scale_x_continuous(breaks = seq(0,160,20),limits = c(0,160))+
          theme(panel.background = element_rect(fill='transparent'),
                text=element_text(face = "bold",size=14),
                legend.position = "top",
                legend.direction = "horizontal",
                legend.title = element_blank(),
                panel.grid.major=element_line(colour = "grey",size = 0.1),
                axis.line = element_line(colour = "black",size = 0.3)
          )+
          geom_vline(xintercept = c(20, 85), linetype = "dashed",col="dark green")
      }
      # ggsave(paste0("~/Documents/Research/MRSEC/cylzcoordist_",moiety,"_.png"),
      #        plot=moietyplt, dpi=600, 
      #        width=5.5*1.5, height=5,units="in")
      print(moietyplt)
    }
    non_PHE <- unique(moietycoordhist$Type)[!unique(moietycoordhist$Type) == "PHE"]
    non_PHE_data <- subset(moietycoordhist, Type %in% non_PHE)
    non_PHE_data_sum <- non_PHE_data %>% 
      group_by(z, Resname) %>% 
      summarize(sum_freq_norm = sum(Freq_norm))
    abs_non_PHE_data_sum <- non_PHE_data_sum %>% 
      mutate(abs_z = abs(z)) %>% 
      group_by(abs_z, Resname) %>% 
      summarize(mean_freq_norm = mean(sum_freq_norm))
    non_PHE_plot <- {
      ggplot() +
        geom_line(data=abs_non_PHE_data_sum, aes(x = abs_z, y = mean_freq_norm/9, col = Resname))+
        #geom_line(data = watzcoordhist, aes(x = z, y = Freq_norm )) +
        xlab("|z| (Å)") +
        ylim(0,0.0020)+
        ylab(expression(paste("Number density (", Å^3, ")")))+
        ggtitle("Non_phe") +
        scale_x_continuous(breaks = seq(0,160,20),limits = c(0,160))+
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )+
        geom_vline(xintercept = c(20, 85), linetype = "dashed",col="dark green")
    }
    ggsave(paste0("~/Documents/Research/MRSEC/cylzcoordist_","nonPHE","_.png"),
           plot=non_PHE_plot, dpi=600, 
           width=5.5*1.5, height=5,units="in")
  }
  if (r) {
    moietyrcoordhist<-read.table(file="/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/moiety_rcoordist_6_to_10_every_25_dat",header = TRUE)
    for (moiety in unique(moietyrcoordhist$Type)) {
      temprcoordhist <- moietyrcoordhist[which(moietyrcoordhist$Type == moiety), ]
      temprcoordhist$r<-as.numeric(temprcoordhist$r)
      temprcoordhist$Freq<-as.numeric(temprcoordhist$Freq)
      
      #/ (((r+binsize/2)^2-(r-binsize/2)^2)*3.1415 * (200) * 320)
      print(
        ggplot() +
          geom_line(data = temprcoordhist, aes(x = r, y = Freq/ (((r+binsize)^2-(r)^2)*3.1415 * (200) * 320) , col = Resname)) +
          #geom_line(data = watzcoordhist, aes(x = z, y = Freq_norm / 20)) +
          xlim(0, 40) +
          labs(title = moiety, x = "z", y = "bulk density")
      )
    }
  }
  
  if (int) {
    # read data file
    moietyint <- read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_3to7_cyl_nottip_every_1_intertypes_evol.dat", header = FALSE)[,-1]
    
    # assign column names
    colnames(moietyint) <- c("PHE-PHE", "AM-AM", "AM1-PHE", "AM2-PHE", "SULF-SULF", "STER")
    
    # calculate column means and standard deviations
    means <- colMeans(moietyint[,1:5])
    stddev <- apply(moietyint[,1:5], 2, sd)
    
    # create data frame for ggplot
    df <- data.frame(Interaction = colnames(moietyint)[1:5], Average = means, StdDev = stddev)
    
    # create bar chart with error bars
    ggplot(df, aes(x = Interaction, y = Average)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = Average - StdDev, ymax = Average + StdDev), width = 0.2) +
      labs(title = "Interaction average in tips", x = "Interaction Type", y = "Average") +
      ylim(0, 0.4)
    
    # read first data file
    moietyint_tips <- read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_3to7_cyl_tip_CSH_every_25_intertypes_evol.dat", header = FALSE)[,-1]
    
    # assign column names
    colnames(moietyint_tips) <- c("PHE-PHE", "AM-AM", "AM1-PHE", "AM2-PHE", "SULF-SULF", "STER")
    
    # calculate column means and standard deviations
    means_tips <- colMeans(moietyint_tips[,1:5])
    stddev_tips <- apply(moietyint_tips[,1:5], 2, sd)
    se_nontips <- stddev_tips / sqrt(nrow(moietyint_tips))
    
    # create data frame for tips
    df_tips <- data.frame(Interaction = colnames(moietyint_tips)[1:5], Average = means_tips, StdDev = stddev_tips, SE =se_nontips,Position = "Tips")
    
    # read second data file
    moietyint_nontips <- read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/mix_27A_3to7_cyl_tip_CSSC_every_25_intertypes_evol.dat", header = FALSE)[,-1]
    
    # assign column names
    colnames(moietyint_nontips) <- c("PHE-PHE", "AM-AM", "AM1-PHE", "AM2-PHE", "SULF-SULF", "STER")
    
    # calculate column means and standard deviations
    means_nontips <- colMeans(moietyint_nontips[,1:5])
    stddev_nontips <- apply(moietyint_nontips[,1:5], 2, sd)
    se_nontips <- stddev_nontips / sqrt(nrow(moietyint_nontips))
    # create data frame for nontips
    df_nontips <- data.frame(Interaction = colnames(moietyint_nontips)[1:5], Average = means_nontips, StdDev = stddev_nontips, SE=se_nontips, Position = "Nontips")
    
    # merge the two data frames
    df <- rbind(df_tips, df_nontips)
    
    # create bar chart with error bars
    intplt<- {
      ggplot(df, aes(x = Interaction, y = Average,fill=Position)) +
        geom_bar(stat = "identity",position=position_dodge(0.9)) +
        geom_errorbar(position=position_dodge(0.9),aes(ymin = Average - StdDev, ymax = Average + StdDev), width = 0.2) +
        labs(#title = "Interaction average in tips",
          x = "Interaction Type", y = "Average") +
        ylim(0, 0.4)+
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = c(0.5,0.9),
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )
      
    }
    ggsave("~/Documents/Research/MRSEC/cylint.png",
           plot=intplt, dpi=600, 
           width=5.5*1.5, height=5,units="in")
    tipint<-read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/mix_27A_3to7_cyl_tip_everyint.dat", header = TRUE)
    
    
    colors <- c("CSH-CSH" = "#D55E00", "CSSC-CSSC" = "#0072B2", "CSH-CSSC" = "Green")
    
    # plot normalized counts by interaction type, colored by pair type
    tipplt<-{
      ggplot(tipint, aes(x = Type, y = Mean, fill = Pair)) +
        geom_bar(stat = "identity", position = position_dodge(0.9)) +
        geom_errorbar(position=position_dodge(0.9),aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
        #scale_fill_manual(values = colors) +
        labs(
          x = "Interaction type in tips", y = "Pervalence") +
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position ="bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )  +
        ylim(0,0.2)
    }
    ggsave("~/Documents/Research/MRSEC/cylint_tip.png",
           plot=tipplt, dpi=600, 
           width=5.5*1.5, height=5,units="in")
    nottipint<-read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/mix_27A_3to7_cyl_nottip_everyint.dat", header = TRUE)
    colors <- c("CSH-CSH" = "#D55E00", "CSSC-CSSC" = "#0072B2", "CSH-CSSC" = "Green")
    
    # plot normalized counts by interaction type, colored by pair type
    interiorplt<-{
      ggplot(nottipint, aes(x = Type, y = Mean, fill = Pair)) +
        geom_bar(stat = "identity", position = position_dodge(0.9)) +
        geom_errorbar(position=position_dodge(0.9),aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
        #scale_fill_manual(values = colors) +
        labs(
          x = "Interaction type in interior", y = "Pervalence") +
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )  +
        ylim(0,0.2)
    }
    ggsave("~/Documents/Research/MRSEC/cylint_interior.png",
           plot=interiorplt, dpi=600, 
           width=5.5*1.5, height=5,units="in")
    Position<-"Tip"
    tipint<-cbind(tipint,Position)
    Position<-"Interior"
    nottipint<-cbind(nottipint,Position)
    totint<-rbind(tipint,nottipint)
    totint <- totint %>%
      group_by(Type,Position) %>%
      summarise(Mean = sum(Mean),SD=sum(SD),SE=sum(SE)) 
    totintplt<-{
      ggplot(totint[-which(totint$Type=="STER"),], aes(x = Type, y = Mean, fill = Position)) +
        geom_bar(stat = "identity", position = position_dodge(0.9)
        ) +
        geom_errorbar(position=position_dodge(0.9),aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
        #scale_fill_manual(values = colors) +
        labs(
          x = "Interaction type", y = "Pervalence") +
        theme(panel.background = element_rect(fill='transparent'),
              text=element_text(face = "bold",size=14),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              panel.grid.major=element_line(colour = "grey",size = 0.1),
              axis.line = element_line(colour = "black",size = 0.3)
        )  +
        ylim(0,0.2)
    }
    ggsave("~/Documents/Research/MRSEC/cylint.png",
           plot=totintplt, dpi=600, 
           width=5.5*1.5, height=5,units="in")
    
  }
}






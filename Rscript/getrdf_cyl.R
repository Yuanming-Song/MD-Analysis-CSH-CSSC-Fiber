sysname="mix_27A"
perbatch<-800
dir="data/"
step=10
firstrun=1
lastrun=8
moietyname <- c("N1", "C1", "O1", "C2", "C3", "S", "N2", "C4", "O2", "C5", "C6", "C7", "C8", "C9", "C10","S1")
moietyname2 <- c("S2", "C11", "C12", "C13", "C14", "N3", "C15", "C16", "C17", "C18", "C19", "C20", "O3", "O4", "N4")
phemoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
ammoietyname <- c("N1", "C1", "O1", "C3", "S", "N2", "C4", "O2", "S1")
thimoietyname <- c("C3", "S","S1")
phemoietyname2 <- c("C12", "C13", "C14", "C16","C18","C19")
ammoietyname2 <- c("S2", "C11", "N3", "C17", "C20", "O3", "O4", "N4")
thimoietyname2 <- c("C11", "S2")

### first frame to be considered
startframe<-1
### for pdb & traj file names
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
#dcdname<-"test.dcd"
pdbname<-"data/mix_27A.pdb"
### how many residues in total
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
totcssc<-maxcssc
totcsh<-maxcsh-maxcssc
### analysis name for output name
analysis="_com_rdf_phe"
### binsize and max for normalization
binsize<-0.1
binmax<-300
binmin<-0
nbins <- ceiling((binmax - binmin) / binsize)
histx<-seq(binmin + binsize / 2, binmax - binsize / 2, by = binsize)
#####  Volume of sphere shell for rdf normalization
nvshell<-4*3.14/3
resnamelist<-c("CSSC","CSH")
### output data file name
outname=paste0("data/",sysname, analysis, ".dat")
outname2=paste0("data/",sysname, analysis, "_normalized.dat")
comoutname<-"COM.rda"
comoutname1<-"COM_phe.rda"
comoutname2<-"COM_am.rda"
comoutname3<-"COM_thiol.rda"
load("COMcutoff.rda")
### load all the library
library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)
library(bigmemory)
### register the core
registerDoParallel(cores = detectCores())
#####   function for calculating distance
euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))
}
#####   function for updating histogram
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
pbcdis <- function(comi, cur, currbox) {
  # Loop through x, y, and z coordinates
  diff <- cur - comi
  for (l in 1:3) {
     # Calculate the difference between the current coordinate and the center of mass coordinate
    absdiff <- abs(diff[l])  # Calculate the absolute value of the difference
    # Update coordinate only when necessary to account for periodic boundary conditions
    if (absdiff > currbox[l]/2) {
      diff[l] <- diff[l] - (diff[l]/abs(diff[l])) * currbox[l] * ceiling((abs(diff[l]) - currbox[l]/2) / currbox[l])
    } 
  }
  # Return updated current coordinates
  sqrt(sum((diff)^2))
}
#####   function for getting COM and update list
getCOM <- function(startresi,finalresi,framenum,moiety1,moiety2) {
  ##### create an empty com matrix
  COMlist<-list()
  counter<-1
  for (i in startresi:finalresi) {
    tempCOMlist<-list()
    if (i<=maxcssc) {
      resname<-"CSSC"
    } else {
      resname<-"CSH"
    }
    #####   select input res as input for getdist()
    refmoietysele<-atom.select(mypdb,resno=i,resid=resname,elety=moiety1,operator = "AND")
    #####   get coordinate
    refmoietyxyz<-mydcd[framenum,refmoietysele$xyz]
    #####   calculate COM of res i
    com1<-as.vector(com.xyz(refmoietyxyz))
    if (resname=="CSSC") {
      #####   select input res as input for getdist()
      refmoietysele<-atom.select(mypdb,resno=i,resid=resname,elety=moiety2,operator = "AND")
      #####   get coordinate
      refmoietyxyz<-mydcd[framenum,refmoietysele$xyz]
      #####   calculate COM of res i
      com2<-as.vector(com.xyz(refmoietyxyz))
      tempCOMlist[[1]]<-com1
      tempCOMlist[[2]]<-com2
    } else {
      tempCOMlist[[1]]<-com1
    }
    COMlist[[counter]]<-tempCOMlist
    ##### update com list
    counter<-counter+1
  }
  COMlist
}
getdist <- function(resname,startresi,finalresi,framenum,resname1,startresk,finalresk,COMlist1,COMlist2=NULL) {
  if (is.null(COMlist2)) {
    COMlist2 <- COMlist1
    ##### k should be updated every time
    newk<-TRUE
  } else {
    ##### k shouldn't be updated every time
    newk<-FALSE
  }
  if (resname == resname1) {
    ##### final i won't have any new pair to be calculated
    finalresi<-finalresi-1
  }
  box<-mycell[framenum,]
  ##### create an empty histogram matrix
  histogram <- matrix(0, nrow = nbins, ncol = 2)
  histogram[, 1] <- histx
  #####   get COM of res i
  for (i in startresi:finalresi) {
    #####   get COM of res i
    comi<-COMlist1[[i]][[1]]
    if (resname == resname1 & newk==TRUE) {
      startresk<-i+1
    }
    for (k in startresk:finalresk) {
      if (i==k) {
        next
      }
      #####   get COM of res k
      comk<-COMlist2[[k]][[1]]
      #####   calculate distance
      dis<-pbcdis(comi,comk,box)
      #####   update matrix
      histogram<-update_histogram(histogram,dis)
      if (k<=maxcssc) {
        #####   calculate COM of res k
        comk<-COMlist2[[k]][[2]]
        #####   calculate distance
        dis<-pbcdis(comi,comk,box)
        #####   update matrix
        histogram<-update_histogram(histogram,dis)
      }
    } 
    if (i<=maxcssc) {
      #####   get COM of res i and second CSX unit
      comi<-COMlist1[[i]][[2]]
      for (k in startresk:finalresk) {
        if (i==k) {
          next
        }
        #####   get COM of res k
        comk<-COMlist2[[k]][[1]]
        #####   calculate distance
        dis<-pbcdis(comi,comk,box)
        #####   update matrix
        histogram<-update_histogram(histogram,dis)
        if (k<=maxcssc) {
          #####   calculate COM of res k
          comk<-COMlist2[[k]][[2]]
          #####   calculate distance
          dis<-pbcdis(comi,comk,box)
          #####   update matrix
          histogram<-update_histogram(histogram,dis)
        }
      } 
    }
  }
  ##### print freqency
  histogram[,2]/sum(histogram[,2])
}
getcutoff<- function(pair1,pair2,COM1,COM2) {
  if (pair1==pair2) {
    cutofflist[[paste0(pair1,"-",pair2)]][COM1,COM2]
  } else {
    cutofflist[["CSH-CSSC"]][COM1,COM2]
  }
}
getnodes<- function(resid,moietyind,COM) {
  if (COM=="AM") {
    corrnum<--2
  } else if (COM=="PHE") {
    corrnum<--1
  } else {
    corrnum<-0
  }
  if (i<maxcssc) {
    if (moietyind==1) {
      corrnum<-corrnum-3
    }
    resid*6+corrnum
  } else {
    maxcssc*3+resid*3+corrnum
  }
}
getContactCom <- function(framenum,COM1,COM2=NULL) {
  if (is.null(COM2)) {
    COM2 <- COM1
    ##### k should be updated every time
    newk<-TRUE
    finalresi<-maxcsh-1
  } else {
    ##### k shouldn't be updated every time
    newk<-FALSE
    finalresi<-maxcsh
  }
  COMlist1<-get(paste0("COMlist",COM1))[[framenum]]
  COMlist2<-get(paste0("COMlist",COM2))[[framenum]]
  box<-mycell[framenum,]
  selflist<-list()
  reflist<-list()
  #####   get COM of res i
  for (i in 1:finalresi) {
    if (newk==TRUE) {
      startresk<-i+1
    } else {
      startresk<-1
    }
    if (i<=maxcssc) {
      pair1<-"CSSC"
      #####   get COM of res i and second CSX unit
      comi<-COMlist1[[i]][[2]]
      cominode<-getnodes(i,2,COM1)
      for (k in startresk:maxcsh) {
        if (i==k) {
          next
        }
        if (k<=maxcssc) {
          contactcutoff<-getcutoff(pair1,"CSSC",COM1,COM2)
          #####   calculate COM of res k
          comk<-COMlist2[[k]][[2]]
          #####   calculate distance
          dis<-pbcdis(comi,comk,box)
          if (dis<=contactcutoff) {
            selflist<-append(selflist,cominode)
            reflist<-append(reflist,getnodes(k,2,COM2))
          }
        } else {
          contactcutoff<-getcutoff(pair1,"CSH",COM1,COM2)
        }
        #####   get COM of res k
        comk<-COMlist2[[k]][[1]]
        #####   calculate distance
        dis<-pbcdis(comi,comk,box)
        if (dis<=contactcutoff) {
          selflist<-append(selflist,cominode)
          reflist<-append(reflist,getnodes(k,1,COM2))
        }
      } 
    } else {
      pair1<-"CSH"
    }
    #####   get COM of res i
    comi<-COMlist1[[i]][[1]]
    cominode<-getnodes(i,1,COM1)
    for (k in startresk:maxcsh) {
      if (i==k) {
        next
      }
      if (k<=maxcssc) {
        contactcutoff<-getcutoff(pair1,"CSSC",COM1,COM2)
        #####   calculate COM of res k
        comk<-COMlist2[[k]][[2]]
        #####   calculate distance
        dis<-pbcdis(comi,comk,box)
        if (dis<=contactcutoff) {
          selflist<-append(selflist,cominode)
          reflist<-append(reflist,getnodes(k,2,COM2))
        }
      } else {
        contactcutoff<-getcutoff(pair1,"CSH",COM1,COM2)
      }
      #####   get COM of res k
      comk<-COMlist2[[k]][[1]]
      #####   calculate distance
      dis<-pbcdis(comi,comk,box)
      if (dis<=contactcutoff) {
        selflist<-append(selflist,cominode)
        reflist<-append(reflist,getnodes(k,1,COM2))
      }
    } 
  }
  out<-c()
  out[[1]]<-selflist
  out[[2]]<-reflist
  out
}
getEdgesCOM <- function(framenum) {
  nodeListI<-c()
  nodeListJ<-c()
  potential_pairs <- c("AM", "PHE", "THI")
  for (COMI in 1:length(potential_pairs)) {
    for (COMJ in 1:length(potential_pairs)) {
      tempcontact<-getContactCom(framenum,COMI,COMJ)
      nodeListI<-append(nodeListI,tempcontact[[1]])
      nodeListJ<-append(nodeListI,tempcontact[[2]])
    }
  }
  weight<-list()
  interactType<-list()
  for (l in 1:length(nodeListI)) {
    Ni<-nodeListI[[l]]
    Nj<-nodeListJ[[l]]
    if (Ni != Nj) {
      pair<-sort(c(Ni,Nj))
      pairindex<-paste(pair[1],pair[2])
      if(is.numeric(weight[[pairindex]])) {
        weight[[pairindex]]<-weight[[pairindex]]+1
      } else {
        weight[[pairindex]]<-1
      }
      if (interactTypeAnal==TRUE) {
        interaction<-paste(network[[pair[1]]]$nname,network[[pair[2]]]$nname)
        if (interaction == "PHE PHE") {
          interactType[[pairindex]] <- "PHE-PHE"
        } else if (interaction %in% c("AM1 AM2", "AM2 AM1", "AM1 AM1", "AM2 AM2")) {
          interactType[[pairindex]] <- "AM-AM"
        } else if (interaction == "CYS CYS") {
          interactType[[pairindex]] <- "DISU"
        } else if (interaction %in% c("AM1 PHE", "PHE AM1")) {
          interactType[[pairindex]] <- "AM1-PHE"
        } else if (interaction %in% c("AM2 PHE", "PHE AM2")) {
          interactType[[pairindex]] <- "AM2-PHE"
        } else {
          interactType[[pairindex]] <- "STER"
        }
      }
    }
  }
  # Create an empty matrix
  weight_matrix <- big.matrix(init=0, nrow = length(weight), ncol = 3)
  # Loop through the weight list
  for (i in 1:length(weight)) {
    # Extract the name and numbers
    name <- names(weight[i])
    numbers <- as.numeric(unlist(strsplit(name, " ")))
    # Set the corresponding row in the matrix
    weight_matrix[i, ] <- as.numeric(c(numbers, weight[[i]]))
  }
  # Print the matrix
  out[[1]]<-as.matrix(weight_matrix)
  if (interactTypeAnal==TRUE) {
    # Create an empty matrix
    interactType_matrix <- big.matrix(init="", nrow = length(interactType), ncol = 3)
    # Loop through the interactType list
    for (i in 1:length(interactType)) {
      # Extract the name and numbers
      name <- names(interactType[i])
      numbers <- unlist(strsplit(name, " "))
      # Set the corresponding row in the matrix
      interactType_matrix[i, ] <- c(numbers, interactType[[i]])
    }
    out[[2]]<-interactType_matrix
  } else {
    out[[2]]<-c()
  }
  out
}
setnetwork<-function() {
  for (count in 1:length(moietynamelist)) {
    moid<<-moid+1
    nname<-moietynamelist[count]
    ntype<-typenameslist[count]
    network[[moid]]<<-list(
      resname = resname,
      resid = resid,
      moid = moid,
      nname = nname,
      ntype = ntype
    )
  }
}
network<-list()
moid<-0
moietynames <- c("AM", "PHE", "THI")
typenames<-c("DIP","NOP","THIO")
resname<-"CSSC"
moietynamelist<-c(moietynames,moietynames)
typenameslist<-c(typenames,typenames)
for (resid in mincssc:maxcssc) {
  setnetwork()
}
resname<-"CSH"
moietynamelist<-moietynames
typenameslist<-typenames
for (i in mincsh:maxcsh) {
  setnetwork()
}
#####  function for normalize rdf
normrdf<- function(dislist,type,endframe) {
  #####  cut data according to bin size & max
  rdf<-as.data.frame(dislist)
  colnames(rdf)<-c("V1","V2")
  ##### (unit: A^-3)
  globalrho<-endframe/box
  #####  normalize according to shell volume and overall concentration  (unit: 1/ A^-3 * A^3 )
  rdf<-rdf %>% mutate(g=V2/(((V1+binsize)^3-(V1)^3)*nvshell*globalrho))
  #####  add Type column
  rdf<-rdf %>% mutate(Type=type)
  #####  delete [x,y] column
  rdf
}
#####   read dcd files for cell dimension 
mycell<-read.dcd(dcdname,cell=T)
box<-mean(mycell[,1])*mean(mycell[,2])*mean(mycell[,3])
print(paste(mycell[,1],mycell[,2],mycell[,3]))
nframes<-nrow(mycell)
if (1) {
  #####	for testing...
  if (0) {
    nframes<-1
    startframe<-nframes
    print(nframes)
  } 
  #####   read pdb files
  if (file.exists(comoutname1)) {
    load(comoutname1)
  } else {
    mypdb<-read.pdb(pdbname)
    mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
    COMlistPHE<-foreach (framenum=startframe:nframes) %dopar% {
      getCOM(mincssc,maxcsh,framenum,phemoietyname,phemoietyname2) 
    }  
    save(COMlistPHE,file = comoutname1)
  }
  if (file.exists(comoutname2)) {
    load(comoutname2)
  } else {
    mypdb<-read.pdb(pdbname)
    mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
    COMlistAM<-foreach (framenum=startframe:nframes) %dopar% {
      getCOM(mincssc,maxcsh,framenum,ammoietyname,ammoietyname2) 
    }  
    save(COMlistAM,file = comoutname2)
  }
  if (file.exists(comoutname3)) {
    load(comoutname3)
  } else {
    mypdb<-read.pdb(pdbname)
    mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
    COMlistTHI<-foreach (framenum=startframe:nframes) %dopar% {
      getCOM(mincssc,maxcsh,framenum,thimoietyname,thimoietyname2) 
    }  
    save(COMlistTHI,file = comoutname3)
  }
  outcssc<-0
  outcsh<-0
  outcshcssc<-0
  for (batch in 1:ceiling((nframes)/perbatch)) {
    begframe<-(batch-1)*perbatch+1
    if (batch!=ceiling(nframes/perbatch)) {
      endframe<-batch*perbatch
    } else {
      endframe<-nframes
    }
    start_time <- Sys.time()
    print(paste("start",batch,start_time))
    #####   run getdistbig for each frame in parallel
    outcssctemp<-foreach (framenum=begframe:endframe,.combine="+") %dopar% {
      getdist("CSSC",mincssc,maxcssc,framenum,"CSSC",mincssc,maxcssc,COMlistPHE[[framenum]],COMlistPHE[[framenum]]) 
    }  
    print("cssc done")
    outcshtemp<-foreach (framenum=begframe:endframe,.combine="+") %dopar% {
      getdist("CSH",mincsh,maxcsh,framenum,"CSH",mincsh,maxcsh,COMlistPHE[[framenum]],COMlistPHE[[framenum]]) 
    }
    print("csh done")
    outcshcssctemp<-foreach (framenum=begframe:endframe,.combine="+") %dopar% {
      getdist("CSSC",mincssc,maxcssc,framenum,"CSH",mincsh,maxcsh,COMlistPHE[[framenum]],COMlistPHE[[framenum]])
    }
    print("csh cssc done")
    outcssc<-outcssc+outcssctemp
    outcsh<-outcsh+outcshtemp
    outcshcssc<-outcshcssc+outcshcssctemp
    outcsscm<-normrdf(cbind(histx,outcssc),"CSSC-CSSC",endframe)
    outcshm<-normrdf(cbind(histx,outcsh),"CSH-CSH",endframe)
    outcshcsscm<-normrdf(cbind(histx,outcshcssc),"CSH-CSSC",endframe)
    outm<-rbind(outcsscm,outcshm,outcshcsscm)
    out<-cbind(histx,outcssc,outcsh,outcshcssc)
    #####   write combined getdistbig() outpu to outname
    write.table(outm, file=outname2,col.names = F,row.names=F,quote = F)
    write.table(out, file=outname,col.names = F,row.names=F,quote = F)
    # Stop tracking time
    end_time <- Sys.time()
    # Calculate the elapsed time
    elapsed_time <- end_time - start_time
    print(paste(begframe,endframe,batch,elapsed_time))
  }
}else { 
  outm<-read.table(file = outname)
  print(paste("read",outname))
  #####  pair matrix for different combination
  analysism<-matrix(nrow=3,ncol=2)
  analysism[1,]<-c("CSH","CSH")
  analysism[2,]<-c("CSSC","CSSC")
  analysism[3,]<-c("CSSC","CSH")
  rdfnorm<-c()
  #####  loop through each pair
  totalpair<-0
  for (i in 1:nrow(analysism)) {
    #####  select dis value corresponding to pair indicated only
    subm<-outm[which(outm[,3]==analysism[i,][1]&outm[,5]==analysism[i,][2]),]
    #####  select dis value corresponding to pair indicated only
    dislist1<-as.numeric(subm[,1])
    if (analysism[i,][1]==analysism[i,][2]){
      #####  largest resid - smallest resid for total number of residue in the system
      subtotalres<-max(as.numeric(subm[,4]))-min(as.numeric(subm[,2]))+1
      #####  Total pairs available
      print(subtotalres)
      normfactor<-(subtotalres-1)*subtotalres/2
      if (analysism[i,][1]=="CSSC") {
        normfactor<-normfactor*2
      }
    } else {
      #####  Total pairs available
      normfactor<-2*(max((as.numeric(subm[,2]))-min(as.numeric(subm[,2])))+1)*(max((as.numeric(subm[,4]))-min(as.numeric(subm[,4])))+1)
      print(max((as.numeric(subm[,2]))-min(as.numeric(subm[,2])))+1)
      print(max((as.numeric(subm[,4]))-min(as.numeric(subm[,4])))+1)
    }
    #####  type indicator for plotting
    type<-paste(analysism[i,][1],analysism[i,][2],sep="-")
    #####  run normrdf
    rdfnorm<-rbind(rdfnorm,normrdf(dislist1,normfactor,type))
    #####  maximum rdf tracking
    print(paste(i, type, length(dislist1),max(as.numeric(rdfnorm[which(rdfnorm[,4]==paste(analysism[i,][1],analysism[i,][2],sep="-")),3]))))
    print(normfactor)
    totalpair<-totalpair+normfactor
  }
  #####  get a list of distance only and run normrdf()
  rdfnorm<-rbind(rdfnorm,normrdf(as.numeric(outm[,1]),totalpair,"All"))
  print(totalpair)
  
  #####   write normalized rdf to outname2
  write.table(rdfnorm, file=outname2,col.names = F,row.names=F,quote = F)
  
}
if (0) {
  COMlist1<-COMlistAM
  COMlist2<-COMlistPHE
  lowcutoff<-2.3
  highcutoff<-3.0
  for (j in 1:5) {
    for (i in mincssc:(maxcsh-1)) {
      #####   get COM of res i
      comi<-COMlist1[[j]][[i]][[1]]
      for (k in 1:maxcssc) {
        if (k==i) {
          next
        }
        #####   get COM of res k
        comk<-COMlist2[[j]][[k]][[1]]
        #####   calculate distance
        dis<-pbcdis(comi,comk,mycell[j,])
        if (dis<highcutoff & dis>lowcutoff) {
          print(paste(j,i,k,dis))
          #break
        }
        if (k<=maxcssc) {
          comk<-COMlist2[[j]][[k]][[2]]
          #####   calculate distance
          dis<-pbcdis(comi,comk,mycell[j,])
          if (dis<highcutoff & dis>lowcutoff) {
            print(paste(j,i,k,dis))
            #break
          }
        }
      }
      if (dis<highcutoff & dis>lowcutoff) {
        #break
      }
      if (i<=maxcssc) {
        comi<-COMlist1[[j]][[i]][[2]]
        for (k in 1:maxcssc) {
          if (k==i) {
            next
          }
          #####   get COM of res k
          comk<-COMlist2[[j]][[k]][[1]]
          #####   calculate distance
          dis<-pbcdis(comi,comk,mycell[j,])
          if (dis<highcutoff & dis>lowcutoff) {
            print(paste(j,i,k,dis))
            #break
          }
          if (k<=maxcssc) {
            comk<-COMlist2[[j]][[k]][[2]]
            #####   calculate distance
            dis<-pbcdis(comi,comk,mycell[j,])
            if (dis<highcutoff & dis>lowcutoff) {
              print(paste(j,i,k,dis))
              #break
            }
          }
        }
      }
    }
    if (dis<highcutoff & dis>lowcutoff) {
      #break
    }
    print(j)
  }
  pherdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_phe_normalized.dat")
  pherdfcylplt<-{
    ggplot(data = pherdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "Phe COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  nonpherdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_nonphe_normalized.dat")
  nonpherdfcylplt<-{
    ggplot(data = nonpherdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "Non Phe COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  amrdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_am_normalized.dat")
  amrdfcylplt<-{
    ggplot(data = amrdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "AM COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  rdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_normalized.dat")
  rdfcylplt<-{
    ggplot(data = rdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "CSX COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  phenonpherdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_phe_nonphe_normalized.dat")
  phenonpherdfcylplt<-{
    ggplot(data = phenonpherdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "Non Phe COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  amthirdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_am_thi_normalized.dat")
  amthirdfcylplt<-{
    ggplot(data = amthirdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "AM-THI COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  phethirdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_phe_thi_normalized.dat")
  phethirdfcylplt<-{
    ggplot(data = phethirdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "PHE-THI COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  thirdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_thi_normalized.dat")
  thirdfcylplt<-{
    ggplot(data = thirdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "THI COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  pheamrdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_phe_am_normalized.dat")
  pheamrdfcylplt<-{
    ggplot(data = pheamrdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "PHE-AM COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  phethirdfcyl<-read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/mix_27A_com_rdf_phe_thi_normalized.dat")
  phethirdfcylplt<-{
    ggplot(data = phethirdfcyl)+
      geom_line(aes(x=V1,y=V3,col=V4))+
      coord_cartesian(
        xlim = c(0,25)
        ,
        ylim = c(0,25)
      )+
      labs(y = "g", x = "PHE-THI COM Distance (Å)",col="Type")+
      plttheme+
      theme(legend.direction = "vertical",
            legend.position = c(0.8,0.8))
  }
  # ggsave("/Users/song/Documents/Research/MRSEC/CSH-CSSC/Restrained_Plot/CylRDFcom.png",
  #       plot=rdfcylplt, dpi=1100,width=4, height=4,units="in")
  #ggsave("/Users/song/Documents/Research/MRSEC/CSH-CSSC/Restrained_Plot/CylRDFcom_phe.png",
  #        plot=pherdfcylplt, dpi=1100,width=4, height=4,units="in")
  #ggsave("/Users/song/Documents/Research/MRSEC/CSH-CSSC/Restrained_Plot/CylRDFcom_nonphe.png",
  #      plot=nonpherdfcylplt, dpi=1100,width=4, height=4,units="in")
  ggplotly(nonpherdfcylplt)
  ggplotly(rdfcylplt)
  ggplotly(phenonpherdfcylplt)
  
  ggplotly(pherdfcylplt)
  ggplotly(amrdfcylplt)
  ggplotly(thirdfcylplt)
  
  ggplotly(amthirdfcylplt)
  ggplotly(phethirdfcylplt)
  ggplotly(pheamrdfcylplt)
  # List of matrix types
  matrix_types <-unique(amthirdfcyl$V4)
  # List of potential pairs
  potential_pairs <- c("am", "phe", "thi")
  # Initialize the cutofflist as an empty list
  cutofflist <- list()
  # Loop through matrix types
  for (matrix_type in matrix_types) {
    # Initialize the symmetric matrix for the current matrix type
    matrix_data <- matrix(0, nrow = 3, ncol = 3, dimnames = list(potential_pairs, potential_pairs))
    for (i in 1:length(potential_pairs)) {
      for (j in 1:length(potential_pairs)) {
        pair1 <- potential_pairs[i]
        pair2 <- potential_pairs[j]
        if (i==j) {
          dataframename<-paste0(pair1,"rdfcyl")
        } else {
          dataframename<-paste0(pair1,pair2,"rdfcyl")
        }
        if(exists(dataframename)) {
          temdata<-get(dataframename)
          value<-temdata$V1[which(temdata$V3==max(temdata$V3[which(temdata$V4==matrix_type)]))]
          matrix_data[pair1, pair2] <- value
          if (i!=j) {
            matrix_data[pair2, pair1] <- value
          }
        }
      }
    }
    colnames(matrix_data) <- toupper(colnames(matrix_data))
    rownames(matrix_data) <-  colnames(matrix_data)
    cutofflist[[matrix_type]]<-matrix_data
    
  }
}

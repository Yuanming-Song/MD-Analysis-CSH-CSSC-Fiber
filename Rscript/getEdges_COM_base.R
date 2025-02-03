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
  if (resid<=maxcssc) {
    if (moietyind==1) {
      corrnum<-corrnum-3
    }
    resid*6+corrnum
  } else {
    maxcssc*3+resid*3+corrnum
  }
}
getContactCom <- function(framenum,COM1,COM2=NULL) {
  if (is.null(COM2)|COM1==COM2) {
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
getEdgesStackCOM <- function(framenum) {
  out<-list()
  nodeListI<-c()
  nodeListJ<-c()
  for (COMI in 1:length(potential_pairs)) {
    for (COMJ in COMI:length(potential_pairs)) {
      tempcontact<-getContactCom(framenum,potential_pairs[COMI],potential_pairs[COMJ])
      nodeListI<-append(nodeListI,tempcontact[[1]])
      nodeListJ<-append(nodeListJ,tempcontact[[2]])
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
        } else if (interaction %in% c("AM AM")) {
          interactType[[pairindex]] <- "AM-AM"
        } else if (interaction == "THI THI") {
          interactType[[pairindex]] <- "DISU"
        } else if (interaction %in% c("AM PHE", "PHE AM")) {
          interactType[[pairindex]] <- "AM1-PHE"
        } else {
          interactType[[pairindex]] <- "THI-PHE"
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
  weight_matrix<-as.matrix(weight_matrix)
  weight_matrix<-rbind(weight_matrix,weight_matrix[,c(2,1,3)])
  weight_matrix<-weight_matrix[!duplicated(weight_matrix),]
  attr(weight_matrix,"n")<-length(weight)
  # Print the matrix
  out[[1]]<-weight_matrix
  if (interactTypeAnal==TRUE) {
    # Create an empty matrix
    interactType_matrix <- matrix(nrow = length(interactType), ncol = 3)
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
for (resid in mincsh:maxcsh) {
  setnetwork()
}
treatedgelist<-function(framenum) {
  edgelist<-edges[[framenum]][[1]]
  edgelist<-rbind(edgelist,edgelist[,c(2,1,3)])
  edgelist<-edgelist[!duplicated(edgelist),]
  attr(edgelist,"n")<-length(network)
  edgelist
}
library(bio3d)
library(doParallel)
registerDoParallel(cores = detectCores())
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
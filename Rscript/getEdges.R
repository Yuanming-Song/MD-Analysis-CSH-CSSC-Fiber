test=FALSE
interactTypeAnal=TRUE
# set myNodeList data/CSH_CSSC_1to1_50mM_water_ions.psf.nodes
myNodeList<-"data/CSH_CSSC_mix_27A_cyl_lattice_redo.psf.nodes"
# if two atoms are further than avoidcutoff, then avoid further computation of thse two residues 
avoidcutoff<-40
maxcssc<-411
totres<-763
step<-10
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
#dcdname<-"test.dcd"
pdbname<-"data/mix_27A.pdb"
begframes<-1
perbatch<-40
resname<-c("CSH","CSSC")
outname<-"mix_27A_edges.rda"
# #------------------------------------------------------------------------

# # Cutoff Parameters 
# set carbonCutoff 5.4
carbonCutoff<-5.4
# set Cutoff 4.6
Cutoff<-4.6
# set sulfurCutoff 6.3
sulfurCutoff<-6.3
# 

# #---------------------END of USER INTERFACE ---------------------------------------------------
# #
# #
# #------------------------------------------------------------------------
# 
# extract resid
getResid <- function(ind) {
  maxcssc_ind <- max(residlist[[maxcssc]])
  if (ind <= maxcssc_ind) {
    resid <- ceiling(ind / 52)
  } else {
    resid <- ceiling((ind - maxcssc_ind) / 27) + maxcssc
  }
  return(resid)
}
myMeasureContacts <- function(sel, ref, cutoff, box, frame, avoidlist = NULL) {
  # proc myMeasureContacts {sel ref cutoff box } {
  # Check if avoidlist is provided
  if (is.null(avoidlist)) {
    avoidlist <- list()
    for (i in 1:totres) {
      avoidlist[[i]]<-c(0)
    }
  }
  #   set selpos [$sel get {x y z}]
  selpos<-as.vector(mydcd[frame,sel$xyz])
  #   set refpos [$ref get {x y z}]
  refpos<-as.vector(mydcd[frame,ref$xyz])
  #   set selind [$sel get index]
  selind<-sel$atom
  #   set refind [$ref get index]
  refind<-ref$atom
  #   set selres [$sel get resid]
  selres<-mypdbatom[which(mypdbatom$eleno%in%selind),]$resno
  #   set refres [$ref get resid]
  refres<-mypdbatom[which(mypdbatom$eleno%in%refind),]$resno
  #   set selList {}
  selList<-c()
  #   set refList {}
  refList<-c()
  #foreach vec1 $selpos ind1 $selind res1 $selres {
  for (i in 1:(length(selpos)/3)) {
    #foreach vec2 $refpos ind2 $refind res2 $refres {
    for (j in 1:(length(refpos)/3)) {
      #       set dist {}
      dist<-c()
      #       if {$res1 != $res2} {
      if (selres[i]!=refres[j]) {
        selresid<- getResid(selind[i])
        refresid<- getResid(refind[j])
        if (selresid %in% avoidlist[[refresid]] | refresid %in% avoidlist[[selresid]]) {
          next
        } else {
          #         set dif [vecsub $vec1 $vec2]
          dif<-selpos[seq(i*3-2,i*3,1)]-refpos[seq(j*3-2,j*3,1)]
          #         foreach x $dif  {
          for (k in 1:length(dif)) {
            #           lappend dist [expr {$x - $box*round($x/$box)}]
            dist[k]<-dif[k]-box[k]*round(dif[k]/box[k])
            #         }
          }
          #         if {[veclength $dist] <= $cutoff} {
          if (sqrt(sum(dif^2))<=cutoff) {
            #           lappend selList $ind1
            selList<-append(selList,selind[i])
            #           lappend refList $ind2
            refList<-append(refList,refind[j])
            #         }
          } else if (sqrt(sum(dif^2))>=avoidcutoff) {
            avoidlist[[selresid]]<-append(avoidlist[[selresid]],refresid)
            avoidlist[[refresid]]<-append(avoidlist[[refresid]],refresid)
          }
        }
        #       }
      }
      #     }
    }
    #   }
  }
  #   return [list $selList $refList]
  out<-c()
  out[[1]]<-selList
  out[[2]]<-refList
  out[[3]]<-avoidlist
  out
  # 
  # }
}
myMeasureContactsCOM <- function(sel, ref, cutoff, box, frame, avoidlist = NULL) {
  # proc myMeasureContacts {sel ref cutoff box } {
  # Check if avoidlist is provided
  if (is.null(avoidlist)) {
    avoidlist <- list()
    for (i in 1:totres) {
      avoidlist[[i]]<-c(0)
    }
  }
  #foreach vec1 $selpos ind1 $selind res1 $selres {
  for (i in mincssc:maxcsh) {
    #foreach vec2 $refpos ind2 $refind res2 $refres {
    for (j in inij:maxcsh) {
      #       set dist {}
      dist<-c()
      #       if {$res1 != $res2} {
      if (selres[i]!=refres[j]) {
        selresid<- getResid(selind[i])
        refresid<- getResid(refind[j])
        if (selresid %in% avoidlist[[refresid]] | refresid %in% avoidlist[[selresid]]) {
          next
        } else {
          #         set dif [vecsub $vec1 $vec2]
          dif<-selpos[seq(i*3-2,i*3,1)]-refpos[seq(j*3-2,j*3,1)]
          #         foreach x $dif  {
          for (k in 1:length(dif)) {
            #           lappend dist [expr {$x - $box*round($x/$box)}]
            dist[k]<-dif[k]-box[k]*round(dif[k]/box[k])
            #         }
          }
          #         if {[veclength $dist] <= $cutoff} {
          if (sqrt(sum(dif^2))<=cutoff) {
            #           lappend selList $ind1
            selList<-append(selList,selind[i])
            #           lappend refList $ind2
            refList<-append(refList,refind[j])
            #         }
          } else if (sqrt(sum(dif^2))>=avoidcutoff) {
            avoidlist[[selresid]]<-append(avoidlist[[selresid]],refresid)
            avoidlist[[refresid]]<-append(avoidlist[[refresid]],refresid)
          }
        }
        #       }
      }
      #     }
    }
    #   }
  }
  #   return [list $selList $refList]
  out<-c()
  out[[1]]<-selList
  out[[2]]<-refList
  out[[3]]<-avoidlist
  out
  # 
  # }
}

getEdges <- function(frame) {
  out<-list()
  #     set sulfurContacts [myMeasureContacts $sulfurAtoms $sulfurAtoms $sulfurCutoff $box]
  sulfurContacts<-myMeasureContacts(sulfurAtoms,sulfurAtoms,sulfurCutoff,mycell[frame,],frame)
  #Ai<-sulfurContacts[[1]]
  #Aj<-sulfurContacts[[2]]
  #     set Contacts [myMeasureContacts $nocarbons $nohs $Cutoff $box]
   Contacts<-myMeasureContacts(nocarbons,nohs,Cutoff,mycell[frame,],frame,sulfurContacts[[3]])
  # #     set carbonContacts [myMeasureContacts $carbonAtoms $carbonAtoms $carbonCutoff $box]
   carbonContacts<-myMeasureContacts(carbonAtoms,carbonAtoms,carbonCutoff,mycell[frame,],frame,Contacts[[3]])
  # #     set Ai [concat [lindex $carbonContacts 0] [lindex $Contacts 0] [lindex $sulfurContacts 0]]
   Ai<-append(append(carbonContacts[[1]],Contacts[[1]]),sulfurContacts[[1]])
  # #     set Aj [concat [lindex $carbonContacts 1] [lindex $Contacts 1] [lindex $sulfurContacts 1]]
   Aj<-append(append(carbonContacts[[2]],Contacts[[2]]),sulfurContacts[[2]])
  #     set weight {}
  weight<-list()
  #     set interactType {}
  interactType<-list()
  #     foreach thing1 $Ai thing2 $Aj {
  for (l in 1:length(Ai)) {
    #       #			puts -nonewline "**** atoms are $thing1 $thing2 "
    #       set Ni [dict get $indices $thing1]
    Ni<-indices[[Ai[[l]]]]
    #       set Nj [dict get $indices $thing2]
    Nj<-indices[[Aj[[l]]]]
    #       if {$Ni != $Nj} {
    if (Ni != Nj) {
      #         set pair [lsort -integer [list $Ni $Nj]]
      pair<-sort(c(Ni,Nj))
      #         #				puts "nodes are $pair"
      #         dict incr weight $pair
      pairindex<-paste(pair[1],pair[2])
      if(is.numeric(weight[[pairindex]])) {
        weight[[pairindex]]<-weight[[pairindex]]+1
      } else {
        weight[[pairindex]]<-1
      }
      if (interactTypeAnal==TRUE) {
      #         if {![dict exists $interactType $pair]} {
      #           set interaction [list [dict get $network [lindex $pair 0] nname] [dict get $network [lindex $pair 1] nname]]
      interaction<-paste(network[[pair[1]]]$nname,network[[pair[2]]]$nname)
      #             switch -- $interaction {
      #             {PHE PHE} {
      #               dict set interactType $pair PHE-PHE
      #             }
      if (interaction == "PHE PHE") {
        interactType[[pairindex]] <- "PHE-PHE"
        #             {AM1 AM2} -
        #               {AM2 AM1} -
        #               {AM1 AM1} -
        #               {AM2 AM2} {
        #                 dict set interactType $pair AM-AM
        #               }
      } else if (interaction %in% c("AM1 AM2", "AM2 AM1", "AM1 AM1", "AM2 AM2")) {
        interactType[[pairindex]] <- "AM-AM"
        #             {CYS CYS} {
        #               dict set interactType $pair DISU
        #             }
      } else if (interaction == "CYS CYS") {
        interactType[[pairindex]] <- "DISU"
        #             {AM1 PHE} -
        #               {PHE AM1} {
        #                 dict set interactType $pair AM1-PHE
        #               }
      } else if (interaction %in% c("AM1 PHE", "PHE AM1")) {
        interactType[[pairindex]] <- "AM1-PHE"
        #             {AM2 PHE} -
        #               {PHE AM2} {
        #                 dict set interactType $pair AM2-PHE
        #               }
      } else if (interaction %in% c("AM2 PHE", "PHE AM2")) {
        interactType[[pairindex]] <- "AM2-PHE"
        #             default {
        #               dict set interactType $pair STER
        #             }
      } else {
        interactType[[pairindex]] <- "STER"
      }
      #           }
      }
    }
    #       }
  }
  #     dict for {key val} $weight {
  #       puts $idf "$key $val"
  #     }
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
  #     dict for {key val} $interactType {
  #       puts $idf2 "$key $val"
  #     }
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
  #   }
  # }
  out[[2]]<-interactType_matrix
  } else {
    out[[2]]<-c()
  }
  #out[[3]]<-carbonContacts[[3]]
  #out[[3]]<-sulfurContacts[[3]]
  out
}

# Create a dictionary based on Eric's node list
# Create empty lists or dictionaries
network <- list()
indices <- list()
residlist<-list()
# set idf [open $myNodeList r]
con <- file(myNodeList)
lines<-readLines(con)
for (nindex in 1:length(lines)) {
  line <- unlist(strsplit(lines[[nindex]], " "))
  line <- lapply(line, function(x) {
    if (is.numeric(as.character(x))) {
      as.numeric(x)
    } else {
      x
    }
  })
  #   set nindex [lindex $line 0]
  #nindex <- line[[1]]
  #   dict set network $nindex resname [lindex $line 1]
  resname <- line[[2]]
  #   dict set network $nindex resid [lindex $line 2]
  resid <- line[[3]]
  #   dict set network $nindex moid [lindex $line 3]
  moid <- line[[4]]
  #   dict set network $nindex nname [lindex $line 4]
  nname <- line[[5]]
  #   dict set network $nindex ntype [lindex $line 5]
  ntype <- line[[6]]
  # Set values in the network dictionary
  network[[nindex]] <x- list(
    resname = resname,
    resid = resid,
    moid = moid,
    nname = nname,
    ntype = ntype
  )
  #   foreach thing [lrange $line 6 end] {
  # Iterate over the remaining items in the line
  for (i in 7:length(line)) {
    #     dict set indices $thing $nindex
    thing <- line[[i]]
    # Set values in the indices dictionary
    indices[[thing]] <- nindex
    #   }
  }
  # }
}
library(bio3d)
library(doParallel)
registerDoParallel(cores = detectCores())
registerDoParallel(cores = perbatch)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
mypdb<-read.pdb(pdbname) #load pdb file
mypdbatom<-mypdb$atom
# set nohs [atomselect $molwork "($mySelection) and noh"]
nohs<-atom.select(mypdb,"noh")
# set nocarbons [atomselect $molwork "($mySelection) and not (hydrogen or carbon or sulfur)"]
nocarbonselements <- c("N1", "N2", "N3", "N4", "O1", "O2", "O3", "O4")
nocarbons<-atom.select(mypdb,elety=nocarbonselements,resid=resname)
# set carbonAtoms [atomselect $molwork "($mySelection) and carbon"]
carbonAtomselements <- c("C1", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C2", "C20", "C3", "C4", "C5", "C6", "C7", "C8", "C9")
carbonAtoms<-atom.select(mypdb,elety=carbonAtomselements,resid=resname)
# set sulfurAtoms [atomselect $molwork "($mySelection) and sulfur"]
sulfurAtoms<-atom.select(mypdb,elety=c("S","S1", "S2"),resid=resname)
for (i in 1:totres) {
  sel<-atom.select(mypdb,resno=i)
  residlist[[i]]<-sel$atom
}
mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
nframes<-nrow(mydcd)
mycell<-read.dcd(dcdname,cell=T)
if (file.exists(outname)) {
  load(outname)
  offsetframe<-length(edges)
} else {
  offsetframe<-0
}
print(paste("offset", offsetframe))
if (test==TRUE) {
  batch<-1
  begframe<-(batch-1)*perbatch+1+offsetframe
  start_time <- Sys.time()
  print(paste("start",0,start_time))
  edgestemp<-getEdges(begframe)
  edges<-list()
  edges[[1]]<-edgestemp
  # Stop tracking time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  save(edges,file="edge_test.rda")
  print(paste("elapse",0, elapsed_time))
  begframe<-2
  for (i in 1:5) {
    registerDoParallel(cores = i)
    endframe<- begframe+i
    # Start tracking time
    start_time <- Sys.time()
    print(paste("start",i,start_time))
    #edgestemp<-foreach(frame=begframes:endframe) %dopar% (
    edgestemp<-foreach(frame=seq(begframe,begframe+i-1,1)) %dopar% (
      getEdges(frame)
    )
    begframe<-begframe+i
    # Stop tracking time
    end_time <- Sys.time()
    # Calculate the elapsed time
    elapsed_time <- end_time - start_time
    print(paste("elapse",i, elapsed_time))
    edges<-append(edges,edgestemp)
    save(edges,file="edge_test.rda")
  }
} else {
  for (batch in 1:ceiling((nframes-offsetframe)/perbatch)) {
    begframe<-(batch-1)*perbatch+1+offsetframe
    if (batch!=ceiling(nframes/perbatch)) {
      endframe<-batch*perbatch+offsetframe
    } else {
      endframe<-nframes
    }
    start_time <- Sys.time()
    print(paste("start",batch,start_time))
    edgestemp<-foreach(frame=begframes:endframe) %dopar% (
      getEdges(frame)
    )
    # Stop tracking time
    end_time <- Sys.time()
    # Calculate the elapsed time
    elapsed_time <- end_time - start_time
    if (!exists("edges")) {
      edges<-edgestemp
    } else {
      edges<-append(edges,edgestemp)
    }
    save(edges,file=outname)
    print(paste(begframe,endframe,length(edges),elapsed_time))
  }
}


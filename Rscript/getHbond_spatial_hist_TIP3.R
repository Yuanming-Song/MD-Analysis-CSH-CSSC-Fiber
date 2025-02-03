source("Rscript/package.R")
source("Rscript/base_fun.R")
source("Rscript/getHbond_spatialcutoff_perframe_TIP3.R")
outputdir<-"/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/cyl1/mix_27A/analysis/data/"
NewAna<-FALSE
midfix=""
if (FALSE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  outputpre<-"mix_27A_fibre_hbond_spati_dist_TIP3"
  uniseltextlist<-c("CSSC","CSH") #unique resname in sim
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (FALSE) {
    midfix="cylcssc."
    outputpre<-"mix_27A_cylcssc_fibre_hbond_spati_dist_TIP3"
  }
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  uniseltextlist<-c("CSSC") #unique resname in sim
  outputpre<-"pure_27A_fibre_hbond_spati_dist_TIP3"
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)
Oindex<-atom.select(mypdb,elety=c("O1","O2","O3","O4"),resid=c("CSH","CSSC"))$atom
Nindex<-atom.select(mypdb,elety=c("N1","N2","N3","N4"),resid=c("CSH","CSSC"))$atom
unique(Nindex-Oindex)
AvoidList<-list()
for (resname in uniseltextlist) {
  for (i in 1:get(paste0("max",resname))) {
    AvoidList[[resname]][[i]]<-atom.select(mypdb,resno=i,resid=resname)$atom
  }
}

dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-50
if (NewAna) {
  pattern <- paste0(outputpre, "_every", step, "_", firstrun, "to(\\d+)")   # Your existing file naming pattern
  files <- list.files(outputdir)   # Get a list of files in the directory
  matching_files <- grep(pattern, files, value = TRUE)   # Filter files based on the pattern
  extractnum<-c()
  for (file_name in matching_files) {
    match_result <- regmatches(file_name, regexpr(pattern, file_name)) #get name before .rda
    extracted_number <- as.numeric(gsub(pattern, "\\1", match_result)) # get the number
    extractnum<-c(extractnum,extracted_number)
  }
  extractnum<-extractnum[which(extractnum==max(extractnum))]
  outname<-paste0(outputdir,outputpre,"_every",step,"_",firstrun,"to",extractnum,".rda")
}
dcdlist<-getdcdlist(firstrun,midfix=midfix)
# Initialize an empty list to store pairs
uniqueResnamePair <- list()
#for tip3 water
TIP3Oindex<-atom.select(mypdb,elety=c("OH2"),resid="TIP3")$atom
TIP3Donorindex<-TIP3Oindex

donorlist <- list(c("OH2","H1", "H2"))
# Create an empty TIP3donorlist with 3 elements
TIP3donorlist1 <- vector("list", 1)
# Loop over each element of CSHdonorlist
for (n in 1:length(TIP3donorlist1)) {
  # Create an empty list to store the clusters
  TIP3donorlist1[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the TIP3donorlist1
  TIP3donorlist1[[n]][[1]] <- donorlist[[n]][1]
  # Add the rest of the elements of the input_list to the TIP3donorlist1
  TIP3donorlist1[[n]][[2]] <- donorlist[[n]][-1]
}
TIP3acceptorlist1<-c("OH2")

# Generate pairs
for (i in seq_along(uniseltextlist)) {
  for (j in seq_along(uniseltextlist)) {
    pair <- paste(sort(c(uniseltextlist[i], uniseltextlist[j])), collapse = "")
    uniqueResnamePair[[pair]] <- TRUE
  }
}

# Extract unique pairs from the list
uniqueResnamePair <- names(uniqueResnamePair)
donorlist1 <- list(c("N1", "H8", "H9"), c("N4", "H21", "H22"))
donorlist2<-list( c("N2", "H7"),c("N3", "H16"))
# Create an empty donorlist with 4 elements
CSSCdonorlist1 <- vector("list", 2)
# Loop over each element of CSSCdonorlist1
for (n in 1:length(CSSCdonorlist1)) {
  # Create an empty list to store the clusters
  CSSCdonorlist1[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSSCdonorlist1
  CSSCdonorlist1[[n]][[1]] <- donorlist1[[n]][1]
  # Add the rest of the elements of the input_list to the CSSCdonorlist1
  CSSCdonorlist1[[n]][[2]] <- donorlist1[[n]][-1]
}
CSSCdonorlist2 <- vector("list", 2)
# Loop over each element of CSSCdonorlist2
for (n in 1:length(CSSCdonorlist2)) {
  # Create an empty list to store the clusters
  CSSCdonorlist2[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSSCdonorlist2
  CSSCdonorlist2[[n]][[1]] <- donorlist2[[n]][1]
  # Add the rest of the elements of the input_list to the CSSCdonorlist2
  CSSCdonorlist2[[n]][[2]] <- donorlist2[[n]][-1]
}
CSSCacceptorlist1<-c("O1", "O4")
CSSCacceptorlist2<-c("O2", "O3")
# Input data
donotlist1 <- list(c("N1", "H9", "H10"))
# Create an empty CSHdonorlist with 4 elements
CSHdonorlist1 <- vector("list", 1)
# Loop over each element of CSHdonorlist
for (n in 1:length(CSHdonorlist1)) {
  # Create an empty list to store the clusters
  CSHdonorlist1[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSHdonorlist1
  CSHdonorlist1[[n]][[1]] <- donotlist1[[n]][1]
  # Add the rest of the elements of the input_list to the CSHdonorlist1
  CSHdonorlist1[[n]][[2]] <- donotlist1[[n]][-1]
}
donotlist2 <- list( c("N2", "H8"))
# Create an empty CSHdonorlist with 4 elements
CSHdonorlist2 <- vector("list", 1)
# Loop over each element of CSHdonorlist2
for (n in 1:length(CSHdonorlist2)) {
  # Create an empty list to store the clusters
  CSHdonorlist2[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSHdonorlist2
  CSHdonorlist2[[n]][[1]] <- donotlist2[[n]][1]
  # Add the rest of the elements of the input_list to the CSHdonorlist2
  CSHdonorlist2[[n]][[2]] <- donotlist2[[n]][-1]
}
CSHacceptorlist1<-c("O1")
CSHacceptorlist2<-c("O2")
#mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno <- with(mypdb$atom[which(mypdb$atom$resid=="CSSC"),], rep(1:ceiling(max(eleno)/26), each = 26)[1:length(eleno)])

processedframe<-0
for (resname in c(uniseltextlist,"TIP3")) {
  for (hbondrole in c("acceptor","donor")) {
    for (l in 1:2) {
      tempcorlistname<-paste0(resname,hbondrole,"list",l,"cor")
      templistname<-paste0(resname,hbondrole,"list",l)
      if (exists(templistname)) {
        templist<-list()
        for (i in 1:length(get(templistname))) {
          templist[[i]]<-list()
          for (j in 1:length(get(templistname)[[i]])) {
            templist[[i]][[j]]<-list()
            for (k in 1:length(get(templistname)[[i]][[j]])) {
              tempsel<-atom.select(mypdb,resid=resname,elety=get(templistname)[[i]][[j]][[k]])#select AM1 donor 
              if (resname=="CSSC") {
                templist[[i]][[j]][[k]]<-tempsel$atom
              } else {
                templist[[i]][[j]][[k]]<-tempsel$atom     
              }
            }
          }
        }
        assign(tempcorlistname,templist)
      }
    }
  }
}
fibresel<-atom.select(mypdb,resid=uniseltextlist) #select heavy atoms only
for (dcdname in dcdlist) {
  print(c(dcdname,format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) #print time
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big = TRUE) #read dcd file
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE) #read box size
  totframe<-nrow(mydcdcell) #how many frames
  print(c(totframe,format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) #print time again
  if (NewAna) {
    load(outname) #add new data to old data
    processedframe<-length(unique(RawHbondTot[,ncol(RawHbondTot)]))
    NewAna<-FALSE
  }
  if (length(seq(1,totframe,step))>1) { #parallel compute with multiple frames
    RawHbond<-foreach(frame=seq(1,totframe,step)) %dopar% (
      getHbond_spatialcutoff_perframe_TIP3(frame,getHist=TRUE)
    )
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) { #only one frame but not last file
      RawHbond<-list()
      RawHbond[[1]]<-getHbond_spatialcutoff_perframe_TIP3(1,getHist=TRUE)
    } else {
      break #ignore last file since it's still running 
    }
  }
  rm(mydcd) #clear space for memory
  rm(mydcdcell) #same
  if (exists("outname")) {
    load(outname) #add new data to old data
    file.remove(outname) #delete old data output
  } else {
    HbondDistHist<-list()
    for (pair in names(RawHbond[[1]])) {
      HbondDistHist[[pair]]<-list()
      for (AMname in names(RawHbond[[1]][[pair]])) {
        HbondDistHist[[pair]][[AMname]]<-list()
        for (coortype in c("Z","R","ZR")) {
          HbondDistHist[[pair]][[AMname]][[coortype]]<-0
        } 
      }
    }
    HbondDistHist[["Frame"]]<-0
    HbondDistHist[["ZBreaks"]]<-seq(-159.5, 159.5,1)
    HbondDistHist[["RBreaks"]]<-seq(0.5, 29.5,1)
    HbondDistHist[["dt"]]<-c()
  }
  for (pair in names(RawHbond[[1]])) {
    for (AMname in names(RawHbond[[1]][[pair]])) {
      for (coortype in names(RawHbond[[1]][[pair]][[AMname]])) {
        for (frame in 1:length(RawHbond)) {
          if (coortype=="Count") {
            HbondDistHist[["dt"]]<-rbind(HbondDistHist[["dt"]],c(RawHbond[[frame]][["Frame"]],sum(RawHbond[[frame]][[pair]][[AMname]][[coortype]]),pair,AMname))
          } else {
            HbondDistHist[[pair]][[AMname]][[coortype]]<-HbondDistHist[[pair]][[AMname]][[coortype]]+RawHbond[[frame]][[pair]][[AMname]][[coortype]]
            
          }
        }
      }
    }
  }
  HbondDistHist[["Frame"]]<-length(RawHbond)+HbondDistHist[["Frame"]]
  outname<-paste0(outputdir,outputpre,"_every",step,"_",firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(HbondDistHist,file = outname) #clear space 
  rm(HbondDistHist) #clear space
  processedframe<-processedframe+totframe
  print(c(dcdname,processedframe))
  rm(RawHbond)
}



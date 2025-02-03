source("Rscript/package.R")
source("Rscript/base_fun.R")
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
NewAna<-FALSE
midfix=""
if (TRUE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  outputpre<-"mix_27A_fibre_hbond_spati"
  uniseltextlist<-c("CSSC","CSH") #unique resname in sim
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (TRUE) {
    midfix="cylcssc."
    outputpre<-"mix_27A_cylcssc_fibre_hbond_spati"
  }
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  uniseltextlist<-c("CSSC") #unique resname in sim
  outputpre<-"pure_27A_fibre_hbond_spati"
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)

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
CSSCacceptorlist1<-c("O2", "O3")
CSSCacceptorlist2<-c("O1", "O4")
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
CSHacceptorlist1<-c("O2")
CSHacceptorlist2<-c("O1")
mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno <- with(mypdb$atom[which(mypdb$atom$resid=="CSSC"),], rep(1:ceiling(max(eleno)/26), each = 26)[1:length(eleno)])

processedframe<-0
for (resname in uniseltextlist) {
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
  print(c(dcdname,Sys.time())) #print time
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big = TRUE) #read dcd file
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE) #read box size
  totframe<-nrow(mydcdcell) #how many frames
  print(c(totframe,Sys.time())) #print time again
  if (NewAna) {
    load(outname) #add new data to old data
    processedframe<-length(unique(RawHbondTot[,ncol(RawHbondTot)]))
  }
  if (length(seq(1,totframe,step))>1) { #parallel compute with multiple frames
    RawHbond<-foreach(frame=seq(1,totframe,step),.combine=rbind) %dopar% (
      getHbond_spatialcutoff_perframe(frame)
    )
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) { #only one frame but not last file
      RawHbond<-getHbond_spatialcutoff_perframe(1)
    } else {
      break #ignore last file since it's still running 
    }
  }
  rm(mydcd) #clear space for memory
  rm(mydcdcell) #same
  if (exists("outname")) {
    load(outname) #add new data to old data
  }
  if (exists("RawHbondTot")) { #with old data present
    RawHbond<-rbind(RawHbondTot,RawHbond)
  } else { #first run
    RawHbondTot<-RawHbond
  }
  rm(RawHbond)
  if (exists("outname")) {
    file.remove(outname) #delete old data output
  }
  outname<-paste0(outputdir,outputpre,"_every",step,"_",firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(RawHbondTot,file = outname) #clear space 
  rm(RawHbondTot) #clear space
  processedframe<-processedframe+length(seq(1,totframe,step))
  print(c(dcdname,processedframe))
}



step=20
maxcssc=411
library(bio3d)
library(sna)
library(doParallel)
library(dplyr)
library(bigmemory)
corenum=detectCores()
print(corenum)
registerDoParallel(cores = corenum)
mypdb<-read.pdb("../setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb")
dcd_prefix <- "/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/nvt"
dcd_files <- c()
# Loop over files and check if they exist
for (i in 1:99) {
  if (i < 10) {
    dcd_name <- paste0(dcd_prefix, "0", i, ".dcd")
  } else {
    dcd_name <- paste0(dcd_prefix, i, ".dcd")
  }
  if (file.exists(dcd_name)) {
    dcd_files <- c(dcd_files, dcd_name)
  } else {
    break
  }
}
dcd_files <- c("data/mix_27A_sum1to8_every_20.dcd")
namelist<-c("C3", "S1", "S2", "C11")
getdihed <- function(j,totframe) {
  myresids<-as.numeric(strsplit(idlists[j],split=" ")[[1]])
  if (length(myresids)>0) {
    xyzlist<-c()
    for (k in 1:length(myresids)) {
      tempres<-myresids[k]
      for (l in 1:4) {
        sel<-atom.select(mypdb,resno = tempres,resid=c("CSSC"),elety=namelist[l],operator="AND")
        xyzlist<-c(xyzlist,mydcd[j,sel$xyz])
      }
    }
    torsion.xyz(xyzlist, atm.inc = 4)
    dihedtemp<-torsion.xyz(xyzlist, atm.inc = 4)
    dihedtemp<-as.data.frame(table(cut(dihedtemp,breaks =seq(-180,180,dihedbinsize))))
    dihedtemp<-dihedtemp %>% mutate(dihed=seq(-180,180-dihedbinsize,dihedbinsize))
    dihedtemp[,2]
  } else {
    c()
  }
}
dihedbinsize<-5
surfacelist<-c("surface","inside")
dihedidstout<-c()
for (dcd_file in dcd_files) {
  mydcd<-read.dcd(dcd_file)
  for (surfacetype in surfacelist) {
    #idlists<-file(paste0("~/Documents/Research/HPC/dfs2/mrsec/cg/analysis/csh_cssc_1to1_27mM_cg/",surfacetype,"id.dat"),"r")
    idlists<-file(paste0("data/dihed_",surfacetype,".ind"),"r")
    idlists<-readLines(idlists)
    totframe<-min(length(idlists),nrow(mydcd))
    print(totframe)
    print(length(idlists))
    print(nrow(mydcd))
    #totframe<-40
    dihed<-foreach (j = 1:totframe,.combine="+") %dopar% {
      getdihed(j,totframe)
    }
    dihedidst<-as.data.frame(cbind(as.numeric(dihed),seq(-180,180-dihedbinsize,dihedbinsize)))
    dihedidst<-dihedidst %>% mutate(normdihed=dihed/sum(dihed),type=surfacetype)
    dihedidstout<-rbind(dihedidst,dihedidstout)
  }
}

colnames(dihedidstout)<-c("Unnorm","dihed","Freq","Type")  
write.table(dihedidstout,file=paste0("data/dihed.dat"))
print("done")

if (0) {
  # Read in the data
  data <- read.table("~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/dihed.dat", header = TRUE)
  
  # Create the plot
  ggplot(data, aes(x = dihed, y = Freq, group = Type, color = Type)) +
    geom_line() +
    labs(x = "Dihed", y = "Freq", color = "Type")
} 







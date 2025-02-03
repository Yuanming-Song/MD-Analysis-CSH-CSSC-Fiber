sysname="mix_27A"
analysis="tip"
dir="data/"
steps=20
firstrun=1
lastrun=8
firstCSH=21373
moietyname <- c("N1", "C1", "O1", "C2", "C3", "S", "N2", "C4", "O2", "C5", "C6", "C7", "C8", "C9", "C10")
moietyname2 <- c("S2", "C11", "C12", "C13", "C14", "N3", "C15", "C16", "C17", "C18", "C19", "C20", "O3", "O4", "N4")
### first frame to be considered
startframe<-1
minz1<--110
maxz1<--90
minz2<-100
maxz2<-130
### for pdb & traj file names
pdbname=paste0(dir,sysname,".pdb")
#pdbname=paste(dir,sysname,"_resonly.pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
### how many residues in total
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
### analysis name for output name
analysis="_com_rdf"
### binsize and max for normalization
binsize<-0.5
binmax<-150
binmin<-0
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
get_counts <- function(tips, firstCSH) {
  # initialize an empty data frame to store the results
  df <- data.frame(frame = integer(), count = numeric(), type = character(), stringsAsFactors = FALSE)
  # loop through each frame
  for (i in seq_along(tips)) {
    # count the number of CSH and CSSC
    csh_count <- sum(firstCSH < tips[[i]]) # number of CSH in this frame
    cssc_count <- length(tips[[i]]) - csh_count # number of CSSC in this frame
    
    # calculate the percentage of CSH and CSSC
    total_count <- csh_count + cssc_count
    csh_percent <- csh_count / total_count 
    cssc_percent <- cssc_count / total_count 
    
    # add the results to the data frame
    df <- rbind(df, data.frame(frame = i, count = csh_percent, type = "CSH", stringsAsFactors = FALSE))
    df <- rbind(df, data.frame(frame = i, count = cssc_percent, type = "CSSC", stringsAsFactors = FALSE))
  }
  
  return(df)
}
#####   function for getting index
get_indices <- function(row_num, third_cols, minrange1, maxrange1, minrange2, maxrange2,within_range = FALSE) {
  # get the relevant columns for the given row
  row <- third_cols[row_num, ]
  # initialize an empty vector to store the results
  matches <- c()
  # loop through each column
  for (i in seq_along(row)) {
    val <- row[i]
    # check if the value is within either range
    if (val >= minrange1 & val <= maxrange1) {
      matches <- c(matches, i)
    } else if (val >= minrange2 & val <= maxrange2) {
      matches <- c(matches, i)
    }
  }
  if (within_range) {
    # only return indices within the ranges
    return(matches)
  } else {
    # return indices outside the ranges
    return(seq_along(row)[-matches])    
    
  }
}
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
#####   read dcd files
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
if (1) {
  #####	for testing...
  if (0) {
    nframes<-1
    startframe<-nframes
    print(nframes)
  }     
  zcoord <- mydcd[, seq(3, ncol(mydcd), 3)]
  tips<-foreach (framenum=startframe:nframes) %dopar% {
    get_indices(framenum, zcoord,minz1,maxz1,minz2,maxz2,TRUE)
  }
  nottips<-foreach (framenum=startframe:nframes) %dopar% {
    get_indices(framenum, zcoord, minz1,maxz1,minz2,maxz2)
  }
}
tipscount<- get_counts(tips, firstCSH)
nottipscount <- get_counts(nottips, firstCSH)
# Add "Tip" and "not Tip" classifiers to count tables
tipscount <- cbind(tipscount, "Tip")
nottipscount <- cbind(nottipscount, "not Tip")
colnames(tipscount)<- c("Frame", "CSH Count", "CSSC Count", "Tip/Not Tip")
colnames(nottipscount)<- c("Frame", "CSH Count", "CSSC Count", "Tip/Not Tip")

# Merge the two count tables
mergedcount <- rbind(tipscount, nottipscount)
# Write the merged count table into a file
write.table(mergedcount, file = "tip_compo.dat", sep = "\t", col.names = c("Frame", "CSH Count", "CSSC Count", "Tip/Not Tip"))
save(tips,file="tips_index.rda")
save(nottips,file="nottips_index.rda")
write_tips(tips, "tips_index.txt")

if (0) {
  # read the table
  data <- read.table("/Users/song/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/tip_compo.dat", header = TRUE)
  
  # create separate data frames for Tip and not Tip
  tip_data <- subset(data, data$`Tip.Not.Tip` == "Tip")
  nottip_data <- subset(data, data$`Tip.Not.Tip` == "not Tip")
  
  # create plot for Tip
  ggplot(tip_data, aes(x = Frame, y = `CSH.Count`, color = `CSSC.Count`)) +
    geom_point() +
    labs(title = "Tip", x = "Frame", y = "Compo")
  
  # create plot for not Tip
  ggplot(nottip_data, aes(x = Frame, y = `CSH.Count`, color = `CSSC.Count`)) +
    geom_point() +
    labs(title = "Not Tip", x = "Frame", y = "Compo")
  
}



# how many residues are there?
totCSSC<-411
totCSH<-352
totalres<-totCSSC*2+totCSH
step<-100
AcceptorDonorCutoff<-3.5
hbondAnglecutoff<-33
# Set the distance bin size max and min for distance
zbinsize <- 0.1
zbinmax <- 200
zbinmin <- -200
# Set the distance bin size max and min for distance
rbinsize <- 0.5
rbinmax <- 35
rbinmin <- 0
### for pdb & traj file names
dir<-"/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/"
pdbname=paste0(dir,"setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb")
# Which runs to consider
firstrun=3
lastrun=19
dcdnamelist<-c()
for (i in firstrun:lastrun) {
  dcdname<-paste0(dir,"nvt",sprintf("%02d", i),".dcd")
  if (file.exists(dcdname)) {
    dcdnamelist<-c(dcdnamelist,dcdname)
  } else {
    break
  }
}
#dcdnamelist<-"npt01.dcd"
# Calculate the number of distance bins needed 
nzbins <- ceiling((zbinmax - zbinmin) / zbinsize)
# Create an array of midpoint values for each distance 
zhistx <- seq(zbinmin + zbinsize / 2, zbinmax - zbinsize / 2, by = zbinsize)
# Create an array of distance bin 
zhisbreak <- c(zhistx - zbinsize / 2, zbinmax)
# Calculate the number of angle bins 
rbins <- ceiling((rbinmax - rbinmin) / rbinsize)
# Create an array of midpoint values 
rhistx <- seq(rbinmin + rbinsize / 2, rbinmax - rbinsize / 2, by = rbinsize)
# Create an array of r bin 
rhisbreak<-c(rhistx-rbinsize/2,rbinmax)
#moietyname for hbond
# Input data
input_list <- list(c("N1", "H8", "H9"), c("N2", "H7"), c("N3", "H16"), c("N4", "H21", "H22"))
# Create an empty donorlist with 4 elements
CSSCdonorlist <- vector("list", 4)
# Loop over each element of CSSCdonorlist
for (n in 1:length(CSSCdonorlist)) {
  # Create an empty list to store the clusters
  CSSCdonorlist[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSSCdonorlist
  CSSCdonorlist[[n]][[1]] <- input_list[[n]][1]
  # Add the rest of the elements of the input_list to the CSSCdonorlist
  CSSCdonorlist[[n]][[2]] <- input_list[[n]][-1]
}
CSSCacceptorlist<-c("O1", "O2", "O3", "O4","N1","N2","N3","N4")
# Input data
input_list <- list(c("N1", "H9", "H10"), c("N2", "H8"))
# Create an empty CSHdonorlist with 4 elements
CSHdonorlist <- vector("list", 2)
# Loop over each element of CSHdonorlist
for (n in 1:length(CSHdonorlist)) {
  # Create an empty list to store the clusters
  CSHdonorlist[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the CSHdonorlist
  CSHdonorlist[[n]][[1]] <- input_list[[n]][1]
  # Add the rest of the elements of the input_list to the CSHdonorlist
  CSHdonorlist[[n]][[2]] <- input_list[[n]][-1]
}
CSHacceptorlist<-c("O1", "O2","N1","N2")
input_list <- list(c("OH2", "H1", "H2"))
# Create an empty TIP3donorlist with 4 elements
TIP3donorlist <- vector("list", 1)
# Loop over each element of TIP3donorlist
for (n in 1:length(TIP3donorlist)) {
  # Create an empty list to store the clusters
  TIP3donorlist[[n]] <- vector("list", 2)
  # Add the first element of the input_list to the TIP3donorlist
  TIP3donorlist[[n]][[1]] <- input_list[[n]][1]
  # Add the rest of the elements of the input_list to the TIP3donorlist
  TIP3donorlist[[n]][[2]] <- input_list[[n]][-1]
}
TIP3acceptorlist<-c("OH2")
library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)
registerDoParallel(cores = detectCores())
# calculate unit normal vector 
getnorm <- function(v1, v2, v3) {
  va <- v2 - v1   # Calculate difference vector between v1 and v2
  vb <- v3 - v1   # Calculate difference vector between v1 and v3
  vc <- cross(va, vb)  # Calculate cross product of va and vb
  vc / sqrt(sum(vc^2))  # Normalize vc and return it
}
# create 2D histogram table
cut2Dhis <- function(coord, breakseq1, breakseq2) {
  # Use the cut() function to bin each variable
  var1_binned <- cut(coord[,1], breaks=breakseq1, include.lowest=TRUE)
  var2_binned <- cut(coord[,2], breaks=breakseq2, include.lowest=TRUE)
  # Use the table() function to create the 2D histogram table
  hist_table <- table(var1_binned, var2_binned)
  # Print the table
  hist_table
}
# calculate the angle between two vectors
getcostheta <- function(ref, comp) {
  sum(ref * comp) / (sqrt(sum(ref^2)) * sqrt(sum(comp^2)))  # Calculate dot product and normalize the vectors to obtain the cosine of the angle
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
#for given reference moiety, create reference info, get raw count value for distance and angle
Get2DHisPerRefHbondCoord<-function(refresid,refresname,framenum,compresname) {
  counter<-1 #counter for which line in matrix to update
  acceptorlist<-get(paste0(refresname,"acceptorlist"))
  donorlist<-get(paste0(compresname,"donorlist"))
  fincomp<-get(paste0("tot",compresname))
  if (refresname==compresname) {
    inicomp<-refresid+1
  } else {
    inicomp<-1
  }
  rawdistm<-matrix(NA,nrow =length(acceptorlist)*length(donorlist)*( totalres),ncol =3 ) #set up a matrix for cut2Dhis function 
  for (oxygen in acceptorlist) {
    hbondrefzsel<-atom.select(mypdb,elety = oxygen,resno =refresid,resid=refresname) # Select atoms belonging to the given moiety
    hbondrefcom<-mydcd[framenum,hbondrefzsel$xyz]   # Obtain the x, y, and z coordinates of the selected atoms
    breakOuterLoop <- FALSE
    for (compresid in sample(inicomp:fincomp) ) { 
      for (donorcount in sample(1:length(donorlist))) {
        hbondNsel<-atom.select(mypdb,elety = donorlist[[donorcount]][[1]],resno =compresid,resid=compresname) # Select atoms belonging to the given moiety
        hbondNxyz<-mydcd[framenum,hbondNsel$xyz] # Get Nitrogen coordinate
        hbondNxyz<-pbc(hbondrefcom,hbondNxyz,mycell[framenum,]) # Correct for PBC
        ONdis<-euclidean(hbondrefcom,hbondNxyz)       # Get NO distance
        if (ONdis<AcceptorDonorCutoff) {
          OHdis<-mycell[framenum,1]*2 # Set an impossible OH distance 
          hydrogenlist<-donorlist[[donorcount]][[2]]
          for (hydrogen in hydrogenlist) {
            hbondOsel<-atom.select(mypdb,elety = hydrogen,resno =compresid,resid=compresname) # Select Oxygen
            hbondOxyz<-mydcd[framenum,hbondOsel$xyz] # Get Oxygen coordinate
            hbondOxyz<-pbc(hbondrefcom,hbondOxyz,mycell[framenum,]) # Correct for PBC
            tepOHdis<-euclidean(hbondrefcom,hbondOxyz) # Get OH distance 
            if (tepOHdis<OHdis) {
              OHdis<-tepOHdis # Compare to set minimal OH distance
              hbondOxyzclose<-hbondOxyz # Save O coord for later
            }
          }
          # Get ONH angle
          hbondAngle<-getcostheta(hbondrefcom-hbondNxyz,hbondOxyzclose-hbondNxyz)
          hbondAngle<-acos(hbondAngle)*180/3.14
          if (hbondAngle<hbondAnglecutoff) {
            breakOuterLoop <- TRUE
            rawdistm<-(hbondNxyz+hbondrefcom)/2
            counter<-1+counter
            break
          }
        }
      }
      if(breakOuterLoop){
        break
      }
    }
  }
  rawdistm<-na.omit(rawdistm)  
  transformedMatrix <- matrix(NA, nrow = nrow(rawdistm), ncol = 2)  # Initialize a matrix to store the transformed data
  transformedMatrix[, 1] <- sqrt(rowSums(rawdistm[, 1:2]^2))  # Compute the square root of the sum of squares of the first two columns of rawdistm and store it in the first column of the new matrix
  transformedMatrix[, 2] <- rawdistm[, 3]  # Copy the third column of rawdistm to the second column of the new matrix
  transformedMatrix
}
Get2DHisPerRefHbondCoordMatrix<-function(refresid,refresname,framenum,compresname) {
  counter<-1 #counter for which line in matrix to update
  acceptorlist<-get(paste0(refresname,"acceptorlist"))
  donorlist<-get(paste0(compresname,"donorlist"))
  rawdistm<-matrix(NA,nrow =length(acceptorlist)*length(donorlist)*(totalres),ncol =3 ) #set up a matrix for cut2Dhis function 
  for (Acceptor in acceptorlist) {
    hbondAcceptorsel<-atom.select(mypdb,elety = Acceptor,resno =refresid,resid=refresname) # Select atoms belonging to the given moiety
    hbondAcceptorXyz<-mydcd[framenum,hbondAcceptorsel$xyz]   # Obtain the x, y, and z coordinates of the selected atoms
    for (donorcount in sample(1:length(donorlist))) {
      hbondNsel<-atom.select(mypdb,elety = donorlist[[donorcount]][[1]],resid=compresname) # Select atoms belonging to the given moiety
      coorddonor<-mydcd[framenum,hbondNsel$xyz] # Get Donor coordinate
      hydrogenlist<-donorlist[[donorcount]][[2]]
      # Calculate the number of atoms
      n_atomsDonor <- length(coorddonor)/3
      # Loop through each atom
      for (i in sample(1:n_atomsDonor)) {
        atomindex<-(3*i-2):(3*i)
        # Get the current atom's coordinates
        hbondDonorXyz <-  pbc(hbondAcceptorXyz,coorddonor[atomindex],mycell[framenum,])
        if (euclidean(hbondAcceptorXyz,hbondDonorXyz)<AcceptorDonorCutoff) {
          AcceptorHdis<-mycell[framenum,1]*2 # Set an impossible acceptor-H distance 
          for (hydrogen in hydrogenlist) {
            hbondHSel<-atom.select(mypdb,elety = hydrogen,resid=compresname) # Select Oxygen
            hbondHXyz<-mydcd[framenum,hbondHSel$xyz][atomindex] # Get Oxygen coordinate
            hbondHXyz<-pbc(hbondAcceptorXyz,hbondHXyz,mycell[framenum,]) # Correct for PBC
            tepdis<-euclidean(hbondAcceptorXyz,hbondHXyz) # Get OH distance 
            if (tepdis<AcceptorHdis) {
              AcceptorHdis<-tepdis # Compare to set minimal OH distance
              hbondHXyzClose<-hbondHXyz # Save O coord for later
            }
          }
          hbondAngle<-getcostheta(hbondAcceptorXyz-hbondDonorXyz,hbondHXyzClose-hbondDonorXyz)
          hbondAngle<-acos(hbondAngle)*180/3.14
          if (hbondAngle<hbondAnglecutoff) {
            breakOuterLoop <- TRUE
            rawdistm[counter,] <- (hbondDonorXyz+hbondAcceptorXyz)/2
            counter<-1+counter
            break
          }
        }
      }
    }
  }
  if (compresname=="TIP3") {
    acceptorlist<-get(paste0(compresname,"acceptorlist"))
    donorlist<-get(paste0(refresname,"donorlist"))
    for (donorcount in sample(1:length(donorlist))) {
      hbondDonorSel<-atom.select(mypdb,elety = donorlist[[donorcount]][[1]],resid=refresname,resno =refresid) 
      hbondDonorXyz<-mydcd[framenum,hbondDonorSel$xyz]
      hydrogenlist<-donorlist[[donorcount]][[2]]
    for (Acceptor in acceptorlist) {
      hbondAcceptorsel<-atom.select(mypdb,elety = Acceptor,resid=compresname) # Select atoms belonging to the given moiety
        coordAcceptor<-mydcd[framenum,hbondAcceptorsel$xyz] # Get Donor coordinate
        # Calculate the number of atoms
        n_atomsAcceptor <- length(coordAcceptor)/3
        # Loop through each atom
        for (i in sample(1:n_atomsAcceptor)) {
          atomindex<-(3*i-2):(3*i)
          # Get the current atom's coordinates
          hbondAcceptorXyz<-  pbc(hbondDonorXyz,coordAcceptor[atomindex],mycell[framenum,])
          if (euclidean(hbondAcceptorXyz,hbondDonorXyz)<AcceptorDonorCutoff) {
            AcceptorHdis<-mycell[framenum,1]*2 # Set an impossible acceptor-H distance 
            for (hydrogen in hydrogenlist) {
              hbondHSel<-atom.select(mypdb,elety = hydrogen,resid=refresname,resno =refresid) # Select Oxygen
              hbondHXyz<-mydcd[framenum,hbondHSel$xyz]# Get Oxygen coordinate
              hbondHXyz<-pbc(hbondDonorXyz,hbondHXyz,mycell[framenum,]) # Correct for PBC
              tepdis<-euclidean(hbondAcceptorXyz,hbondHXyz) # Get Acceptor H distance 
              if (tepdis<AcceptorHdis) {
                AcceptorHdis<-tepdis # Compare to set minimal OH distance
                hbondHXyzClose<-hbondHXyz # Save O coord for later
              }
            }
            hbondAngle<-getcostheta(hbondAcceptorXyz-hbondDonorXyz,hbondHXyzClose-hbondDonorXyz)
            hbondAngle<-acos(hbondAngle)*180/3.14
            if (hbondAngle<hbondAnglecutoff) {
              breakOuterLoop <- TRUE
              rawdistm[counter,]<-(hbondDonorXyz+hbondAcceptorXyz)/2
              counter<-1+counter
              break
            }
          }
        }
      }
    }
    
    
  }
  rawdistm <- na.omit(rawdistm)
  transformedMatrix <- matrix(NA, nrow = nrow(rawdistm), ncol = 2)
  transformedMatrix[, 1] <- sqrt(rawdistm[, 1]^2 + rawdistm[, 2]^2)
  transformedMatrix[, 2] <- rawdistm[, 3]
  transformedMatrix
}
#for given frame, get raw count value for distance and angle
merge2refeachframe<-function(framenum,refresname,compresname) {
  temphis<-c()
  if (refresname==compresname) {
    resdiff<-1
  } else {
    resdiff<-0
  }
  fini<-get(paste0("tot",refresname))-resdiff
  # loop through resid then update 2d historgram
  for (i in 1:fini){
    temphis<-rbind(temphis,Get2DHisPerRefHbondCoordMatrix(i,refresname,framenum,compresname))
  }
  #temphis
  cut2Dhis(temphis,rhisbreak,zhisbreak)
}
tothistgramcshcssc<-0 #initial empty histogram
tothistgramcsh<-0 #initial empty histogram
tothistgramcssc<-0 #initial empty histogram
tothistgramCSHTIP3<-0 #initial empty histogram
tothistgramCSSCTIP3<-0 #initial empty histogram
mypdb<-read.pdb(pdbname,hex=TRUE) #load pdb file
dcdcounter<-firstrun #for tracking 
totalframes<-0 #for tracking 
outlist<-c("CSHCSSC","CSSC","CSH")
outlist<-c("CSHTIP3","CSSCTIP3")

for (dcdname in dcdnamelist) {
  #####   read dcd files and get cell size
  mydcd<-read.dcd(dcdname,verbose=FALSE)
  nframes<-nrow(mydcd)
  mycell<-read.dcd(dcdname,cell=T,verbose=FALSE)
  totalframes<-length(seq(1,nframes,step))+totalframes
  print(paste(dcdname,totalframes))
	if (0){
    #in parallel, get 2d histogram for each frame then add up
  histogramperdcd<-foreach (framenum = 1:nframes, .combine = "+") %dopar% {
    merge2refeachframe(framenum,"CSSC","CSSC")
  }
  #update total histogram
  tothistgramcssc<-tothistgramcssc+histogramperdcd
  #in parallel, get 2d histogram for each frame then add up
  histogramperdcd<-foreach (framenum = 1:nframes, .combine = "+") %dopar% {
    merge2refeachframe(framenum,"CSH","CSH")
  }
  #update total histogram
  tothistgramcsh<-tothistgramcsh+histogramperdcd
  #in parallel, get 2d histogram for each frame then add up
  histogramperdcd<-foreach (framenum = 1:nframes, .combine = "+") %dopar% {
    merge2refeachframe(framenum,"CSSC","CSH")
  }
  #update total histogram
  tothistgramcshcssc<-tothistgramcshcssc+histogramperdcd
  #write current total histogram
  }
  #in parallel, get 2d histogram for each frame then add up
  histogramperdcd<-foreach (framenum = seq(1,nframes,step), .combine = "+") %dopar% {
    merge2refeachframe(framenum,"CSSC","TIP3")
  }
  tothistgramCSSCTIP3<-tothistgramCSSCTIP3+histogramperdcd
  histogramperdcd<-foreach (framenum = seq(1,nframes,step), .combine = "+") %dopar% {
    merge2refeachframe(framenum,"CSH","TIP3")
  }
  tothistgramCSHTIP3<-tothistgramCSHTIP3+histogramperdcd
  tothistgramout<-c()
  #convert it to a 3 column table for easy plotting 
  for (sufix in outlist) {
    dataname<-paste0("tothistgram",sufix)
    type<-sufix
    for (j in 1:nrow(get(dataname))) {
      #each row is for a unique distance, so get first column is always same, and 2nd & 3rd columns are angle bin value and raw count
      distance<-rhistx[j]
      tothistgramout<-rbind(tothistgramout,cbind(distance,zhistx,as.numeric(get(dataname)[j,]),type))
    }  
    outfile2dname<-paste0("outfile2d",sufix)
    if (dcdcounter>firstrun) {
      file.remove(get(outfile2dname))
    }
    assign(outfile2dname,paste0("data/",firstrun,"upto_npt",dcdcounter,"_coordis_hbond_",sufix,".dat"))
    write.table(get(dataname), file=get(outfile2dname),quote = F)

  }
  #add column name
  tothistgramout<-as.data.frame(tothistgramout)
  colnames(tothistgramout)<-c("r","z","Freq","Type")
  #delete and write file
  if (dcdcounter>firstrun) {
    file.remove(outfile)
  }
  outfile<-paste0("data/",firstrun,"upto_npt",dcdcounter,"_coordis_hbond", ".dat")
  write.table(tothistgramout, file=outfile,quote = F)
  # print dcdname, current total frames, and how many pairs were analyzed for current dcd file in each frame
  print(paste(dcdname,totalframes,sum(histogramperdcd)/nframes))
  #update trackers
  dcdcounter<-dcdcounter+1
}
#print final tracker info
print(totalframes)
#for plotting on local R
if (0) {  # Set the directory to the current directory
  dircylmix <- "~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/data/"
  # Loop over values of i from 5 to 20
  for (i in 108:1) {
    # Construct the name of the output file for this value of i
    outfile <- paste0("3upto_npt", i, "_coordis_hbond", ".dat")
    # Check if the output file exists in the current directory
    if (outfile %in% list.files(dircylmix)) {
      # Print the name of the output file
      outfile<-paste0(dircylmix,outfile)
      outfile2d <- paste0(dircylmix,"88upto_npt", i, "_rdf_hbond_2Dhis", ".dat")
      break
    }
  }
  #read table
  hbond2dhiscylmix<-read.table(file=outfile,header=TRUE,row.names = 1)
  hbond2dhiscylmixz <- aggregate(Freq ~ z + Type, data = hbond2dhiscylmix, sum)
  hbond2dhiscylmixznorm <- hbond2dhiscylmixz %>%
    group_by(Type) %>%
    mutate(Freq_norm = Freq/((50)*(27.5^2)*3.1415*zbinsize)) %>%
    ungroup()
  abshbond2dhiscylmixnorm <- hbond2dhiscylmixznorm %>% 
    mutate(abs_z = abs(z)) %>% 
    group_by(abs_z, Type)   
  abshbond2dhiscylmixnorm <- aggregate(Freq_norm ~ abs_z + Type, data = abshbond2dhiscylmixnorm, mean)
  cylzcoordistplt<-{
    #ggplot(hbond2dhiscylmixznorm, aes(x = z, y = Freq_norm, col = Type)) +
    ggplot(abshbond2dhiscylmixnorm, aes(x = abs_z, y = Freq_norm, col = Type)) +
      geom_line() +
      ggtitle("All noh")+
      #xlim(0, 160) +
      ylim(0,0.007)+
      xlab("z (Å)") +
      # geom_vline(xintercept = c(20, 85), linetype = "dashed",col="dark green")+
      # geom_rect(xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, fill = "antiquewhite", alpha = 0.02,color=NA) +
      # geom_rect(xmin = 85, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "cyan", alpha = 0.002,color=NA) +
      #ylab(paste0("Number density (",expression(Å^3),")"))+
      ylab(expression(paste("Number density (", Å^3, ")")))+
      #scale_color_manual(values = c("CSSC" = "blue", "CSH" = "red", "Water" = "black")) +
      #scale_x_continuous(breaks = seq(0,160,20),limits = c(0,160))+
      theme(panel.background = element_rect(fill='transparent'),
            text=element_text(face = "bold",size=14),
            legend.position = c(0.5,0.9),
            legend.direction = "horizontal",
            legend.title = element_blank(),
            panel.grid.major=element_line(colour = "grey",size = 0.1),
            axis.line = element_line(colour = "black",size = 0.3)
      )    
  }
  # ggsave("~/Documents/Research/MRSEC/abscylzcoordist.png",
  #        plot=cylzcoordistplt, dpi=300, 
  #        width=5.5*1.5, height=5,units="in")
  hbond2dhiscylmixr <- aggregate(Freq ~ r + Type, data = hbond2dhiscylmix, sum)
  hbond2dhiscylmixrnorm <- hbond2dhiscylmixr %>%
    group_by(Type) %>%
    mutate(Freq_norm = Freq/(((r+rbinsize)^2-(r)^2)*3.1415 * (50) * 215)) %>%
    ungroup()
  cylrcoordistplt<-{
    ggplot(hbond2dhiscylmixrnorm, aes(x = r, y = Freq_norm, col = Type)) +
      geom_line() +
      ggtitle("All noh")+
      #xlim(0, 160) +
      ylim(0,0.007)+
      xlab("r (Å)") +
      # geom_vline(xintercept = c(20, 85), linetype = "dashed",col="dark green")+
      # geom_rect(xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, fill = "antiquewhite", alpha = 0.02,color=NA) +
      # geom_rect(xmin = 85, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "cyan", alpha = 0.002,color=NA) +
      #ylab(paste0("Number density (",expression(Å^3),")"))+
      ylab(expression(paste("Number density (", Å^3, ")")))+
      #scale_color_manual(values = c("CSSC" = "blue", "CSH" = "red", "Water" = "black")) +
      #scale_x_continuous(breaks = seq(0,160,20),limits = c(0,160))+
      theme(panel.background = element_rect(fill='transparent'),
            text=element_text(face = "bold",size=14),
            legend.position = c(0.5,0.9),
            legend.direction = "horizontal",
            legend.title = element_blank(),
            panel.grid.major=element_line(colour = "grey",size = 0.1),
            axis.line = element_line(colour = "black",size = 0.3)
      )    
  }
  
  #box size / (sphere volume factor) * (total pairs x total frames)
  outlist<-unique(hbond2dhiscylmix$Type)
  hbondrdf1to1<-c()
  hbond2dhiscylmixnorm<-c()
  for (Type in outlist) {
    normfactorname<-paste0(Type,"normfactor")
    tempdata<-hbond2dhiscylmix[which(hbond2dhiscylmix$Type==Type),]
    print(paste(Type,sum(tempdata[which(tempdata$D<3.5 & tempdata$Angle<33),]$Freq)/(24602*16)))
    assign(normfactorname,100^3/((4*3.14/3)*sum(tempdata$Freq)))
    temprdf<-cbind(tempdata$D,tempdata$Freq)
    hbondrdftempsum <- aggregate(V2 ~ V1, data =  temprdf, sum)
    hbondrdftempsum <- hbondrdftempsum %>% mutate(g=get(normfactorname)*V2/(((V1+zbinsize/2)^3-(V1-zbinsize/2)^3)),Type=Type)
    hbondrdf1to1<-rbind(hbondrdftempsum,hbondrdf1to1)
    tempdata <- tempdata %>% mutate(g=get(normfactorname)*Freq/(((D+zbinsize/2)^3-(D-zbinsize/2)^3)))
    # print({
    #   ggplot(tempdata[which(tempdata$D<5),], aes(x = Angle, y = D, fill = Freq/(sum(Freq)))) +
    #     geom_tile() +
    #     scale_fill_viridis_c(
    #       #limits = c(0, max(phe2dhiscssc$NormFreq_D))
    #     )+
    #     scale_x_continuous(breaks = seq(0,50,2),limits = c(0,50))+
    #     scale_y_continuous(breaks = seq(0,20,0.2),limits = c(2.6,5))+
    #     labs(x = "θ", y = paste(Type,"PHE COM r"))
    # })
  }
  hbondrdfplt<-{
    ggplot() +
      geom_line(data = hbondrdf1to1, aes(x=V1, y=g, col=Type)) +
      #scale_color_manual(values=c("blue", "red"))+
      scale_x_continuous(breaks = seq(0,30,5),limits = c(0,50))+
      scale_y_continuous(breaks = seq(0,160,20),limits = c(0,160))+
      xlab("r") +
      ylab("N-O g(r)")
  }
  normfactor1to1<-100^3/((4*3.14/3)*sum(hbond2dhiscylmix$Freq))
  hbondrdf1to1sum <- aggregate(hbond2dhiscylmix$Freq ~ hbond2dhiscylmix$D, data =  hbond2dhiscylmix, sum)
  colnames(hbondrdf1to1sum)<- c("r","g")
  hbondrdf1to1sum <- hbondrdf1to1sum %>% mutate(g=normfactor1to1*g/(((r+zbinsize/2)^3-(r-zbinsize/2)^3)))
  hbondrdfplt<-{
    ggplot() +
      geom_line(data = hbondrdfcsscsum, aes(x=V1, y=V2*csscnormfactor/(((V1+zbinsize/2)^3-(V1-zbinsize/2)^3)), color="CSSC")) +
      geom_line(data = hbondrdfcshsum, aes(x=V1, y=V2*cshnormfactor/(((V1+zbinsize/2)^3-(V1-zbinsize/2)^3)), color="CSH")) +
      geom_line(data = hbondrdf1to1sum, aes(x=r, y=g,color="50%")) +
      xlab("r") +
      ylab("N-O g(r)") +
      scale_color_manual(values=c("black","red", "blue"))+
      scale_x_continuous(breaks = seq(0,30,5),limits = c(0,50))+
      scale_y_continuous(breaks = seq(0,160,20),limits = c(0,160)
      )
  }
  
  hbond2dhiscsscm<-read.table(file=outfile2d,header=TRUE,row.names = 1)
  for (i in 1:nrow(hbond2dhiscsscm)) {
    hbond2dhiscsscm[i,]<-hbond2dhiscsscm[i,] * csscnormfactor / ((dishistx[i] + zbinsize/2)^3 - (dishistx[i] - zbinsize/2)^3)
  }
  hbond3dplt<-{plot_ly(y = dishistx, x = anglehistx, z = as.matrix(hbond2dhiscsscm))%>% 
    add_surface(colorscale = 'inferno', cmin = 1, cmax = max(hbond2dhiscsscm)) %>%
    layout(scene = list(
      yaxis = list(range = c(0, 20),title="hbond COM difference (Å)"),
      xaxis = list(title = "cos(angle between normal vector)"),
      zaxis = list(title = "g(r)"#,range = c(0, 15)
                   )
    )
    )
  }
  sum(hbond2dhiscsscm[which(hbond2dhiscsscm$D<7),]$Freq)/sum(hbond2dhiscsscm[which(hbond2dhiscsscm$D>7),]$Freq)
  sum(hbond2dhiscsscm$Freq)/2025
  sum(hbond2dhiscylmix[which(hbond2dhiscylmix$D<3.5 & hbond2dhiscylmix$Angle<33),]$Freq)/(24602*32)
  
}


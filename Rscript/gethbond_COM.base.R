### how many residues in total
mincssc<-1
maxcssc<-411 #411
maxcsh<-763
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)
library(dplyr)
perbatch<-detectCores()
mincsh<-maxcssc+1
totcssc<-maxcssc
totcsh<-maxcsh-maxcssc
donorlist1 <- list(c("N1", "H8", "H9"), c("N4", "H21", "H22"))
donorlist2<-list(c("N3", "H16"), c("N2", "H7"))
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
amtypelist<-c("AM1","AM2")
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
registerDoParallel(cores = detectCores())
# calculate unit normal vector 
getnorm <- function(v1, v2, v3) {
  va <- v2 - v1   # Calculate difference vector between v1 and v2
  vb <- v3 - v1   # Calculate difference vector between v1 and v3
  vc <- cross(va, vb)  # Calculate cross product of va and vb
  vc / sqrt(sum(vc^2))  # Normalize vc and return it
}
# create 2D histogram table
cuthis <- function(coord, breakseq1, breakseq2) {
  # Use the cut() function to bin each variable
  var1_binned <- cut(as.numeric(coord[,1]), breaks=breakseq1, include.lowest=TRUE)
  var2_binned <- cut(as.numeric(coord[,2]), breaks=breakseq2, include.lowest=TRUE)
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
pbc <- function(comi, cur, currcell) {
  # Loop through x, y, and z coordinates
  for (l in 1:3) {
    diff <- cur[l] - comi[l]  # Calculate the difference between the current coordinate and the center of mass coordinate
    absdiff <- abs(diff)  # Calculate the absolute value of the difference
    currbox<-currcell[l]
    # Update coordinate only when necessary to account for periodic boundary conditions
    if (absdiff > currbox/2) {
      cur[l] <- cur[l] - (diff/absdiff) * currbox * ceiling((absdiff- currbox/2) / currbox)
    } 
  }
  # Return updated current coordinates
  cur
}
#base function file
#get hbond function
gethbond<-function(temptempnodeA,temptempnodeB,temptempacceptorlist,temptempdonorlist,framenum) {
  output<-c()
  for (acceptor in temptempacceptorlist) {
    hbondAzsel<-atom.select(mypdb,elety = acceptor,resno =network[[temptempnodeA]]$resid,resid=network[[temptempnodeA]]$resname) # Select atoms belonging to the given moiety
    hbondAcom<-mydcd[framenum,hbondAzsel$xyz]   # Obtain the x, y, and z coordinates of the selected atoms
    hydrogenlist<-temptempdonorlist[[2]]
    OHdis<-mycell[framenum,1]*2 # Set an impossible OH distance 
    for (hydrogen in hydrogenlist) {
      hbondHsel<-atom.select(mypdb,elety = hydrogen,resno =network[[temptempnodeB]]$resid,resid=network[[temptempnodeB]]$resname) # Select acceptor
      hbondHxyz<-mydcd[framenum,hbondHsel$xyz] # Get acceptor coordinate
      hbondHxyz<-pbc(hbondAcom,hbondHxyz,mycell[framenum,]) # Correct for PBC
      tepOHdis<-euclidean(hbondAcom,hbondHxyz) # Get OH distance 
      if (tepOHdis<OHdis) {
        OHdis<-tepOHdis # Compare to set minimal OH distance
        hbondHxyzclose<-hbondHxyz # Save O coord for later
      }
    }
    hbondDsel<-atom.select(mypdb,elety = temptempdonorlist[[1]],resno =network[[temptempnodeB]]$resid,resid=network[[temptempnodeB]]$resname) # Select atoms belonging to the given moiety
    hbondDxyz<-mydcd[framenum,hbondDsel$xyz] # Get Nitrogen coordinate
    hbondDxyz<-pbc(hbondAcom,hbondDxyz,mycell[framenum,]) # Correct for PBC
    ADdis<-euclidean(hbondAcom,hbondDxyz)       # Get NO distance
    # Get OHN angle
    AHDang<-getcostheta(hbondAcom-hbondHxyzclose,hbondDxyz-hbondHxyzclose)
    acosAHDang<-acos(AHDang)*180/3.14
    output<-rbind(output,c(ADdis,AHDang,acosAHDang,pbc(c(0,0,0),(hbondDxyz+hbondAcom)/2,mycell[framenum,]),hbondAcom-hbondDxyz)) # Update matrix (distance, angle)
  }
  output
}
getrawhbond<-function(frame) {
  #distll am-am edge out by checking which node in the network is AM type
  tempnode<-edges[[frame]][[1]]
  amedgelist<-tempnode[which((tempnode[,1]+2)%%3==0 & (tempnode[,2]+2)%%3==0 & tempnode[,2]>tempnode[,1]),]
  #AM1-AM1 AM2-AM2 AM1-AM2
  rawhbond<-c()
  for (tempnodeline in 1:nrow(amedgelist)) {
    tempnodeA<-amedgelist[tempnodeline,1]
    tempnodeB<-amedgelist[tempnodeline,2]
    for (nodecount in c("A","B")) {
      nodecountind<-ifelse(nodecount=="A",1,2)
      tempnode<-amedgelist[tempnodeline,nodecountind]
      assign(paste0("tempnode",nodecount),tempnode)
      if (tempnode<=maxcssc*6){
        tempprefix<-"CSSC"
        if (((tempnode+2)/3)%%2) {
          tempindex<-2
        } else {
          tempindex<-1
        }
      } else {
        tempprefix<-"CSH"
        tempindex<-1
      }
      for (number in 1:2) {  # Assuming the range of numbers
        for (temptype in c("acceptor","donor")) {
          templist<-get(paste0(tempprefix,temptype,"list",number))
          assign(paste0(temptype,"list",nodecount, number), templist[[tempindex]])
          # print(paste0(tempprefix,temptype,"list",number))
          # print(paste0(temptype,"list",nodecount, number))
        }
      }
    }
    for (i in c("A","B")) {
      i1<-i
      i2<-ifelse(i == "A", "B", "A")
      for (j in 1:2) {
        for (k in 1:2) {
          temprawhbond<-gethbond(get(paste0("tempnode", i1)),
                                 get(paste0("tempnode", i2)),
                                 get(paste0("acceptorlist",i1,j)),
                                 get(paste0("donorlist",i2,k)),
                                 frame
          )
          amtype<-paste(amtypelist[j],amtypelist[k],sep="-")
          temprawhbond<-c(temprawhbond,amtype)
          rawhbond<-rbind(rawhbond,temprawhbond)
        }
      }
    }
  }
  row.names(rawhbond)<-NULL
  rawhbond
}
seinitbin<-function() {
  # Create an array of midpoint values for each distance 
  dishistx <<- seq(disbinmin + disbinsize / 2, disbinmax - disbinsize / 2, by = disbinsize)
  # Create an array of distance bin 
  dishisbreak <<- c(dishistx - disbinsize / 2, disbinmax)
  # Create an array of midpoint values 
  anglehistx <<- seq(anglebinmin + anglebinsize / 2, anglebinmax - anglebinsize / 2, by = anglebinsize)
  # Create an array of angle bin 
  anglehisbreak<<-c(anglehistx-anglebinsize/2,anglebinmax)
}
seinitbin()
get2dhis<-function(framenum) {
  rawhbondtemp <- totrawhbond[[framenum]]
  rawhbondtemp<-rawhbondtemp[,c(1,2)]
  cuthis(rawhbondtemp,dishisbreak,anglehisbreak)
}
screenhbond<-function(framenum,hang,hdis){
  rawhbondtemp <- totrawhbond[[framenum]]
  rawhbondtemp<-rawhbondtemp[which(as.numeric(rawhbondtemp[,3])>hang & as.numeric(rawhbondtemp[,1])<hdis),]
  rawhbondtemp
}
gethbondcoordhisr<-function(framenum,rbin,rmax=FALSE) {
  rawhbondtemp <- as.data.frame(tothbond[[framenum]])
  rawhbondtemp <- rawhbondtemp %>% mutate (r=sqrt(as.numeric(V4)^2+as.numeric(V5)^2))
  if (rmax==FALSE) {
    rmax=max(rawhbondtemp$r)
  }
  table(cut(rawhbondtemp$r, breaks=seq(0,rmax,rbin), include.lowest=TRUE))
}
gethbondcoordhisz<-function(framenum,zbin,zmax=FALSE) {
  rawhbondtemp <- as.data.frame(tothbond[[framenum]])
  rawhbondtemp <- abs(as.numeric(rawhbondtemp[,6]))
  if (zmax==FALSE) {
    zmax=max(rawhbondtemp)
  }
  table(cut(rawhbondtemp, breaks=seq(0,zmax,zbin), include.lowest=TRUE))
}
gethbondtypehis<-function(framenum){
  temptype<-tothbond[[frame]][,7]
  # Now 'temptype' has "AM1-AM2" instead of "AM2-AM1"
  #temptype <- lapply(temptype, function(x) gsub("AM2-AM1", "AM1-AM2", x))
  temptypehis<-factor(temptype, levels = c("AM1-AM1", "AM2-AM2","AM1-AM2","AM2-AM1"))
  table(temptypehis)
}
getradialdir<-function(framenum,hbondtypeinp=NULL,dim=1) {
  out<-c()
  if (is.null(hbondtypeinp)) {
    hbondtemp<-tothbond[[framenum]]
  } else {
    hbondtemp<-tothbond[[framenum]][which(tothbond[[framenum]][,10]==hbondtypeinp),]
  }
  if (length(hbondtemp)>0) {
    if (is.null(nrow(hbondtemp))) {
      out<-getcostheta(c(0,0,1),as.numeric(hbondtemp[c(7,8,9)]))
    } else {
      for (i in 1:nrow(hbondtemp)) {
        HbondAng<-getcostheta(c(0,0,1),as.numeric(hbondtemp[i,c(7,8,9)]))
        #acosHbondAngg<-acos(HbondAng)*180/3.14
        out<-c(out,HbondAng)
      }
    }
    out<-cut(as.numeric(abs(out)), breaks=seq(0,1,orientbin), include.lowest=TRUE)
    if (dim==2){
      out1<-cut(as.numeric(hbondtemp[,2]), breaks=seq(-1,1,orientbin), include.lowest=TRUE)
      out<-table(out,out1)
    } else {
      out<-table(out)
    }
    out
  } else {
    0
  }
}



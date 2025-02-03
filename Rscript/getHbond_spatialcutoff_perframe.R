source("Rscript/getHbond_spatialcutoff.R")

#get donor acceptor distance, and donor hydrogen acceptor angle per frame
getHbond_spatialcutoff_perframe<-function(frame,getHist=FALSE,zlim=100) {
  if (getHist) {
    ONcutoff=3.5
  } else {
    ONcutoff=4
  }
  #loop through donor
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  output<-c()
  #loop through unique resname
  for (resname in uniseltextlist) {
    #since there are 2 AM group
    for (AMindex in 1:2) {
      #get the index for extracting coordinates
      DonorListCor<-get(paste0(resname,"donorlist",AMindex,"cor"))
      #get the specific donor from that 
      DonorList<-get(paste0(resname,"donorlist",AMindex))
      for (residindex in 1:length(DonorList)) {
        outputtemp<-getHbond_spatialcutoff(AMindex,particle,residindex,resname,DonorListCor,DonorList,frame,ONcutoff=ONcutoff)
        output<-rbind(output,outputtemp)
        # filtered<-as.data.frame(outputtemp[which(outputtemp[,1]<3.5 & outputtemp[,3]>147),])
        # print(c(nrow(filtered),nrow(filtered)/(2*1174),resname,residindex))
      }
    }
  }
  if (getHist) {
    output<-cbind(output,apply(output[, c(12, 13)], 1, function(row) {
      paste(sort(row), collapse = "")
    }))
    HisOut<-list()
    filtered<-as.data.frame(output[which(output[,1]<3.5 & output[,3]>147),])
    filtered[,10][which(filtered[,10]==21)]<-12
    for(i in 1:11) {
      filtered[,i]<-as.numeric(filtered[,i])
    }
    for (pair in unique(filtered[,10])) {
      pairname<-paste0("P",pair)
      HisOut[[pairname]]<-list()
      HbondTemp<-as.data.frame(filtered[which(filtered[,10]==pair),])
      for(i in 1:11) {
        HbondTemp[,i]<-as.numeric(HbondTemp[,i])
      }
      
      HisOut[[pairname]][["Z"]]<-table(cut(HbondTemp[,6],breaks = seq(-160, 160, by = 1), include.lowest = TRUE))
      
      ResPairHis<-table(HbondTemp[,14])
      ResPairHis[setdiff(uniqueResnamePair,rownames(ResPairHis))] <- 0
      ResPairHis<-ResPairHis[order(names(ResPairHis))]
      HisOut[[pairname]][["ResPair"]] <- ResPairHis
      HbondTemp<-HbondTemp[which(abs(HbondTemp[,6])<zlim),]
      HisOut[[pairname]][["R"]]<-table(cut(sqrt(HbondTemp[,4]^2+HbondTemp[,5]^2),breaks = seq(0, 30, by = 1), include.lowest = TRUE))
    }
    HisOut
  } else {
    output
  }
}


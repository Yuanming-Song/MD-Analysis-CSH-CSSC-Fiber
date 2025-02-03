# Define the moiety name vectors
PHEmoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
AMmoietyname <- c("N1", "C1", "O1", "N2", "C4", "O2", "C2")
THImoietyname <- c("C3", "S", "S1")
PHEmoietyname2 <- c("C12", "C13", "C14", "C16", "C18", "C19")
AMmoietyname2 <- c("N3", "C17", "C20", "O3", "O4", "N4","C15")
THImoietyname2 <- c("C11", "S2")

for (moiety in potential_pairs) {
  
  tempXYZlist<-list()
  for (i in 1:1000) {
    #####   select input res as input for getdist()
    moietysel<-atom.select(mypdb,resno=i,elety=get(paste0(moiety,"moietyname")),operator = "AND")
    if (length(moietysel$xyz)==0) {
      break
    }
    tempXYZlist[[i]]<-list()
    #####   store to a list
    tempXYZlist[[i]][[1]]<-moietysel$xyz
    
    #####   select input res as input for getdist()
    moietysel<-atom.select(mypdb,resno=i,elety=get(paste0(moiety,"moietyname2")),operator = "AND")
    if (length(moietysel$xyz)==0) {
      next
    }
    #####   store to a list
    tempXYZlist[[i]][[2]]<-moietysel$xyz
  }
  assign(paste0(moiety,"XYZ"),tempXYZlist)
  rm(tempXYZlist)
}


getCOM_per_moiety <- function(tempdcd) {
  moietyout<-list()    
  for(moiety in potential_pairs) {
    tempCOMlist<-list()
    tempXYZlist<-get(paste0(moiety,"XYZ"))
    for (i in 1:length(tempXYZlist)) {
      tempCOMlist[[i]]<-list()
      for (j in 1:length(tempXYZlist[[i]])) {
        TempMoietyXYZ<-tempXYZlist[[i]][[j]]
        #####   get coordinate
        tempmoietyxyz<-tempdcd[TempMoietyXYZ]
        #####   calculate COM of res i
        tempCOMlist[[i]][[j]]<-as.vector(com.xyz(tempmoietyxyz))
      }
    }
    moietyout[[moiety]]<-tempCOMlist
  }
  moietyout
}


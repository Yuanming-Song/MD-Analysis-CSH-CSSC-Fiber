comoutname<-"data/COM.rda"
comoutname1<-"data/COM_phe.rda"
comoutname2<-"data/COM_am.rda"
comoutname3<-"data/COM_thiol.rda"
load("data/COMcutoff.rda")
potential_pairs <- c("AM", "PHE", "THI")

getrho_fibre_network <- function(frame) {
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  #edges<-getedge_cutoff(frame)
  edges<-getedge_cutoff_general(frame,tempdcd)
  attr(edges,"n")<-max(mypdb$atom$resno)
  compdist<-component.dist(edges)
  Resid_in_net<-which(compdist$membership==which(compdist$csize==max(compdist$csize)))
  nCSSC<-length(which(Resid_in_net>maxcssc))
  fibresel<-atom.select(mypdb,resno=Resid_in_net) #select heavy atoms only
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  c(frame,length(fibresel$atom)/convhulln(particle,"FA")$vol,length(Resid_in_net),nCSSC) #calculate heavy atom density 
}


# Define a function to convert COMlist entries into a data frame

getedge_cutoff<-function(framenum,firstpasscutoff=10){
  # Initialize an empty list to store the results
  tempedgelist <- c()
  
  #loop through each moiety name combine
  for (COM1ind in 1:length(potential_pairs) ) {
    COM1<-potential_pairs[COM1ind]
    for (COM2ind in COM1ind:length(potential_pairs) ) {
      COM2<-potential_pairs[COM2ind]
      
      # Define the cutoff based on the moiety indices (use current_row$moietyind and retained_rows$moietyind)
      
      COMlist1<-get(paste0("COMlist",COM1))[[framenum]]
      COMlist2<-get(paste0("COMlist",COM2))[[framenum]]
      
      # Convert COMlist1 and COMlist2 into data frames
      COMlist1_df <- convert_COMlist_to_df(COMlist1)
      COMlist2_df <- convert_COMlist_to_df(COMlist2)
      
      
      
      
      # Loop through each row in COMlist1_df
      for (i in seq_len(nrow(COMlist1_df))) {
        current_row <- COMlist1_df[i, ]
        
        # Retain rows in COMlist2_df with resid larger than the current row's resid
        retained_rows <- COMlist2_df %>%
          filter(resid > current_row$resid)
        
        # Subtract the current coordinates from all retained rows
        retained_rows <- retained_rows %>%
          mutate(
            x_diff = x - current_row$x,
            y_diff = y - current_row$y,
            z_diff = z - current_row$z
          )
        
        # Retain only rows where the absolute values of all first 3 columns are smaller than firstpasscutoff
        retained_rows <- retained_rows %>%
          filter(abs(x_diff) < firstpasscutoff, abs(y_diff) < firstpasscutoff, abs(z_diff) < firstpasscutoff)
        
        # Compute the Euclidean distance for each retained row
        retained_rows <- retained_rows %>%
          mutate(
            distance = sqrt(x_diff^2 + y_diff^2 + z_diff^2)
          )
        
        
        # For rows where the distance is smaller than cutoff, append to tempedgelist
        for (j in seq_len(nrow(retained_rows))) {
          pair1 <- ifelse(current_row$resid <= maxcssc , "CSSC", "CSH")
          pair2 <- ifelse(retained_rows$resid[j] <= maxcssc, "CSSC", "CSH")
          contactcutoff<-getcutoff(pair1,pair2,COM1,COM2)
          
          if (retained_rows$distance[j] < contactcutoff) {
            tempedgelist <- rbind(tempedgelist, c(current_row$resid, retained_rows$resid[j]))
          }
        }
      }
    }
  }
  # Convert tempedgelist to a data frame if it's not already one
  tempedgelist_df <- as.data.frame(tempedgelist)
  
  # Retain only unique rows
  tempedgelist_unique <- distinct(tempedgelist_df)
  tempedgelist_unique<-as.matrix(cbind(tempedgelist_unique,1))
  tempedgelist_unique<-rbind(tempedgelist_unique,tempedgelist_unique[,c(2,1,3)])
  tempedgelist_unique
}

getedge_cutoff_general<-function(framenum,tempdcd,firstpasscutoff=10){
  # Initialize an empty list to store the results
  tempedgelist <- c()
  COMlist<-getCOM_per_moiety(tempdcd)
  #loop through each moiety name combine
  for (COM1ind in 1:length(potential_pairs) ) {
    COM1<-potential_pairs[COM1ind]
    for (COM2ind in COM1ind:length(potential_pairs) ) {
      COM2<-potential_pairs[COM2ind]
      
      # Define the cutoff based on the moiety indices (use current_row$moietyind and retained_rows$moietyind)
      
      COMlist1<-COMlist[[COM1]]
      COMlist2<-COMlist[[COM2]]
      
      # Convert COMlist1 and COMlist2 into data frames
      COMlist1_df <- convert_COMlist_to_df(COMlist1)
      COMlist2_df <- convert_COMlist_to_df(COMlist2)
      
      # Loop through each row in COMlist1_df
      for (i in seq_len(nrow(COMlist1_df))) {
        current_row <- COMlist1_df[i, ]
        
        # Retain rows in COMlist2_df with resid larger than the current row's resid
        retained_rows <- COMlist2_df %>%
          filter(resid > current_row$resid)
        
        # Subtract the current coordinates from all retained rows
        retained_rows <- retained_rows %>%
          mutate(
            x_diff = x - current_row$x,
            y_diff = y - current_row$y,
            z_diff = z - current_row$z
          )
        
        # Retain only rows where the absolute values of all first 3 columns are smaller than firstpasscutoff
        retained_rows <- retained_rows %>%
          filter(abs(x_diff) < firstpasscutoff, abs(y_diff) < firstpasscutoff, abs(z_diff) < firstpasscutoff)
        
        # Compute the Euclidean distance for each retained row
        retained_rows <- retained_rows %>%
          mutate(
            distance = sqrt(x_diff^2 + y_diff^2 + z_diff^2)
          )
        
        
        # For rows where the distance is smaller than cutoff, append to tempedgelist
        for (j in seq_len(nrow(retained_rows))) {
          pair1 <- ifelse(current_row$resid <= maxcssc , "CSSC", "CSH")
          pair2 <- ifelse(retained_rows$resid[j] <= maxcssc, "CSSC", "CSH")
          contactcutoff<-getcutoff(pair1,pair2,COM1,COM2)
          
          if (retained_rows$distance[j] < contactcutoff) {
            tempedgelist <- rbind(tempedgelist, c(current_row$resid, retained_rows$resid[j]))
          }
        }
      }
    }
  }
  # Convert tempedgelist to a data frame if it's not already one
  tempedgelist_df <- as.data.frame(tempedgelist)
  
  # Retain only unique rows
  tempedgelist_unique <- distinct(tempedgelist_df)
  tempedgelist_unique<-as.matrix(cbind(tempedgelist_unique,1))
  tempedgelist_unique<-rbind(tempedgelist_unique,tempedgelist_unique[,c(2,1,3)])
  tempedgelist_unique
}

#for the initially wrapped dcd file 

#slice it at each y value

#pick a radom CSH/CSSC, wrap only that section of the cell

#check if the new com is within quarter of box origin

getrho_fibre_network_v4 <- function(frame) {
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  #edges<-getedge_cutoff(frame)
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  edges<-record_contact_edges(particle, cutoff_matrix, AtomIndexLib, AtomTypeLib)
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  compdist<-component.dist(edges)
  Resid_in_net<-which(compdist$membership==which(compdist$csize==max(compdist$csize)))
  nCSSC<-length(which(Resid_in_net>maxcssc))
  fibresel<-atom.select(mypdb,resno=Resid_in_net) #select heavy atoms only
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  c(frame,length(fibresel$atom)/convhulln(particle,"FA")$vol,length(Resid_in_net),nCSSC) #calculate heavy atom density 
}





getrdf_COM<-function(framenum,firstpasscutoff=20){
  # Initialize an empty list to store the results
  rdfout<-list()
  rdfout[["CSSC"]]<-list()
  rdfout[["CSH"]]<-list()
  rdfout[["MIX"]]<-list()
  
  #loop through each moiety name combine
  for (COM1ind in 1:length(potential_pairs) ) {
    COM1<-potential_pairs[COM1ind]
    for (COM2ind in COM1ind:length(potential_pairs) ) {
      COM2<-potential_pairs[COM2ind]
      rdfout[["CSSC"]][[paste(COM1,COM1,sep="-")]]<-c()
      rdfout[["CSH"]][[paste(COM1,COM1,sep="-")]]<-c()
      rdfout[["MIX"]][[paste(COM1,COM1,sep="-")]]<-c()
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
          if(pair1==pair2) {
            rdfout[[pair1]][[paste(COM1,COM1,sep="-")]]<-c(rdfout[[pair1]][[paste(COM1,COM1,sep="-")]],retained_rows$resid[j])
          } else {
            rdfout[["MIX"]][[paste(COM1,COM1,sep="-")]]<-c(rdfout[["MIX"]][[paste(COM1,COM1,sep="-")]],retained_rows$resid[j])
          }
        }
      }
    }
  }
  for (name in names(rdfout)) {
    for (pair in names(rdfout[[name]])) {
      rdfout[[name]][[pair]]<-table(cut(rdfout[[name]][[pair]],breaks=seq(0,30,0.5)))
    }
  }
  rdfout
}

# step<-25
# # Set the directory where your files are stored
# dataDir <- "/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/edge/"
# 
# # Set the prefix and extension of your files
# infile <- "mix_27A_3to7_cyl_nottip_every"
# ext <- ".edgeatt"
# 
# # Set the number of frames to loop through
# totframes <- 200
# 
# # Create a list of file names
# file_list <- c()
# 
# for (frame in 1:totframes) {
#   framenum <- sprintf("%03d", frame)
#   file_name <- paste0(dataDir, infile, "_", step,"_",framenum, ext)
#   file_list[[frame]] <- file_name
# }
# # define maxcssc number
# mincssc <- 1409-1
# 
# # initialize histogram table
# hist_table <- data.frame(Pair = character(),
#                          Type = character(),
#                          Count = numeric(),
#                          stringsAsFactors = FALSE)
# 
# # get list of files in directory
# 
# # loop over files
# for (file_name in file_list) {
#   # read table
#   table_data <- read.table(file_name, header = FALSE, col.names = c("Num1", "Num2", "Type"))
#   # iterate over rows in table
#   for (i in 1:nrow(table_data)) {
#     # determine pair type based on maxcssc number
#     if (table_data[i, "Num1"] > mincssc & table_data[i, "Num2"] > mincssc) {
#       pair_type <- "CSSC-CSSC"
#     } else if (table_data[i, "Num1"] > mincssc | table_data[i, "Num2"] > mincssc) {
#       pair_type <- "CSH-CSSC"
#     } else {
#       pair_type <- "CSH-CSH"
#     }
#     # update histogram table
#     hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] <- hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] + 1
#     if (nrow(hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], ]) == 0) {
#       hist_table <- rbind(hist_table, data.frame(Pair = pair_type, Type = table_data[i, "Type"], Count = 1, stringsAsFactors = FALSE))
#     }
#   }
# }
# write.table(hist_table,file = paste0(infile,"int.dat"),row.names=FALSE,quote = F)

step <- 25
# Set the directory where your files are stored
dataDir <- "/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/edge/"

# Set the prefix and extension of your files
infile <- "mix_27A_3to7_cyl_tip_every"
ext <- ".edgeatt"
mincssc <- 1409-1

# Set the number of frames to loop through
totframes <- 200
tablelist<-c()
for (frame in 1:totframes) {
  framenum <- sprintf("%03d", frame)
  file_name <- paste0(dataDir, infile, "_", step, "_", framenum, ext)
  tablelist[[frame]] <- read.table(file_name, header = FALSE, col.names = c("Num1", "Num2", "Type"))
}

getframehis<-function(table_data) {
  # Create an empty histogram table for this frame
  hist_table <- data.frame(Pair = character(),
                           Type = character(),
                           Count = numeric(),
                           stringsAsFactors = FALSE)
  
  # read table
  
  for (i in 1:nrow(table_data)) {
    # determine pair type based on maxcssc number
    if (table_data[i, "Num1"] > mincssc & table_data[i, "Num2"] > mincssc) {
      pair_type <- "CSSC-CSSC"
    } else if (table_data[i, "Num1"] > mincssc | table_data[i, "Num2"] > mincssc) {
      pair_type <- "CSH-CSSC"
    } else {
      pair_type <- "CSH-CSH"
    }
    # update histogram table
    hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] <- hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] + 1
    if (nrow(hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], ]) == 0) {
      hist_table <- rbind(hist_table, data.frame(Pair = pair_type, Type = table_data[i, "Type"], Count = 1, stringsAsFactors = FALSE))
    }
  }
  hist_table$Count<-hist_table$Count/sum(hist_table$Count)
  # add histogram for this frame to the list
  hist_table
}
getframehis2<-function(table_data) {
  # Create an empty histogram table for this frame
  hist_table <- data.frame(Pair = character(),
                           Type = character(),
                           Count = numeric(),
                           stringsAsFactors = FALSE)
  
  # read table
  table_data<-as.data.frame(table_data)
  table_data[,1]<-as.numeric(table_data[,1])
  table_data[,2]<-as.numeric(table_data[,2])
  
  colnames(table_data) <- c("Num1", "Num2", "Type")
  
  for (i in 1:nrow(table_data)) {
    # determine pair type based on maxcssc number
    if (table_data[i, "Num1"] <= maxcssc*6 & table_data[i, "Num2"] <= maxcssc*6) {
      pair_type <- "CSSC-CSSC"
    } else if (table_data[i, "Num1"] > maxcssc*6 & table_data[i, "Num2"] > maxcssc*6) {
      pair_type <- "CSH-CSH"
    } else {
      pair_type <- "CSH-CSSC"
    }
    # update histogram table
    hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] <- hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], "Count"] + 1
    if (nrow(hist_table[hist_table$Pair == pair_type & hist_table$Type == table_data[i, "Type"], ]) == 0) {
      hist_table <- rbind(hist_table, data.frame(Pair = pair_type, Type = table_data[i, "Type"], Count = 1, stringsAsFactors = FALSE))
    }
  }
  hist_table$Count<-hist_table$Count/sum(hist_table$Count)
  # add histogram for this frame to the list
  hist_table
}

library(doParallel)
registerDoParallel(cores = detectCores())
for (j in 1:length(edges)) {
  a<-edges[[j]][[1]][seq(1,nrow(edges[[j]][[1]])/2,1),]
  interactType_matrix <- matrix(nrow = nrow(a), ncol = 3)
  for (i in 1:nrow(a)) {
    pair<-sort(c(a[i,c(1,2)]))
    pairindex<-paste(pair[1],pair[2])
    interaction<-paste(network[[pair[1]]]$nname,network[[pair[2]]]$nname)
    if (interaction == "PHE PHE") {
      interactType <- "PHE-PHE"
    } else if (interaction %in% c("AM AM")) {
      interactType <- "AM-AM"
    } else if (interaction == "THI THI") {
      interactType<- "DISU"
    } else if (interaction %in% c("AM PHE", "PHE AM")) {
      interactType <- "AM-PHE"
    } else if (interaction %in% c("THI PHE", "PHE THI")) {
      interactType <- "THI-PHE"
    } else if (interaction %in% c("THI AM", "AM THI")) {
      interactType <- "THI-AM"
    } else {
      print(c(pair,interaction))
      break
    }
    interactType_matrix[i,]<-c(pair,interactType)
  }
  edges[[j]][[2]]<-interactType_matrix
  if (nrow(a)!=nrow(interactType_matrix)) {
    print("error")
    break
  }
}

# Create a list to store histograms for each frame
hist_frames  <- foreach (frame = 1:length(gs),.combine=rbind) %dopar% {
  #getframehis(tablelist[[frame]])
  getframehis2(edges[[frame]][[2]])
}

Mean <- aggregate(hist_frames$Count, by = list(Pair = hist_frames$Pair, Type = hist_frames$Type), FUN = function(x) c(mean(x)))
SD <- aggregate(hist_frames$Count, by = list(Pair = hist_frames$Pair, Type = hist_frames$Type), FUN = function(x) c(sd(x)))[,3]
SE <- aggregate(hist_frames$Count, by = list(Pair = hist_frames$Pair, Type = hist_frames$Type), FUN = function(x) c(sd(x)/sqrt(length(x))))[,3]
hist_frames_out<-cbind(Mean, SD,SE)
colnames(hist_frames_out) <- c("Pair", "Type", "Mean", "SD", "SE")
hist_frames_out$Type <- ifelse(hist_frames_out$Type == "DISU", "SULF-SULF", hist_frames_out$Type)

# write table to file
write.table(hist_frames_out, file = paste0(infile, "int.dat"), row.names = FALSE, quote = FALSE)


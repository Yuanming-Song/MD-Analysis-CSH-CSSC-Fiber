source("Rscript/getHbond_spatialcutoff_visual.R")

getHbond_spatialcutoff_perframe_visual<-function(frame,getHist=TRUE) {
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
  #loop through unique resname
getHbond_spatialcutoff_visual(particle,frame,ONcutoff=ONcutoff)

}

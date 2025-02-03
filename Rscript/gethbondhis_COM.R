# Define parameters for hbond
hdis <- 3.5
hang <- 180 - 33
rbin <- 0.5
rmax <- 30
zbin <- 3
zmax <- 90
# Set the distance bin size max and min for distance
disbinsize <- 0.1
disbinmax <- 40# sqrt(2 * 100^2) / 2
disbinmin <- 0
# Set the distance bin size max and min for distance
anglebinsize <- 0.05
anglebinmax <- 1
anglebinmin <- -1
# File names
rawhbonddisname <- "hbonddis_COM.rda"
rawhbondname <- "rawhbond_COM.rda"
tothbondname <- "hbond_COM.rda"
hbondcoordhisname <- "hbondcoordis_COM.rda"
hbondtypehisname<-"hbondtypedis_COM.rda"
hbondradialbytypename<-"hbondradial_COM_type.rda"
hbondradialname<-"hbondradial_COM.rda"

# Load necessary functions
source("Rscript/gethbond_COM.base.R")
# Load rawhbond data
load(rawhbondname)
nframes <- length(totrawhbond)
# Calculate 2D histogram of distances
rawdistmout <- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
  get2dhis(framenum)
}
save(rawdistmout, file = rawhbonddisname)
# Calculate screened hydrogen bonds
tothbond <- foreach(framenum = 1:nframes) %dopar% {
  screenhbond(framenum, hang, hdis)
}
save(tothbond, file = tothbondname)

# Calculate radial distribution function (r) and axial distribution function (z)
rawhbondcoordhisr <- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
  gethbondcoordhisr(framenum, rbin, rmax)
}
rawhbondcoordhisr <- as.data.frame(cbind(seq(0 + rbin, rmax, rbin), rawhbondcoordhisr))
colnames(rawhbondcoordhisr) <- c("r", "Freq")
hbondcoordhisr <- rawhbondcoordhisr %>% mutate(Density = Freq / (sum(Freq) * pi * (r^2 - (r - rbin)^2)))

rawhbondcoordhisz <- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
  gethbondcoordhisz(framenum, zbin, zmax)
}
rawhbondcoordhisz <- as.data.frame(cbind(seq(0 + zbin, zmax, zbin), rawhbondcoordhisz))
colnames(rawhbondcoordhisz) <- c("z", "Freq")
hbondcoordhisz <- rawhbondcoordhisz %>% mutate(Density = Freq / (sum(Freq)))

# Combine radial and axial distributions into a list
hbondcoordhis <- list()
hbondcoordhis[[1]] <- hbondcoordhisr
hbondcoordhis[[2]] <- hbondcoordhisz
save(hbondcoordhis, file = hbondcoordhisname)

# Calculate hydrogen bond type distribution
rawhbondtypehis <- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
  gethbondtypehis(framenum)
}
rawhbondtypehis <- as.data.frame(rawhbondtypehis)
colnames(rawhbondtypehis) <- c("Type", "Freq")
hbondtypehis <- rawhbondtypehis %>% mutate(Density = Freq / (sum(Freq)))
save(hbondtypehis,file = hbondtypehisname)
orientbin<-0.1
hbondradialbytypetotal<-c()
for (hbondtype in unique(tothbond[[1]][,10])) {
  rawhbondradial <- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
    getradialdir(framenum,hbondtype,1)
  }
  rawhbondradial<-as.data.frame(cbind(seq(0+orientbin/2,1-orientbin/2,orientbin),rawhbondradial,hbondtype))
  colnames(rawhbondradial)<-c("Angle","Freq","Type")
  rawhbondradial$Freq<-as.numeric(rawhbondradial$Freq)
  rawhbondradial<-rawhbondradial %>% mutate (Density=Freq/(sum(Freq)))
  hbondradialbytypetotal<-rbind(hbondradialbytypetotal,rawhbondradial)
}
save(hbondradialbytypetotal,file=hbondradialbytypename)
hbondradialtotal<- foreach(framenum = 1:nframes, .combine = "+") %dopar% {
  getradialdir(framenum,dim=2)
}
save(hbondradialtotal,file=hbondradialname)


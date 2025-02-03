systemlist<-c("")
#systemlist<-c("_CSSC","_CSH")
### how many residues in total
steps<-10
numclus<-5
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
totcssc<-maxcssc
totcsh<-maxcsh-maxcssc
infilepre<-"mix_27A_COM"
outPre<-paste("data/", "component_resid_CSSC_", infilepre, "_every",steps,sep="")
library(sna)
library(doParallel)
corenum=24
corenum=detectCores()
registerDoParallel(cores = corenum)
setnodeslist<-function() {
  for (count in 1:length(moietynamelist)) {
    moid<<-moid+1
    nname<-moietynamelist[count]
    ntype<-typenameslist[count]
    nodeslist[[moid]]<<-c(
       resname,
       resid,
      moid,
       nname,
      ntype
    )
  }
}
nodeslist<-list()
moid<-0
moietynames <- c("AM", "PHE", "THI")
typenames<-c("DIP","NOP","THIO")
resname<-"CSSC"
moietynamelist<-c(moietynames,moietynames)
typenameslist<-c(typenames,typenames)
for (resid in mincssc:maxcssc) {
  setnodeslist()
}
resname<-"CSH"
moietynamelist<-moietynames
typenameslist<-typenames
for (resid in mincsh:maxcsh) {
  setnodeslist()
}
nodes<-c()
for (i in 1:length(nodeslist)) {
  nodes<-rbind(nodes,c(i,nodeslist[[i]]))
}
colnames(nodes)<-c("id","resname","resid","molid","name","type")
getResid<-function(g,nodes,firstclu,lastclu){
  g<-g[which(g[,1]<=maxcssc*6 & g[,2]<=maxcssc*6),c(1,2)]
  attr(g,"n")<-length(nodes)
  cd<-component.dist(g)
  myclus<-order(cd$csize,decreasing = T)[firstclu:lastclu]
  myresids<-lapply(myclus,function(z,cd,nodes){
    mynodes=which(cd$membership==z)
    return(unique(nodes[mynodes,"resid"]))
  },nodes=nodes,cd=cd)
  return(myresids)
}
getturnover<-function(i,a,b){
  cssca<-a[which(a<=maxcssc)]
  csha<-a[which(a>maxcssc)]
  csscb<-b[which(b<=maxcssc)]
  cshb<-b[which(b>maxcssc)]
  if (length(which(cssca %in% csscb))==0) {
    csscout<-cssca 
  } else {
    csscout<-cssca[-which(cssca %in% csscb)]
  }
  if (length(which(csha %in% cshb))==0) {
    cshout<-csha
  } else {
    cshout<-csha[-which(csha %in% cshb)]
  }
  csscin<-csscb[-which(csscb %in% cssca)]
  cshin<-cshb[-which(cshb %in% csha)]
  if (length(csha)==0) {
    cshturnover<-0
  } else {
    cshturnover<-length(cshout)/length(csha)
  }
  if (length(cssca)==0) {
    csscturnover<-0
  } else {
    csscturnover<-length(csscout)/length(cssca)
  }
  rbind(c(i,cshturnover, "CSH"), c(i,csscturnover, "CSSC"))
}
for (system in systemlist) { 
  infile<-paste0(infilepre,system,".edgel.stack.rda")
  load(infile)
  gs<-edges
  rm(edges)
  gsize<-3*(maxcssc+maxcsh)#attr(gs[[1]],"n")
  print(length(gs))
  resids<-c()
  emptylist<-c()
  for (frame in seq(1,length(gs),corenum)) {
    lastframe=frame+corenum-1
    if (lastframe>length(gs)) {
      lastframe=length(gs)
    }
    residtemp<-foreach(j=frame:lastframe) %dopar% (
      #resids<-foreach(j=1:10) %dopar% (
      getResid(gs[[j]][[1]],nodes=nodes,1,numclus)
    )
    counter<-1
    for (i in frame:lastframe) {
      if (length(residtemp[[counter]])==0) {
        emptylist<-c(emptylist,i)
      }
      resids[[i]]<-residtemp[[counter]]
      counter<-counter+1
    }
    save(resids,file=paste0(outPre,system,"1-",numclus,".rda"))
    print(paste("Just finished frame",lastframe,"and emptylist is"))
    print(emptylist)
  }
  #load("component_resid_csh_cssc_1to1_50mM_cg_big2_every10.rda")
  for (j in 1:length(resids)) {
    #for (j in 1:10) {
    if (length(resids[[j]])==0) {
      emptylist<-c(emptylist,j)
      #resids[[j]]<-getResid(gs[[j]],nodes=nodes,4,10)
    }
  }
  if (length(emptylist)>1) {
    residtemp<-foreach(j=1:length(emptylist)) %dopar% (
      getResid(gs[[emptylist[j]]],nodes=nodes,1,numclus)
    )	
    for (j in 1:length(emptylist)) {
      resids[[emptylist[j]]]<-residtemp[[j]]
    }
    save(resids,file=paste0(outPre,system,"1-",numclus,"2.rda"))
  }	
  if (system=="") {
    totclus<-length(resids[[1]])
    lastframe<-length(resids)
    replist<-matrix(0, nrow = length(resids),ncol = totclus)
    tempreplist<-which(replist[lastframe,]==0)
    tempcluslist<-sapply(tempreplist,function(x){
      #length(resids[[i]][[which(resids[[i]][[x]]<=maxcssc)]])
      length(which(resids[[lastframe]][[x]]<=maxcssc))
    })
    clusid<-length(tempcluslist[which(tempcluslist>2)])
    totclus<-length(resids[[1]])
    replist<-matrix(0, nrow = length(resids),ncol = totclus)
    tempreplist<-which(replist[lastframe,]==0)
    tempcluslist<-sapply(tempreplist,function(x){
      #length(resids[[i]][[which(resids[[i]][[x]]<=maxcssc)]])
      #length(which(resids[[lastframe]][[x]]<=maxcssc))
      length(resids[[lastframe]][[x]])
    })
    clusid<-length(tempcluslist[which(tempcluslist>4)])
    clusidbig<-clusid
    for (k in 1:clusid) {
      if (k==1) {
        residsnew<-c()
        compo<-c()
        missframe<-c()
        missclu<-c()
      }
      componame<-paste0("compo",k)
      assign(componame,c())
      for (i in seq(lastframe,2,-1)) {
        if (i==lastframe) {
          residsnew[[lastframe]][[k]]<-resids[[lastframe]][[k]]
          compo<-rbind(compo,c(lastframe,length(resids[[lastframe]][[k]]),paste("cluster",k)))
          assign(componame,rbind(get(componame),c(lastframe,length(resids[[lastframe]][[k]]),paste("cluster",k))))
          replist[lastframe,k]<-k
        }
        residlist<-residsnew[[i]][[k]]
        turnoverlist<-c()
        for (j in 1:totclus) {
          if (replist[i-1,j]==0 & length(resids[[i-1]][[j]])>10) {
            turnovertep<-getturnover(i,resids[[i-1]][[j]],residlist)[1,2]
          } else {
            turnovertep<-100
          }
          turnoverlist<-c(turnoverlist,turnovertep)
        }
        if (turnoverlist[which.min(turnoverlist)]==1) {
          missframe<-c(missframe,i-1)
          missclu<-c(missclu,k)
          break
        } else {
          replist[i-1,which.min(turnoverlist)]<-k
          residsnew[[i-1]][[k]]<-resids[[i-1]][[which.min(turnoverlist)]]
          compo<-rbind(compo,c(i-1,length(resids[[i-1]][[which.min(turnoverlist)]]),paste("cluster",k)))
          assign(componame,rbind(get(componame),c(i-1,length(resids[[i-1]][[which.min(turnoverlist)]]),paste("cluster",k))))
        }
      }
    }
    for (i in seq(length(resids),50,-1)) {
      while (TRUE) {
        tempreplist<-which(replist[i,]==0)
        tempcluslist<-sapply(tempreplist,function(x){
          #length(resids[[i]][[which(resids[[i]][[x]]<=maxcssc)]])
          #length(which(resids[[i]][[x]]<=maxcssc))
          length(resids[[i]][[x]])

        })
        k<-tempreplist[which.max(tempcluslist)]
        if (tempcluslist[which.max(tempcluslist)]>10) {
          lastmissframe<-i
          clusid<-clusid+1
          residsnew[[lastmissframe]][[clusid]]<-resids[[lastmissframe]][[k]]
          replist[lastmissframe,k]<-clusid
          tempilist<-c(lastmissframe)
          tempklist<-c(k)
          residlist<-resids[[lastmissframe]][[k]]
          turnoverlist<-c()
          for (j in 1:clusidbig) {
            turnovertep<-getturnover(lastmissframe+1,residlist,resids[[lastmissframe+1]][[j]])[2,2]
            
            turnoverlist<-c(turnoverlist,turnovertep)
          }
          parentcluid<-which.min(turnoverlist)
          compotemp<-c(lastmissframe,length(resids[[lastmissframe]][[k]]),paste("cluster",clusid,"sub",parentcluid))
          componame<-paste0("compo",parentcluid)
          
          for (i in seq(lastmissframe,2,-1)) {
            residlist<-residsnew[[i]][[clusid]]
            turnoverlist<-c()
            for (j in 1:totclus) {
              if (replist[i-1,j]==0 & length(resids[[i-1]][[j]])>10) {
                turnovertep<-getturnover(i,resids[[i-1]][[j]],residlist)[2,2]
              } else {
                turnovertep<-100
              }
              turnoverlist<-c(turnoverlist,turnovertep)
            }
            tempk<-which.min(turnoverlist)
            if (turnoverlist[tempk]==1) {
              missframe<-c(missframe,i-1)
              missclu<-c(missclu,k)
              break
            } else {
              replist[i-1,tempk]<-clusid
              residsnew[[i-1]][[clusid]]<-resids[[i-1]][[tempk]]
              compotemp<-rbind(compotemp,c(i-1,length(resids[[i-1]][[tempk]]),paste("cluster",clusid,"sub",parentcluid)))
              tempilist<-c(tempilist,i-1)
              tempklist<-c(tempklist,tempk)
            }
          }
          if (length(nrow(compotemp))==0 || length(compotemp[,2])<50) {
            for (l in 1:length(tempilist)){
              tempk<- tempklist[[l]]
              tempi<- tempilist[[l]]
              residsnew[[tempi]][[clusid]]<-c()
              replist[tempi,tempk]<-0
            }
            clusid<-clusid-1
            
          }else {
            compo<-rbind(compo,compotemp)
            assign(componame,rbind(get(componame),compotemp))
            print(paste(clusid,i,tempk))
            #break
          }
        } else {
          break
        }
      }
    }
    turnover<-c()
    for (i in 1:clusidbig) {
      oldres<-resids[[1]][[i]]
      for (frame in 2:lastframe) {
        newres<-resids[[frame]][[j]]
        turnover<-rbind(turnover,getturnover(i,oldres,newres))
        oldres<-newres
      }
      write.table(format(turnover, digits=3),file=paste0("turnover",j,".dat"),col.names=F,row.names=F,quote=F)
      componame<-paste0("compo",i)
      compotemp<-as.data.frame(get(componame))
      compotemp$V1<-as.numeric(compotemp$V1)
      compotemp$V2<-as.numeric(compotemp$V2)
      write.table(compotemp,file=paste0("compo",i,".dat"))
      outname<-paste0("resids",i,".dat",col.names=F,row.names=F,quote=F)
      if (file.exists(outname)==TRUE) {
        file.remove(outname)
      }
      for (j in seq(1,length(resids),1)) {
        myresids<-resids[[j]][[i]]
        #myresids<-myresids[which(myresids<maxcssc & myresids%%2!=0)]
        #if (length(myresids)>0) {
        #write.table(t(j),file=outname,append=T,col.names=F,row.names=F,quote=F)
        write.table(t(resids[[j]][[i]]),file=outname,append=T,col.names=F,row.names=F,quote=F)
        #}
      } 
    }
    save(resids,file=paste0("component_resid_csh_cssc_1to1_50mM_cg_big2_every",steps,"1-",totclus,"sorted.rda"))
  }
}



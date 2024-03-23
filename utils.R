
# equiv to lapply with progress bar
apply.to.all <- function(my.list, func){
  pb <- txtProgressBar(
    min=1, 
    max = length(my.list)+1, 
    style=3)

  i <- 1
  setTxtProgressBar(pb, i)
  
  time.start <- Sys.time()
  result <- c()  
  
  for(el in my.list){
    resel <- func(el)
    if(!is.na(resel)){
      result <- append(result, resel)
    }
    if(i%%2 & Sys.time()-time.start > 40*i){ # display if each takes more than 40 secs
      time.left <- (length(my.list)/i-1)*as.numeric(Sys.time()-time.start)
      print(paste("Temps restant estimé :", time.left %/% 60, "minutes et", time.left %/% 60, "secondes ..."))
    }
    i <- i+1
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(result)
}

create.checkpoint <- function(checkpoint.num, seurobj_vec){
  i<-1
  for(seurOBJ in seurobj_vec){
    SaveSeuratRds(seurOBJ, paste("./checkpoint",num,"/sobj",i,".rds",sep=""))
    i <- i+1
  }
}


# use : seurobj_vec <- resurrect(amount)
resurrect <- function(checkpoint.num, number.of.rds.files, seurobj_vec){
  # resurrect
  seurobj_vec <- c()
  for(i in 1:number.of.rds.files){
    seurobj_dir <- paste("./checkpoint",checkpoint.num,"/sobj",i,".rds", sep="")
    
    seurobj <- LoadSeuratRds(seurobj_dir)
    seurobj_vec <- append(seurobj_vec, seurobj)
  }
  return(seurobj_vec)
}



save.plot.png <- function(plot, dir, width=600, height=600){
  # No switch of directories to save in global working dir
  png(dir, width, height)
  plot(plot)
  dev.off()
}

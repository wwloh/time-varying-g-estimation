rm(list=ls())
library("data.table")

## prep data ##################################################################
## save the below data file into a subfolder "osfstorage-archive"
## data from: https://osf.io/cv2u9/files/osfstorage/60747466f6585f051b61b85e

rawData <- readRDS("osfstorage-archive/processed_data.rds")
rawData <- data.table(rawData)
setkey(rawData)

setcolorder(rawData,c("id","wave",
                      names(rawData)[!(names(rawData) %in% c("id","wave"))]))
setkey(rawData)

# treatment
rawData[,mean(books_used,na.rm=TRUE),by=wave]
# outcome
rawData[,mean(life_satisfaction,na.rm=TRUE),by=wave]

# covariates
Lnames <- c("gender","age",paste0("well_being_",1:2),
            # add other media usage as possible confounders
            "music_time","films_time","tv_time",
            "games_time","magazines_time","audiobooks_time")

# analysis ####################################################################
source("../time_varying-helper_funs.R")

allData <- rawData
allData[, t := wave-1]
allData <- allData[,c("id","t",Lnames,
                      "books_used","life_satisfaction"), with=FALSE]
setkey(allData)

# inspect a single individual at random
some.id <- sample(allData$id,1)
allData[id==some.id]

# grand mean-center age
allData[, age := age - allData[, unique(age), by=id][,mean(V1,na.rm=TRUE)]]

# well-being 2 at each time point t>0 as difference from baseline (t=0)
allData[t==0, wb2_t0 := well_being_2]
allData[, well_being_2 := well_being_2-.SD[t==0, well_being_2], by=id]
allData[t==0, well_being_2 := wb2_t0]
allData[, wb2_t0 := NULL]

allData[id==some.id]
rm(some.id)

str(allData)
summary(allData)

## consider outcome at each wave as end of study outcome in turn ##############
res.list <- NULL
for (T_ in 5:1) {
  # drop individuals with missing data
  drop_id <- allData[t<=T_,any(is.na(.SD)),by=list(id,t),
                     .SDcols=c(Lnames,"books_used","life_satisfaction")][
                       V1==TRUE,id]
  drop_id <- unique(unlist(drop_id))
  Data <- data.table(allData[t<=T_ & !(id %in% drop_id)])
  setkey(Data)
  all_id <- Data[,unique(id)]
  N <- length(all_id)
  cat("T =",T_, "; N =",allData[,length(unique(id))], 
      "; dropped =",length(drop_id), "; complete =",N, "\n")

  # check treatment prevalence at each time point
  print(Data[,mean(books_used,na.rm=TRUE),by=t])
  
  res <- OneEst(Data,
                T_=T_,Lnames=Lnames,
                Treat_name="books_used",Outcome_name="life_satisfaction",
                Znames=c("age","well_being_2"),
                use.DR=TRUE)
  cat("T =", T_, "\n")
  print(lapply(res$est,round,2))
  
  boot.res <- lapply(1:10000, function(bb) {
    bootOK <- FALSE
    while(!bootOK) {
      boot_id <- sort(unique(sample(all_id,N,replace=TRUE)))
      boot_D <- Data[id %in% boot_id]
      setkey(boot_D)
      stm.boot <- tryCatch(system.time(
        boot_res <- OneEst(boot_D,
                           T_=T_,Lnames=Lnames,
                           Treat_name="books_used",
                           Outcome_name="life_satisfaction",
                           Znames=c("age","well_being_2"),
                           use.DR=TRUE)
      )[3], error=function(cond) return(NA))
      if (!is.na(stm.boot)) {
        bootOK <- TRUE
        return(unlist(boot_res$est))
      }
    }
  })
  boot_res <- do.call(rbind,boot.res)
  print(round(t(apply(boot_res,2,quantile,probs=c(.025,.975), na.rm=TRUE)),2))
  
  res.list[[T_]] <- list("res"=res,"boot"=boot_res)
}

save.image("data_prep_analysis.Rdata")
q()

# load results ################################################################
rm(list=ls())
load("data_prep_analysis.Rdata")
library("data.table")
library("xtable")

res.table <- t(rbindlist(lapply(res.list, function(x) {
  x.boot <- x$boot
  x.df <- rbind(
    "est"=unlist(x$res$est)*100,
    "se"=apply(x.boot,2,sd,na.rm=TRUE)*100,
    apply(x.boot,2,quantile,probs=c(.025,.975), na.rm=TRUE)*100)
  data.frame(apply(x.df,2,function(x.res) {
    # fix at two decimal places
    x.res.vals <- format(round(x.res,1), nsmall = 1)
    # include brackets and comma for CI
    c(x.res.vals[1:2], paste0("(",x.res.vals[3]),",",paste0(x.res.vals[4],")"))
  }))
}),fill=TRUE))
xtable(res.table)

# print summary tables of fitted outcome regression models
lapply(res.list,function(one.res) 
  lapply(one.res$res$fitted, function(one.res.fitted)
    summary(one.res.fitted$Outcome)$coef))

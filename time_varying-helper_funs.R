OneEst <- function(Data, # long format
                   T_,  # time point for end-of-study outcome
                   Lnames, # covariate names
                   Treat_name,
                   Outcome_name,
                   Znames=NULL, # effect modifiers in SNMM
                   use.DR=TRUE,
                   fit.ps="mle"
                   ) {
  ## names cannot be "P" or "Ry" or contain "_."
  
  # prepare dataset for fitting ###############################################
  D <- data.table(Data)
  setkey(D)
  
  # create rows with NAs for observations with missing waves
  if (D[,.N,by=id][,any(N<T_+1)]) {
    for (i in D[,.N,by=id][N<T_+1,id]) {
      # create new rows
      D <- rbind(D,data.table("id"=i,"t"=(0:T_)[!((0:T_) %in% D[id==i,t])]),
                 fill=TRUE)
      setkey(D)
    }
  }
  
  # set repeated covariates to value at baseline and to NA at other times
  Lbaseline <- apply(D[,lapply(.SD, function(x) 
    length(unique(x[!is.na(x)]))==1L), by="id"][
      ,!(colnames(D)%in%c("id","t")),with=FALSE],2,all)
  Lbaseline <- Lnames[Lnames %in% names(Lbaseline[Lbaseline])]
  D[t>0,(Lbaseline) := NA]
  rm(Lbaseline)
  
  # reshape to wide format
  timesep <- "_." # used to indicate time point for a variable
  D <- data.table::dcast(data = D, 
                         formula = id ~ t, 
                         value.var = names(D)[!(names(D) %in% c("id","t"))],
                         sep = timesep,
                         fill = NA)
  setkey(D)
  
  ## ignore variables (besides outcome) observed at or after time T
  var_T <- grep(paste0(timesep,T_),colnames(D),value=TRUE)
  D[, var_T[var_T!=paste0(Outcome_name,timesep,T_)]] <- NULL
  rm(var_T)
  
  ## remove columns with all NAs
  D <- D[, colSums(is.na(D))<nrow(D), with=FALSE]
  setkey(D)
  
  res.coef <- res.formulae <- res.fit <- NULL
  for (tt in (T_-1L):0) {
    # treatment name
    A_t <- paste0(Treat_name,timesep,tt)
    # PS name
    P_t <- paste0("P",timesep,tt)
    # outcome name
    if (tt==T_-1L) {
      Rt <- paste0(Outcome_name,timesep,T_)
    } else {
      Rt <- paste0("Ry",timesep,tt+1L)
    }
    
    # variable history
    tt.predictors <- c(
      sapply(Lnames,paste,0:tt,sep=timesep),        # covariate history
      paste0(Outcome_name,timesep,0:max(0,(tt-1))), # outcome history
      paste0(Treat_name,timesep,0:max(0,(tt-1)))    # treatment history
    )
    names(tt.predictors) <- NULL
    tt.predictors <- tt.predictors[(tt.predictors %in% colnames(D)) & 
                                     (tt.predictors != A_t)]
    tt.predictors <- sort(tt.predictors)
    
    tt.ps <- as.formula(paste0(A_t,"~",paste(tt.predictors,collapse="+")))
    
    # fit PS model
    ps.fit <- NULL
    if (use.DR==TRUE) {
      if (fit.ps=="mle") {
        ps.fit <- glm(formula=tt.ps,data=D,family=binomial("logit"))
        ps.hat <- ps.fit$fitted.values
      } else if (fit.ps=="cbps") {
        ps.fit <- CBPS::CBPS(formula=tt.ps,data=D,ATT=0,iterations=10000)
        ps.hat <- ps.fit$fitted.values
      } else if (fit.ps=="gbm") {
        ps.fit <- twang::ps(formula=tt.ps, data=as.data.frame(D),
                            estimand="ATE", stop.method="es.max",verbose=FALSE)
        ps.hat <- ps.fit$ps[,1]
      }
      if (nrow(D)>length(ps.hat)) {
        ps.hat.all <- rep(NA,nrow(D))
        ps.hat.all[as.integer(names(ps.hat))] <- ps.hat
        ps.hat <- ps.hat.all
        rm(ps.hat.all)
      }
      D[, paste0("P",timesep,tt) := ps.hat]
    }
    
    ## predictors from SNMM
    if(is.null(Znames) || !("1" %in% Znames)) {
      Znames <- c("1",Znames)
    }
    
    tt.ZP <- unlist(lapply(Znames, function(z) {
      if (use.DR==TRUE) {
        z_ap <- c(A_t,P_t)
      } else {
        z_ap <- A_t
      }
      if (z=="1") {
        # constant unconditional effect
        return(z_ap)
      } else {
        # interaction with all previous occurrences of the covariate
        z_t <- grep(z,colnames(D),value=TRUE)
        ## drop any future occurrences
        z_t <- z_t[as.integer(lapply(
          strsplit(z_t,split=paste0(z,timesep)),"[",2)) <= tt]
        z_t <- sapply(z_t, function(zz) paste0(zz,":",z_ap))
        colnames(z_t) <- NULL
        return(z_t)
      }
    }))
    
    # fit outcome model
    tt.out <- as.formula(paste0(Rt,"~",paste(c(tt.predictors,tt.ZP),collapse="+")))
    out.fit <- lm(tt.out,data=D)
    psi.hat <- coef(out.fit)[grepl(A_t,names(coef(out.fit)))]
    res.coef[[tt+1]] <- psi.hat
    
    # blip down outcome
    if (tt>0) {
      blip_down.t <- (model.matrix(out.fit)[,names(psi.hat),drop=FALSE] %*% 
                        psi.hat)[,1]
      if (nrow(D)>length(blip_down.t)) {
        blip_down.t.all <- rep(NA,nrow(D))
        blip_down.t.all[as.integer(names(blip_down.t))] <- blip_down.t
        blip_down.t <- blip_down.t.all
        rm(blip_down.t.all)
      }
      Rt_minus1 <- unlist(D[,Rt,with=FALSE]) - blip_down.t
      D[, paste0("Ry",timesep,tt) := Rt_minus1]
      rm(Rt_minus1,blip_down.t)
    }
    
    res.formulae[[tt+1]] <- list("PS"=tt.ps,"Outcome"=tt.out)
    res.fit[[tt+1]] <- list("PS"=ps.fit,"Outcome"=out.fit)
    rm(ps.fit,out.fit,psi.hat,tt.ps,tt.out)
  }
  return(list("est"=res.coef,"formulae"=res.formulae,"fitted"=res.fit))
}

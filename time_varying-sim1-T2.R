rm(list=ls())
library("data.table")
library("CBPS")
library("twang")
source("time_varying-helper_funs.R")

# simulation settings
simsettings <- expand.grid("n"=1e2*c(1,4,16), # sample sizes
                           "kc"=c(TRUE,FALSE) # effect of (hidden) covariates
)
simsettings

for (seed in 1:nrow(simsettings)) {

# generate data ###############################################################
simsettings[seed,]
(N <- simsettings[seed,"n"])   # sample size
(kc_A <- 0.7*simsettings[seed,"kc"]) # coefficients for unmeasured confounding
(kc_LY <- exp(0.7)*simsettings[seed,"kc"])

(nuL <- (1:6)/8)
T_ <- 2 # time point for end-of-study outcome
REP <- 1000L # number of sims

OneData <- function(n=10) {
  U <- rnorm(n)
  V <- rnorm(n)
  Data <- NULL
  for (tt in 0:T_) {
    if (tt==0) {
      # true covariates used to generate treatment and outcome
      L0 <- matrix(rnorm(n*6,mean=U),ncol=6)
      
      At.ast <- (cbind(V,L0) %*% c(kc_A,nuL))[,1]
      At.ast <- exp(At.ast)/(1+exp(At.ast))
      At <- rbinom(n,size=1,prob=At.ast)
      
      Yt.ast <- (cbind(U,L0) %*% c(kc_LY,exp(nuL)))[,1]
      Yt <- Yt.ast + rnorm(n)
      Yt <- scale(Yt)
      
      Lt <- L0
    } else {
      # previous time points
      Lt <- matrix(NA,nrow=n,ncol=6)
      Lt[,1:4] <- L0[,1:4]
      Lt[,5:6] <- as.matrix(Data[[tt]][,paste0("L.",5:6)])
      A_tminus1 <- Data[[tt]][,"A"]
      Y0 <- Data[[1]][,"Y"]
      
      if (tt<T_) {
        # true covariates used to generate treatment and outcome
        L5.ast <- (cbind(U,U*Y0,Lt) %*% c(kc_LY,kc_LY,exp(nuL)))[,1]
        L6.ast <- (cbind(U,U*Y0,Lt) %*% c(kc_LY,kc_LY,exp(nuL)))[,1]
        L5.ast <- L5.ast + rnorm(n)
        L6.ast <- L6.ast + rnorm(n)
        Lt[,5] <- scale(L5.ast) + 0.35*A_tminus1
        Lt[,6] <- scale(L6.ast) + 0.70*A_tminus1
        
        At.ast <- (cbind(V,Lt,Y0) %*% c(kc_A,nuL,0.7))[,1]
        At.ast <- exp(At.ast)/(1+exp(At.ast))
        At <- rbinom(n,size=1,prob=At.ast)
        
        Yt <- rep(NA,n)
      } else {
        At <- rep(NA,n)
        
        Yt.ast <- (cbind(U,Lt,Y0*exp(U/2),
                         ## arbitrary functions of the baseline covariates
                         sqrt(abs(apply(Lt[,1:4],1,prod))),
                         pmin((Lt[,4]*Y0)^2,qchisq(.95,df=1)),
                         ## direct effect
                         A_tminus1) %*%
                     c(kc_LY,exp(nuL[1:4]),1,0.5,kc_LY,1,1,0.7))[,1]
        Yt <- Yt.ast + rnorm(n)
        
        Lt <- matrix(NA,nrow=n,ncol=6)
      }
    }
    
    Data[[tt+1]] <- data.frame("id"=1:n,
                               "t"=tt,
                               "L"=Lt,
                               "A"=At,
                               "Y"=Yt)
    rm(Lt,At,Yt)
  }
  Data <- rbindlist(Data)
  setkey(Data)
  return(Data)
}

simres <- NULL
meths <- c("dr.mle","outcome.only","single.reg")
ptm <- proc.time()[3]
for (i in 1:REP) {
  set.seed(9000+i)
  
  Data <- OneData(n=N)
  
  simres.i <- matrix(NA,nrow=length(meths),ncol=T_)
  row.names(simres.i) <- meths
  
  res.mle <- OneEst(Data,T_,
                    Lnames=paste0("L.",1:6),Treat_name="A",Outcome_name="Y",
                    use.DR=TRUE,fit.ps="mle",return.data=TRUE)
  simres.i["dr.mle",] <- unlist(res.mle$est)
  
  res <- OneEst(Data,T_,Lnames=paste0("L.",1:6),Treat_name="A",Outcome_name="Y",
                use.DR=FALSE,return.data=TRUE)
  simres.i["outcome.only",] <- unlist(res$est)
  rm(res)
  
  # remove all PS from the single regression model for the outcome
  form.y <- res.mle$formulae[[2]]$Outcome
  form.y <- paste0(as.character(form.y)[2],"~",paste(
    grep("P_.",strsplit(as.character(form.y)[3],split=" [+] ")[[1]],
         value=TRUE,invert=TRUE),collapse="+"))
  fitY <- coef(lm(as.formula(form.y),data=res.mle$data.wide))
  simres.i["single.reg",] <- fitY[sort(grep("A_.",names(fitY),value=TRUE))]
  rm(form.y,fitY)
  
  simres[[i]] <- data.frame("meths"=row.names(simres.i),simres.i)
  
  cat(i,"|",round((proc.time()[3]-ptm)/60),"mins \n")
}

simres.dt <- rbindlist(simres)
simres.dt <- cbind(data.table(simsettings[seed,]),simres.dt)
setkey(simres.dt)
print(simres.dt[,lapply(.SD,mean),by=eval(key(simres.dt)[1:3])])

save(simres.dt, file=paste0("time_varying-sim1-T2_",seed,".Rdata"))
}
q()

# load results ################################################################
rm(list=ls())
library("data.table")
library("xtable")
files_to_load <- list.files()[grep(".Rdata",list.files())]
simres <- NULL
for (seed in 1:length(files_to_load)) {
  load(file=paste0("time_varying-sim1-T2_",seed,".Rdata"))
  simres[[seed]] <- simres.dt
  rm(simres.dt)
}
simres.dt <- rbindlist(simres)
setkey(simres.dt)

# number of simulations
simres.dt[, .N, by=eval(key(simres.dt)[1:3])][,unique(N)]

# check number of sims that converged
simres.dt[,lapply(.SD, function(x) sum(!is.na(x))),
          by=eval(key(simres.dt)[1:3]),.SDcols="X1"]

# difference with true value
simres.dt[, bias.X1 := X1 - 0.7]
simres.dt[, bias.X2 := X2 - 0.7]
setkey(simres.dt)

res.table <- NULL
res.table[["bias"]] <- simres.dt[,lapply(.SD, function(x) mean(x,na.rm=TRUE)),
                                 by=eval(key(simres.dt)[1:3]),
                                 .SDcols=c("bias.X1","bias.X2")]
res.table[["ese"]] <- simres.dt[,lapply(.SD, function(x) sd(x,na.rm=TRUE)),
                                by=eval(key(simres.dt)[1:3]),
                                .SDcols=c("X1","X2")]
setnames(res.table[["ese"]],c("X1","X2"),c("ese.X1","ese.X2"))
res.table[["mse"]] <- simres.dt[,lapply(.SD, function(x) sqrt(mean(x*x,na.rm=TRUE))),
                                by=eval(key(simres.dt)[1:3]),
                                .SDcols=c("bias.X1","bias.X2")]
setnames(res.table[["mse"]],c("bias.X1","bias.X2"),c("mse.X1","mse.X2"))
(res.table <- Reduce(function(...) merge(..., all = TRUE), res.table))

# relabel entries
res.table[, kc := kc!=0]
res.table[meths=="dr.mle", meths := "DR-MLE"]
res.table[meths=="outcome.only", meths := "Outcome only"]
res.table[meths=="single.reg", meths := "Single regression"]
setkey(res.table,kc)
print(xtable(res.table,digits=c(0,0,1,rep(2,7))),include.rownames=FALSE)


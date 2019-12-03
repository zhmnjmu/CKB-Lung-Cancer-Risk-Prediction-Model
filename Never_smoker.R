#!/usr/bin/Rscript

library("flexsurv")
library("rstpm2")
library("foreign")
load("model_for_neversmoker.rdata")
load("model_for_mor.rdata")

dat <- read.csv("data_for_neversmoker.csv",h=T)
dat$smoke <- as.factor(dat$smoke)

stpm2cif_nev <- function(model_nev,model_mor,newdata,X) {
    stpm2if <- function(X) {
        ## setup up newdata
        newdat <- data.frame(as.data.frame(lapply(newdata, rep, 21)))
		newdat$time=X
		newdat$timesurvival=X
		
		## calculate cause-specific survival
        Sk1 <- predict(model_nev, newdata=newdat)
		Sk2 <- predict(model_mor, newdata=newdat)
        S <- apply(cbind(Sk1,Sk2),1,prod) #±©Â¶·çÏÕ=Éú´æ·çÏÕ*Î´·¢²¡·çÏÕ#
        ## calculate cause-specific hazard
        S*predict(model_nev,newdata=newdat,type="hazard") #±©Â¶·çÏÕ*·¢²¡·çÏÕ#
    }
    sapply(X, function(xi)
           ifelse(xi==0,0,
                  integrate(stpm2if,0,xi)$value))
}

results_nev <- as.data.frame(matrix(NA,ncol=2,nrow=1))
colnames(results_nev) <- c("studyid","cumhaz")


results_nev$studyid="ID"
results_nev$cumhaz=stpm2cif_nev(model_nev,model_mor,dat[1,],10)	

write.csv(results_nev,"Post_estimation_of_neversmoker.csv",row.names=F)

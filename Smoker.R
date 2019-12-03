#!/usr/bin/Rscript

library("flexsurv")
library("rstpm2")
library("foreign")


load("model_for_smoker.rdata")
load("model_for_mor.rdata")

dat <- read.csv("data_for_smoker.csv",h=T)
dat$smoke <- as.factor(dat$smoke)

stpm2cif_sm <- function(model_sm,model_mor,newdata,X) {
    stpm2if <- function(X) {
        ## setup up newdata
        newdat <- data.frame(as.data.frame(lapply(newdata, rep, 21)))
		newdat$time=X
		newdat$timesurvival=X
		
		## calculate cause-specific survival
        Sk1 <- predict(model_sm, newdata=newdat)
		Sk2 <- predict(model_mor, newdata=newdat)
        S <- apply(cbind(Sk1,Sk2),1,prod) #±©Â¶·çÏÕ=Éú´æ·çÏÕ*Î´·¢²¡·çÏÕ#
        ## calculate cause-specific hazard
        S*predict(model_sm,newdata=newdat,type="hazard") #±©Â¶·çÏÕ*·¢²¡·çÏÕ#
    }
    sapply(X, function(xi)
           ifelse(xi==0,0,
                  integrate(stpm2if,0,xi)$value))
}

results_sm <- as.data.frame(matrix(NA,ncol=2,nrow=1))
colnames(results_sm) <- c("studyid","cumhaz")


results_sm$studyid="ID"
results_sm$cumhaz=stpm2cif_sm(model_sm,model_mor,dat[1,],10)	


write.csv(results_sm,"Post_estimation_of_smoker.csv",row.names=F)

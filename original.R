#This is the version of the HMM model that works with the constant paramters.
devMULTIEVENT <- function(b,data,eff,e, nh,km1){ 
	
		
# data= encounter histories, eff= counts/frequency
# e vector of dates of first captures
# garb vector of initial states 
# km1 number of recapture occasions 
# nh number individuals


# OBSERVATIONS (+1)
# 0 = non-detected
# 1 = not resighted/recovered, captured
# 2 = resighted not recovered, not captured
# 3 = resighted not recovered, captured
 # 4 = resighted and recovered (Note that these are assummed not possible in Barker's model)
 # 5 = not resighted, recovered
 
 # OBSERVATIONS for last time period
 # 6 = non-detected
 # 7 = resighted not recovered
 # 8 = not resighted, recovered
   
   
 # STATES
 # 1 = alive 
 # 2 = newly dead
 # 3 = dead
 
 # PARAMETERS
 # phi  survival prob
 # p    detection prob
 # r    recovery prob
 # R    resight prob given alive at next time (see Barker 1997, parameter R)
 # Rr   resight (and not recovered) prob given not alive at next time (see Barker 1997, parameter R')
 
 # logit link for all parameters
 # note: below, we decompose the state and obs process in two steps composed of binomial events, 
 # which makes the use of the logit link appealing; 
 # if not, a multinomial (aka generalised) logit link should be used
 phi <- 1/(1+exp(-b[1]))
 p <- 1/(1+exp(-b[2]))
 R <- 1/(1+exp(-b[3]))
 Rr <- 1/(1+exp(-b[4]))
 r <- 1/(1+exp(-b[5]))
 
 
 # prob of obs (rows) cond on states (col)
 P0 = matrix(c((1-R)*(1-p),0,0,0,(1-Rr)*(1-r),0,0,0,1),nrow=3,ncol=3,byrow=T)
 P1 = matrix(c((1-R)*p,0,0,0,0,0,0,0,0),nrow=3,ncol=3,byrow=T)
 P2 = matrix(c(R*(1-p),0,0,0,Rr*(1-r),0,0,0,0),nrow=3,ncol=3,byrow=T)
 P3 = matrix(c(R*p,0,0,0,0,0,0,0,0),nrow=3,ncol=3,byrow=T)
 P4 = matrix(c(0,0,0,0,Rr*r,0,0,0,0),nrow=3,ncol=3,byrow=T)
 P5 = matrix(c(0,0,0,0,(1-Rr)*r,0,0,0,0),nrow=3,ncol=3,byrow=T)
 P6 = matrix(c((1-R),0,0,0,(1-Rr)*(1-r),0,0,0,1),nrow=3,ncol=3,byrow=T)
 P7 = matrix(c(0,0,0,0,(1-Rr)*r,0,0,0,0),nrow=3,ncol=3,byrow=T)
 P8 = matrix(c(R,0,0,0,Rr*(1-r),0,0,0,0),nrow=3,ncol=3,byrow=T)
 
 
 P =array(rep(0,3*3*9), dim=c(3,3,9))
 P[,,1]=P0
 P[,,2]=P1
 P[,,3]=P2
 P[,,4]=P3
 P[,,5]=P4
 P[,,6]=P5
 P[,,7]=P6
 P[,,8]=P7
 P[,,9]=P8
 
 # prob of states at t+1 given states at t
 A <- matrix(c(phi,1-phi,0,0,0,1,0,0,1),nrow=3,ncol=3,byrow=T)
 
 # init states (all animals are alive at release and we condition on release)
 # c(prob(alive), prob(newly dead), prob(dead))
 PI <- c(1,0,0)
 
 # likelihood
    l <- 0
    for (i in 1:nh) # loop on ind
    {
       ei <- e[i] # date of first det
       evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
       ALPHA=PI
       # cond on first capture, might have to change km1 depending on how structure the histories, also might have to change from ei:km1 for loop.
       for (j in 1:km1) 
       {
  #     	if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!
         ALPHA <- (ALPHA %*% A)%*%P[,,evennt[j]]
       }
       l <- l + logprot(sum(ALPHA))*eff[i]
    }
     l <- -l
     l
 }
 
 # avoid explosion of log(v) for small values of v
 logprot <- function(v){
 eps <- 2.2204e-016
 u <- log(eps) * (1+vector(length=length(v)))
 index <- (v>eps)
 u[index] <- log(v[index])
 u
 }
 
 # read in data
 
 # Read in data as character string to be converted to a history matrix later
 oyster.inp=read.table(file="oystercatchers.txt", colClasses="character", header=T)
 sex<- as.numeric(unlist(oyster.inp$Male)) #sex where Male=1, Female=0
 
 # Read characters into history matrix where "times" is the frequency of that particular history
 LD.history <- matrix(as.numeric(unlist(strsplit(rep(oyster.inp$history,times=1), ""))),ncol=nchar(oyster.inp$history), byrow=TRUE)  
 
 nh= dim(LD.history)[1] #number of animals observed in experiment
 
 nsample= dim(LD.history)[2]/2 #number of sample times used in experiment
 nsample

 
 
 #The "first" matrix specifies the first sample time at which each animal was observed 
 first<- NULL
 for(i in 1:nobs){
    first<- c(first,min(which(LD.history[i,]!=0)))
 }
   
 # Create event history matrix 
 history=matrix(rep(0,nobs*nsample), nrow=nobs, ncol=nsample)
 nspot=nsample*2 #number of spots in an LiDi type history
 
 # Convert LiDi type event history with x=0,1,2,3,4,5,6,7,8 type history, conditional on first release. 
 # Thus actually looking at the DiLi+1 events and mapping these to x.
 # For example LiDi history of 10 12 01 would map to event history of 128  
 for(i in 1:nobs){
 	index=(first[i]+1)/2
 	j<-first[i]+1
 	while(j<nspot){
 		if(LD.history[i,j]==0 & LD.history[i,j+1]==0){history[i,index]=0 #not observed
 		}else if(LD.history[i,j]==0 & LD.history[i,j+1]==1){history[i,index]=1 #not resighted, recapture
 		}else if(LD.history[i,j]==2 & LD.history[i,j+1]==0){history[i,index]=2 #resight, not recaptured
         }else if(LD.history[i,j]==2 & LD.history[i,j+1]==1){history[i,index]=3 #resight, recapture
         }else if(LD.history[i,j]==1 & LD.history[i,j+1]==0){history[i,index]=5 #recovered
         }
         j<-j+2
         index=index+1 
         #print(c(j, index))
          	
 	 }#endwhile
 	 if(j==nspot){#Final time period only has resight and recovery possibility, thus x=6,7, or 8
 			if(LD.history[i,j]==0){history[i,index]<-6 #no observation
 				}else if(LD.history[i,j]==1){history[i,index]<- 7 #recovery
 					} else if(LD.history[i,j]==2){history[i,index]<-8} #resight
 		#print(c(i,j,index,LD.history[i,j], history[i,index]))
 		}
 }




 # define various quantities#
 km1 <- nsample
 
 # counts
 eff <- rep(1,nh)
 
 # compute the date of first capture fc, and state at initial capture init.state
 fc <- first
 
 # initial values
 binit <- runif(5)# where 5 is the number of parameters phi, p, R, Rr, r
 
 # transpose data
 data <- t(history)



 # fit model 
 deb=Sys.time()
 tmpmin <- optim(binit,devMULTIEVENT,NULL,hessian=TRUE,data,eff,fc,nh,km1,method="BFGS",control=list(trace=1, REPORT=1))

 fin=Sys.time()
 fin-deb 

 
 # get estimates and back-transform
 x <- tmpmin$par
 phi <- 1/(1+exp(-x[1]))
 p <- 1/(1+exp(-x[2]))
 R <- 1/(1+exp(-x[3]))
 Rr <- 1/(1+exp(-x[4]))
 r <- 1/(1+exp(-x[5]))
 
 
 
 # Get standard errors
 H=tmpmin$hessian
 
 # Calulate the variance via the delta method
 VAR<- function(design_matrix,Hess,stdparms)
 { library(MASS)
   inv_Hess=ginv(Hess,tol=sqrt(.Machine$double.eps));
   Cov_XB=design_matrix%*%inv_Hess%*%t(design_matrix);
   std.vect=stdparms*(1-stdparms);
   COV=diag(std.vect)%*%Cov_XB%*%diag(std.vect);
   VAR.vect=vector(mode="double");
   i=1;
   for(r in 1:nrow(design_matrix)){ VAR.vect[i]=COV[r,r]; i=i+1;}
   VAR.vect
 }
 # Create a design matrix, here it is the identity matrix for the constant model
 X = matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),nrow=5,ncol=5,byrow=T)
 stdparms=1/(1+exp(-x))
 std.VAR=VAR(X,H,stdparms)
 std.se=sqrt(std.VAR)
 
cat("Parameter estimates and estimated standard errors","\n", 
"phi=", round(phi, digits=3), "se=", round(std.se[1], digits=3),"\n",
"p  =", round(p, digits=3), "se=", round(std.se[2], digits=3),"\n",
"R  =", round(R, digits=3), "se=", round(std.se[3], digits=3),"\n",
"R' =", round(Rr, digits=3), "se=", round(std.se[4], digits=3),"\n",
"r  =", round(r, digits=3), "se=", round(std.se[5], digits=3),"\n")
Parameter estimates and estimated standard errors 
 phi= 0.942 se= 0.006 
 p  = 0.418 se= 0.012 
 R  = 0.088 se= 0.007 
 R' = 0.065 se= 0.043 
 r  = 0.196 se= 0.036 
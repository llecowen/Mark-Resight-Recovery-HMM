#' Simulate data from a mark-resight-recapture model for use with HMM

#' Data are in the form of observations from 0 to 5 and represent between and at sample time observations.  
#0 not resighted nor recovered between $t$ and $t+1$, not captured at sample time $t+1$ 
#1 not resighted nor recovered between $t$ and $t+1$, captured at sample time $t+1$  
#2 resighted not recovered between $t$ and $t+1$, not captured at sample time $t+1$
#3 resighted not recovered between $t$ and $t+1$, captured at sample time $t+1$
#4 resighted and recovered between $t$ and $t+1$
#5 not resighted but recovered between $t$ and $t+1$ 



# define various quantities
K = 3     # number of sample times
N = 1000       # number of individuals
p = 0.9   # detection prob
phi = 0.8 # survival prob
R = 0.5   # resight given survives prob
Rr = 0.5  # resight given dies prob
r  = 0.9  # recovery prob


# pre-allocate memory

alive <- array(dim = c(N, K)) # Latent alive states
y <- array(NA, dim = c(N, K-1)) # Detection histories
muZ=NA
capture<- array(dim=c(N,K)) # observed captures
resight<- array(dim=c(N,K)) # observed resights
recovery<- array(dim=c(N,K)) # observed recovery


# define state process
# first year/season
alive[,1] <- rbinom(N, 1, 1) # Initial occupancy state, all are alive at time 1

# subsequent years
for(i in 1:N){ # Loop over individual
	for(j in 2:K){ # Loop over years
		muZ <- alive[i, j-1]*phi  # Prob for alive		
		alive[i,j] <- rbinom(1, 1, muZ)
		}
}

# define observation processes

for(i in 1:N){
	capture[i,1]<-1
	recovery[i,1]<-0
	resight[i,1]<-0
	for(k in 2:K){
		prob <- alive[i,k] * p
		capture[i,k]<- rbinom(1,1,prob)
		prob2<- alive[i,k-1]*(1-alive[i,k])*r
		recovery[i,k] <- rbinom(1,1,prob2)
		prob3<-alive[i,k-1]*alive[i,k]*R + alive[i,k-1]*(1-alive[i,k])*Rr
		resight[i,k]<- rbinom(1,1,prob3)
		}
}



for(i in 1:N){
	for(k in 1:(K-1)){
		if (capture[i,k+1]==1){
			if(resight[i,k+1]==1){
				y[i,k]<-3
			} else { y[i,k]<-1
				}#endelse
		} else if (capture[i,k+1]==0){
			if(resight[i,k+1]==1){
				if(recovery[i,k+1]==1){
					y[i,k]<-4
				} else {
					y[i,k]<-2
					}#endelse
			} else if(resight[i,k+1]==0){ 
				if(recovery[i,k+1]==1){
					y[i,k]<-5
				} else {
					y[i,k]<-0
					}#endelse
				}#endresight
			}#endcapture
		}#endfor
}#endfor

#print out observations
data=y
#data=as.data.frame(y)

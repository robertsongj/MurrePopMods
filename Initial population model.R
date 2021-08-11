# projection model for murre to assess impacts of current [2015] harvest

# 18 August 2015

# approach
#       define parameters and parameterize
#       project population
#       assess impact of harvest


################DEFINE PARAMETERS########################

# Pre-breeding matrix, as in Wiese et al. 2004
# reduce the number of age classes to 6 just to make things a little more manageable
# start with dummy numbers just to get things started

N0 <- 62420

# build matrix
sj <- 0.55
sa <- 0.90

pb <- c(0,0,0,0.367,0.7,0.985)

m <- c(0,0,0,0.5, 0.6, 0.7)

F <- pb*m*0.5

F

# build matrix

# build sub-matrix of survival elements

S <- diag(5) * sa

S <- cbind(S, S[,5])
S

#and stick them together

A  <- rbind(F, S)
A

# basic projection tools in popbio

library(popbio)

pop.projection(A, rep(10, 6), 25)
reproductive.value(A)/sum(reproductive.value(A))
generation.time(A)
damping.ratio(A)
eigen.analysis(A)
matplot2(A)

# spread the population over N

N <- N0 / sum(stable.stage(A)[4:6]) * stable.stage(A)
N

A2<-splitA(A)

A$F
A$T


# define the oil and hunt mortality

# remember to 1/2 for females only

Noil <- 3000/2 * stable.stage(A)

hunt.vul <- c(5, 3, 2, 1,1,1)


hunt.vec <- (hunt.vul*stable.stage(A))/ sum(hunt.vul * stable.stage(A))

Nhunt <- 20000/2 * hunt.vec

#define the length of the projection

years <- 10

Nout <- rbind(N, matrix(0, nrow = years, ncol = dim(A2$T)[1]))

#guts of the projection


for (i in  1:years) {
  
  # project the fall flight
  
     Nhold <-  (A2$F*sqrt(sj) + sqrt(A2$T)) %*% Nout[i,]

     #kill step
     
     Nhold <- Nhold-Noil-Nhunt
     
  # transition the population until spring   
     
     Nout[i+1,] <-  sqrt(diag(c(sj, rep(sa, 5)))) %*% Nhold
}

Nout


#Heyde-Cohen for calculating lambda

lam <- exp((log(Nout[years+1, 6]) - log(Nout[1,6]))/ (years - 1))
lam







# double check projection is working


Nout <- rbind(N, matrix(0, nrow = years, ncol = dim(A2$T)[1]))




for (i in  1:years) {
  
  Nout[i+1,] <- (A2$F*sj + A2$T) %*% Nout[i,]
}

Nout

pop.projection(A, rep(10,6), 25)

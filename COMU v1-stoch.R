#########################
# COMU population model v1.0
# 3 regions - Lab, Northern Newfoundland and southern Newfoundland
# 
#
#########################




###############################################

# build projection matrix

sj <- 0.56      # based on Harris et al. 2007
sa <- 0.95      # based on Robertson et al. 2006 - mean of m and f

sesj <- 0.10 # approx based on Harris et al. 2007
sesa <- 0.02 # based on Robertson et al. 2006

# annual estimates from Harris et al. 2007, distribution quite flat
# and not happy with an se draw
sjsamples <- c(0.4,0.31,0.52,0.91,0.3,0.55,0.35,0.33,0.85,0.86,0.67,0.77,0.52,0.73,0.75,0.67,0.67,0.45,0.38)

pb <- c(0,0,0.025,0.367,0.7,0.985) # based on Wiese et al. 2004 for TBMU

#use m.adj to change breeding conditions
m.adj <- 0.5

m.point <- 0.643 * m.adj # based on Anne's data Great/Gull Island (1997-2014)
m <- c(0,0,0,rep(m.point,3))

sem <- 0.144 # from Anne's data for Great/Gull Island (1997-2014)

F <- pb*m*0.5*sj
F


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
#matplot2(A)

#build specific matrices for Af (fall flight) and Aw (rest of winter, S only)

Ff <- pb*m*0.5*sqrt(sj)

Af <- rbind(Ff, sqrt(S))
Af

Aw <- sqrt(diag(c(sj, rep(sa, 5))))
Aw

##############################################

#aportioning the harvest

###############################################

# vectors run north to south

# change units from individuals to females /2

# assume even sex ratio

Nbreed <-c(73354, 855048, 530000)/2
breed_locs <-c("Lab","N_NL","S_NL")
names(Nbreed) <- breed_locs

Nbreed

# harvest statistics use 2010-2014 with annual species-specific corrections

total_harvest <- c(6993, 16401, 7829, 47119, 44213)/2
xTH <- mean(total_harvest)
summary(total_harvest)

# allocate harvest to breeding regions

## two approaches 

## 1 - relative use of nearshore habitats from McFarlane Tranquilla

regional_vul_habitats <- c(0.24, 0.04, 0.11)/sum(c(0.24, 0.04, 0.11))
regional_vul_habitats

## 2 - relative direct recovery rates, data from 2005-2013

regional_vul_f<- c(0.0396, 0.0086,0.0071)/sum(c(0.0396, 0.0086,0.0071))
regional_vul_f

regional_vul <- rbind(regional_vul_habitats, regional_vul_f)


# allocate harvest to age classes

I_hunt <- c(2, 25, 24, 81, 87)
A_hunt <- c(12, 27, 33, 60, 120)
year <- seq(2010:2014)

prop_juv <- I_hunt/(I_hunt + A_hunt)
prop_juv

bin.m <- glm(cbind(I_hunt, A_hunt) ~ as.factor(year), family = binomial)
summary(bin.m)
bin.null <- glm(cbind(I_hunt,A_hunt) ~ 1, family = binomial)
summary(bin.null)
anova(bin.m, bin.null, test = "Chisq")
AIC(bin.m, bin.null) 
confint(bin.null)
exp(confint(bin.null))/(1 +exp(confint(bin.null)))
coef(bin.null)
mean_prop_juv <- exp(coef(bin.null))/(1+exp(coef(bin.null)))
sqrt(vcov(bin.null))

# grabbed from glm, could be simpler like mean of prop_juv vector

mean_prop_juv

# bring in a stable stage distribution

stable.stage(A)

get the total population size

Ntot <- Nbreed / sum(stable.stage(A)[4:6])

Ntot

Nsummer <- stable.stage(A)%o%Ntot
Nsummer

#create vector of regional vulnerabilities based on colony size

vul_breed <- (regional_vul*(Nbreed/(sum(Nbreed))))/(sum(regional_vul*(Nbreed/(sum(Nbreed)))))
vul_breed

#create matrix of age vulnerabilities, aportion prop.Juv to stage 1, rest to 
#based on age structure 

age_vul <- (stable.stage(A) / sum(stable.stage(A)[2:6])) * (1-mean_prop_juv)
age_vul[1] <- mean_prop_juv

# combine the age and colony specific vulner

tot_vul <- age_vul%o%vul_breed  
tot_vul
sum(tot_vul)

# assign the hunt to the vulnerabilities

Nhunt <- tot_vul * xTH
Nhunt

# now sort out the oil

total_oil <- c(2262, 3855, 3666, 931)/2
mean_oil <- mean(total_oil)
summary(total_oil)

# oiling is only driven by relative population size and stable.stage

Noil <- (stable.stage(A)%o%(Nbreed/sum(Nbreed))) * mean_oil
Noil

# wrapper to run X iterations of the randomization

maxiter <- 500

# set up lam to hold lamdba values

lam <- array(0, dim=c(maxiter, length(Nbreed)))

for (iter in 1:maxiter) {

#define the length of the projection

years <- 20

# set up Nout to handle projections

Nout <- array(0, dim = c(years+1,dim(A)[1], length(Nbreed)))
Nout[1,,1:3] <- Nsummer

#guts of the projection


for (i in  1:years) {
  
  # add some error elements of Af
  
  
  #########simple, just random normal error
  #Af.draw <- matrix(rnorm(Af, Af, sAf), nrow=6, ncol=6)
  #Af.draw
  
  ############more complex, bring in elements from proper draws
  
  m.draw <- c(rep(0,3), rep(rnorm(1, m.point, sem),3))
  m.draw 
  
  #sj.draw <- beta.draws(sj,sesj,1)
  sj.draw <- sample(sjsamples, 1)
  
  sa.draw <- beta.draws(sa,sesa,1)
  
  Ff <- pb*m.draw*0.5*sqrt(sj.draw)
  
  Af.draw <- rbind(Ff, sqrt(cbind(diag(5)*sa.draw, c(rep(0,4), sa.draw))))
  
  # project the fall flight
  
  Nfall <-  Af.draw %*% Nout[i,,1:3]
  
  # draw either the habitat or recovery regional_vul vector
  
  pick <- sample(1:2,1) 
  vul_breed <- (regional_vul[pick,]*(Nbreed/(sum(Nbreed))))/(sum(regional_vul[pick,]*(Nbreed/(sum(Nbreed)))))
  
  # draw proportion juvenile from data
  
  pick2 <- sample(2:5,1)
  age_vul <- (stable.stage(A) / sum(stable.stage(A)[2:6])) * (1-prop_juv[pick2])
  age_vul[1] <- prop_juv[pick2]
  age_vul
  
  # gather up the vulnerabilities
  
  tot_vul <- age_vul%o%vul_breed  
  
  # draw a random value from total_harvest
  
  Nhunt <- tot_vul * sample(total_harvest,1)
  
  #draw a random value from Total_oil
  
  Noil <- (stable.stage(A)%o%(Nbreed/sum(Nbreed))) * sample(total_oil,1)
  
  #kill step
  
  Nwinter <- Nfall - Noil - Nhunt
  
  # transition the population until spring   
  
  #########simple, just random normal error
  #Aw.draw <- matrix(rnorm(Aw, Aw, Aw*0.01), nrow=6, ncol=6)
  
  ############ more complex, bring in elements from proper draws
  
  
  #sj.draw <- beta.draws(sj,sesj,1)
  sj.draw <- sample(sjsamples, 1)
  
  sa.draw <- beta.draws(sa,sesa,1)
    
  Aw.draw <- sqrt(diag(c(sj.draw, rep(sa.draw, 5))))  
  Aw.draw
  
  Nout[i+1,,1:3] <-  Aw.draw %*% Nwinter
}

#Nout

#Heyde-Cohen for calculating lambda

lam[iter,] <- exp((log(Nout[years+1, 6,1:3]) - log(Nout[1,6,1:3]))/ (years - 1))

}

lam
colnames(lam) <- breed_locs

library(lattice)


par(mfcol=c(3,1))
hist(lam[,1])
hist(lam[,2])
hist(lam[,3])


setwd("C:\\Users\\robertsong\\Documents\\seabirds\\murres\\harvest\\2015 assessment")

tiff("COMU-poor.tiff", width = 3.2, height = 4.8, units = 'in', res = 300, pointsize = 8)

par(mai=c(1.1,0.6,0.2,0.1))

par(mfcol=c(1,1))
ticklabels <- c("Labrador","Northern NL", "Southern NL")


boxplot(lam[,1:3], xlab = "", ylab = 
          expression(lambda), names = ticklabels, line = 0, cex.lab = 1.5, 
        cex.axis = 1.5, las =3)
abline(h=1, lty =2 )
abline(h=eigen.analysis(A)$lambda1, lty =3 )

dev.off()

apply(lam, 2, mean)
apply(lam, 2, sd)      



stats <- apply(lam, 2, quantile)

quants <- array(0, dim = c(7,3))
for (n  in 1:3) {
  x <- quantile(lam[,n], c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  quants[,n] <- quantile(lam[,n], c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
}

rownames(quants) <- names(x)
colnames(quants) <- breed_locs
quants

# function to transform mean and SD into beta shape parameters and return a vector of randomized draws

beta.draws <- function(p, sdp, draws) {
  
  varp <- sdp^2
  alpha <- p*(((p*(1-p))/varp) - 1)
  beta <- (1-p)*(((p*(1-p))/varp) - 1)
  rbeta(draws,alpha,beta)
  
}


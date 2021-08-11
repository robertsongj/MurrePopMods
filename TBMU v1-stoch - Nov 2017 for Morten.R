#########################
# TBMU population model v1.0
# 7 regions - Arctic Canada, Atlantic Canada, Hubson Bay, NW Greenland
# SW Greenland, Iceland, Spitsbergen
# 
# not included, Bjornoya - only 304 birds of almost 700K - assumed no impact
# August 2015
# Gregory J Robertson - Wildlife Research Division, ECCC
#########################

#########################
# additional annotation made November 2017 in advance of sending code to Morten
#
# general approach follows Wiese et al. (2004 - Biol Cons) and Gilliland et al. (2009 - Wild. Biol.)
#
# basically populations are set up in breeding areas and are projected to mid-winter
# at which point numbers are removed based on their representation in various regions with
# specific kills of oiling (Canada only) and harvest (in Canada and Greenland).
# The populations are then projected back to the spring with a transition matrix
# (seasonal survival rates from winter to spring along the main diagonal)
#
# Wiese et al. 2004 used 11 stages/ages in projection matrix, reduced to 6 for this exercise
#
# based on a pre-breeding matrix projection
#
# values used in projection matrix are quite high, and deliberately so.
# Basically represents a population of murres producing at R-naught (around 1.07). 
# So no additional impacts, and no density-dependence. 
# Harvest and oiling kills removed directly, while R-naught can be ramped
# down via adjustments (e.g. m.adj) or additional data from colonies not producing very well
#
########################

########################
#
# common abbreviations
# s - generally survival rates
# m - fecundity
# A - projection matrices
# se - proceeding variable name, generally is a standard error
# f - following a variable name, generally means fall
# w - following a variable name, generally means winter
# j - following a variable name, generally means juvenile
# a - following a variable name, generally means adult
# pb - proportion breeding
# m - reproductive output (fledglings/pair)
# F - Fecundity in the sense of the first row of A
# S - survival elements of projection matrix A

# other variable names, mainly vectors and matrices related to numbers killed and vulnerabilities
# to that kill are self-explanatory (I hope)

# in general values for Canada do not have further subscripts
# Greenland specific values are annotated _Green
# 
########################

########################
# little function to transform mean and SD into beta shape parameters and return a vector of randomized draws
# needed for the randomization process in the projections

beta.draws <- function(p, sdp, draws) {
  
  varp <- sdp^2
  alpha <- p*(((p*(1-p))/varp) - 1)
  beta <- (1-p)*(((p*(1-p))/varp) - 1)
  rbeta(draws,alpha,beta)
  
}
########################

# build projection matrix


# annual survival rates
sj <- 0.56      # based on Harris et al. 2007
sa <- 0.95      # based on Wiese et al. 2004 - sa with oil and hunt removed

sesj <- 0.10 # approx based on Harris et al. 2007 - not used in main runs
sesa <- 0.033 # based on Smith and Gaston 2012 - 0.033 is process variance on 0.90

# annual survival estimates from Harris et al. 2007, distribution quite flat
# and not happy with an se draw
# of course these are common murre data, so be great to use something better if we have it for TBMU

sjsamples <- c(0.4,0.31,0.52,0.91,0.3,0.55,0.35,0.33,0.85,0.86,0.67,0.77,0.52,0.73,0.75,0.67,0.67,0.45,0.38)

#proportion breeding
pb <- c(0,0,0.025,0.367,0.7,0.985) # based on Wiese et al. 2004 for TBMU

# use m.adj to change breeding conditions, set to 1 for average conditions
# right now the fecundity values are quite high, use this gadget to lower production
# can replace this with real values from Europe

m.adj <- 1.0

# fecundity
m.point <- 0.69 * m.adj # from Smith and Gaston (2012)

# vectorize m, needed to create transition matrix A
m <- c(0,0,0,rep(m.point,3))

sem <- 0.088 # from Smith and Gaston (2012)

# create first row of transition matrix, the 0.5 is to limit to females
F <- pb*m*0.5*sj
F


# build sub-matrix of survival elements
S <- diag(5) * sa

S <- cbind(S, S[,5])
S

#and stick them together

A  <- rbind(F, S)
A

# projection matrix A complete


# basic projection tools in popbio - just to have a look at A

library(popbio)

pop.projection(A, rep(10, 6), 25)
reproductive.value(A)/sum(reproductive.value(A))
generation.time(A)
damping.ratio(A)
eigen.analysis(A)
#matplot2(A)

#build specific matrices for Af (fall flight) and Aw (rest of winter, S only)

# in the absense of more specific seasonal data, s is simply split by taking the sqrt

Ff <- pb*m*0.5*sqrt(sj)

Af <- rbind(Ff, sqrt(S))
Af

Aw <- sqrt(diag(c(sj, rep(sa, 5))))
Aw

##############################################

#apportioning the harvest

###############################################

# vectors run northwest to southeast - sort of

# change units from individuals to females /2

# assume even sex ratio

Nbreed <-c(1080000, 16352, 2000000, 625026, 27448, 660000, 1470000)/2
breed_locs <-c("Arctic_Can","Atl_Can","HB", "NW_Green", "SW_Green", "Iceland", "Spits")
names(Nbreed) <- breed_locs

Nbreed

# Canadian harvest statistics - use 2010-2014 with annual species-specific corrections

total_harvest <- c(45588, 49467, 44362, 66970, 63101)/2
xTH <- mean(total_harvest)
summary(total_harvest)

# harvest statistics for Greeland - should be double checked by Morten 
total_harvest_Green <- c(89305, 84409, 64660, 62843, 64468, 66935)/2

# allocate harvest and oiling to breeding regions

## relative use of Newfoundland shelf waters from 
## Frederiksen et al. 2016

regional_vul <- c(176938, 8229, 274662, 189046, 1751, 0, 40902)/sum(c(176938, 8229, 274662, 189046, 1751, 0, 40902))
regional_vul

## and same for Greenland

regional_vul_Green <- c(86229+2173, 0, 8124+387, 47019+1212, 3673+6612, 218810+42844, 156172+238596)/sum(c(86229+2173, 0, 8124+387, 47019+1212, 3673+6612, 218810+42844, 156172+238596))
regional_vul_Green

#simple calculation to look at how many birds from each breeding region are killed by Canada

(total_harvest%*%t(regional_vul))*2

#simple calculation to look at how many birds from each breeding region are killed by Greenland

(total_harvest_Green%*%t(regional_vul_Green))*2


# allocate harvest to age classes - based on Canadian Wing bee data, 2010-2014
# I in this case is immature (could be a J for juvenile - wing bee data separates first year from
# all else - called 'adults')

I_hunt <- c(70, 97, 155, 133, 187)
A_hunt <- c(21, 60, 168, 67,108)
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

# now some juvenile ratio data for Greenland - this will need input from Morten et al.!

prop_juv_Green <- 0.6
se_prop_juv_Green <- 0.05

# bring in a stable stage distribution

stable.stage(A)

# to get the total population size. We only have breeding population size, so need to re-construct 
# the younger pre-breeding cohorts - uses the stable stage distribution in the absense of better 
# information

Ntot <- Nbreed / sum(stable.stage(A)[4:6])

Ntot

#set-up starting population vector - this is the population going into the breeding season

Nsummer <- stable.stage(A)%o%Ntot
Nsummer

#create vector of regional vulnerabilities based on colony size - taken from COMU coding

#vul_breed <- (regional_vul*(Nbreed/(sum(Nbreed))))/(sum(regional_vul*(Nbreed/(sum(Nbreed)))))

# but no need for heroics here for TBMU, work done by Morten, and inputted directly from Appendix F

vul_breed <- regional_vul
vul_breed

#create matrix of age-specific harvest vulnerabilities, apportion prop.Juv to stage 1, rest allocated 
#based on age structure
# this is an oversimplification to a degree, we have some data showing that 2nd and 3rd years appears
# more vulnerable to harvest - but poorly quantified

age_vul <- (stable.stage(A) / sum(stable.stage(A)[2:6])) * (1-mean_prop_juv)
age_vul[1] <- mean_prop_juv

# combine the age and colony specific vulnerability vectors to create a total vulnerability matrix

tot_vul <- age_vul%o%vul_breed  
tot_vul
sum(tot_vul) # should sum to 1

# assign the hunt to the vulnerabilities

Nhunt <- tot_vul * xTH
Nhunt

# now sort out the oil - from Robertson et al. 2014 - 37th AMOP proceedings

total_oil <- c(74719, 60131, 7695, 15124)/2
mean_oil <- mean(total_oil)
summary(total_oil)

# COMU - oiling is only driven by relative population size and stable.stage
# TBMU - oiling is driven by the same vulnerabilities as the hunt

Noil <- (stable.stage(A)%o%(regional_vul)) * mean_oil
Noil

# look at the oiling numbers per region

colSums(Noil)*2

# wrapper to run X iterations of the randomization

maxiter <- 500

# set up lam to hold lamdba values

lam <- array(0, dim=c(maxiter, length(Nbreed)))

# begin iterations of the projections

for (iter in 1:maxiter) {

#define the length of each projection

years <- 20

# set up Nout to handle projections
# Nout stores the age-, breeding population- and annual-specific populations in a 3D array
# where rows are years of projection, columns are the 6 age classes, and depth are the 7
# breeding regions

Nout <- array(0, dim = c(years+1,dim(A)[1], length(Nbreed)))

# set the first row of each matrix to the initial population size

Nout[1,,1:length(Nbreed)] <- Nsummer

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
  
  Nfall <-  Af.draw %*% Nout[i,,1:length(Nbreed)]
  
  # establish regional_vul vector
   
  vul_breed <- regional_vul
  
  # draw proportion juvenile from data - for Canada one of 5 years
  
  pick2 <- sample(1:5,1)
  age_vul <- (stable.stage(A) / sum(stable.stage(A)[2:6])) * (1-prop_juv[pick2])
  age_vul[1] <- prop_juv[pick2]
  age_vul
  
  # gather up the vulnerabilities
  
  tot_vul <- age_vul%o%vul_breed  
  
  # draw a random value from total_harvest
  
  Nhunt <- tot_vul * sample(total_harvest,1)
  
  # now do it all again for Greenland
  
  vul_breed_Green <- regional_vul_Green
  
  prop_juv_Green_draw <- beta.draws(prop_juv_Green, se_prop_juv_Green, 1)
  
  age_vul_Green <- (stable.stage(A) / sum(stable.stage(A)[2:6])) * (1-prop_juv_Green_draw)
  age_vul_Green[1] <- prop_juv_Green_draw
  age_vul_Green

  tot_vul_Green <- age_vul_Green%o%vul_breed_Green
  
  Nhunt_Green <- tot_vul_Green * sample(total_harvest_Green,1)
  
  #draw a random value from Total_oil
  
  Noil <- (stable.stage(A)%o%(regional_vul)) * sample(total_oil,1)
  
  #kill step - simple subtraction - can remove one or more of these to assess impacts of each
  
  Nwinter <- Nfall   - Nhunt - Noil - Nhunt_Green 
  
  # transition the population until spring   
  
  #########simple, just random normal error
  #Aw.draw <- matrix(rnorm(Aw, Aw, Aw*0.01), nrow=6, ncol=6)
  
  ############ more complex, bring in elements from proper draws
  
  #sj.draw <- beta.draws(sj,sesj,1)
  sj.draw <- sample(sjsamples, 1)
  
  sa.draw <- beta.draws(sa,sesa,1)
    
  Aw.draw <- sqrt(diag(c(sj.draw, rep(sa.draw, 5))))  
  Aw.draw
  
  # write the new population size for each year to the storage array Nout
  Nout[i+1,,1:length(Nbreed)] <-  Aw.draw %*% Nwinter
}

Nout

#Heyde-Cohen for calculating lambda

lam[iter,] <- exp((log(Nout[years+1, 6,1:length(Nbreed)]) - log(Nout[1,6,1:length(Nbreed)]))/ (years - 1))

}

# lam contains lambda values for each breeding region, 
# with as many rows as simluations requested in maxiter

colnames(lam) <- breed_locs
lam

# the rest that follows are basic summarization tools for presenting the stochastic lambdas.


library(lattice)


par(mfcol=c(4,2))
hist(lam[,1])
hist(lam[,2])
hist(lam[,3])
hist(lam[,4])
hist(lam[,5])
hist(lam[,6])
hist(lam[,7])

setwd("C:\\Users\\robertsong\\Documents\\seabirds\\murres\\harvest\\2015 assessment")

tiff("TBMU-poor.tiff", width = 4.8, height = 4.8, units = 'in', res = 300, pointsize = 8)


par(mfcol=c(1,1))

ticklabels <- c("Arctic Canada","Atlantic Canada", "Hudson Bay", "NW Greenland", 
                 "SW Greenland", "Iceland", "Spitsbergen")

par(mai=c(1.35,0.6,0.2,0.1))


boxplot(lam[,1:length(Nbreed)], ylab = 
          expression(lambda), names = ticklabels, line = 0, 
          cex.axis = 1.5, las =3, cex.lab = 1.5,
          pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))
abline(h=1, lty =2 )
abline(h=eigen.analysis(A)$lambda1, lty =3 )


dev.off()


apply(lam, 2, mean)
apply(lam, 2, sd)      
stats <- apply(lam, 2, quantile)

quants <- array(0, dim = c(7,length(Nbreed)))
for (n  in 1:length(Nbreed)) {
  x <- quantile(lam[,n], c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  quants[,n] <- quantile(lam[,n], c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
}

rownames(quants) <- names(x)
colnames(quants) <- breed_locs
quants



#########################
# routine to aportion harvest and oil to common murres
# 19 August 2015
#########################


# vectors run north to south

# change units from individuals to females /2

# assume even sex ratio

Nbreed <-c(73354, 855048, 530000)/2
names(Nbreed) <- c("Lab","N_NL","S_NL")

Nbreed

# harvest statistics use 2010-2014 with annual species-specific corrections

total_harvest <- c(6993, 16401, 7829, 47119, 44213)/2
xTH <- mean(total_harvest)
summary(total_harvest)

# allocate harvest to breeding regions

## two approaches 

## 1 - relative use of nearshore habitats from McFarlane Tranquilla

regional_vul <- c(0.24, 0.04, 0.11)/sum(c(0.24, 0.04, 0.11))
regional_vul

## 2 - relative direct recovery rates, data from 2005-2013

regional_vul <- c(0.0396, 0.0086,0.0071)/sum(c(0.0396, 0.0086,0.0071))
regional_vul


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

Nfall <- stable.stage(A)%o%Ntot
Nfall

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

# oiling is only driven by population size and stable.stage

Noil <- (stable.stage(A)%o%(Nbreed/sum(Nbreed))) * mean_oil
Noil

Nwinter <- Nfall - Nhunt - Noil
Nwinter



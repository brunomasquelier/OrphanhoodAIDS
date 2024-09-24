# R script for generating microsimulations and developing new adjustments for AIDS-related biases
# Bruno Masquelier - 2024
rm(list = ls()) 
set.seed(3)

# Packages required
library(gdata) 
library(RColorBrewer, warn.conflicts = F)
library(boot)
library(RJSONIO)
library(WDI)
library(survival)
library(zoo)

# Basic colors
palette_G = rep(brewer.pal(9,"Greens")[3:9], 10)
palette_B = rep(brewer.pal(9,"Blues")[3:9], 10)
palette_R = rep(brewer.pal(9,"Reds")[3:9], 10)
palette_Gr = rep(brewer.pal(9,"Greys")[3:9], 10)
palette_Pp = rep(brewer.pal(9,"Purples")[3:9], 10)
palette_O = rep(brewer.pal(9,"Oranges")[3:9], 10)

# Basic parameters
max_age = 100
start_seg =  c(1800, seq(2000, 2045, by = 5)); length_seg = c(diff(start_seg), (start_seg[length(start_seg)] - start_seg[length(start_seg) - 1])); mid_seg = start_seg + length_seg/2
end_seg = start_seg + length_seg; nb_segments = length(start_seg)
age_ref = c(0,1:4, seq(5,120, by = 1))*12
age_ref5 = c(0,1, seq(5,120, by = 5))*12

#________________________________________________________________________________________________________     
# 1 - PREPARATION OF THE SIMULATIONS
#________________________________________________________________________________________________________     

# Mortality parameters (Brass logits, General Standard)
alpha_mort = c(0, -0.4, -0.8)# c(0.2, -0.2, -0.6, -1.0)
beta_mort =  c(0.8, 1.1)
n_mort = length(alpha_mort)*length(beta_mort)

mort_param <- expand.grid(alpha_mort=alpha_mort, beta_mort=beta_mort)
mort_param$n_mort = 1:n_mort

# Fertility parameters (Booth standard)
alpha_fert =  c(-0.35, 0.25)#c(-0.5, -0.2,0.1,0.4)
beta_fert =  c(0.85,1.15) #c( 0.7, 1.0, 1.3    )
n_fert = length(alpha_fert)*length(beta_fert)

# Growth rate
r = 0.02

fert_param <- expand.grid(alpha_fert=alpha_fert, beta_fert=beta_fert)
fert_param$n_fert = 1:n_fert

# HIV parameters
alpha_HIV = c(5,7)
beta_HIV =  c(3, 5)
H_HIV =  c(0.08, 0.13)
n_HIV = length(alpha_HIV)*length(beta_HIV)*length(H_HIV)

HIV_param <- expand.grid(alpha_HIV=alpha_HIV, beta_HIV=beta_HIV, H_HIV = H_HIV)
HIV_param$n_HIV = 1:n_HIV

# ART/PMTCT parameters
coverage_ART =  c(0, 1, 2) # no coverage, low coverage (73% in the end), high coverage (86% at the end)
n_ART = length(coverage_ART)
ART_param <- expand.grid(coverage_ART=coverage_ART)
ART_param$n_ART = 1:n_ART

all_param <- expand.grid(n_fert=1:n_fert, n_mort=1:n_mort, n_HIV=1:n_HIV, n_ART= 1:n_ART)
all_param = merge(all_param, fert_param, by = "n_fert")
all_param = merge(all_param, mort_param, by = "n_mort")
all_param = merge(all_param, HIV_param, by = "n_HIV")
all_param = merge(all_param, ART_param, by = "n_ART")
all_param = all_param[order(all_param$n_mort, all_param$n_fert, all_param$n_HIV, all_param$n_ART),]
all_param$no_sim = 1:nrow(all_param)
nrow(all_param)

print_param = rbind(
  cbind('Background mortality', ""),  
  cbind('Alpha', paste(alpha_mort, collapse = ', ')),
  cbind('Beta', paste(beta_mort, collapse = ', ')),
  cbind('Fertility (seronegative women)', ""),  
  cbind('Alpha', paste(alpha_fert, collapse = ', ')),
  cbind('Beta', paste(beta_fert, collapse = ', ')),
  cbind('Incidence curve', ""),
  cbind('Alpha', paste(alpha_HIV, collapse = ', ')),
  cbind('Beta', paste(beta_HIV, collapse = ', ')),
  cbind('H', paste(H_HIV, collapse = ', ')),
  cbind('ART', ""),
  cbind('Max coverage', paste(coverage_ART, collapse = ', ')))

    library(xtable)
    param.latex <- print(xtable(print_param, caption = "Parameters used to set up the microsimulations", label = "param.latex"), include.rownames = F, include.colnames = F, hline.after = c(1,4, 7))

# number of simulations 
n_sim = n_mort * n_fert * n_HIV * n_ART

# Mortality and fertility age patterns
logitsBrass <- read.xls(paste0(getwd(), "\\Data\\Brass General Standard to SOCSIM.xls"), sheet = 1)
Boothstandard <- read.xls(paste0(getwd(), "\\Data\\Booth Standard to SOCSIM.xls"), sheet = 1)

F_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(F_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 
M_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(M_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 

# Ages
F_life_tables[,13,,] = c(0:(dim(F_life_tables)[1] - 1)) 
M_life_tables[,13,,] = c(0:(dim(M_life_tables)[1] - 1)) 

# lx
F_life_tables[1,6,,] = 10000; M_life_tables[1,6,,] = 10000 
for(i in 2:(nrow(F_life_tables))) {
  for(j in 1:dim(F_life_tables)[4]) {
    F_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,2]))))[i-1]*10000
    M_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,2]))))[i-1]*10000
  }}

# dx
for(i in 1:(nrow(F_life_tables)-1)) {
  F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
  M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
}

F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]

# px
F_life_tables[,1,,] = (F_life_tables[,6,,] - F_life_tables[,7,,]) / F_life_tables[,6,,]
M_life_tables[,1,,] = (M_life_tables[,6,,] - M_life_tables[,7,,]) / M_life_tables[,6,,]

# Smoothing at older ages
for(i in 1:dim(F_life_tables)[3]) {
  for(j in 1:dim(F_life_tables)[4]) {
    F_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(F_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
    M_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(M_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
  }}

# qx
F_life_tables[,2,,] = 1-F_life_tables[,1,,]
M_life_tables[,2,,] = 1-M_life_tables[,1,,]

# revise table based on smooth probs at older ages
M_life_tables[,6,i,j] = M_life_tables[1,6,i,j] * cumprod(c(1,1-M_life_tables[,2,i,j]))[1:101]
F_life_tables[,6,i,j] = F_life_tables[1,6,i,j] * cumprod(c(1,1-F_life_tables[,2,i,j]))[1:101]

for(i in 1:(nrow(F_life_tables)-1)) {
  F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
  M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
}

F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]

# Age groups
F_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)
M_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)

# nax
for(j in 1:n_mort){
  for(i in 1:101){
    F_life_tables[i,3,,j] = sum(0.5:11.5/12 * -diff(c(1, (1-F_life_tables[i,2,1,j])^(1:12/12))))/sum(-diff(c(1, (1-F_life_tables[i,2,1,j])^(1:12/12))))
  }
}
M_life_tables[,3,,] = F_life_tables[,3,,]
stopifnot(unlist(M_life_tables[,2,,] == F_life_tables[,2,,]))

# n
F_life_tables[,4,,] = 1; M_life_tables[,4,,] = 1

# Lx
for(i in 1:(nrow(F_life_tables)-1)) {
  F_life_tables[i,8,,] = F_life_tables[i,4,,]*F_life_tables[i+1,6,,] + F_life_tables[i,3,,]*F_life_tables[i,7,,]
  M_life_tables[i,8,,] = M_life_tables[i,4,,]*M_life_tables[i+1,6,,] + M_life_tables[i,3,,]*M_life_tables[i,7,,]
}

F_life_tables[nrow(F_life_tables),8,,] = F_life_tables[nrow(F_life_tables),3,,]*F_life_tables[nrow(F_life_tables),7,,]
M_life_tables[nrow(M_life_tables),8,,] = M_life_tables[nrow(M_life_tables),3,,]*M_life_tables[nrow(M_life_tables),7,,]

# Tx
for(i in 1:dim(F_life_tables)[3]) {
  for(j in 1:dim(F_life_tables)[4]) {
    F_life_tables[1,9,i,j] = sum(F_life_tables[,8,i,j]) 
    M_life_tables[1,9,i,j] = sum(M_life_tables[,8,i,j]) 
  }}

for(i in 1:dim(F_life_tables)[3]) {
  for(j in 1:dim(F_life_tables)[4]) {
    F_life_tables[2,9,i,j] = sum(F_life_tables[,8,i,j]) - F_life_tables[1,8,i,j]
    M_life_tables[2,9,i,j] = sum(M_life_tables[,8,i,j]) - M_life_tables[1,8,i,j]
  }}

for(i in 3:(nrow(F_life_tables))) {
  for(j in 1:dim(F_life_tables)[4]) {
    F_life_tables[i,9,,j] = apply(F_life_tables[,8,,j], c(2), sum) - apply(F_life_tables[1:(i-1),8,,j], c(2), sum)			
    M_life_tables[i,9,,j] = apply(M_life_tables[,8,,j], c(2), sum) - apply(M_life_tables[1:(i-1),8,,j], c(2), sum)
  }}

# ex
F_life_tables[,10,,] = F_life_tables[,9,,]/F_life_tables[,6,,]
M_life_tables[,10,,] = M_life_tables[,9,,]/M_life_tables[,6,,]

# e0
round(range(M_life_tables[1,10,1,]),1); round(mean(M_life_tables[1,10,1,]))
# e15
round(range(M_life_tables[16,10,1,])); round(mean(M_life_tables[16,10,1,]), 1)

# F/TF = exp(-exp(-Y)) with Y = alpha + beta*Ys
for(i in 1:dim(F_life_tables)[3]) {
  F_life_tables[12:50,15,i,]  =  exp(-exp(-(fert_param[i,1] + fert_param[i,2] * Boothstandard[,2])))
}
F_life_tables[12:49,15,,]  =  F_life_tables[13:50,15,,1] - F_life_tables[12:49,15,,1]
F_life_tables[50,15,,] = 1 - F_life_tables[50,15,,]

# Mean age at childbearing
for(i in 1:n_fert) {		
  fert_param$MAC[i] = weighted.mean(F_life_tables[,13,i,1]+0.5, F_life_tables[,15,i,1])
}

F_life_tablesplot = F_life_tables

# Crude death rate
both = rep(NA, length(1:dim(F_life_tables)[4]))
for(i in 1:dim(F_life_tables)[4]) {
  both[i] = 1/sum(exp(-r * (F_life_tables[,13,1,i] + 0.5)) * (F_life_tables[,8,1,i] + F_life_tables[,8,1,i])/(F_life_tables[1,6,1,i] + F_life_tables[1,6,1,i]))
}

for(k in 1:dim(F_life_tables)[4]) {
  F_life_tables[,12,,k] = (both[k]) * exp(-r * (F_life_tables[,13,,k] + 0.5)) * F_life_tables[,8,,k]/F_life_tables[1,6,,k]
  M_life_tables[,12,,k] = (both[k]) * exp(-r * (M_life_tables[,13,,k] + 0.5)) * M_life_tables[,8,,k]/M_life_tables[1,6,,k]
}

# Number of births (when ISF == 1)
for(i in 1:dim(F_life_tables)[4]) {
  F_life_tables[,16,,i] =  F_life_tables[,15,,i] *F_life_tables[,12,,i]
}

for(i in 1:n_sim) {		
  all_param$M[i] = weighted.mean(F_life_tables[,13,all_param$n_fert[i],all_param$n_mort[i]]+0.5, F_life_tables[,16,all_param$n_fert[i],all_param$n_mort[i]])
}

round(range( all_param$M), 1); round(mean( all_param$M), 1)				

# Proportions of orphans during the stable period calculated from stable population theory
start_group = seq(5, 45, 5)
num_to_int = as.list(1:n_mort)
denom_to_int = as.list(1:n_mort)
Sx = as.list(1:n_mort)

for(x in 1:n_mort) {
  num_to_int[[x]] = as.list(1:n_fert)
  denom_to_int[[x]] = as.list(1:n_fert)
  Sx[[x]] = as.list(1:n_fert)
  
  for(y in 1:n_fert) {
    
    num_to_int[[x]][[y]] = rep(NA, length(start_group))
    denom_to_int[[x]][[y]] = rep(NA, length(start_group))
    for(k in 1:9) {	
      num_to_int[[x]][[y]][k] = integrate(function(t) { 
        sapply(t, function(t) {
          exp(-r*t)* (F_life_tables[t + 1,6,y,x]/10000)*integrate(function(a) (exp(-r*(a+0.5)) * ((F_life_tables[a + 1 + t,6,y,x] + F_life_tables[a + 2 + t,6,y,x])/20000) * F_life_tables[a + 1,15,y,x]), 11, 49, subdivisions = 1000, rel.tol = 0.01)$value
        })}, start_group[k], start_group[k]+5, rel.tol = 0.01)$value
      denom_to_int[[x]][[y]][k] = integrate(function(t) { 
        sapply(t, function(t) {
          exp(-r*t)* (F_life_tables[t + 1,6,y,x]/10000)*integrate(function(a) (exp(-r*(a+0.5)) * ((F_life_tables[a + 1    ,6,y,x] + F_life_tables[a + 2    ,6,y,x])/20000) * F_life_tables[a + 1,15,y,x]), 11, 49, subdivisions = 1000, rel.tol = 0.01)$value
        })}, start_group[k], start_group[k] +5, rel.tol = 0.01)$value	   
    }}}

for(x in 1:n_mort) {
  for(y in 1:n_fert) {
    Sx[[x]][[y]] = data.frame(num_to_int[[x]][[y]]/denom_to_int[[x]][[y]], x, y, start_group+5)
    names(Sx[[x]][[y]]) = c("Sn_5", "n_mort", "n_fert", "n")
    
    Sx[[x]][[y]]$l2 = F_life_tables[3,6,y,x]/F_life_tables[1,6,y,x]
    # n = start_group + 5
    for(z in 1:9) {
      Sx[[x]][[y]]$np25[z] = F_life_tables[26+5+ start_group[z],6,y,x]/F_life_tables[26,6,y,x]
      Sx[[x]][[y]]$l25n[z] = F_life_tables[26+5+ start_group[z],6,y,x]/F_life_tables[1,6,y,x]
      
    }}
  Sx[[x]] = do.call(rbind, Sx[[x]])
}

Sx = data.frame(do.call(rbind, Sx))

# Underlying rates - life exp, 5q0, 35q15 and age-specific nqx
for(i in 1:n_sim) {		
  all_param$e0[i] = F_life_tables["0",'ex',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q5[i] = 1-F_life_tables["5",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["0",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q35_15[i] = 1-F_life_tables["50",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["15",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  
  all_param$q10_25[i] = 1-F_life_tables["35",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q15_25[i] = 1-F_life_tables["40",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q20_25[i] = 1-F_life_tables["45",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q25_25[i] = 1-F_life_tables["50",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q30_25[i] = 1-F_life_tables["55",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q35_25[i] = 1-F_life_tables["60",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q40_25[i] = 1-F_life_tables["65",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q45_25[i] = 1-F_life_tables["70",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
  all_param$q50_25[i] = 1-F_life_tables["75",'lx',all_param$n_fert[i],all_param$n_mort[i]]/F_life_tables["25",'lx',all_param$n_fert[i],all_param$n_mort[i]]
   }

Sx = merge(Sx, all_param, by = c("n_fert", "n_mort"))

# Coefficients developed by Timaeus in 1992
coef_TM92 = rbind(c(- 0.2894, 0.00125,  1.2559),c(- 0.1718, 0.00222,  1.1123),	c(- 0.1513, 0.00372,  1.0525),		c(- 0.1808, 0.00586,  1.0267),		c(- 0.2511, 0.00885,  1.0219),
                  c(- 0.3644, 0.01287,  1.0380),		c(- 0.5181, 0.01795,  1.0753),		c(- 0.6880,  0.02342,  1.1276),		c(- 0.8054, 0.02721,  1.1678))
rownames(coef_TM92) = seq(10, 50, 5)
# compute the orphanhood-based nqx
for(k in 1:9){
  n = seq(10, 50, 5)[k]
  Sx$qn_25orph[Sx$n == n] = 1-(coef_TM92[k,1]+coef_TM92[k,2]*Sx$M[Sx$n == n]+coef_TM92[k,3]*Sx$Sn_5[Sx$n == n])
}

# Compute ASFRs - ISF is defined by the mortality rates and the age structure of fertility
# Stability over 150 years
for(k in 1:dim(F_life_tables)[4]) {
  stopifnot(F_life_tables[,"1ca",,k] == (both[k]) * exp(-r * (F_life_tables[,"age",,k] + 0.5)) * F_life_tables[,"Lx",,k]/F_life_tables[1,"lx",,k])
}

for(k in 1:dim(F_life_tables)[4]) {
  for(j in 1:dim(F_life_tables)[3]) {
    up = both[k]/sum(F_life_tables[,"fx",j,k]*F_life_tables[,"1ca",j,k]/2)
    F_life_tables[,"fx",j,k] =  F_life_tables[,"fx",j,k]*up
  }}

check =matrix(NA,  dim(F_life_tables)[3], dim(F_life_tables)[4])
for(k in 1:dim(F_life_tables)[4]) {
  for(j in 1:dim(F_life_tables)[3]) {
    check[j,k] =  sum(exp(-r * (F_life_tables[12:50,"age",j,k] + 0.5)) * (F_life_tables[12:50,"Lx",j,k]/F_life_tables[1,"lx",j,k])  *  F_life_tables[12:50,"fx",j,k]/2)
  }}
stopifnot(unique(round(check, 9) == 1))

# Generate incidence curves, ART and PMTCT trends
incidence_curve = matrix(NA, ncol = length(start_seg)-1, nrow = n_sim)
rownames(incidence_curve) = 1:n_sim
colnames(incidence_curve) = start_seg[2:length(start_seg)]

plot(200,200, type = 'l', ylim = c(0, 6), xlim = c(200, 250), ylab = 'Incidence (15-49)', xlab = 'Years')
for(i in 1:n_sim){#which(!duplicated(all_param$n_HIV))){
  gf = function(t) (((t^(all_param$alpha_HIV[i]-1))*exp(-(t)/all_param$beta_HIV[i]))/(factorial(all_param$alpha_HIV[i]-1)*all_param$beta_HIV[i]^all_param$alpha_HIV[i]))
  intgf = function(y) integrate(gf, lower = y-5, upper = y)$value
  points(mid_seg[2:length(mid_seg)]-1800, (1-(exp(-sapply(seq(5, 50, 5), FUN = intgf)*all_param$H_HIV[i])))*100, type = 'l', col = palette_G[all_param$n_HIV[i]], lwd = 2)
  #abline(v = ((all_param$alpha_HIV[i] - 1)*all_param$beta_HIV[i])+2000)
  incidence_curve[i,] = (1-(exp(-sapply(seq(5, 50, 5), FUN = intgf)*all_param$H_HIV[i])))
}

# Trends in ART and PMTCT from UNAIDS (both sexes)
art = WDI(indicator='SH.HIV.ARTC.ZS', start=1990, end=2020, extra = TRUE); art = art[order(art$year),]
pmtct = WDI(indicator='SH.HIV.PMTC.ZS', start=1990, end=2020, extra = TRUE); pmtct = pmtct[order(pmtct$year),]
incidence = WDI(indicator='SH.HIV.INCD.ZS', start=1990, end=2020, extra = TRUE); incidence = incidence[order(incidence$year),]

country = c('Botswana', 'Cameroon', 'Central African Republic', "Cote d'Ivoire", "Eswatini", "Kenya", 'Lesotho', 'Malawi', 'Mozambique', 'Namibia', 'Rwanda', 'South Africa', 'Tanzania', 'Uganda', 'Zambia', "Zimbabwe")
for(i in country[country %in% "Mozambique" == FALSE]){
par(mfrow = c(1,2))
plot(incidence[incidence$country == i,]$year, incidence[incidence$country == i,]$SH.HIV.INCD.ZS, type = 'l', col = "grey50", xlab = 'Years', ylab = "Incidence", xlim = c(1990, 2025), main = paste("Incidence - ", i))
plot(art[art$country == i,]$year, art[art$country == i,]$SH.HIV.ARTC.ZS/100, type = 'l', col = "grey50", ylim = c(0,1), xlab = 'Years', ylab = "ART coverage", xlim = c(1990, 2025), main = "ART")
}

library(car)
par(mfrow = c(1,2))
i = country[1]
plot(art[art$country == i,]$year, art[art$country == i,]$SH.HIV.ARTC.ZS/100, type = 'l', col = "grey50", ylim = c(0,1), xlab = 'Years', ylab = "ART coverage", xlim = c(2000, 2025), main = "ART")
countrylow = art[art$country %in% country & art$year == 2020 & art$SH.HIV.ARTC.ZS < 83 & !is.na(art$SH.HIV.ARTC.ZS),]$country
countryhigh = art[art$country %in% country & art$year == 2020 & art$SH.HIV.ARTC.ZS >= 83 & !is.na(art$SH.HIV.ARTC.ZS),]$country
for(i in countrylow){
  points(art[art$country == i,]$year, art[art$country == i,]$SH.HIV.ARTC.ZS/100, type = 'l', col = "#e34a33")
}
for(i in countryhigh){
  points(art[art$country == i,]$year, art[art$country == i,]$SH.HIV.ARTC.ZS/100, type = 'l', col = "#31a354")
}
datalow = art[art$country %in% countrylow & art$year >= 2000,]
initcoef = coef(lm(car::logit(SH.HIV.ARTC.ZS/100)~year,data=datalow))
nlslow<-nls(SH.HIV.ARTC.ZS~phi1/(1+exp(-(phi2+phi3*year))),
            start=list(phi1=100,phi2=initcoef[1],phi3=initcoef[2]),data=datalow,trace=TRUE)
phi1<-coef(nlslow)[1]
phi2<-coef(nlslow)[2]
phi3<-coef(nlslow)[3]
x<-c(2000:2025) #construct a range of x values bounded by the datalow
y<-phi1/(1+exp(-(phi2+phi3*x))) #predicted mass
predictlow<-data.frame(x,y)
points(predictlow$x, predictlow$y/100, type = 'l', lwd = 3, col = "#e34a33")
artpredictlow = predictlow

datahigh = art[art$country %in% countryhigh  & art$year >= 2000,]
initcoef = coef(lm(car::logit(SH.HIV.ARTC.ZS/100)~year,data=datahigh))
nlshigh<-nls(SH.HIV.ARTC.ZS~phi1/(1+exp(-(phi2+phi3*year))),
             start=list(phi1=100,phi2=initcoef[1],phi3=initcoef[2]),data=datahigh,trace=TRUE)
phi1<-coef(nlshigh)[1]
phi2<-coef(nlshigh)[2]
phi3<-coef(nlshigh)[3]
x<-c(2000:2025) #construct a range of x values bounded by the datahigh
y<-phi1/(1+exp(-(phi2+phi3*x))) #predicted mass
predicthigh<-data.frame(x,y)
points(predicthigh$x, predicthigh$y/100, type = 'l', lwd = 3, col = "#31a354")
artpredicthigh = predicthigh

legend('topleft', c("Rapid scale-up scenario", "Slow scale-up scenario"), col = c("#31a354", "#e34a33"), lty = 1, lwd = 2, horiz = F, bty = "n")

#pmtct
i = country[1]
plot(pmtct[pmtct$country == i,]$year, pmtct[pmtct$country == i,]$SH.HIV.PMTC.ZS/100, type = 'l', col = "grey50", ylim = c(0,1), xlab = 'Years', ylab = "PMTCT coverage", xlim = c(2000, 2025), main = "PMTCT")
for(i in countrylow){
  points(pmtct[pmtct$country == i,]$year, pmtct[pmtct$country == i,]$SH.HIV.PMTC.ZS/100, type = 'l', col = "#e34a33")
}
for(i in countryhigh){
  points(pmtct[pmtct$country == i,]$year, pmtct[pmtct$country == i,]$SH.HIV.PMTC.ZS/100, type = 'l', col = "#31a354")
}
datalow = pmtct[pmtct$country %in% countrylow,]
initcoef = coef(lm(car::logit(SH.HIV.PMTC.ZS/100)~year,data=datalow))
nlslow<-nls(SH.HIV.PMTC.ZS~phi1/(1+exp(-(phi2+phi3*year))),
            start=list(phi1=100,phi2=initcoef[1],phi3=initcoef[2]),data=datalow,trace=TRUE)
phi1<-coef(nlslow)[1]
phi2<-coef(nlslow)[2]
phi3<-coef(nlslow)[3]
x<-c(2000:2025) #construct a range of x values bounded by the datalow
y<-phi1/(1+exp(-(phi2+phi3*x))) #predicted mass
predictlow<-data.frame(x,y)
points(predictlow$x, predictlow$y/100, type = 'l', lwd = 3, col = "#e34a33")
pmtctpredictlow = predictlow

datahigh = pmtct[pmtct$country %in% countryhigh,]
initcoef = coef(lm(car::logit(SH.HIV.PMTC.ZS/100)~year,data=datahigh))
nlshigh<-nls(SH.HIV.PMTC.ZS~phi1/(1+exp(-(phi2+phi3*year))),
             start=list(phi1=100,phi2=initcoef[1],phi3=initcoef[2]),data=datahigh,trace=TRUE)
phi1<-coef(nlshigh)[1]
phi2<-coef(nlshigh)[2]
phi3<-coef(nlshigh)[3]
x<-c(2000:2025) #construct a range of x values bounded by the datahigh
y<-phi1/(1+exp(-(phi2+phi3*x))) #predicted mass
predicthigh<-data.frame(x,y)
points(predicthigh$x, predicthigh$y/100, type = 'l', lwd = 3, col = "#31a354")
pmtctpredicthigh = predicthigh
legend('bottomright', c("Rapid scale-up scenario", "Slow scale-up scenario"), col = c("#31a354", "#e34a33"), lty = 1, lwd = 2, horiz = F, bty = "n")

art_curve = matrix(0, ncol = length(start_seg)-1, nrow = n_sim)
rownames(art_curve) = 1:n_sim
colnames(art_curve) = start_seg[2:length(start_seg)]
for(i in 1:n_sim){
  if(all_param$coverage_ART[i] == 1){
  art_curve[i,6:10] = round(approx(artpredictlow$x, artpredictlow$y, xout = c(2002.5+seq(0, 20,5)))$y)
  }
  if(all_param$coverage_ART[i] == 2){
    art_curve[i,6:10] = round(approx(artpredicthigh$x, artpredicthigh$y, xout = c(2002.5+seq(0, 20,5)))$y)
  }
}

pmtct_curve = matrix(0, ncol = length(start_seg)-1, nrow = n_sim)
rownames(pmtct_curve) = 1:n_sim
colnames(pmtct_curve) = start_seg[2:length(start_seg)]
for(i in 1:n_sim){
  if(all_param$coverage_ART[i] == 1){
    pmtct_curve[i,6:10] = round(approx(pmtctpredictlow$x, pmtctpredictlow$y, xout = c(2002.5+seq(0, 20,5)))$y)
  }
  if(all_param$coverage_ART[i] == 2){
    pmtct_curve[i,6:10] = round(approx(pmtctpredicthigh$x, pmtctpredicthigh$y, xout = c(2002.5+seq(0, 20,5)))$y)
  }
}

# Prepare stable populations for SOCSIM
# Stable population at the start: 458 to have 50000*2 at the end 
starting_pop =  round((50000/exp(r*200))/2, 0)

pop = as.list(1:nrow(all_param))

for(i in 1:nrow(all_param)) {
  pop[[i]] = cbind(append(round(starting_pop * M_life_tables[,12,all_param$n_fert[i],all_param$n_mort[i]], 0), round(starting_pop * F_life_tables[,12,all_param$n_fert[i],all_param$n_mort[i]], 0)), c(rep(0,101), rep(1,101)), c(rep(c(1), 202)), rep(3,202), rep(0:100,2), rep(99,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(1,202), rep(0,202), rep(0,202))
  pop[[i]][,6] = (100 - pop[[i]][,5])*12
  pop[[i]] = pop[[i]][,-c(5)]
  pop[[i]] = as.data.frame(cbind(c(1:sum(pop[[i]][,1])),pop[[i]][rep(1:dim(pop[[i]])[1], times = pop[[i]][,1]),c(-1)][order(pop[[i]][rep(1:dim(pop[[i]])[1], times = pop[[i]][,1]),c(-1)][,4]),]))
  rownames(pop[[i]]) = 1:nrow(pop[[i]])
  colnames(pop[[i]]) = c("person_id", "sex", "group", "next", "month_d_o_b", "mother", "father", "e_s_mom", "e_s_dad", "lborn", "last_marr", "mstatus", "month_d_o_d", "fmult")
} 

# Prepare HIV-free life tables for SOCSIM
taux_TM = cbind(
  c(c(0, 1), 2:100),
  c(c(1,0), rep(0, 99)))
birth_f = cbind(c(9:51, 100), rep(0, 44))

# Number of groups
nb_groups = 10
# 1= Sero-negatives
# 2= HIV-positive in incubation phase
# 3 = full blown aids
# 4= terminal phase (instant death)
# 5= Vertically-infected children (5:8), with 5 being perinatal and 6-8 postnatal infection
# 9= Under ART
# 10= Vertical infection avoided (instantly sent to 1)

# Mortality rates (by sex) for seronegatives: RATES_M_M and RATES_M_F are combined in MORT_SEG
RATES_M_M = as.list(1:nrow(all_param)); RATES_M_F = as.list(1:nrow(all_param)); MORT_SEG = as.list(1:nrow(all_param))

for(i in 1:nrow(all_param)) {
  RATES_M_M[[i]] = cbind(taux_TM[2:nrow(taux_TM),], format(rbind(matrix(data = c(#all_param$nmrFq[i], 1-((1-all_param$pnmrFq[i])^(1/11)), 
    round(1-((1-M_life_tables[1:99,"nqx",all_param$n_fert[i],all_param$n_mort[i]])^(1/12)),12), 1-(0.01^(1/12))), nrow = (max_age), ncol = 1))[,1], scientific = F))
  RATES_M_F[[i]] = cbind(taux_TM[2:nrow(taux_TM),], format(rbind(matrix(data = c(#all_param$nmrFq[i], 1-((1-all_param$pnmrFq[i])^(1/11)), 
    round(1-((1-F_life_tables[1:99,"nqx",all_param$n_fert[i],all_param$n_mort[i]])^(1/12)),12), 1-(0.01^(1/12))), nrow = (max_age), ncol = 1))[,1], scientific = F))
  
  MORT_SEG[[i]]= as.list(1:nb_groups)
  # HIV-free life tables
  for(j in c(1,2,3, 5:8,9, 10)) {	# including those under ART (9)
    MORT_SEG[[i]][[j]] =rbind(
      cbind("death ", j, " M single"),
      RATES_M_M[[i]],
      cbind("death ", j, " F single"),
      RATES_M_F[[i]])
  }
  # Terminal phase: instant death
  MORT_SEG[[i]][[4]] =rbind(
    cbind("death ", 4, " M single"),
    matrix(c(10, 100, 0,  0, 0.99, 0.99), ncol = 3),
    cbind("death ", 4, " F single"),
    matrix(c(10, 100, 0,  0, 0.99, 0.99), ncol = 3)
  )
}		

# Fertility rates
# With reduction among seropositive mothers
FERT_SEGHIV = as.list(1:nrow(all_param))		
for(i in 1:nrow(all_param)) {
  FERT_SEGHIV[[i]]= as.list(1:8)
  FERT_SEGHIV[[i]][[1]] = rbind(c("birth   ", 1, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEGHIV[[i]][[2]] = rbind(c("birth   ", 2, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                                c("birth   ", 2, paste("F single", 1, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                                c("birth   ", 2, paste("F single", 2, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                                c("birth   ", 2, paste("F single", 3, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                                c("birth   ", 2, paste("F single", 4, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEGHIV[[i]][[3]] = rbind(c("birth   ", 3, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1.26, 12), rep(0.76, 5), rep(0.67, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEGHIV[[i]][[4]] = rbind(c("birth   ", 4, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                                c("birth   ", 4, paste("F single", 1, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                                c("birth   ", 4, paste("F single", 2, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                                c("birth   ", 4, paste("F single", 3, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                                c("birth   ", 4, paste("F single", 4, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEGHIV[[i]][[5]] = rbind(c("birth   ", 5, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEGHIV[[i]][[6]] = rbind(c("birth   ", 6, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEGHIV[[i]][[7]] = rbind(c("birth   ", 7, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEGHIV[[i]][[8]] = rbind(c("birth   ", 8, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEGHIV[[i]][[9]] = rbind(c("birth   ", 9, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEGHIV[[i]][[10]] = rbind(c("birth   ", 10, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
}
# Revised fertility : 1.26 for the 15-19 years-old, 0.76 for the 20-24 years-old, 0.67 later on (Lewis, 2004)

# Without reduction among seropositive mothers
FERT_SEG = as.list(1:nrow(all_param))		
for(i in 1:nrow(all_param)) {
  FERT_SEG[[i]]= as.list(1:8)
  FERT_SEG[[i]][[1]] = rbind(c("birth   ", 1, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEG[[i]][[2]] = rbind(c("birth   ", 2, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                             c("birth   ", 2, paste("F single", 1, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                             c("birth   ", 2, paste("F single", 2, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                             c("birth   ", 2, paste("F single", 3, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)),
                             c("birth   ", 2, paste("F single", 4, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEG[[i]][[3]] = rbind(c("birth   ", 3, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(c(rep(1, 12), rep(1, 5), rep(1, 26))*F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEG[[i]][[4]] = rbind(c("birth   ", 4, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                             c("birth   ", 4, paste("F single", 1, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                             c("birth   ", 4, paste("F single", 2, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                             c("birth   ", 4, paste("F single", 3, sep = ' ')), birth_f[c(1,44),c(1,2,2)],
                             c("birth   ", 4, paste("F single", 4, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEG[[i]][[5]] = rbind(c("birth   ", 5, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEG[[i]][[6]] = rbind(c("birth   ", 6, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEG[[i]][[7]] = rbind(c("birth   ", 7, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEG[[i]][[8]] = rbind(c("birth   ", 8, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
  FERT_SEG[[i]][[9]] = rbind(c("birth   ", 9, paste("F single", 0, sep = ' ')), cbind(birth_f, format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F))) 
  FERT_SEG[[i]][[10]] = rbind(c("birth   ", 10, paste("F single", 0, sep = ' ')), birth_f[c(1,44),c(1,2,2)])
}

# No marriage rates for any group at any time
NUPT_SEG =      rbind(
  cbind("marriage ", 1, " M single"),
  cbind(taux_TM[c(11,101),1:2], c( 0, 0)),
  cbind("marriage ", 1, " F single"),
  cbind(taux_TM[c(11,101),1:2], c(0, 0)))

# Prepare transition rates

# INCUBATION PERIOD (2 to 3)
# Duration specific, with weibull distribution
# Spectrum 2009 = @TECHREPORT{Stover2009,  author = {Stover, J},  title = {AIM: A Computer Program for Making HIV/AIDS Projections and Examining the Demographic and Social Impacts of AIDS}, institution = {Futures Group International, Health Policy Initiative, Task Order 1.}, year = {2009}, address = {Washington, DC}, owner = {masquelier}, timestamp = {2009.12.10}}
# African pattern (CD4 count < 200)
incub = array(0, dim = c(100,15,2))
dimnames(incub)[[2]] = c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", as.character(1)) 
dimnames(incub)[[3]] = c('males', 'females')
incub[,'lx','males'] = 1- pweibull(0:(99), shape  = 1.71, scale = 10.2)
incub[,'lx','females'] = 1- pweibull(0:(99), shape  = 1.67, scale = 10.65)
incub[,'age',] = c(0:(dim(incub)[1] - 1)) 
for(i in 1:(dim(incub)[1]-1)) { incub[i,'dx',] = incub[i,'lx',] - incub[i+1,'lx',]}
incub[dim(incub)[1],'dx',] = incub[dim(incub)[1],'lx',]
incub[,'npx',] = (incub[,'lx',] - incub[,'dx',]) / incub[,'lx',] 
incub[,'nqx',] = 1-incub[,'npx',]
# median :  10.2*(log(2))^(1/1.71)
# Males median should be around 8
# median :  10.65*(log(2))^(1/1.67)

INCUB = rbind(
  # From 2 (infected) to 3 (full blown aids)
  cbind("duration_specific transit ", 2, paste(" M single", 3)),
  cbind("transit ", 2, paste(" M single", 3)),
  cbind(taux_TM[c(2:35, 101),1:2], round(c(1-incub[1:34,'npx','males']^(1/12),  1-(0.01^(1/12))), 6)),
  cbind("duration_specific transit ", 2, paste(" F single", 3)),
  cbind("transit ", 2, paste(" F single", 3)),
  cbind(taux_TM[c(2:35, 101),1:2], round(c(1-incub[1:34,'npx','females']^(1/12),  1-(0.01^(1/12))), 6)))

# FULL BLOWN AIDS (3 to 4)
# Duration specific, with weibull distribution
# Spectrum 2009 = @TECHREPORT{Stover2009,  author = {Stover, J},  title = {AIM: A Computer Program for Making HIV/AIDS Projections and Examining the Demographic and Social Impacts of AIDS}, institution = {Futures Group International, Health Policy Initiative, Task Order 1.}, year = {2009}, address = {Washington, DC}, owner = {masquelier}, timestamp = {2009.12.10}}
fullblownaids = array(0, dim = c(100,15,2))
dimnames(fullblownaids)[[2]] = c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", as.character(1)) 
dimnames(fullblownaids)[[3]] = c('males', 'females')
fullblownaids[,'lx','males'] = 1- pweibull(0:(99), shape  = 1, scale = 2.16)
fullblownaids[,'lx','females'] = 1- pweibull(0:(99), shape  = 1, scale = 2.82)
fullblownaids[,'age',] = c(0:(dim(fullblownaids)[1] - 1)) 
for(i in 1:(dim(fullblownaids)[1]-1)) { fullblownaids[i,'dx',] = fullblownaids[i,'lx',] - fullblownaids[i+1,'lx',]}
fullblownaids[dim(fullblownaids)[1],'dx',] = fullblownaids[dim(fullblownaids)[1],'lx',]
fullblownaids[,'npx',] = (fullblownaids[,'lx',] - fullblownaids[,'dx',]) / fullblownaids[,'lx',] 
fullblownaids[,'nqx',] = 1-fullblownaids[,'npx',]
# median :  2.82*(log(2))^(1/1)

AIDS = rbind( 
  # From 3 (full blown aids) to 4 (death)
  cbind("duration_specific transit ", 3, paste(" M single", 4)),
  cbind("transit ", 3, paste(" M single", 4)),
  cbind(taux_TM[c(2:35, 101),1:2], round(c(1-fullblownaids[c(1:34),'npx','males']^(1/12),  1-(0.01^(1/12))), 6)),
  cbind("duration_specific transit ", 3, paste(" F single", 4)),
  cbind("transit ", 3, paste(" F single", 4)),
  cbind(taux_TM[c(2:35, 101),1:2], round(c(1-fullblownaids[c(1:34),'npx','females']^(1/12),  1-(0.01^(1/12))), 6)))

# Survival under ART
  surv_art = array(0, dim = c(100,15,1))
  dimnames(surv_art)[[2]] = c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", as.character(1)) 
  
  surv_art[,'lx',1] = c(1, cumprod(c(0.85, rep(0.95, 98))))
  surv_art[,'age',] = c(0:(dim(surv_art)[1] - 1)) 
  for(i in 1:(dim(surv_art)[1]-1)) { surv_art[i,'dx',] = surv_art[i,'lx',] - surv_art[i+1,'lx',]}
  surv_art[dim(surv_art)[1],'dx',] = surv_art[dim(surv_art)[1],'lx',]
  surv_art[,'npx',] = (surv_art[,'lx',] - surv_art[,'dx',]) / surv_art[,'lx',]; surv_art[,'nqx',] = 1-surv_art[,'npx',]
  for(i in 2:(dim(surv_art)[1]-1)) { surv_art[i,'nax',] = ((-1/24)*surv_art[i-1,'dx',] + 0.5*surv_art[i,'dx',] + (1/24)*surv_art[i+1,'dx',])/surv_art[i,'dx',]}
  for(i in c(1,(dim(surv_art)[1]-10):dim(surv_art)[1])) {surv_art[i,'nax',] = 1/2	}; surv_art[,'nx',] = 1
  surv_art[,'tx',] = surv_art[,'nqx',] / (surv_art[,'nx',]-surv_art[,'nqx',]*(surv_art[,'nx',]-surv_art[,'nax',]))
  surv_art[,'tx',][is.na(surv_art[,'tx',]) | surv_art[,'tx',] >= 1] = 0.99 
  
  par(mfrow = c(1,1))
  plot(0:30,  surv_art[1:31,'lx',], type="l", ylim = c(0, 1), xlim = c(0, 31), lwd = 3, col = palette_G[3], ylab = "S(x)", xlab = "Time spent under ART", cex.lab = 2, cex.axis = 2) ## standard exponential distribution
  
  ART_TO_DEATH = rbind(c("duration_specific ", "transit 9", "M single", 4), c("transit", 9, "F single", 4), cbind(surv_art[,'age',]+1,0, round(surv_art[,'tx',1]/12, 7), ""))
  ART_TO_DEATH = ART_TO_DEATH[rep(1:nrow( ART_TO_DEATH), 6),]
  ART_TO_DEATH[grep("single", ART_TO_DEATH[,3]),][,3] =  rep(paste(rep(c("M", "F"), each = 3), rep(c("single", "married", "widowed"))), each = 2)
  ART_TO_DEATH[,3] = paste(ART_TO_DEATH[,3], ART_TO_DEATH[,4], sep = ' '); ART_TO_DEATH = ART_TO_DEATH[,1:3]
  # ART_TO_DEATH[ART_TO_DEATH[,1] %in% c("transit","duration_specific "),]

# Parameters for no vertical transmission
INFCHILDNOVERT = rbind(     
  cbind("transit ", 2, paste(" M single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 2, paste(" F single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 4, paste(" M single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 4, paste(" F single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3),
  cbind("transit ", 9, paste(" F single", 1)),
  matrix(c(10, 100,  0 , 0, 0.99,   0), ncol = 3))

# Morality of children vertically infected (aids only) ("Updates to the Spectrum ... "- Sex Trans Infect 2012)
survchild = matrix(c(0.646, 1.336, 1.062, 0.058, 2.195, 0.44, 1.015, 1.484, 0.058, 2.2, 0.248, 1.241, 2.11, 0.058, 2.2, 0.048, 1.873, 1.708, 0.058, 2.2), ncol = 4)
rownames(survchild) = c('pi', 'lambda1', 'mu1', 'lambda2', 'mu2')
colnames(survchild) = c('peri', 'post0_180', 'post181_365', 'post365')

childHIV = array(0, dim = c(25,15,4))
dimnames(childHIV)[[2]] = c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", as.character(1)) 
dimnames(childHIV)[[3]] =   colnames(survchild)
for(k in 1:4){
  X = -1*(survchild['lambda1',k]*0:24)^(survchild['mu1',k])
  Y = -1*(survchild['lambda2',k]*0:24)^(survchild['mu2',k])
  childHIV[,'lx',colnames(survchild)[k]] = (survchild['pi',k]*exp(X) + (1-survchild['pi',k])*exp(Y))
}
childHIV[,'age',] = c(0:(dim(childHIV)[1] - 1)) 
for(i in 1:(dim(childHIV)[1]-1)) { childHIV[i,'dx',] = childHIV[i,'lx',] - childHIV[i+1,'lx',]}
childHIV[dim(childHIV)[1],'dx',] = childHIV[dim(childHIV)[1],'lx',]
childHIV[,'npx',] = (childHIV[,'lx',] - childHIV[,'dx',]) / childHIV[,'lx',] 
childHIV[,'nqx',] = 1-childHIV[,'npx',]
#plot(0:24, childHIV[,'lx',1], type = 'l')
#for(i in 2:4){points(0:24, childHIV[,'lx',i], type = 'l')}

SURVCHILD = rbind( 
  # From 5 (perinatal) to 4 (instant death)
  cbind("duration_specific transit ", 5, paste(" M single", 4)),
  cbind("transit ", 5, paste(" M single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','peri'])^(1/12), 6), 0.99)), ncol = 3),
  cbind("duration_specific transit ", 5, paste(" F single", 4)),
  cbind("transit ", 5, paste(" F single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','peri'])^(1/12), 6), 0.99)), ncol = 3),
  # From 6 (Postnatal 0-180) to 4 (death)
  cbind("duration_specific transit ", 6, paste(" M single", 4)),
  cbind("transit ", 6, paste(" M single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post0_180'])^(1/12), 6), 0.99)), ncol = 3),
  cbind("duration_specific transit ", 6, paste(" F single", 4)),
  cbind("transit ", 6, paste(" F single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post0_180'])^(1/12), 6), 0.99)), ncol = 3),
  # From 7 (Postnatal 181-365) to 4 (death)
  cbind("duration_specific transit ", 7, paste(" M single", 4)),
  cbind("transit ", 7, paste(" M single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post181_365'])^(1/12), 6), 0.99)), ncol = 3),
  cbind("duration_specific transit ", 7, paste(" F single", 4)),
  cbind("transit ", 7, paste(" F single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post181_365'])^(1/12), 6), 0.99)), ncol = 3),
  # From 7 (Postnatal 365+) to 4 (death)
  cbind("duration_specific transit ", 8, paste(" M single", 4)),
  cbind("transit ", 8, paste(" M single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post365'])^(1/12), 6), 0.99)), ncol = 3),
  cbind("duration_specific transit ", 8, paste(" F single", 4)),
  cbind("transit ", 8, paste(" F single", 4)),
  matrix(c(c(1:24, 100), rep(0, 25), c(round(1-(1-childHIV[1:24,'nqx','post365'])^(1/12), 6), 0.99)), ncol = 3))

# INFECTIONS VERTICALES
INFCHILDnoPMTCT = rbind(     
  # The probability of HIV infection at birth for a child born to an HIV-positive mother is assumed to be 20% in the absence of prophylaxis
  # From born from infected mother (2, 3 or 4) to infected child (5): 20% chance in utero -> fast progressors
  cbind("transit ", 2, paste(" M single", 5)),
  matrix(c(0, 100, 1,  0, 0.2/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
  cbind("transit ", 2, paste(" F single", 5)),
  matrix(c(0, 100, 1,  0, 0.2/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 5)),
  matrix(c(0, 100, 1,  0, 0.2/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 5)),
  matrix(c(0, 100, 1,  0, 0.2/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
  # From born from infected mother (2, 3 or 4) to non infected (1): everyone is transferred back to 1 after age 2
  cbind("transit ", 2, paste(" M single", 1)),
  matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
  cbind("transit ", 2, paste(" F single", 1)),
  matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 1)),
  matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 1)),
  matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
  # born to mother under ART - immediately sent back to 1, assumed to be under PMTCT
  cbind("transit ", 9, paste(" M single", 1)),
  matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
  cbind("transit ", 9, paste(" F single", 1)),
  matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
  # The probability of infection through breastfeeding is assumed to be 1.5% per month for mixed feeding during the first 6 months, 0.75% per month for exclusive feeding during the first 6 months, 0.75% per month for months 7 and later 
  # (Sex Trans Inf 2010 - The Spectrum projection package: improvements in estimating incidence by age and sex, mother-to-child transmission, HIV progression in children and double orphans)
  # To post-natal 0-180
  cbind("transit ", 2, paste(" M single", 6)),
  matrix(c(0,  100, 6, 0,  (0.679)*(0.015)+(0.313)*(0.0075), 0), ncol = 3),
  cbind("transit ", 2, paste(" M single", 6)),
  matrix(c(0,  100, 6, 0,  (0.679)*(0.015)+(0.313)*(0.0075), 0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 6)),
  matrix(c(0,  100, 6, 0,  (0.679)*(0.015)+(0.313)*(0.0075), 0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 6)),
  matrix(c(0,  100, 6, 0,  (0.679)*(0.015)+(0.313)*(0.0075), 0), ncol = 3),
  # To post-natal 181-365
  cbind("transit ", 2, paste(" M single", 7)),
  matrix(c(0, 1,  100, 6, 0, 0, 0, (1-0.064)*(0.0075),  0), ncol = 3),
  cbind("transit ", 2, paste(" F single", 7)),
  matrix(c(0, 1,  100, 6, 0, 0, 0, (1-0.064)*(0.0075),  0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 7)),
  matrix(c(0, 1,  100, 6, 0, 0, 0, (1-0.064)*(0.0075),  0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 7)),
  matrix(c(0, 1,  100, 6, 0, 0, 0, (1-0.064)*(0.0075),  0), ncol = 3),
  # To post-natal 365 +
  cbind("transit ", 2, paste(" M single", 8)),
  matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-0.234)*(0.0075), 0), ncol = 3),
  cbind("transit ", 2, paste(" F single", 8)),
  matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-0.234)*(0.0075), 0), ncol = 3),
  cbind("transit ", 3, paste(" M single", 8)),
  matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-0.234)*(0.0075), 0), ncol = 3),
  cbind("transit ", 3, paste(" F single", 8)),
  matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-0.234)*(0.0075), 0), ncol = 3))

# Create the directories
dir.create("TM92HIV")
dir.create("90days")

Rates = as.list(1:nrow(all_param))

# Specify if vertical infection and lower fertility are allowed
vertical = TRUE 

  
#________________________________________________________________________________________________________     
# 2 - CREATE MICROSIMULATED POPULATIONS
#________________________________________________________________________________________________________     
  
    for(i in 1:n_sim) {
    for(vertical in c(TRUE, FALSE)){
    keeptimestart = Sys.time()
    # Generate the starting populations
    cat(paste("write.table(as.data.frame(pop[[", i, "]]), file = \"", "TM92HIV", "/", "Pop_START_", paste("sim", i, sep = "_"), ".opop\", quote = F, row.names=F, col.names=F); cat(NULL, file = \"", "TM92HIV", "/", "Pop_START_", paste("sim", i, sep = "_"), ".opox\"); cat(NULL, file = \"", "TM92HIV", "/", "Pop_START_", paste("sim", i, sep = "_"), ".omar\")",  sep = ""), file = paste("TM92HIV", "/", "Generate_Pop_START_", paste("sim", i, sep = "_"),".R", sep = ""))				
    source(paste(getwd(), "/TM92HIV/Generate_Pop_START_", "sim_", i,".R", sep = ""))
    unlink(paste(getwd(), "/TM92HIV/Generate_Pop_START_", "sim_", i,".R", sep = ""))
    
    # Create .sup file
    nb_seg = length(c(1800, seq(2000, 2045, 5)))
    
    SEGMENT = as.list(1:nb_seg)
    SEGMENT[[1]]= c(paste("\nduration", 200*12, sep = " "),
                    "\nbint 0", 
                    paste("\nsex_ratio", 0.5, sep = ' '),
                    "\nhetfert 0", 
                    "\nalpha  0",
                    "\nbeta 1",
                    paste("\ninclude ",  "TM92HIV", "/", "Rates_TM92_", paste("sim", i, sep = "_"), sep = ""),
                    "\nrun")
    cat("segments", 1,
        paste("\ninput_file ",  "TM92HIV",  "/", "Pop_START_", paste("sim", i, sep = "_"), sep = ""),
        paste("\noutput_file ", "90days", "/", "Pop_1_", paste("sim", i, sep = "_"), sep = ""),
        unlist(SEGMENT[[1]]),
        file = paste("TM92HIV", "/", "TM92_", paste("sim", i, sep = "_"), "_", 1, ".sup", sep = ""))
    
    # Create rates for first segment
    if(vertical == TRUE){
      Rates[[i]] = rbind(c("*" , "SIM", i), 
                         do.call(rbind,MORT_SEG[[i]]), 
                         do.call(rbind,FERT_SEGHIV[[i]]),
                         unlist(NUPT_SEG)
                         #,unlist(TRANS_SEG)
      )}
    if(vertical == FALSE){
      Rates[[i]] = rbind(c("*" , "SIM", i), 
                         do.call(rbind,MORT_SEG[[i]]), 
                         do.call(rbind,FERT_SEG[[i]]),
                         unlist(NUPT_SEG)
                         #,unlist(TRANS_SEG)
      )}
    
    cat(paste("write.table(as.data.frame(Rates[[", i, "]]), file = \"", "TM92HIV", "/", "Rates_TM92_", paste("sim", i, sep = "_"), "\", quote = F, row.names=F, col.names=F)",  sep = ""), 
        file = paste("TM92HIV", "/", "GenerateRate_TM92_", paste("sim", i, sep = "_"), ".R", sep = ""))		
    source(paste("TM92HIV", "/", "GenerateRate_TM92_", paste("sim", i, sep = "_"), ".R", sep = ""))
    unlink(paste("TM92HIV", "/", "GenerateRate_TM92_", paste("sim", i, sep = "_"), ".R", sep = ""))
    
    # Fromdos
    system(paste('fromdos ', "TM92HIV", '/', 'TM92_', paste("sim", i, sep = "_"),"_", 1, '.sup ', sep = ''), wait = T) 
    system(paste('fromdos ', "TM92HIV", '/', 'Rates_TM92_', paste("sim", i, sep = "_"), sep = ''), wait = T) 
   
    # Now run the first segment
    system(paste('socsim ', "TM92HIV", '/', 'TM92_', paste("sim", i, sep = "_"), "_", 1, '.sup ', i, sep = ''), wait = T) 
    
    # Distribution of female infections in Spectrum (2013)
    propinf = data.frame(seq(15, 65, 5), c(0.288, 0.239,  0.182, 0.105, 0.057, 0.086, 0.029, 0.002, 0.002, 0.005, 0.005))
    names(propinf) = c("age", 'prop')
    # De 15-19 ? 45-49
    propinf$propc = propinf$prop/sum(propinf$prop[1:7])
    # plot(seq(15, 65, 5), propinf$propc, type = 'l')
    
    # Following segments suivants
    end_seg = 1201+c(200, seq(205, 250, 5))*12
    
    for(k in 2:nb_seg){
      final_pop <- data.frame(read.table(file = paste("90days", '/Pop_', k-1, '_sim_', i, ".opop", sep = ''))); colnames(final_pop) = c("pid", "sex", "group", "nextevent", "month_dob", "mother", "father", "nesibm", "e_s_dad", "lborn", "last_marr", "mstatus", "month_dod", "fmult")
      if(k == 2){cat("\nGrowth rate (first seg.):", 100*(log(nrow(final_pop[final_pop$month_dod == 0,])/nrow(final_pop[final_pop$mother == 0,])))/200)}
      
      final_pop$age = trunc(((end_seg[k-1] - final_pop$month_dob)/12))
      final_pop$age5 = 5*trunc(final_pop$age/5)
      
      # allocate infections based on the incidence curve
      final_pop$atrisk = 0
      final_pop$atrisk[final_pop$month_dod == 0 & final_pop$age >= 15 & final_pop$age < 50] = 1
      propinf$inc = NA
      for(z in 1:7){
        propinf$inc[z] = propinf$propc[z]*incidence_curve[i,k-1]/(nrow(final_pop[final_pop$group == 1 & final_pop$atrisk == 1 & final_pop$sex == 0 & final_pop$age5 == seq(15, 45, 5)[z],])/
                                                                    nrow(final_pop[final_pop$atrisk == 1 & final_pop$sex == 0 & final_pop$group == 1,])  )
      }
      INC_SEG = rbind(
        # From 1 (non infected) to 2 (infected)
        cbind("transit ", 1, paste(" M single", 2)),
        cbind(taux_TM[c(seq(16, 51, 5), 101),1:2], format(round(c( 0,  1-(1-propinf$inc[1:7])^(1/12), 0), 6), sc= F)),
        cbind("transit ", 1, paste(" F single", 2)),
        cbind(taux_TM[c(seq(16, 51, 5), 101),1:2], format(round(c( 0,  1-(1-propinf$inc[1:7])^(1/12), 0), 6), sc= F)))
      
      #by default, no reduction in vertical transmission (no PTMCT)
      INFCHILD = INFCHILDnoPMTCT
      
      # Treatment - after year 2025 and for scenarios with ART
      if(all_param$coverage_ART[i]>0 & k > 7){
       # Transition to ART (women only)
       if (length(which(final_pop$group %in% c(2,3) & final_pop$atrisk == 1 & final_pop$sex == 1)) > 1){
         needingART = length(which(final_pop$group %in% c(2,3) & final_pop$atrisk == 1 & final_pop$sex == 1))
         underART = length(which(final_pop$group %in% c(9) & final_pop$atrisk == 1 & final_pop$sex == 1))
         estimatedcov = underART/(needingART+underART)
         nbnew = round(((art_curve[i,k-1]/100)-estimatedcov)*(needingART+underART))
         if(nbnew > 0){
         # first, allocate treatments to women aged 15-49 who have reached stage 3 and will give birth soon    
         set1 = which(final_pop$group %in% c(3) & final_pop$atrisk == 1 & final_pop$nextevent == 0) 
         # second, allocate treatments to women aged 15-49 who have reached stage 3    
         set2 = which(final_pop$group %in% c(3) & final_pop$atrisk == 1 & final_pop$sex == 1 & final_pop$nextevent != 0) 
         # third, allocate treatments to HIV-positive women aged 15-49 who will give birth soon
         set3 = which(final_pop$group %in% c(2) & final_pop$atrisk == 1 & final_pop$nextevent == 0) 
         # fourth, allocate treatments to HIV-positive women aged 15-49
         set4 = which(final_pop$group %in% c(2) & final_pop$atrisk == 1 & final_pop$sex == 1 & final_pop$nextevent != 0) 
         to_art = c(set1, set2, set3, set4)[1:nbnew]
         fromwhere =  final_pop$group[to_art]
         final_pop$group[to_art] = 9
         write.table(final_pop, file = paste(paste("90days", '/Pop_', k-1, '_sim_', i, ".opop", sep = "")), quote = F, row.names=F, col.names=F)
         transit_pop <- data.frame(read.table(file = paste("90days", '/Pop_', k-1, '_sim_', i, ".otx", sep = ''))); colnames(transit_pop) = c("pid", "date", "fromg", "tog", "sequence")
         transit_pop$sequence[transit_pop$pid %in% to_art] = transit_pop$sequence[transit_pop$pid %in% to_art] - 1 
         transit_pop = as.data.frame(rbind(as.matrix(transit_pop), cbind(to_art, rep(end_seg[k-1], nbnew), fromwhere, rep(9, nbnew), 0)))
         transit_pop = transit_pop[order(transit_pop[,1], -transit_pop[,5]),]
         write.table(transit_pop, file = paste(paste("90days", '/Pop_', k-1, '_sim_', i, ".otx", sep = "")), quote = F, row.names=F, col.names=F)
         }
        # PMTCT in addition to ART
        INFCHILDPMTCT = rbind(     
           # The probability of HIV infection at birth for a child born to an HIV-positive mother is assumed to be 20% in the absence of prophylaxis
           
           # From born from infected mother (2, 3 or 4) to infected child (5): 20% chance in utero -> fast progressors
           cbind("transit ", 2, paste(" M single", 5)),
           matrix(c(0, 100, 1,  0, ((1-pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
           cbind("transit ", 2, paste(" F single", 5)),
           matrix(c(0, 100, 1,  0, ((1-pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
            cbind("transit ", 3, paste(" M single", 5)),
           matrix(c(0, 100, 1,  0, ((1-pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 5)),
           matrix(c(0, 100, 1,  0, ((1-pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), 0), ncol = 3),
           # From born from infected mother (2, 3 or 4) to non infected (1): everyone is transferred back to 1 after age 2
           cbind("transit ", 2, paste(" M single", 1)),
           matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
           cbind("transit ", 2, paste(" F single", 1)),
           matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
           cbind("transit ", 3, paste(" M single", 1)),
           matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 1)),
           matrix(c(2, 10, 100, 0,  0, 0, 0, 0.95, 0), ncol = 3),
           # born to mother under ART - immediately sent back to 1, assumed to be under PMTCT
           cbind("transit ", 9, paste(" M single", 1)),
           matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
           cbind("transit ", 9, paste(" F single", 1)),
           matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
           # born to mother under PMTCT - immediately sent back to 1
           cbind("transit ", 10, paste(" M single", 1)),
           matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
           cbind("transit ", 10, paste(" F single", 1)),
           matrix(c(10, 100,   0, 0,  0.95, 0), ncol = 3),
           # The probability of infection through breastfeeding is assumed to be 1.5% per month for mixed feeding during the first 6 months, 0.75% per month for exclusive feeding during the first 6 months, 0.75% per month for months 7 and later 
           # (Sex Trans Inf 2010 - The Spectrum projection package: improvements in estimating incidence by age and sex, mother-to-child transmission, HIV progression in children and double orphans)
           # To post-natal 0-180
           cbind("transit ", 2, paste(" M single", 6)),
           matrix(c(0,  100, 6, 0,  (1-pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075)), 0), ncol = 3),
           cbind("transit ", 2, paste(" M single", 6)),
           matrix(c(0,  100, 6, 0,  (1-pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075)), 0), ncol = 3),
           cbind("transit ", 3, paste(" M single", 6)),
           matrix(c(0,  100, 6, 0,  (1-pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075)), 0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 6)),
           matrix(c(0,  100, 6, 0,  (1-pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075)), 0), ncol = 3),
           # To post-natal 181-365
           cbind("transit ", 2, paste(" M single", 7)),
           matrix(c(0, 1,  100, 6, 0, 0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075),  0), ncol = 3),
           cbind("transit ", 2, paste(" F single", 7)),
           matrix(c(0, 1,  100, 6, 0, 0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075),  0), ncol = 3),
           cbind("transit ", 3, paste(" M single", 7)),
           matrix(c(0, 1,  100, 6, 0, 0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075),  0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 7)),
           matrix(c(0, 1,  100, 6, 0, 0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075),  0), ncol = 3),
           # To post-natal 365 +
           cbind("transit ", 2, paste(" M single", 8)),
           matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           cbind("transit ", 2, paste(" F single", 8)),
           matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           cbind("transit ", 3, paste(" M single", 8)),
           matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 8)),
           matrix(c( 1, 2, 100, 0, 0,  0, 0, (1-pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           # Vertical transmission prevented through PMTCT,
           cbind("transit ", 2, paste(" M single", 10)),
           matrix(c(0, 0, 100, 1, 2, 0, 0, ((pmtct_curve[i,k-1]/100)), 0), ncol = 3),
           cbind("transit ", 2, paste(" F single", 10)),
           matrix(c(0, 0, 100, 1, 2, 0, 0, ((pmtct_curve[i,k-1]/100)), 0), ncol = 3),
           cbind("transit ", 3, paste(" M single", 10)),
           matrix(c(0, 0, 100, 1, 2, 0, 0, ((pmtct_curve[i,k-1]/100)), 0), ncol = 3),
           cbind("transit ", 3, paste(" F single", 10)),
           matrix(c(0, 0, 100, 1, 2, 0, 0, ((pmtct_curve[i,k-1]/100)), 0), ncol = 3)
           # cbind("transit ", 2, paste(" M single", 10)),
           # matrix(c(0, 0, 0, 1, 2, 100, 1, 2, 6, 0, 0, 0, 0, ((pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), (pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075))^(6/4),  (pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075), (pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           # cbind("transit ", 2, paste(" F single", 10)),
           # matrix(c(0, 0, 0, 1, 2, 100, 1, 2, 6, 0, 0, 0, 0, ((pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), (pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075))^(6/4),  (pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075), (pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           # cbind("transit ", 3, paste(" M single", 10)),
           # matrix(c(0, 0, 0, 1, 2, 100, 1, 2, 6, 0, 0, 0, 0, ((pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), (pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075))^(6/4),  (pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075), (pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3),
           # cbind("transit ", 3, paste(" F single", 10)),
           # matrix(c(0, 0, 0, 1, 2, 100, 1, 2, 6, 0, 0, 0, 0, ((pmtct_curve[i,k-1]/100)*0.2)/F_life_tables[1,"npx",all_param$n_fert[i],all_param$n_mort[i]]^(1/12), (pmtct_curve[i,k-1]/100)*((0.679)*(0.015)+(0.313)*(0.0075))^(6/4),  (pmtct_curve[i,k-1]/100)*(1-0.064)*(0.0075), (pmtct_curve[i,k-1]/100)*(1-0.234)*(0.0075), 0), ncol = 3)
           )
        INFCHILD = INFCHILDPMTCT
       }
        
       }
       
      if(vertical == TRUE){
        SEGMENT[[k]]= c(paste("\nduration", 5*12, sep = " "),
                        "\nbint 0", 
                        paste("\nsex_ratio", 0.5, sep = ' '),
                        "\nhetfert 0", 
                        "\nalpha  0",
                        "\nbeta 1",
                        paste("\ninclude ",  "TM92HIV", "/", "Rates_TM92_", paste("sim", i, sep = "_"), sep = ""),
                        paste('\n', apply(INCUB, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(AIDS, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(INFCHILD, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(SURVCHILD, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(INC_SEG, 1, function(x) paste(x, collapse = " "))),
                        #paste('\n', apply(TO_ART, 1, function(x) paste(x, collapse = " "))),
                        "\n*birthtarget  1 958", 
                        "\nrun")
      }
      if(vertical == FALSE){
        SEGMENT[[k]]= c(paste("\nduration", 5*12, sep = " "),
                        "\nbint 0", 
                        paste("\nsex_ratio", 0.5, sep = ' '),
                        "\nhetfert 0", 
                        "\nalpha  0",
                        "\nbeta 1",
                        paste("\ninclude ",  "TM92HIV", "/", "Rates_TM92_", paste("sim", i, sep = "_"), sep = ""),
                        paste('\n', apply(INCUB, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(AIDS, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(INFCHILDNOVERT, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(SURVCHILD, 1, function(x) paste(x, collapse = " "))),
                        paste('\n', apply(INC_SEG, 1, function(x) paste(x, collapse = " "))),
                        #paste('\n', apply(TO_ART, 1, function(x) paste(x, collapse = " "))),
                        "\n*birthtarget  1 958", 
                        "\nrun")
      }
      
      cat("segments", 1,
          paste("\ninput_file ",  "90days", "/", "Pop_", k-1, "_", paste("sim", i, sep = "_"), sep = ""),
          paste("\noutput_file ", "90days", "/", "Pop_", k, "_", paste("sim", i, sep = "_"), sep = ""),
          unlist(SEGMENT[[k]]),
          file = paste("TM92HIV", "/", "TM92_", paste("sim", i, sep = "_"), "_", k, ".sup", sep = ""))
      
      system(paste('fromdos ', "TM92HIV", '/', 'TM92_', paste("sim", i, sep = "_"),"_", k, '.sup ', sep = ''), wait = T) 
      system(paste('socsim ', "TM92HIV", '/', 'TM92_', paste("sim", i, sep = "_"), "_", k, '.sup ', i, sep = ''), wait = T) 
      
      if(k == nb_seg & vertical == TRUE){
        final_pop <- data.frame(read.table(file = paste("90days", '/Pop_', k, '_sim_', i, ".opop", sep = ''))); colnames(final_pop) = c("pid", "sex", "group", "nextevent", "month_dob", "mother", "father", "nesibm", "e_s_dad", "lborn", "last_marr", "mstatus", "month_dod", "fmult")
        save(final_pop, file = paste(getwd(), '/90days/Pop_TM92_sim_', i, "_out.rda", sep = ''), ascii = F)
        final_trans <- data.frame(read.table(file = paste("90days", '/Pop_', k, '_sim_', i, ".otx", sep = '')))
        colnames(final_trans) = c("pid", "date", "fromg", "tog", "sequence")
        save(final_trans, file = paste(getwd(), "/90days", '/Trans_TM92_sim_', i, "_out.rda", sep = ''), ascii = F)
      }
      if(k == nb_seg & vertical == FALSE){
        final_pop <- data.frame(read.table(file = paste("90days", '/Pop_', k, '_sim_', i, ".opop", sep = ''))); colnames(final_pop) = c("pid", "sex", "group", "nextevent", "month_dob", "mother", "father", "nesibm", "e_s_dad", "lborn", "last_marr", "mstatus", "month_dod", "fmult")
        save(final_pop, file = paste(getwd(), '/90days', '/Pop_TM92_sim_', i, "novert_out.rda", sep = ''), ascii = F)
        final_trans <- data.frame(read.table(file = paste("90days", '/Pop_', k, '_sim_', i, ".otx", sep = '')))
        colnames(final_trans) = c("pid", "date", "fromg", "tog", "sequence")
        save(final_trans, file = paste(getwd(), "/90days", '/Trans_TM92_sim_', i, "novert_out.rda", sep = ''), ascii = F)
      } 
      keeptimeend = Sys.time()
      
   # }
 
  }
    system(paste("rm ", "90days", "/*.opop", sep = ''))
    system(paste("rm ", "90days", "/*.opox", sep = ''))
    system(paste("rm ", "90days", "/*.omar", sep = ''))
    system(paste("rm ", "90days", "/*.pyr", sep = ''))
    system(paste("rm ", "90days", "/*.otx", sep = ''))
  # End of simulation

#________________________________________________________________________________________________________     
# 2 - RECALCULATE UNDERLYING MORTALITY AND COMPUTE ORPHAN PREVALENCE
#________________________________________________________________________________________________________     
  
  for(i in 1:n_sim){
    
    if(vertical == TRUE){
      stopifnot(file.exists(paste(getwd(), '/90days/Pop_TM92_sim_', i, "_out.rda", sep = '')))
    }
    if(vertical == FALSE){
      stopifnot(file.exists(paste(getwd(), '/90days/Pop_TM92_sim_', i, "novert_out.rda", sep = '')))
    }
    }

for(i in 1:n_sim){
    cat("\n", "** Start -", paste("sim", i, sep = "_"), "**", "\n")
   
    results_sim = vector("list", 1)
    names(results_sim) = c("MACB")
    
    #___________________________________________________________________________
    # 1. Load file, compute d_o_d and d_o_b
    if(vertical == TRUE){
      load(paste(getwd(), '/90days/Pop_TM92_sim_', i, "_out.rda", sep = ''))
    }
    if(vertical == FALSE){
      load(paste(getwd(), "/90days\\", "Pop_TM92_sim_", i, "novert_out.rda", sep = ""))
    }
    
    file = final_pop;
    names(file) = sub("e_s_mom", "nesibm", names(file)); names(file) = sub("d_o_b", "dob", names(file)); names(file) = sub("d_o_d", "dod", names(file)); names(file) = sub("person_id", "pid", names(file))
    names(file) = sub("last_marr", "union_id", names(file))
    
    file$month_dod[file$month_dod == 0] = NA
    
    start_sim = max(file$month_dob[file$mother == 0]+1, na.rm = T); end_sim = max(c(file$month_dod, file$month_dob),na.rm = TRUE) + 1; length_sim = end_sim - start_sim
    time.cut = seq(start_sim, end_sim + 12, by = 12); time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    # Growth on the first 200 years
    cat(100*(log(nrow(file[c(is.na(file$month_dod) | file$month_dod >= start_sim+12*200 )& file$month_dob <= start_sim+12*200,])/nrow(file[file$mother == 0,])))/200)
    
    w = file[,c(grep("pid|dob|mother|father|dod|sex|union_id|nesibm|group", colnames(file)))]
    
    # 1. DOB and DOD
    w$journ = trunc(runif(nrow(w), min =0 , max = 31)); w$journ[w$month_dob ==  start_sim] = trunc(runif(nrow(w[w$month_dob ==  start_sim,]), min = 1 , max = 31))
    w$dob = w$month_dob + w$journ/31#; w$dob[w$dob %in%  ((start_seg - 1800)*12 + 1)] = w$dob[w$dob %in%  ((start_seg - 1800)*12 + 1)] - 1/31 
    w$month_dod[w$month_dod == 0] = NA
    w$jourd = trunc(runif(nrow(w), min =0 , max = 31)); w$jourd[!is.na(w$month_dod) & w$month_dod ==  start_sim] = trunc(runif(nrow(w[!is.na(w$month_dod) & w$month_dod ==  start_sim,]), min = 1 , max = 31))
    w$jourd[which(w$month_dod == w$month_dob)] = trunc(runif(nrow(w[which(w$month_dod == w$month_dob),]), min =w$journ[which(w$month_dod == w$month_dob)]+1, max = 31))
    w$dod = w$month_dod + w$jourd/31
    w$dob[!is.na(w$dod) & w$dob == w$dod] = w$dob[!is.na(w$dod) & w$dob == w$dod] - 1/31
    
    w = w[order(w$pid),]; stopifnot(which(rownames(w) != w$pid) == 0)
    
    # Slight change if age of mothers at birth is a multiple of 12 months 
    w$dob_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dob[w$mother[(length(w$mother[w$mother == 0]) +1):nrow(w)]])
    w$dob[which(abs(floor((w$dob - w$dob_m)/12) - (w$dob - w$dob_m)/12) == 0)] = w$dob[which(abs(floor((w$dob - w$dob_m)/12) - (w$dob - w$dob_m)/12) == 0)] + 1/31 
    
    # Check that children are not born after their mother's death and get mother's dob
    w$dod_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dod[w$mother[(length(w$mother[w$mother == 0]) +1):nrow(w)]])
    w$dob[!is.na(w$dod_m) & w$dob >= w$dod_m] = w$dob[!is.na(w$dod_m) & w$dob >= w$dod_m] - 1
    w$dob_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dob[w$mother[length(w$mother[w$mother == 0]):nrow(w)]])
    
    # Date of birth of next eldest sibling
    zna<-function(x){return(ifelse(x==0,NA,x))} 
    w$dob_nesibm = w$dob[zna(w$nesibm)]
    while( nrow(w[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm,]) != 0) {
      w$dob[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm] = w$dob_nesibm[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm] + 1/31
      w$dob_nesibm = w$dob[zna(w$nesibm)]
    }
    w$dod[!is.na(w$dod) & w$dob >= w$dod] = w$dob[!is.na(w$dod) & w$dob >= w$dod] + 1/31
    w$dob[w$mother > 0 & w$dob < start_sim] = start_sim + 0.01
    w$dob[w$dob > end_sim] = end_sim - 0.01
    
    #___________________________________________________________________________
    # Mean age at chilbearing
    temp = w[w$mother > 0,]
    temp$MACB = trunc((temp$dob - temp$dob_m)/12)
    temp$yeardob<- trunc((temp$dob-start_sim)/12 + 1800)
    temp = xtabs(~temp$MACB + temp$yeardob)
    
    results_sim$MACB = apply(temp, 2, function(x) weighted.mean(as.numeric(rownames(temp)), w = x))+0.5
    
    #___________________________________________________________________________
    # Estimation of mortality rates (overall)
    
    bio.mort = w[,-c(grep("month|jour|dod_m|dob_m|bornfromHIV|intrapartum", colnames(w)))]
    
    # Follow-up time pour les personnes decedees et les survivants
    bio.mort$fu.time = (bio.mort$dod - bio.mort$dob)
    bio.mort$fu.time[is.na(bio.mort$fu.time)] = (end_sim - bio.mort$dob[is.na(bio.mort$fu.time)])
    
    # Folluw-up time pour les personnes decedees nees avant start_sim
    bio.mort$fu.time[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] = (bio.mort$dod[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] - start_sim)
    stopifnot(nrow(bio.mort[is.na(bio.mort$fu.time),]) == 0)
    
    # Event
    bio.mort$event = as.numeric(!is.na(bio.mort$dod))
    
    # Survival object	
    bio.mort_Surv = with(bio.mort, Surv(time = bio.mort$fu.time, event = bio.mort$event, type = "right"))
    
    bio.mort$startage = NA
    bio.mort$startage[which(bio.mort$dob >= start_sim)] = 0 
    bio.mort$startage[which(bio.mort$dob < start_sim)] = (start_sim - bio.mort$dob[which(bio.mort$dob < start_sim)])/12
    stopifnot(nrow(bio.mort[is.na(bio.mort$startage),]) == 0)
    
    age <- tcut(bio.mort$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
    
    time.cut = seq(start_sim, end_sim + 12, by = 12)
    time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    annee <- tcut(bio.mort$dob, time.cut, labels=time.cut.labels)
    annee[bio.mort$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
    
    bio.mort_pp <- pyears(bio.mort_Surv ~ annee + age + sex , bio.mort, scale = 12, data.frame=TRUE)
    stopifnot(round(sum(bio.mort$fu.time)/12 - sum(bio.mort_pp$data$pyears), 4) == 0)
    mort_pp <- bio.mort_pp$data; stopifnot(nrow(mort_pp[mort_pp$event > 1 & mort_pp$pyears == 0,]) == 0)
    mort_pp$age <- as.numeric(as.character(mort_pp$age))/12
    mort_pp$age[mort_pp$age > 80] = 80
    mort_pp$agegrp = trunc(mort_pp$age/5)*5
    mort_pp$annee<- (as.numeric(as.character(mort_pp$annee))-start_sim)/12 + 1800
    
    mort_pp$sim.id = i
    
    results_sim$mort_pp = mort_pp
    
    bio.mort$agedece = (bio.mort$dod-bio.mort$dob)/12
    bio.mort$agedecetr = trunc((bio.mort$dod-bio.mort$dob)/12)
    bio.mort$agedecetr[bio.mort$agedecetr > 80] = 80
    
    bio.mort$anneedec5<- trunc(((bio.mort$dod-start_sim)/12 + 1800)/5)*5
    
    nax = tapply(bio.mort$agedece-bio.mort$agedecetr, list(bio.mort$agedecetr, bio.mort$anneedec5), mean)
    
    results_sim$nax = nax
    #___________________________________________________________________________
    # Exposure time of children who are not orphaned
    bio.nonorph_m = w[w$mother != 0,-c(grep("month|mother|jour|dob_m|bornfromHIV|intrapartum", colnames(w)))]
    
    # Exposure time for the deceased and surviving
    bio.nonorph_m$fu.time[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod >= bio.nonorph_m$dod_m] = (bio.nonorph_m$dod_m[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod >= bio.nonorph_m$dod_m] - bio.nonorph_m$dob[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod >= bio.nonorph_m$dod_m])
    bio.nonorph_m$fu.time[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod < bio.nonorph_m$dod_m] = (bio.nonorph_m$dod[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod < bio.nonorph_m$dod_m] - bio.nonorph_m$dob[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & bio.nonorph_m$dod < bio.nonorph_m$dod_m])
    
    bio.nonorph_m$fu.time[is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m)] = bio.nonorph_m$dod_m[is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m)] - bio.nonorph_m$dob[is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m)]
    bio.nonorph_m$fu.time[!is.na(bio.nonorph_m$dod) & is.na(bio.nonorph_m$dod_m)] = bio.nonorph_m$dod[!is.na(bio.nonorph_m$dod) & is.na(bio.nonorph_m$dod_m)] - bio.nonorph_m$dob[!is.na(bio.nonorph_m$dod) & is.na(bio.nonorph_m$dod_m)]
    
    bio.nonorph_m$fu.time[is.na(bio.nonorph_m$fu.time)] = (end_sim - bio.nonorph_m$dob[is.na(bio.nonorph_m$fu.time)])
    
    # Exposure time for those who died and were born before start_sim
    bio.nonorph_m$fu.time[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & (bio.nonorph_m$dob < start_sim) & bio.nonorph_m$dod > bio.nonorph_m$dod_m] = (bio.nonorph_m$dod[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & (bio.nonorph_m$dob < start_sim) & bio.nonorph_m$dod > bio.nonorph_m$dod_m] - start_sim)
    bio.nonorph_m$fu.time[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & (bio.nonorph_m$dob < start_sim) & bio.nonorph_m$dod < bio.nonorph_m$dod_m] = (bio.nonorph_m$dod_m[!is.na(bio.nonorph_m$dod) & !is.na(bio.nonorph_m$dod_m) & (bio.nonorph_m$dob < start_sim) & bio.nonorph_m$dod < bio.nonorph_m$dod_m] - start_sim)
    
    stopifnot(nrow(bio.nonorph_m[is.na(bio.nonorph_m$fu.time),]) == 0)
    
    # Event
    bio.nonorph_m$event = as.numeric(!is.na(bio.nonorph_m$dod_m))
    
    # Survival object	
    bio.nonorph_m_Surv = with(bio.nonorph_m, Surv(time = bio.nonorph_m$fu.time, event = bio.nonorph_m$event, type = "right"))
    
   
    bio.nonorph_m$startage = NA
    bio.nonorph_m$startage[which(bio.nonorph_m$dob >= start_sim)] = 0 
    bio.nonorph_m$startage[which(bio.nonorph_m$dob < start_sim)] = (start_sim - bio.nonorph_m$dob[which(bio.nonorph_m$dob < start_sim)])/12
    stopifnot(nrow(bio.nonorph_m[is.na(bio.nonorph_m$startage),]) == 0)
    
    age <- tcut(bio.nonorph_m$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
    
    time.cut = seq(start_sim, end_sim + 12, by = 12)
    time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    annee <- tcut(bio.nonorph_m$dob, time.cut, labels=time.cut.labels)
    annee[bio.nonorph_m$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
    
    bio.nonorph_m_pp <- pyears(bio.nonorph_m_Surv ~ annee + age + sex, bio.nonorph_m, scale = 12, data.frame=TRUE)
    
    # Check that offtable  = 0
    stopifnot(round(sum(bio.nonorph_m$fu.time)/12 - sum(bio.nonorph_m_pp$data$pyears), 4) == 0)
    
    nonorph_m_pp <- bio.nonorph_m_pp$data; stopifnot(nrow(nonorph_m_pp[nonorph_m_pp$event > 1 & nonorph_m_pp$pyears == 0,]) == 0)
    nonorph_m_pp$age <- as.numeric(as.character(nonorph_m_pp$age))/12
    nonorph_m_pp$age[nonorph_m_pp$age > 80] = 80
    nonorph_m_pp$agegrp = trunc(nonorph_m_pp$age/5)*5
    
    nonorph_m_pp$annee<- (as.numeric(as.character(nonorph_m_pp$annee))-start_sim)/12 + 1800
    
    nonorph_m_pp$sim.id = i
    
    results_sim$nonorph_m_pp = nonorph_m_pp
    
    #___________________________________________________________________________
    # Exposure time of children who are orphaned
    
    bio.orph_m = w[w$mother != 0 & !is.na(w$dod_m) & (is.na(w$dod) | (!is.na(w$dod) & w$dod > w$dod_m)),-c(grep("month|mother|jour|dob_m|bornfromHIV|intrapartum", colnames(w)))]
    
    # Exposure time for the deceased and surviving
    bio.orph_m$fu.time = bio.orph_m$dod - bio.orph_m$dod_m
    bio.orph_m$fu.time[is.na(bio.orph_m$fu.time)] = end_sim - bio.orph_m$dod_m[is.na(bio.orph_m$fu.time)]
    stopifnot(nrow(bio.orph_m[is.na(bio.orph_m$fu.time),]) == 0)
    
    stopifnot(round(sum(bio.nonorph_m$fu.time) + sum(bio.orph_m$fu.time) - sum(bio.mort$fu.time[bio.mort$mother > 0])) < 5)
    
    # Event
    bio.orph_m$event = as.numeric(!is.na(bio.orph_m$dod))
    
    # Survival object	
    bio.orph_m_Surv = with(bio.orph_m, Surv(time = bio.orph_m$fu.time, event = bio.orph_m$event, type = "right"))
    
    bio.orph_m$startage = NA
    bio.orph_m$startage[which(bio.orph_m$dob >= start_sim)] = (bio.orph_m$dod_m[which(bio.orph_m$dob >= start_sim)] - bio.orph_m$dob[which(bio.orph_m$dob >= start_sim)])/12 
    bio.orph_m$startage[which(bio.orph_m$dob < start_sim)] = (bio.orph_m$dod_m[which(bio.orph_m$dob >= start_sim)] - start_sim)/12
    stopifnot(nrow(bio.orph_m[is.na(bio.orph_m$startage),]) == 0)
    
    age <- tcut(bio.orph_m$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
    
    time.cut = seq(start_sim, end_sim + 12, by = 12)
    time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    annee <- tcut(bio.orph_m$dod_m, time.cut, labels=time.cut.labels)
    
    #bio.orph_m$cohort = (trunc(bio.orph_m$dob/(12*5))*5 + 1800)
    
    bio.orph_m_pp <- pyears(bio.orph_m_Surv ~ annee + age + sex, bio.orph_m, scale = 12, data.frame=TRUE)
    
    # Check that offtable  = 0
    stopifnot(round(sum(bio.orph_m$fu.time)/12 - sum(bio.orph_m_pp$data$pyears), 4) == 0)
    
    orph_m_pp <- bio.orph_m_pp$data; stopifnot(nrow(orph_m_pp[orph_m_pp$event > 1 & orph_m_pp$pyears == 0,]) == 0)
    orph_m_pp$age <- as.numeric(as.character(orph_m_pp$age))/12
    orph_m_pp$age[orph_m_pp$age > 80] = 80
    orph_m_pp$agegrp = trunc(orph_m_pp$age/5)*5
    
    orph_m_pp$annee<- (as.numeric(as.character(orph_m_pp$annee))-start_sim)/12 + 1800
    
    orph_m_pp$sim.id = i
    
    results_sim$orph_m_pp = orph_m_pp
    
    #___________________________________________________________________________
    # Incidence & prevalence
    bio.inf = w[,-c(grep("month|jour|dod_m|dob_m|bornfromHIV|intrapartum", colnames(w)))]
    if(vertical == TRUE){
      load(file = paste("90days", '/Trans_TM92_sim_', i, "_out.rda", sep = ''))
    }
    if(vertical == FALSE){
      load(file = paste("90days", '/Trans_TM92_sim_', i, "novert_out.rda", sep = ''))
    }
    
    # Proportion of children born from HIV-infected mothers who aquire HIV (either in utero, during childbirth or through breastfeeding)
    # Children born from HIV-infected mothers, with vertical transmission
    MTCHIV = final_trans[final_trans$fromg %in% c(2:3) & final_trans$tog %in% c(5:8), grep('pid|date', names(final_trans))]
    # children born from HIV-infected mothers, without vertical transmission
    MTCnoHIV = final_trans[final_trans$fromg %in% c(2:3) & final_trans$tog %in% c(1), grep('pid|date', names(final_trans))]
    # children born from HIV-infected mothers, no vertical transmission prevented thanks to PMTCT
    MTCPMTCT = final_trans[final_trans$fromg %in% c(2:3) & final_trans$tog %in% c(10), grep('pid|date', names(final_trans))]
    # children born from mothers on ART, vertical transmission averted thanks to ART
    MTCART = final_trans[final_trans$fromg %in% c(9) & final_trans$tog %in% c(1), grep('pid|date', names(final_trans))]
    stopifnot(unique(MTCHIV$pid %in% MTCnoHIV$pid) == FALSE)
    stopifnot(unique(MTCnoHIV$pid %in% MTCHIV$pid) == FALSE)
    stopifnot(unique(MTCPMTCT$pid %in% MTCHIV$pid) == FALSE)
    stopifnot(unique(MTCART$pid %in% MTCHIV$pid) == FALSE)
    stopifnot(unique(!duplicated(MTCHIV$pid)))
    stopifnot(unique(!duplicated(MTCnoHIV$pid)))
    
    # All children born to HIV-positive mother
    w$bornfromHIV = 0
    w$bornfromHIV[w$pid %in% c(MTCHIV$pid, MTCnoHIV$pid, MTCPMTCT$pid, MTCART$pid)] = 1
    
    # Vertically infected
    w$vertinf = 0
    w$vertinf[w$pid %in% c(MTCHIV$pid)] = 1
    
    # Vertical transmission prevented through PMTCT or ART
    w$pmtctprev = 0
    if(length(MTCART$pid) > 0 | length(MTCPMTCT$pid) > 0){
    w$pmtctprev[w$pid %in% c(MTCART$pid, MTCPMTCT$pid)] = 1
    }
    
    w$year_dob = trunc((w$month_dob-start_sim)/12 + 1800)
    
    # Total number of children born (1) from all mothers, (2) from HIV-infected mothers, (3) vertically infected
    w$count = 1
    results_sim$births = rbind(tapply(w$count, w$year_dob, sum),
                               tapply(w$bornfromHIV, w$year_dob, sum),
                               tapply(w$vertinf, w$year_dob, sum))
    
    results_sim$pmtct = NA
    if(all_param$coverage_ART[i] > 0){
      results_sim$pmtc_cov = tapply(w$pmtctprev, w$year_dob, sum)/tapply(w$bornfromHIV, w$year_dob, sum)
    }
    if(all_param$coverage_ART[i] == 0){
      results_sim$pmtc_cov = c(rep(0,200), rep(0,50))
    }  
    
    infections = final_trans[final_trans$fromg == 1 & final_trans$tog == 2, grep('pid|date', names(final_trans))]
    stopifnot(unique(duplicated(infections$pid)) == FALSE)
    bio.inf = merge(bio.inf, infections, by = 'pid', all.x= T, all.y = T)
    # Exposure time for the infected, deceased and surviving
    bio.inf$fu.time = (bio.inf$date - bio.inf$dob) # infect?es
    bio.inf$fu.time[is.na(bio.inf$fu.time)] = (bio.inf$dod[is.na(bio.inf$fu.time)] - bio.inf$dob[is.na(bio.inf$fu.time)])# d?c?d?es
    bio.inf$fu.time[is.na(bio.inf$fu.time)] = (end_sim - bio.inf$dob[is.na(bio.inf$fu.time)]) # survivants
    
    # Exposure time for the deceased born before start_sim
    bio.inf$fu.time[!is.na(bio.inf$dod) &  (bio.inf$dob < start_sim)] = (bio.inf$dod[!is.na(bio.inf$dod) &  (bio.inf$dob < start_sim)] - start_sim)
    stopifnot(nrow(bio.inf[is.na(bio.inf$fu.time),]) == 0)
    
    # No folluw-up time for infected children
    bio.inf$fu.time[bio.inf$group == 5] =0
    
    # Event
    bio.inf$event = as.numeric(!is.na(bio.inf$date))
    stopifnot(sum(bio.inf$event[bio.inf$group == 5]) < 10)
    
    # Survival object	
    bio.inf_Surv = with(bio.inf, Surv(time = bio.inf$fu.time, event = bio.inf$event, type = "right"))
    
    bio.inf$startage = NA
    bio.inf$startage[which(bio.inf$dob >= start_sim)] = 0 
    bio.inf$startage[which(bio.inf$dob < start_sim)] = (start_sim - bio.inf$dob[which(bio.inf$dob < start_sim)])/12
    stopifnot(nrow(bio.inf[is.na(bio.inf$startage),]) == 0)
    
    age <- tcut(bio.inf$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
    
    time.cut = seq(start_sim, end_sim + 12, by = 12)
    time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    annee <- tcut(bio.inf$dob, time.cut, labels=time.cut.labels)
    annee[bio.inf$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
    
    bio.inf_pp <- pyears(bio.inf_Surv ~ annee + age + sex , bio.inf, scale = 12, data.frame=TRUE)
    stopifnot(round(sum(bio.inf$fu.time)/12 - sum(bio.inf_pp$data$pyears), 4) == 0)
    inf_pp <- bio.inf_pp$data; stopifnot(nrow(inf_pp[inf_pp$event > 1 & inf_pp$pyears == 0,]) == 0)
    inf_pp$age <- as.numeric(as.character(inf_pp$age))/12
    inf_pp$age[inf_pp$age > 80] = 80
    inf_pp$agegrp = trunc(inf_pp$age/5)*5
    inf_pp$annee<- (as.numeric(as.character(inf_pp$annee))-start_sim)/12 + 1800
    
    inf_pp = inf_pp[inf_pp$age >= 15 & inf_pp$age < 50,]
    inf_pp = inf_pp[inf_pp$annee > 1999,]
    
   
    inf_pp$sim.id = i
    results_sim$inf_pp = inf_pp
    
    # pmtct coverage: births born to ART mother/all births born to HIV+ mother
      
    #___________________________________________________________________________
    # Population untreated with ART
    if(all_param$coverage_ART[i] > 0){
    bio.art = w[,-c(grep("month|jour|dod_m|dob_m|bornfromHIV|intrapartum", colnames(w)))]
    if(vertical == TRUE){
      load(file = paste("90days", '/Trans_TM92_sim_', i, "_out.rda", sep = ''))
    }
    if(vertical == FALSE){
      load(file = paste("90days", '/Trans_TM92_sim_', i, "novert_out.rda", sep = ''))
    }
    
    initiationART = final_trans[final_trans$fromg %in% c(2,3) & final_trans$tog == 9, grep('pid|date', names(final_trans))]
    stopifnot(unique(duplicated(initiationART$pid)) == FALSE)
    bio.art = merge(bio.art, initiationART, by = 'pid', all.x= T, all.y = T)
    # Exposure time under treatment
    bio.art$fu.time = (bio.art$date - bio.art$dob) # mis sous traitement
    bio.art$fu.time[is.na(bio.art$fu.time)] = (bio.art$dod[is.na(bio.art$fu.time)] - bio.art$dob[is.na(bio.art$fu.time)])# d?c?d?es
    bio.art$fu.time[is.na(bio.art$fu.time)] = (end_sim - bio.art$dob[is.na(bio.art$fu.time)]) # survivants
    
    # Exposure time for the deceased born before start_sim
    bio.art$fu.time[!is.na(bio.art$dod) &  (bio.art$dob < start_sim)] = (bio.art$dod[!is.na(bio.art$dod) &  (bio.art$dob < start_sim)] - start_sim)
    stopifnot(nrow(bio.art[is.na(bio.art$fu.time),]) == 0)
    
    # Event
    bio.art$event = as.numeric(!is.na(bio.art$date))
    
    # Survival object	
    bio.art_Surv = with(bio.art, Surv(time = bio.art$fu.time, event = bio.art$event, type = "right"))
    
    bio.art$startage = NA
    bio.art$startage[which(bio.art$dob >= start_sim)] = 0 
    bio.art$startage[which(bio.art$dob < start_sim)] = (start_sim - bio.art$dob[which(bio.art$dob < start_sim)])/12
    stopifnot(nrow(bio.art[is.na(bio.art$startage),]) == 0)
    
    age <- tcut(bio.art$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
    
    time.cut = seq(start_sim, end_sim + 12, by = 12)
    time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
    
    annee <- tcut(bio.art$dob, time.cut, labels=time.cut.labels)
    annee[bio.art$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
    
    bio.art_pp <- pyears(bio.art_Surv ~ annee + age + sex , bio.art, scale = 12, data.frame=TRUE)
    stopifnot(round(sum(bio.art$fu.time)/12 - sum(bio.art_pp$data$pyears), 4) == 0)
    art_pp <- bio.art_pp$data; stopifnot(nrow(art_pp[art_pp$event > 1 & art_pp$pyears == 0,]) == 0)
    art_pp$age <- as.numeric(as.character(art_pp$age))/12
    art_pp$age[art_pp$age > 80] = 80
    art_pp$agegrp = trunc(art_pp$age/5)*5
    art_pp$annee<- (as.numeric(as.character(art_pp$annee))-start_sim)/12 + 1800
    
    art_pp = art_pp[art_pp$age >= 15 & art_pp$age < 50,]
    art_pp = art_pp[art_pp$annee > 1999,]
    
    art_pp$sim.id = i
    results_sim$art_pp = art_pp
    tabart =  (tapply(mort_pp$pyears[mort_pp$annee > 1999 & mort_pp$age >= 15 & mort_pp$age < 50 & mort_pp$sex == 1], mort_pp$annee[mort_pp$annee > 1999 & mort_pp$age >= 15 & mort_pp$age < 50 & mort_pp$sex == 1], sum)-tapply(art_pp$pyears[art_pp$annee > 1999 & art_pp$sex == 1], art_pp$annee[art_pp$annee > 1999 & art_pp$sex == 1], sum))/
      (tapply(mort_pp$pyears[mort_pp$annee > 1999 & mort_pp$age >= 15 & mort_pp$age < 50 & mort_pp$sex == 1], mort_pp$annee[mort_pp$annee > 1999 & mort_pp$age >= 15 & mort_pp$age < 50 & mort_pp$sex == 1], sum)-tapply(inf_pp$pyears[inf_pp$annee > 1999 & inf_pp$sex == 1], inf_pp$annee[inf_pp$annee > 1999 & inf_pp$sex == 1], sum))
    results_sim$art_cov = c(rep(0,200), approx(as.numeric(names(tabart))+0.5, tabart, xout = 2000:2049 + 0.5)$y)
    }
    if(all_param$coverage_ART[i] == 0){
      results_sim$art_pp = i
      results_sim$art_cov = c(rep(0,200), rep(0,50))
    }
    
    results_sim$param = all_param[i,]
     if(vertical == TRUE){
       save(results_sim, file = paste(getwd(), "/90days/", "results_",  paste("sim", i, sep = "_"),  ".rda", sep = ""), ascii = FALSE)
     }
    if(vertical == FALSE){
      save(results_sim, file = paste(getwd(), "/90days/", "results_",  paste("sim", i, sep = "_"),  "novert.rda", sep = ""), ascii = FALSE)
    }
    }

# check that all results files have been produced
for(i in c(1:n_sim)){
  if(vertical == TRUE){
  stopifnot(file.exists(paste(getwd(), "/90days/", "results_",  paste("sim", i, sep = "_"), ".rda",  sep = "")))
  stopifnot(file.exists(paste(getwd(), "/90days/", "Pop_TM92_sim_",  i, "_out.rda", sep = "")))
  }
  if(vertical == FALSE){
    stopifnot(file.exists(paste(getwd(), "/90days/", "results_",  paste("sim", i, sep = "_"),  "novert.rda",  sep = "")))
    stopifnot(file.exists(paste(getwd(), "/90days/", "Pop_TM92_sim_",  i, "novert_out.rda", sep = "")))
  }
}

# Life table function
lt.mx = function (nmx, sex = "female", age = c(0, 1, seq(5, 110, 5)), 
                  nax = NULL) {
  if (is.null(nax)) {
    nax <- rep(2.5, length(age))
    if (sex == "male") {
      if (nmx[1] >= 0.107) {
        nax[1] <- 0.33
        nax[2] <- 1.352
      }
      else {
        nax[1] <- 0.045 + 2.684 * nmx[1]
        nax[2] <- 1.651 - 2.816 * nmx[1]
      }
    }
    if (sex == "female") {
      if (nmx[1] >= 0.107) {
        nax[1] <- 0.35
        nax[2] <- 1.361
      }
      else {
        nax[1] <- 0.053 + 2.8 * nmx[1]
        nax[2] <- 1.522 - 1.518 * nmx[1]
      }
    }
  }
  else {
    nax = nax
  }
  n <- c(diff(age), 999)
  nqx <- (n * nmx)/(1 + (n - nax) * nmx)
  nqx <- c(nqx[-(length(nqx))], 1)
  for (i in 1:length(nqx)) {
    if (nqx[i] > 1) 
      nqx[i] <- 1
  }
  nage <- length(age)
  nqx <- round(nqx, 4)
  npx <- 1 - nqx
  l0 = 1e+06
  lx <- round(cumprod(c(l0, npx)))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLx <- n * lxpn + ndx * nax
  Tx <- rev(cumsum(rev(nLx)))
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  lt <- cbind(Age = age, nax = round(nax, 3), nmx = round(nmx, 
                                                          4), nqx = round(nqx, 4), npx = round(npx, 4), ndx = ndx, 
              lx = lx, nLx = round(nLx), Tx = round(Tx), ex = round(ex, 
                                                                    2))
  #lt <- lt[lt[, 6] != 0, ]
  e0 <- lt[1, 10]
  lt.45q15 <- 1 - (lx[age == 60]/lx[age == 15])
  lt.5q0 <- 1 - (lx[age == 5]/lx[age == 0])
  return(list(e0 = e0, lt.5q0 = lt.5q0, lt.45q15 = lt.45q15, 
              lt = lt))
}

# With vertical transmission and reduced fertility
mort = vector("list", n_sim) 
MACB = matrix(NA, nrow = n_sim, ncol = 250)
orph = array(NA, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
PY = array(NA, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
orph10 = array(NA, dim =  c(length(seq(0, 45, 5)), length(1800:2049), n_sim))

for(i in 1:n_sim) {
  
  load(file = paste("90days\\", "results_",  paste("sim", i, sep = "_"),  '.rda', sep = ""))
  temp = results_sim$mort_pp[results_sim$mort_pp$sex == 1,]
  temp$annee5 = trunc(temp$annee/5)*5 # Create 5-year periods
  taux_F = (tapply(temp$event, list(temp$age, temp$annee5), sum)/tapply(temp$pyears, list(temp$age, temp$annee5), sum))
  taux_F[is.na(taux_F)]= 0
  nax = results_sim$nax
  
  nax = 1+(1/taux_F) - 1/(1-exp(-taux_F))# as per Preston, p.62
  nax[is.na(nax)] = mean(F_life_tables[2:99,3,all_param$n_fert[i],all_param$n_mort[i]])
  
  mort[[i]] = as.list(1:length(seq(1800, 2045, 5))) # for periods starting in 1800, 1805, etc. until 2045
  for(j in 1:length(seq(1800, 2045, 5))){
    mort[[i]][[j]] = lt.mx(taux_F[,j], sex="female", age=c(0:80), nax = nax[,j])
  }
  
  temporph = results_sim$orph_m_pp
  temporph$annee = factor(x = temporph$annee, 1800:2050, labels = 1800:2050)
  tempnonorph = results_sim$nonorph_m_pp
  tempnonorph$annee = factor(x = tempnonorph$annee, 1800:2050, labels = 1800:2050)
  
  # Trends in orphan prevalence
  pyearsorph = tapply(temporph$pyears, list(temporph$agegrp, temporph$annee), sum)[1:11,]
  pyearsnonorph = tapply(tempnonorph$pyears, list(tempnonorph$agegrp, tempnonorph$annee), sum)[1:11,]
  # Annual estimates
  orph[,,i] = pyearsorph[1:10,1:250]/(pyearsorph[1:10,1:250]+pyearsnonorph[1:10,1:250])
  orph10[,,i] = (pyearsorph[1:10,1:250]+pyearsorph[2:11,1:250] )/((pyearsorph[1:10,1:250]+pyearsorph[2:11,1:250] )+(pyearsnonorph[1:10,1:250]+pyearsnonorph[2:11,1:250] ))
  PY[,,i] = (pyearsorph[1:10,1:250]+pyearsnonorph[1:10,1:250])
  
  MACB[i,] = results_sim$MACB[1:250]
  print(i)
} 
dimnames(orph)[[2]] = 1800:2049 + 0.5
dimnames(orph)[[3]] = 1:n_sim

dimnames(PY)[[2]] = 1800:2049 + 0.5
dimnames(PY)[[3]] = 1:n_sim

save(mort, file = "90days\\mort.rda", ascii = FALSE)
save(orph, file = "90days\\orph.rda", ascii = FALSE)
save(PY, file = "90days\\PY.rda", ascii = FALSE)
save(orph10, file = "90days\\orph10.rda", ascii = FALSE)
save(MACB, file = "90days\\MACB.rda", ascii = FALSE)

load(file = "90days\\mort.rda")
load(file = "90days\\orph.rda")
load(file = "90days\\orph10.rda")
load(file = "90days\\MACB.rda")
load(file = "90days\\PY.rda")

# Adult and child mortality in DHS and simulations
library(rdhs)
indicators <- dhs_indicators()
U5MR <- dhs_data(indicatorIds = "CM_ECMT_C_U5M")
q3515 <- dhs_data(indicatorIds = "MM_AMPB_W_AMP")
HIVDHS <- dhs_data(indicatorIds = "HA_HIVP_W_HIV")

# variations in mortality according to parameters
tapply((1-sapply(mort, function(x) x[[50]]$lt['50',7])/sapply(mort, function(x) x[[50]]$lt['15',7]))*1000,
all_param$coverage_ART, mean)
tapply((1-sapply(mort, function(x) x[[50]]$lt['50',7])/sapply(mort, function(x) x[[50]]$lt['15',7]))*1000,
       all_param$alpha_mort, mean)
tapply((1-sapply(mort, function(x) x[[50]]$lt['50',7])/sapply(mort, function(x) x[[50]]$lt['15',7]))*1000,
       all_param$alpha_HIV, mean)

# Number of survivors

  nbsurv = as.vector(rep(NA, n_sim))
  for(i in 1:n_sim) {
    load(paste("90days\\", "Pop_TM92_sim_", i, "_out.rda", sep = ""))
    file = final_pop;
    names(file) = sub("e_s_mom", "nesibm", names(file)); names(file) = sub("d_o_b", "dob", names(file)); names(file) = sub("d_o_d", "dod", names(file)); names(file) = sub("person_id", "pid", names(file))
    names(file) = sub("last_marr", "union_id", names(file))
    nbsurv[i] = nrow(file[file$month_dod == 0,])
    print(i)
  }
  save(nbsurv, file = "90days\\nbsurv.rda", ascii = FALSE)
     
load(file = "90days\\nbsurv.rda")

# Incidence and prevalence trends
#_______________________________________________________

prevHIV = matrix(NA, nrow = n_sim, ncol =250)
prevHIVagegrp = array(NA, dim = c(7,length(1800:2049), n_sim))
incHIV = matrix(NA, nrow = n_sim, ncol =250)
ARTcov = matrix(NA, nrow = n_sim, ncol =250)

for(i in 1:n_sim) {
  load(file = paste("90days\\", "results_sim_",  paste("sim", i, sep = "_"),   sep = ""))
  inf_pp = results_sim$inf_pp[results_sim$inf_pp$sex == 1,]
  mort_pp = results_sim$mort_pp[results_sim$mort_pp$sex == 1,]
  
  tabinc = tapply(inf_pp$pyears, list(inf_pp$agegrp, inf_pp$annee), sum)
  tabmort =  tapply(mort_pp$pyears, list(mort_pp$agegrp, mort_pp$annee), sum)
  tabmort = tabmort[rownames(tabmort) %in% rownames(tabinc), colnames(tabmort) %in% colnames(tabinc)]
  
  # trends in HIv prevalence
  prevHIV[i,] =  c(rep(0,200), approx(as.numeric(colnames(tabinc))+0.5, apply(tabmort - tabinc, 2, sum)/apply(tabmort, 2, sum), xout = 2000:2049 + 0.5)$y)
  for(k in 1:length(seq(15, 45, 5))){
    prevHIVagegrp[k,,i] = c(rep(0, 200), approx(as.numeric(colnames(tabinc))+0.5, (tabmort[k,] - tabinc[k,])/tabmort[k,], xout = 2000:2049 + 0.5)$y)
  }
  
  # Trends in HIV incidence
  incHIV[i,] =  c(rep(0,200), approx(as.numeric(colnames(tabinc))+0.5, 
                                       tapply(inf_pp$event, list(inf_pp$annee), sum)/tapply(inf_pp$pyears, list(inf_pp$annee), sum)
                                      , xout = 2000:2049 + 0.5)$y)
  if(all_param$coverage_ART[i] != 0){
     ARTcov[i,] =  results_sim$art_cov
  }
  print(i)
}

# Trends in PMTCT coverage
PMTCTcov = matrix(NA, nrow = n_sim, ncol =250)
for(i in which(all_param$coverage_ART > 0)) {
  load(file = paste("90days\\", "results_sim_",  paste("sim", i, sep = "_"),   sep = ""))
   PMTCTcov[i,] =  results_sim$pmtc_cov[names(results_sim$pmtc_cov) %in% 1800:2049]
}  
save(prevHIV, file = "90days\\prevHIV.rda", ascii = FALSE)
save(prevHIVagegrp, file = "90days\\prevHIVagegrp.rda", ascii = FALSE)
save(incHIV, file = "90days\\incHIV.rda", ascii = FALSE)
save(ARTcov, file = "90days\\ARTcov.rda", ascii = FALSE)
save(PMTCTcov, file = "90days\\PMTCTcov.rda", ascii = FALSE)

load(file = "90days\\prevHIV.rda")
load(file = "90days\\prevHIVagegrp.rda")
load(file = "90days\\incHIV.rda")
load(file = "90days\\ARTcov.rda")
load(file = "90days\\PMTCTcov.rda")

# trends in ART
 plot(artpredictlow$x, artpredictlow$y, type = 'l', lwd = 2, ylim = c(0, 90))
 points(artpredicthigh$x, artpredicthigh$y, type = 'l', lwd = 2, col = 2)
 for(i in which(all_param$coverage_ART == 1)){
 points(2025:2050-25, ARTcov[i,225:250]*100, type = 'l', col = 1)
 }
 for(i in which(all_param$coverage_ART == 2)){
   points(2025:2050-25, ARTcov[i,225:250]*100, type = 'l', col = 2)
 }
# trends in PMTCT
 plot(pmtctpredictlow$x, pmtctpredictlow$y, type = 'l', lwd = 2, ylim = c(0, 100))
 points(pmtctpredicthigh$x, pmtctpredicthigh$y, type = 'l', lwd = 2, col = 2)
 for(i in which(all_param$coverage_ART == 1)){
   points(2025:2050-25, PMTCTcov[i,225:250]*100, type = 'l', col = 1)
 }
 for(i in which(all_param$coverage_ART == 2)){
   points(2025:2050-25, PMTCTcov[i,225:250]*100, type = 'l', col = 2)
 }
# Prevalence at the time of birth of children (among women 15-49), annual estimates
prevHIVbirth = array(0, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
for(k in 1:9){
  prevHIVbirth[k,,] = apply(prevHIV, 1, function(x) approx(1800:2049 + 0.5, x, 1800:2049 + 0.5 -seq(2.5, 47.5, 5)[k])$y)
}
# Prevalence at the time of birth of children (among women 15-49), estimates for 5-year periods, so 2047.5 is our last estimate for the period 2045-2050
prevHIVbirth5 = array(0, dim = c(length(seq(0, 45, 5)), length(seq(1800,2045, 5)+2.5), n_sim))
for(k in 1:9){
  prevHIVbirth5[k,,] = apply(prevHIV, 1, function(x) approx(1800:2049 + 0.5, x, seq(1800,2045, 5)+2.5 -seq(2.5, 47.5, 5)[k])$y)
}
# Prevalence at time t (among women 15-49), estimates for 5-year periods
prevHIV5 = array(0, dim = c(length(seq(0, 45, 5)), length(seq(1800,2045, 5)+2.5), n_sim))
for(k in 1:9){
  prevHIV5[k,,] = apply(prevHIV, 1, function(x) approx(1800:2049 + 0.5, x, seq(1800,2045, 5)+2.5)$y)
}
# ART coverage at time of birth, annual estimates
ARTcov[is.na(ARTcov)]= 0
ARTcovbirth = array(0, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
for(k in 1:9){
  ARTcovbirth[k,,] = apply(ARTcov, 1, function(x) approx(1800:2049 + 0.5, x, 1800:2049 + 0.5 -seq(2.5, 47.5, 5)[k])$y)
}
# PMTCT coverage at time of birth; annual estimates
PMTCTcov[is.na(PMTCTcov)]= 0
PMTCTcovbirth = array(0, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
for(k in 1:9){
  PMTCTcovbirth[k,,] = apply(PMTCTcov, 1, function(x) approx(1800:2049 + 0.5, x, 1800:2049 + 0.5 -seq(2.5, 47.5, 5)[k])$y)
}
# Orphan prevalence, estimates for 5-year periods
orphyear5 = array(0, dim = c(length(seq(0, 45, 5)), length(seq(1800,2045, 5)+2.5), n_sim))
for(k in 1:9){
  #orphyear5[k,,] = apply(orph[k,,], 2, function(x) approx(1900:2049 + 0.5, x, seq(1900,2045, 5)+2.5)$y)
  orphyear5[k,,] = apply(orph[k,,], 2, function(x) rollapply(x, width = 5, mean, by=5))
}

# MACB at time of birth, annual estimates
MACBbirth = array(0, dim = c(length(seq(0, 45, 5)), length(1800:2049), n_sim))
for(k in 1:9){
  MACBbirth[k,,] = apply(MACB, 1, function(x) approx(1800:2049 + 0.5, x, c(1800:2049) + 0.5 -seq(2.5, 47.5, 5)[k])$y)
}

#________________________________________________________________________________________________________     
# 3 - TIMAEUS 1992 COEFFICIENTS, NO ADJUSTMENT, ONE CENSUS/SURVEY
#________________________________________________________________________________________________________     

# Underlying life table survivorship, for 5-year periods, centered on mid-period (2047.5 is the period for the last estimate)
np25 = array(NA, dim = c(length(seq(5, 45, 5)), length(seq(1800,2045, 5)), n_sim))
dimnames(np25)[[1]] = seq(10, 50, 5)
for(k in 1:length(seq(10, 50, 5))){
  n = seq(10, 50, 5)[k]
  for(i in 1:n_sim){
    np25[k, ,i] = sapply(mort[[i]], function(x) x$lt[paste(25+n),'lx']/x$lt[paste(25),'lx'])[1:50]
}}

# Brass & Bamgboye time location, reference year for annual estimates  
  Z_location = read.xls("C:\\Projets_en_cours\\ORPHANHOOD_ESTIMATES\\Orphanhood AIDS\\Data\\Z_location.xls")
  tn = array(NA, dim = c(length(seq(10, 50, 5)), length(1800:2049), n_sim))
  dimnames(tn)[[1]] = seq(10, 50, 5)
  for(k in 1:length(seq(10, 50, 5))){
    tn[k,,] = (seq(10, 50, 5)[k] /2)*(1 - (log(1-orph10[k+1,,])/3 + matrix(approx(Z_location$age, y = Z_location$z, xout = t(MACB)+seq(10, 50, 5)[k] , method="linear")$y, ncol = n_sim)  + 0.0037*(27-t(MACB))))
  }

# TM 92 coefficients with uncorrected proportions, interpolated for 5 year periods (2047.5 is our last point)
  nq25TM92 = array(0, dim = dim(np25))
  dimnames(nq25TM92)[[1]] = seq(10, 50, 5)
  for(k in 1:length(seq(10, 45, 5))){
    for(i in 1:n_sim){
      nq25TM92[k, ,i] = approx(c(1800:2049)+0.5 - tn[k,,i], 1-(coef_TM92[k,1]+coef_TM92[k,2]*MACBbirth[k,,i]+coef_TM92[k,3]*(1-orph[k+1,,i])), xout = seq(1800,2045, 5) + 2.5)$y
    }}
  # this is the risk of dying nq25 at the time seq(1800,2045, 5) + 2.5

# ratio of the estimated to true odds of surviving 
storeratioodds = matrix(NA, nrow = 8, ncol = 11)
plot(0,0, type = 'l', ylim = c(0.65, 2),xlim = c(1965, 2050)-1800, col = 'grey100', ylab = 'Ratios', main = "(b) Ratio of estimated to true odds of surviving", xlab = "Simulation years", log = 'y')
for(k in 1:8){ #
  temp = ((1-nq25TM92[k,,])/(nq25TM92[k,,]))* ((1-np25[k,,])/np25[k,,])    #100*( (1-np25[k,,all_param$n_ART == 1]) - nq25TM92[k,,all_param$n_ART == 1])/(1-np25[k,,all_param$n_ART == 1])
  #temp = ((1-nq25TM92[k,,all_param$n_ART == 1])/(nq25TM92[k,,all_param$n_ART == 1]))* ((1-np25[k,,all_param$n_ART == 1])/np25[k,,all_param$n_ART == 1])    #100*( (1-np25[k,,all_param$n_ART == 1]) - nq25TM92[k,,all_param$n_ART == 1])/(1-np25[k,,all_param$n_ART == 1])
  temp[!is.na(temp) & temp == Inf] = NA
  points(seq(1800,2045, 5)+2.5-1800, apply(temp, 1, function(x) median(x, na.rm = T)), type = 'l', col = brewer.pal(9, "Greens")[9:1][k+1], lwd = 2)
  
  storeratioodds[k,] = (apply(temp, 1, function(x) median(x, na.rm = T))[c(seq(1800,2045, 5)+2.5-1800) %in% seq(192.5, 242.5, 5)])  
}
abline(h = 1, lty = 2)
abline(v = 200)
text(197, 1.3, srt = 90, "HIV onset")
legend("topleft", c( expression(paste(""[10], "p"[25])),expression(paste(""[15], "p"[25])),expression(paste(""[20], "p"[25])),expression(paste(""[25], "p"[25])),expression(paste(""[30], "p"[25])),expression(paste(""[35], "p"[25])),expression(paste(""[40], "p"[25])),expression(paste(""[45], "p"[25]))), bty = "n",col = brewer.pal(9, "Greens")[9:2], lty = 1, pch = NA, lwd = 2)

#________________________________________________________________________________________________________     
# 4 - TIMAEUS & NUNN 1997, WITH BRASS & BAMGBOYE TIME LOCATION
#________________________________________________________________________________________________________     
 
# First adjust the proportions
# if F = 0.8 and h = 0.25
# Timaeus 2013: (1-(1-(1/3))*0.75)->0.5
orphadj = array(0, dim(orph))
orphadj[2,,] = 1- (1-orph[2,,])*(1-(0.5/2)*prevHIVbirth[2,,])
orphadj[3,,] = 1-  (1-orph[3,,])*(1-(0.5*3/4)*prevHIVbirth[3,,])
for(k in 4:length(seq(5, 50, 5))){
  orphadj[k,,] =  1-((1-orph[k,,])*(1-(0.5)*prevHIVbirth[k,,]))
}
# Then apply the coefficients
coef_TN97 = data.frame(rbind(c(-0.3611, 0.00125, 1.2974), c(-0.4030, 0.00222, 1.3732),
                             c(-0.212, 0.00372, 1.1342), c(-0.2389, 0.00586, 1.1131),
                             c(-0.2513, 0.00885, 1.0223))); names(coef_TN97) = c("a","b", "c"); rownames(coef_TN97) = seq(10, 30, 5)
# we need tnadj based on new proportions to date the estimates
orph10adj = array(0, dim(orph))
for(i in 1:9){
  orph10adj[i,,] = ((orphadj[i,,]*PY[i,,])+(orphadj[i+1,,]*PY[i+1,,]))/(PY[i,,] + PY[i+1,,])
}
tnadj = array(NA, dim = c(length(seq(10, 50, 5)), length(1800:2049), n_sim))
dimnames(tnadj)[[1]] = seq(10, 50, 5)
for(k in 1:length(seq(10, 50, 5))){
  tnadj[k,,] = (seq(10, 50, 5)[k] /2)*(1 - (log(1-orph10adj[k+1,,])/3 + matrix(approx(Z_location$age, y = Z_location$z, xout = t(MACB)+seq(10, 50, 5)[k] , method="linear")$y, ncol = n_sim)  + 0.0037*(27-t(MACB))))
}
# Now we can compute the estimates obtained if we used Timaeus and Nunn with adjusted proportions, based on a single set
# This is again referring to mid periods, with 2047.5 as our last point
nq25TN97 = array(0, dim = dim(np25))
dimnames(nq25TN97)[[1]] = seq(10, 50, 5)
for(k in 1:length(seq(10, 30, 5))){
  for(i in 1:n_sim){
    nq25TN97[k, ,i] = approx(c(1800:2049)+0.5 - tnadj[k,,i], 1-(coef_TN97[k,1]+coef_TN97[k,2]*MACBbirth[k,,i]+coef_TN97[k,3]*(1-orphadj[k+1,,i])), xout = seq(1800,2045, 5) + 2.5)$y
  }}
# replace with nq25TM92 when HIV prevalence at birth was lower than 5%
stopifnot(dim(prevHIVbirth5[2:10,,]) == dim(nq25TN97))
nq25TN97[prevHIVbirth5[2:10,,] <  0.05 & !is.na(prevHIVbirth5[2:10,,])] = nq25TM92[prevHIVbirth5[2:10,,] <  0.05 & !is.na(prevHIVbirth5[2:10,,])]

# ratio of the estimated to true odds of surviving 
storeratioodds = matrix(NA, nrow = 8, ncol = 11)
#title <- expression(atop("(f) Ratio of estimated to true odds of surviving", "Tim\u00E6us and Nunn (1997) method"))
plot(0,0, type = 'l', ylim = c(0.6, 2),xlim = c(1965, 2050)-1800, col = 'grey100', ylab = 'Ratios', main = "(b) Ratio of estimated to true odds of surviving", xlab = "Simulation years", log = 'y')
for(k in 1:8){ #
  temp = ((1-nq25TN97[k,,])/(nq25TN97[k,,]))* ((1-np25[k,,])/np25[k,,])    #100*( (1-np25[k,,all_param$n_ART == 1]) - nq25TM92[k,,all_param$n_ART == 1])/(1-np25[k,,all_param$n_ART == 1])
  temp[!is.na(temp) & temp == Inf] = NA
  points(seq(1800,2045, 5)+2.5-1800, apply(temp, 1, function(x) median(x, na.rm = T)), type = 'l', col = brewer.pal(9, "Greens")[9:1][k+1], lwd = 2)
  
  storeratioodds[k,] = (apply(temp, 1, function(x) median(x, na.rm = T))[c(seq(1800,2045, 5)+2.5-1800) %in% seq(192.5, 242.5, 5)])  
}
abline(h = 1, lty = 2)
abline(v = 200)
text(197, 1.3, srt = 90, "HIV onset")
legend("topleft", c( expression(paste(""[10], "p"[25])),expression(paste(""[15], "p"[25])),expression(paste(""[20], "p"[25])),expression(paste(""[25], "p"[25])),expression(paste(""[30], "p"[25])),expression(paste(""[35], "p"[25])),expression(paste(""[40], "p"[25])),expression(paste(""[45], "p"[25]))), bty = "n",col = brewer.pal(9, "Greens")[9:2], lty = 1, pch = NA, lwd = 2)

#________________________________________________________________________________________________________     
# 5 - NEW ADJUSTMENT TO PROPORTIONS
#________________________________________________________________________________________________________     

# First test a revised adjustment on proportions, accounting for the effect of PMTCT on h and ART on F
orphadjart = array(0, dim(orph))
orphadjart[2,,] = 1- (1-orph[2,,])*(1-((1-(1-0.33*(1-PMTCTcovbirth[2,,]))*(1-0.25*(1-ARTcovbirth[2,,])))/2)*prevHIVbirth[2,,])
orphadjart[3,,] = 1-  (1-orph[3,,])*(1-((1-(1-0.33*(1-PMTCTcovbirth[3,,]))*(1-0.25*(1-ARTcovbirth[3,,])))*3/4)*prevHIVbirth[3,,])
for(k in 4:length(seq(5, 50, 5))){
  orphadjart[k,,] =  1-((1-orph[k,,])*(1-(1-(1-0.33*(1-PMTCTcovbirth[k,,]))*(1-0.25*(1-ARTcovbirth[k,,])))*prevHIVbirth[k,,]))
}


# Compare all unbiased proportions with the biased proportions 
orphunbiased = array(data = NA, dim = dim(orph))
for(i in 1:n_sim){
cat(i)
load(file = paste(getwd(), "/90days/", "results_",  paste("sim", i, sep = "_"),  "novert.rda", sep = ""))
temporph = results_sim$orph_m_pp
temporph$annee = factor(x = temporph$annee, 1800:2050, labels = 1800:2050)
tempnonorph = results_sim$nonorph_m_pp
tempnonorph$annee = factor(x = tempnonorph$annee, 1800:2050, labels = 1800:2050)
pyearsorph = tapply(temporph$pyears, list(temporph$agegrp, temporph$annee), sum)[1:11,]
pyearsnonorph = tapply(tempnonorph$pyears, list(tempnonorph$agegrp, tempnonorph$annee), sum)[1:11,]
orphunbiased[,,i] = pyearsorph[1:10,1:250]/(pyearsorph[1:10,1:250]+pyearsnonorph[1:10,1:250])
}

library(caret)
# Candidate models
predicterrors = as.list(1:9)
for(k in 1:9){ 
  cat(paste( "\n", k, "-"))
  predicterrors[[k]] = as.data.frame(matrix(NA, nrow = 7, ncol = 3))
  names(predicterrors[[k]]) = c("formula", "RMSE_in", "RMSE_out")
  predicterrors[[k]]$formula = c("bias ~  prevHIVreg",
                                 "bias ~  prevHIVbirthreg",
                                 "bias ~  prevHIVbirthreg + ARTcovreg",
                                 "bias ~  prevHIVbirthreg + PMTCTcovbirthreg",
                                 "bias ~  prevHIVbirthreg + PMTCTcovbirthreg + ARTcovreg",
                                 "bias ~  PMTCTgapbirthreg",
                                 "bias ~ PMTCTgapbirthreg + ARTcovreg")
  predicterrors[[k]]$age = seq(0, 40, 5)[k]
  
  bias = as.vector(unlist((1- orphunbiased[k,seq(1, 250, 5),])/(1-orph[k,seq(1, 250, 5),])))
  prevHIVbirthreg = as.vector(unlist(prevHIVbirth[k,seq(1, 250, 5),])); prevHIVbirthreg[is.na(prevHIVbirthreg)] = 0
  MACBbirthreg = as.vector(unlist(MACBbirth[k,seq(1, 250, 5),]))
  prevHIVreg = as.vector(unlist(t(prevHIV[,seq(1, 250, 5)])))
  PMTCTcovbirthreg = as.vector(unlist(PMTCTcovbirth[k,seq(1, 250, 5),])); PMTCTcovbirth[is.na(PMTCTcovbirth)] = 0
  ARTcovbirthreg = as.vector(unlist(ARTcovbirth[k,seq(1, 250, 5),]))
  ARTcovreg = as.vector(unlist(t(ARTcov[,seq(1, 250, 5)])))
  PMTCTgapbirthreg = prevHIVbirthreg*(1-PMTCTcovbirthreg)
  ARTgapreg = prevHIVreg*(1-ARTcovreg)
  sim_id = as.vector(matrix(all_param$no_sim, nrow = 50, ncol = 576, byrow = TRUE))
  
  temp = data.frame(bias, prevHIVbirthreg, prevHIVreg, PMTCTcovbirthreg, ARTcovbirthreg, ARTcovreg, PMTCTgapbirthreg, ARTgapreg, sim_id)
  temp = temp[!is.na(temp$bias),]
  training = temp[temp$sim_id %in% sample(1:576, 461, replace  = FALSE),]
  test = temp[temp$sim_id %in% training$sim == FALSE,]
  for(j in 1:7){
    cat(j)
    model1= lm(predicterrors[[k]]$formula[j], data = training)
    # in-sample
    training$predict = NA
    training$predict = predict(model1)
    predicterrors[[k]]$RMSE_in[j] = round(RMSE(training$bias, training$predict),5)
    # out-of-sample
    test$predict = NA
    test$predict = predict(model1, newdata = test)
    predicterrors[[k]]$RMSE_out[j] = round(RMSE(test$bias, test$predict), 5)
   }
}
tablepredictors = rbind(cbind(predicterrors[[1]][,c(2:3)], predicterrors[[2]][,c(2:3)], predicterrors[[3]][,c(2:3)]),
                        cbind(predicterrors[[4]][,c(2:3)], predicterrors[[5]][,c(2:3)], predicterrors[[6]][,c(2:3)]),
                        cbind(predicterrors[[7]][,c(2:3)], predicterrors[[8]][,c(2:3)], predicterrors[[9]][,c(2:3)]))
library(xtable)
param.latex <- print(xtable(tablepredictors, digits = 5, caption = "Prediction errors for candidate models for bias in proportions of mothers surviving", label = "param.latex"), include.rownames = F, include.colnames = F, hline.after =c(seq(7, 7*3, 7)))

# Model the bias on all sims
orphadjreg = orphadj
bias_reg = as.list(1:8)
RMSER2 = as.list(1:8)
for(k in 1:8){
print(k)
bias = as.vector(unlist((1- orphunbiased[k,seq(1, 250, 5),])/(1-orph[k,seq(1, 250, 5),])))
prevHIVbirthreg = as.vector(unlist(prevHIVbirth[k,seq(1, 250, 5),])); prevHIVbirthreg[is.na(prevHIVbirthreg)] = 0
print(range(prevHIVbirthreg, na.rm = TRUE))
MACBbirthreg = as.vector(unlist(MACBbirth[k,seq(1, 250, 5),]))
prevHIVreg = as.vector(unlist(t(prevHIV[,seq(1, 250, 5)])))
PMTCTcovbirthreg = as.vector(unlist(PMTCTcovbirth[k,seq(1, 250, 5),])); PMTCTcovbirth[is.na(PMTCTcovbirth)] = 0
ARTcovbirthreg = as.vector(unlist(ARTcovbirth[k,seq(1, 250, 5),]))
ARTcovreg = as.vector(unlist(t(ARTcov[,seq(1, 250, 5)])))
PMTCTgapbirthreg = prevHIVbirthreg*(1-PMTCTcovbirthreg)
ARTgapreg = prevHIVreg*(1-ARTcovreg)
ARTscenarioreg = as.vector(matrix(all_param$n_ART, nrow = 50, ncol = 576, byrow = TRUE))

#summary(lm(bias~prevHIVbirthreg))
#summary(lm(bias[ARTscenarioreg == 1]~prevHIVbirthreg[ARTscenarioreg == 1]))
#summary(lm(bias[ARTscenarioreg >= 1 & ARTcovreg == 0 & PMTCTcovbirthreg==0]~prevHIVbirthreg[ARTscenarioreg >= 1 & ARTcovreg == 0 & PMTCTcovbirthreg==0]))
#summary(lm(bias~prevHIVbirthreg + prevHIVreg))
#summary(lm(bias~ prevHIVreg))
#summary(lm(bias~prevHIVbirthreg + PMTCTcovbirthreg))
#summary(lm(bias~prevHIVbirthreg + ARTcovreg))
#summary(lm(bias~PMTCTgapbirthreg ))
#plot(prevHIVbirth[1,seq(1, 250, 5),]*(1-PMTCTcovbirth[1,seq(1, 250, 5),]), (1- orphunbiased[1,seq(1, 250, 5),])/(1-orph[1,seq(1, 250, 5),]), ylim = c(0.9, 1.10), pch = 19, col = colartcov[1,seq(1, 250, 5),])
#summary(lm(bias~PMTCTgapbirthreg + ARTcovreg))
#summary(lm(bias[ARTcovreg >= 0 & PMTCTcovbirthreg>=0]~PMTCTgapbirthreg[ARTcovreg >= 0 & PMTCTcovbirthreg>=0] + ARTcovreg[ARTcovreg >= 0 & PMTCTcovbirthreg>=0]))
#summary(lm(bias[ARTcovreg == 0 & PMTCTcovbirthreg==0]~PMTCTgapbirthreg[ARTcovreg == 0 & PMTCTcovbirthreg==0] + ARTcovreg[ARTcovreg == 0 & PMTCTcovbirthreg==0]))

bias_reg[[k]] = lm(bias~PMTCTgapbirthreg + ARTcovreg)
RMSER2[[k]] = c(RMSE(bias[!is.na(bias)], predict(bias_reg[[k]])), R2(bias[!is.na(bias)], predict(bias_reg[[k]])))
orphadjreg[k,,] = 1- (1-orph[k,,])*(1+coef(bias_reg[[k]])[2]*(prevHIVbirth[k,,]*(1-PMTCTcovbirth[k,,]))+coef(bias_reg[[k]])[3]*t(ARTcov))

par(mfrow = c(1,3))
plot(prevHIVbirth[k,seq(1, 250, 5),]*(1-PMTCTcovbirth[k,seq(1, 250, 5),]), (1- orphunbiased[k,seq(1, 250, 5),])/(1-orph[k,seq(1, 250, 5),]), ylim = c(0.8, 1.2), pch = 19, col = colartcov[k,seq(1, 250, 5),], main = paste(seq(5, 45, 5)[k], 'Unadjusted'), log = 'y', xlab = "HIV*(1-PMTCT)", ylab = 'Sn(unbiased)/Sn'); abline(h = 1)
text(0.05, 0.85, round(sqrt(sum(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orph[k,seq(1, 250, 5),]))^2, na.rm = TRUE)/sum(is.na(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orph[k,seq(1, 250, 5),]))^2))),5))
plot(prevHIVbirth[k,seq(1, 250, 5),]*(1-PMTCTcovbirth[k,seq(1, 250, 5),]), (1- orphunbiased[k,seq(1, 250, 5),])/(1-orphadj[k,seq(1, 250, 5),]), ylim = c(0.8, 1.2), pch = 19, col = colartcov[k,seq(1, 250, 5),], main = paste(seq(5, 45, 5)[k], 'Adjusted with T&N97'), log = 'y', xlab = "HIV*(1-PMTCT)", ylab = 'Sn(unbiased)/Sn'); abline(h = 1)
text(0.05, 0.85, round(sqrt(sum(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orphadj[k,seq(1, 250, 5),]))^2, na.rm = TRUE)/sum(is.na(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orphadj[k,seq(1, 250, 5),]))^2))),5))
plot(prevHIVbirth[k,seq(1, 250, 5),]*(1-PMTCTcovbirth[k,seq(1, 250, 5),]), (1- orphunbiased[k,seq(1, 250, 5),])/(1-orphadjreg[k,seq(1, 250, 5),]), ylim = c(0.8, 1.2), pch = 19, col = colartcov[k,seq(1, 250, 5),], main =paste(seq(5, 45, 5)[k],  'Adjusted with reg. coef'), log = 'y', xlab = "HIV*(1-PMTCT)", ylab = 'Sn(unbiased)/Sn'); abline(h = 1)
text(0.05, 0.85, round(sqrt(sum(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orphadjreg[k,seq(1, 250, 5),]))^2, na.rm = TRUE)/sum(is.na(((1- orphunbiased[k,seq(1, 250, 5),]) - (1-orphadjreg[k,seq(1, 250, 5),]))^2))),5))
}

coeftomanuscript = cbind(seq(0, 30, 5), round(t(sapply(bias_reg, function(x) x$coef))[1:7,], 4), round(do.call(rbind, RMSER2)[1:7,], 4))
library(xtable)
param.latex <- print(xtable(coeftomanuscript, digits = 4, caption = "Coefficients for adjusting proportions of mothers surviving for HIV-related bias", label = "param.latex"), include.rownames = F, include.colnames = F, hline.after = c(0,7))

# MACB for 5-year periods, centered on midpoint 2002.5, 2007.5, etc.
MACByear5 = array(0, dim = c(length(seq(0, 45, 5)), length(seq(1800,2045, 5)+2.5), n_sim))
MACByear5[1,,] = apply(MACB, 1, function(x) rollapply(x, width = 5, mean, by=5))
for(k in 2:10){MACByear5[k,,] = MACByear5[1,,]}
# MACB at birth for 5-year periods, centered on midpoint 2002.5, 2007.5, etc.
MACBbirthyear5 = array(0, dim = c(length(seq(0, 45, 5)), length(seq(1800,2045, 5)+2.5), n_sim))
for(k in 1:10){
  MACBbirthyear5[k,,] = apply(MACBbirth[k,,], 2, function(x) rollapply(x, width = 5, mean, by=5))
}

#________________________________________________________________________________________________________     
# 6 - DEVELOPMENT OF NEW COEFFICIENTS
#________________________________________________________________________________________________________     

# Preston method for calculating growth-corrected proportions, first without adjustment for HIV-related bias
cutpoints = which(floor((1800:2049)/5) - ((1800:2049)/5) == 0)
ecart = diff(cutpoints)[1]

GR_prop_not_orph = array(0, c(dim(orph[,,1][, cutpoints]), n_sim)-c(0,1,0))
dimnames(GR_prop_not_orph)[[1]] = seq(5, 50, 5)
dimnames(GR_prop_not_orph)[[2]] = as.numeric(dimnames(orph[,,1])[[2]][cutpoints])[1:dim(GR_prop_not_orph)[2]]
dimnames(GR_prop_not_orph)[[3]] = 1:n_sim

GR_prop_orph = array(0, c(dim(orph[,,1][, cutpoints]), n_sim)-c(0,1,0))
dimnames(GR_prop_orph)[[1]] = seq(5, 50, 5)
dimnames(GR_prop_orph)[[2]] = as.numeric(dimnames(orph[,,1])[[2]][cutpoints])[1:dim(GR_prop_orph)[2]]
dimnames(GR_prop_orph)[[3]] = 1:n_sim

for (j in 1:dim(GR_prop_orph)[2]) {
  GR_prop_not_orph[,j,] = log((1-orph[, cutpoints[j+1],])/(1-orph[, cutpoints[j],]) )/ecart
  GR_prop_orph[1,j,] = 1-(sqrt((1-orph[1, cutpoints[j+1],])*(1-orph[1, cutpoints[j],]))*exp(GR_prop_not_orph[1,j,] * 2.5))
  GR_prop_orph[2,j,] = 1-(sqrt((1-orph[2, cutpoints[j+1],])*(1-orph[2, cutpoints[j],]))*exp(GR_prop_not_orph[2,j,] * 2.5 + GR_prop_not_orph[1,j,] * 5))
  
  for(k in 3:length(GR_prop_orph[,j,1])) {
    GR_prop_orph[k,j,] = 1-(sqrt((1-orph[k, cutpoints[j+1],])*(1-orph[k, cutpoints[j],]))*exp(GR_prop_not_orph[k,j,] * 2.5 + apply(GR_prop_not_orph[1:(k-1),j,], 2, sum) * 5))
  }
}
# This provides growth-corrected proportions for mid-periods centered on 1803, 1808...2043
# As we used two sets of proportions referring to, for example 2040.5 and 2044.5
years_GR = (diff(as.numeric(dimnames(GR_prop_orph)[[2]]))/2)[1] + as.numeric(dimnames(GR_prop_orph)[[2]][1:dim(GR_prop_not_orph)[2]]) 

# Preston method based on growth-corrected prop with adjusted proportions   
GR_prop_not_orphadjreg = array(0, c(dim(orphadjreg[,,1][, cutpoints]), n_sim)-c(0,1,0))
dimnames(GR_prop_not_orphadjreg)[[1]] = seq(5, 50, 5)
dimnames(GR_prop_not_orphadjreg)[[2]] = as.numeric(dimnames(orph[,,1])[[2]][cutpoints])[1:dim(GR_prop_not_orphadjreg)[2]]
dimnames(GR_prop_not_orphadjreg)[[3]] = 1:n_sim
      
GR_prop_orphadjreg = array(0, c(dim(orphadjreg[,,1][, cutpoints]), n_sim)-c(0,1,0))
dimnames(GR_prop_orphadjreg)[[1]] = seq(5, 50, 5)
dimnames(GR_prop_orphadjreg)[[2]] = as.numeric(dimnames(orph[,,1])[[2]][cutpoints])[1:dim(GR_prop_orphadjreg)[2]]
dimnames(GR_prop_orphadjreg)[[3]] = 1:n_sim
      
for (j in 1:dim(GR_prop_orphadjreg)[2]) {
    GR_prop_not_orphadjreg[,j,] = log((1-orphadjreg[, cutpoints[j+1],])/(1-orphadjreg[, cutpoints[j],]) )/ecart
    GR_prop_orphadjreg[1,j,] = 1-(sqrt((1-orphadjreg[1, cutpoints[j+1],])*(1-orphadjreg[1, cutpoints[j],]))*exp(GR_prop_not_orphadjreg[1,j,] * 2.5))
    GR_prop_orphadjreg[2,j,] = 1-(sqrt((1-orphadjreg[2, cutpoints[j+1],])*(1-orphadjreg[2, cutpoints[j],]))*exp(GR_prop_not_orphadjreg[2,j,] * 2.5 + GR_prop_not_orphadjreg[1,j,] * 5))
      
    for(k in 3:length(GR_prop_orphadjreg[,j,1])) {
    GR_prop_orphadjreg[k,j,] = 1-(sqrt((1-orphadjreg[k, cutpoints[j+1],])*(1-orphadjreg[k, cutpoints[j],]))*exp(GR_prop_not_orphadjreg[k,j,] * 2.5 + apply(GR_prop_not_orphadjreg[1:(k-1),j,], 2, sum) * 5))
    }
  }

# Dataframe of proportions, mortality rates, HIV prevalence, ART coverage, PMTCT, etc.
dfsim = as.list(1)
for(k in 1:8){
  print(k)
  dfsim[[k]] = as.list(1)
  for(i in all_param$no_sim){
    dfsim[[k]][[i]] = data.frame(sim = i,
                                 time = c(seq(1800,2045, 5)+2.5)[21:50],
                                 np25 = np25[k,21:50,i],
                                 Sn_5 = approx(years_GR,1- GR_prop_orph[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 Sn_5adj = approx(years_GR,1- GR_prop_orphadjreg[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 # MACB at time at birth
                                 MACBbirth = approx(1800:2049+0.5, MACBbirth[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 # prevalence of HIV at time t (mid-period), t-2.5 (first survey/census), t+2.5 (second survey/census) and at birth
                                 prevHIV= approx(1800:2049+0.5, prevHIV[i,], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 prevHIVfirst= approx(1800:2049+0.5, prevHIV[i,], xout = c(seq(1800,2045, 5)+0)[21:50])$y,
                                 prevHIVsecond= approx(1800:2049+0.5, prevHIV[i,], xout = c(seq(1800,2045, 5)+5)[21:50])$y,
                                 prevHIVbirth = approx(1800:2049+0.5, prevHIVbirth[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 # ART coverage at time t and and birth
                                 ARTcov = approx(1800:2049+0.5, ARTcov[i,], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 ARTcovfirst = approx(1800:2049+0.5, ARTcov[i,], xout = c(seq(1800,2045, 5)+0)[21:50])$y,
                                 ARTcovsecond = approx(1800:2049+0.5, ARTcov[i,], xout = c(seq(1800,2045, 5)+5)[21:50])$y,
                                 ARTcovbirth = approx(1800:2049+0.5, ARTcovbirth[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 # PMTCT coverage at time t and at birth
                                 PMTCTcov = approx(1800:2049+0.5, PMTCTcov[i,], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 PMTCTcovbirth = approx(1800:2049+0.5, PMTCTcovbirth[k+1,,i], xout = c(seq(1800,2045, 5)+2.5)[21:50])$y,
                                 # Expanding or receding epidemic?
                                 receding =round(approx(1800:2049+0.5, as.numeric(c(1:250)> which.max(prevHIV[i,])), xout = c(seq(1800,2045, 5)+2.5)[21:50])$y))
    
  }
  dfsim[[k]] = do.call(rbind, dfsim[[k]])
  dfsim[[k]] = dfsim[[k]][!is.na(dfsim[[k]]$Sn_5),]
}

temp = dfsim[[1]]
temp$logitSn_5 = gtools::logit(temp$Sn_5)
temp$logitSn_5adj = gtools::logit(temp$Sn_5adj)
temp$logitnp25 = gtools::logit(temp$np25)
temp$ARTgapbirth = temp$prevHIVbirth*(1-temp$ARTcovbirth)
temp$ARTgap = temp$prevHIV*(1-temp$ARTcov)

# Standard regression, pre-HIV, similar to Timaeus 1992
model1 = lm(np25 ~  MACBbirth + Sn_5adj, data = temp[temp$time < 2000,])
points((seq(0.4, 0.99, 0.01)),-0.297859 +  1.273258*((seq(0.4, 0.99, 0.01))) + 0.000881*25, type = 'l',  col = brewer.pal(5, "Reds")[3], lwd = 2)
predicted = predict(model1, newdata = temp[temp$time < 2000,])
observed = temp[temp$time < 2000,]$np25
plot(predicted, observed, cex = 0.5); abline(0, 1)
#library(caret)
R2(predicted, observed)
RMSE(predicted, observed)

# Standard regression, HIV, no treatment
model2 = lm(np25 ~  MACBbirth + Sn_5adj + prevHIV, data = temp[temp$ARTcov == 0 & temp$PMTCTcov == 0 ,])
summary(model2)
predicted = predict(model2, newdata = temp[ temp$ARTcov == 0 & temp$PMTCTcov == 0 ,])
observed =  temp[ temp$ARTcov == 0 & temp$PMTCTcov == 0 ,]$np25
plot(observed, predicted, cex = 0.5); abline(0, 1)
R2(observed, predicted)
RMSE(observed, predicted)

model3 = lm(np25 ~  MACBbirth + Sn_5adj*ARTgap, data = temp)
summary(model3)
predicted = predict(model3, newdata = temp)
observed = temp$np25
plot(observed, predicted, cex = 0.5); abline(0, 1)
R2(observed, predicted)
RMSE(observed, predicted)

model3 = lm(np25 ~  MACBbirth + Sn_5adj + prevHIV*ARTcov, data = temp)
summary(model3)
predicted = predict(model3, newdata = temp)
observed = temp$np25
plot(observed, predicted, cex = 0.5); abline(0, 1)
R2(observed, predicted)
RMSE(observed, predicted)
sqrt((sum((observed - predicted)^2))/length(observed)) # Max RMSE over 5-year periods
mean((observed - predicted)/predicted) # max mean relative error over 5-year periods
# min median ratio of estimated to true odds of surviving 
# max median ratio of estimated to true odds of surviving 

predicterrors = as.list(1:8)
# First model
for(k in 1:8){ 
  cat(k)
  predicterrors[[k]] = as.data.frame(matrix(NA, nrow = 9, ncol = 7))
  names(predicterrors[[k]]) = c("formula", "RMSE_in", "Max_in", "Min_in", "RMSE_out", "Max_out", "Min_out")
  predicterrors[[k]]$formula = c("np25 ~  MACBbirth + Sn_5adj",
                                "np25 ~  MACBbirth + Sn_5adj + prevHIV + ARTcov",
                                "np25 ~  MACBbirth + Sn_5adj + prevHIV + ARTcov + prevHIVbirth + ARTcovbirth",
                                "np25 ~  MACBbirth + Sn_5adj + prevHIV + ARTcov + prevHIVdelta + ARTcovdelta",
                                "np25 ~  MACBbirth + Sn_5adj + ARTgap",
                                "np25 ~  MACBbirth + Sn_5adj*ARTgap",
                                "np25 ~  MACBbirth + Sn_5adj + ARTgap + ARTgapdelta",
                                "two steps: HIV/ ARTgap",
                                "two steps: HIV + HIVdelta/ ARTgap + ARTdelta")
  predicterrors[[k]]$age = seq(10, 40, 5)[k]
  temp = dfsim[[k]]
  
  temp$ARTgap = temp$prevHIV*(1-temp$ARTcov)
  temp$ARTgapfirst = temp$prevHIVfirst*(1-temp$ARTcovfirst)
  temp$ARTgapsecond = temp$prevHIVsecond*(1-temp$ARTcovsecond)
  temp$ARTgapdelta = temp$ARTgapsecond-temp$ARTgapfirst
  
  temp$prevHIVdelta = temp$prevHIVfirst-temp$prevHIVsecond
  temp$ARTcovdelta = temp$ARTcovfirst-temp$ARTcovsecond
  
  training = temp[temp$sim %in% sample(1:576, 461, replace  = FALSE),]
  test = temp[temp$sim %in% training$sim == FALSE,]
  # One step
  for(j in 1:7){
  cat(j)
  model1= lm(predicterrors[[k]]$formula[j], data = training)
  # in-sample
  training$predict = NA
  training$predict = predict(model1)
  predicterrors[[k]]$RMSE_in[j] = RMSE(training$np25, training$predict)
  training$ratio = ((training$predict)/(1-(training$predict)))*((1-training$np25)/training$np25)
  predicterrors[[k]]$Max_in[j] = max(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_in[j] = min(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  # out-of-sample
  test$predict = NA
  test$predict = predict(model1, newdata = test)
  predicterrors[[k]]$RMSE_out[j] = RMSE(test$np25, test$predict)
  test$ratio = ((test$predict)/(1-(test$predict)))*((1-test$np25)/test$np25)
  predicterrors[[k]]$Max_out[j] = max(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_out[j] = min(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  }
  
  # two-steps
  j = 8
  step1= lm(np25 ~  MACBbirth + Sn_5adj + prevHIV, data = training[training$ARTcov == 0,])
  step2= lm(np25 ~  MACBbirth + Sn_5adj + ARTgap, data = training[training$ARTcov > 0,])
  # in-sample
  training$predict = NA
  training$predict[training$ARTcov == 0] = predict(step1)
  training$predict[training$ARTcov > 0] = predict(step2)
  predicterrors[[k]]$RMSE_in[j] = RMSE(training$np25, training$predict)
  training$ratio = ((training$predict)/(1-(training$predict)))*((1-training$np25)/training$np25)
  predicterrors[[k]]$Max_in[j] = max(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_in[j] = min(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  # out-of-sample
  test$predict = NA
  test$predict[test$ARTcov == 0] = predict(step1, newdata = test[test$ARTcov == 0,])
  test$predict[test$ARTcov > 0] = predict(step2, newdata = test[test$ARTcov > 0,])
  predicterrors[[k]]$RMSE_out[j] = RMSE(test$np25, test$predict)
  test$ratio = ((test$predict)/(1-(test$predict)))*((1-test$np25)/test$np25)
  predicterrors[[k]]$Max_out[j] = max(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_out[j] = min(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  
  # two-steps
  j = 9
  step1= lm(np25 ~  MACBbirth + Sn_5adj + prevHIV + prevHIVdelta, data = training[training$ARTcov == 0,])
  step2= lm(np25 ~  MACBbirth + Sn_5adj + ARTgap + ARTgapdelta, data = training[training$ARTcov > 0,])
  # in-sample
  training$predict = NA
  training$predict[training$ARTcov == 0] = predict(step1)
  training$predict[training$ARTcov > 0] = predict(step2)
  predicterrors[[k]]$RMSE_in[j] = RMSE(training$np25, training$predict)
  training$ratio = ((training$predict)/(1-(training$predict)))*((1-training$np25)/training$np25)
  predicterrors[[k]]$Max_in[j] = max(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_in[j] = min(tapply(training$ratio, training$time , function(x) median(x, na.rm = T)))
  # out-of-sample
  test$predict = NA
  test$predict[test$ARTcov == 0] = predict(step1, newdata = test[test$ARTcov == 0,])
  test$predict[test$ARTcov > 0] = predict(step2, newdata = test[test$ARTcov > 0,])
  predicterrors[[k]]$RMSE_out[j] = RMSE(test$np25, test$predict)
  test$ratio = ((test$predict)/(1-(test$predict)))*((1-test$np25)/test$np25)
  predicterrors[[k]]$Max_out[j] = max(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  predicterrors[[k]]$Min_out[j] = min(tapply(test$ratio, test$time , function(x) median(x, na.rm = T)))
  
  
}
# compare models in terms of
# Root mean square error	Root median square error	Mean relative error	Median relative error
# training (80%) and test (20%)
# 1. macb + Sn_5adj 
# 2. macb + Sn_5adj + HIV + ART
# 3. macb + Sn_5adj + HIV + HIVbirth + ART + ARTbirth
# 4. macb + Sn_5adj + HIV + deltaHIV + ART + deltaART
# 5. macb + Sn_5adj + HIVgap
# 6. macb + Sn_5adj*HIVgap
# 7. macb + Sn_5adj + HIVgap + deltaHIVgap
# 8. two sets of coefficients (pre ART, post ART):  macb + Sn_5adj + HIV \\ macb + Sn_5adj + HIV + ART
# 9. two sets of coefficients (pre ART, post ART):  macb + Sn_5adj + HIV + deltaHIV \\ macb + Sn_5adj + ARTgap + deltaARTgap
# so our last model is to be preferred

tablepredictors = rbind(cbind(predicterrors[[1]][,c(2:3)], predicterrors[[2]][,c(2:3)], predicterrors[[3]][,c(2:3)]),
                        cbind(predicterrors[[4]][,c(2:3)], predicterrors[[5]][,c(2:3)], predicterrors[[6]][,c(2:3)]),
                        cbind(predicterrors[[7]][,c(2:3)], predicterrors[[8]][,c(2:3)], predicterrors[[8]][,c(2:3)]))
library(xtable)
param.latex <- print(xtable(tablepredictors, digits = 5, caption = "Prediction errors for candidate models for bias in proportions of mothers surviving", label = "param.latex"), include.rownames = F, include.colnames = F, hline.after =c(seq(7, 7*3, 7)))


xtable(predicterrors[[1]][,2:7], digits= 4)


# Revised coefficients
coefset1 = as.list(1)
coefset2 = as.list(1)
plot(0,0, type = 'l', ylim = c(0.7, 2),xlim = c(1965, 2050)-1800, col = 'grey100', ylab = 'Ratios', main = expression(paste("(b) Ratio of estimated to true odds of surviving")), xlab = "Simulation years")
abline(h =1)
for(k in 1:8){ 
  temp = dfsim[[k]]
  # Pre-treatment
  temp$prevHIVdelta = temp$prevHIVsecond-temp$prevHIVfirst
  refreg1= lm(np25 ~  MACBbirth + Sn_5adj + prevHIV + prevHIVdelta, data = temp[ temp$ARTcov == 0,])
  coefset1[[k]] = refreg1
  temp$predict = NA
  temp[ temp$ARTcov == 0,]$predict = predict(refreg1)
  # With treatment
  temp$ARTgap = temp$prevHIV*(1-temp$ARTcov)
  temp$ARTgapfirst = temp$prevHIVfirst*(1-temp$ARTcovfirst)
  temp$ARTgapsecond = temp$prevHIVsecond*(1-temp$ARTcovsecond)
  temp$ARTgapdelta = temp$ARTgapsecond-temp$ARTgapfirst
  refreg2= lm(np25 ~  MACBbirth + Sn_5adj + ARTgap + ARTgapdelta, data = temp[ (temp$ARTcov > 0),])
  coefset2[[k]] = refreg2
  temp[ (temp$ARTcov > 0),]$predict = predict(refreg2)
  temp$ratio = ((temp$predict)/(1-(temp$predict)))*((1-temp$np25)/temp$np25)
  points(temp$time[temp$sim == 1]-1800, tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T)), type = 'l', col = brewer.pal(9, "Greens")[9:1][k+1], lwd = 2) 
  
  #round(tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T))[19:29],2)
  
  dfsim[[k]]$prediction = temp$predict
  if(k == 1){
    text(temp$time[temp$sim == 1][which.max(tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T)))]-1800, 1.5, round(max(tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T))), 2))
    text(temp$time[temp$sim == 1][which.min(tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T)))]-1800, 1, round(min(tapply(temp$ratio, temp$time , function(x) median(x, na.rm = T))), 2))
  }
}

preARTcoef = as.data.frame(rbind(coef(coefset1[[1]]), coef(coefset1[[2]]), coef(coefset1[[3]]), coef(coefset1[[4]]), coef(coefset1[[5]]), coef(coefset1[[6]]), coef(coefset1[[7]]), coef(coefset1[[8]])))
names(preARTcoef) = c("intercept1", "macb1", "snadj1", "hiv1", "hivdelta1")
preARTcoef$n = seq(10, 45, 5)
postARTcoef = as.data.frame(rbind(coef(coefset2[[1]]), coef(coefset2[[2]]), coef(coefset2[[3]]), coef(coefset2[[4]]), coef(coefset2[[5]]), coef(coefset2[[6]]), coef(coefset2[[7]]), coef(coefset2[[8]])))
names(postARTcoef) = c("intercept2", "macb2", "snadj2", "artgap2", "artgapdelta2")
postARTcoef$n = seq(10, 45, 5)



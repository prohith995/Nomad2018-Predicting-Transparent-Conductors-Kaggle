###
### This is my working directory.  Your will be different.
###

###
### Extremely Important: You will have to set your directories appropriately!
### Also I have already un-zipped my data.
###

###
### Upload the original data.  These are not the
###  geometry data:
###

train = read.csv("./train.csv/train.csv", header = T)
test  = read.csv("./test.csv/test.csv", header = T)

###
### Save the id's for potential later use:
###

id_train = train$id
id_test  = test$id

train$id = NULL
test$id  = NULL

###
### Transform responses so it matches the evaluation criteria.
### I will add this back later.
###

formation_energy_ev_natom = log(train$formation_energy_ev_natom + 1)
bandgap_energy_ev         = log(train$bandgap_energy_ev + 1)

train$formation_energy_ev_natom = NULL
train$bandgap_energy_ev         = NULL

###
### combine data:
###

combined = rbind(train,test)

rows_for_train = 1:2400
rows_for_test  = 2401:3000

###
### Define two functions that will help with feature extraction:
###
### The first one is geometric mean.
### The second one is just the population standard deviation.
###

gmean  = function(x) { exp(mean(log(x))) }

pop_sd = function(x)
{
   x_mean = mean(x)
   out    = sqrt(mean((x - x_mean)^2))
   out
}

###
### Now I will extract features for the train & test data:
###

###
### Get the ratios for each element:
###

al_ratio = 2*combined$percent_atom_al/5
ga_ratio = 2*combined$percent_atom_ga/5
in_ratio = 2*combined$percent_atom_in/5
o_ratio  = rep(3/5 , nrow(combined))

###
### Electron affinities (eV) for the free atoms calculated using 
###   the local density approximation (LDA):
###

al_ea = -0.2563
ga_ea = -0.1081
in_ea = -0.3125
o_ea  = -0.225633333

mean_ea = (al_ratio * al_ea) + (ga_ratio * ga_ea) + (in_ratio * in_ea) + (o_ratio * o_ea)
sd_ea   =  sqrt( al_ratio * (al_ea - mean_ea)^2 +  ga_ratio * (ga_ea - mean_ea)^2 + 
           in_ratio * (in_ea - mean_ea)^2 + o_ratio * (o_ea - mean_ea)^2 )

combined$mean_ea = mean_ea
combined$sd_ea   = sd_ea

###
### Electronegativities:
###

al_en = 1.61
ga_en = 1.81
in_en = 1.78
o_en  = 3.44

mean_en = (al_ratio * al_en) + (ga_ratio * ga_en) + (in_ratio * in_en) + (o_ratio * o_en)
sd_en   =  sqrt( al_ratio * (al_en - mean_en)^2 +  ga_ratio * (ga_en - mean_en)^2 + 
           in_ratio * (in_en - mean_en)^2 + o_ratio * (o_en - mean_en)^2 )

combined$mean_en = mean_en
combined$sd_en   = sd_en

###
### homo
###

al_homo = -2.697
ga_homo = -2.732
in_homo = -2.784
o_homo  = -2.74

mean_homo = (al_ratio * al_homo) + (ga_ratio * ga_homo) + (in_ratio * in_homo) + (o_ratio * o_homo)
sd_homo   =  sqrt( al_ratio * (al_homo - mean_homo)^2 +  ga_ratio * (ga_homo - mean_homo)^2 + 
             in_ratio * (in_homo - mean_homo)^2 + o_ratio * (o_homo - mean_homo)^2 )

combined$mean_homo = mean_homo
combined$sd_homo   = sd_homo

###
### Ionization potentials (eV) for free atoms calculated using the local density approximation (LDA)
###

al_ionpot = -5.78
ga_ionpot = -5.8182
in_ionpot = -5.5374
o_ionpot  = -5.711866667

mean_ionpot = (al_ratio * al_ionpot) + (ga_ratio * ga_ionpot) + (in_ratio * in_ionpot) + (o_ratio * o_ionpot)
sd_ionpot   =  sqrt( al_ratio * (al_ionpot - mean_ionpot)^2 +  ga_ratio * (ga_ionpot - mean_ionpot)^2 + 
             in_ratio * (in_ionpot - mean_ionpot)^2 + o_ratio * (o_ionpot - mean_ionpot)^2 )

combined$mean_ionpot = mean_ionpot
combined$sd_ionpot   = sd_ionpot

###
### Lowest unoccupied molecular orbitals (eV) for free atoms calculated using the local density approximation (LDA)
###

al_lumo = 0.368
ga_lumo = 0.13
in_lumo = 0.695
o_lumo  = 0.397666667

mean_lumo = (al_ratio * al_lumo) + (ga_ratio * ga_lumo) + (in_ratio * in_lumo) + (o_ratio * o_lumo)
sd_lumo   =  sqrt( al_ratio * (al_lumo - mean_lumo)^2 +  ga_ratio * (ga_lumo - mean_lumo)^2 + 
             in_ratio * (in_lumo - mean_lumo)^2 + o_ratio * (o_lumo - mean_lumo)^2 )

combined$mean_lumo = mean_lumo
combined$sd_lumo   = sd_lumo

###
### Atomic mass
###

al_amu = 26.9815386
ga_amu = 69.723
in_amu = 114.818
o_amu  = 15.9994

mean_amu = (al_ratio * al_amu) + (ga_ratio * ga_amu) + (in_ratio * in_amu) + (o_ratio * o_amu)
sd_amu   =  sqrt( al_ratio * (al_amu - mean_amu)^2 +  ga_ratio * (ga_amu - mean_amu)^2 + 
             in_ratio * (in_amu - mean_amu)^2 + o_ratio * (o_amu - mean_amu)^2 )

combined$mean_amu = mean_amu
combined$sd_amu   = sd_amu

###
### d orbital radius (Angstrom) for free atoms calculated using 
###   the local density approximation (LDA)
###

al_rd_max = 3.11
ga_rd_max = 2.16
in_rd_max = 1.94
o_rd_max  = 2.403333333

mean_rd_max = (al_ratio * al_rd_max) + (ga_ratio * ga_rd_max) + (in_ratio * in_rd_max) + (o_ratio * o_rd_max)
sd_rd_max   =  sqrt( al_ratio * (al_rd_max - mean_rd_max)^2 +  ga_ratio * (ga_rd_max - mean_rd_max)^2 + 
             in_ratio * (in_rd_max - mean_rd_max)^2 + o_ratio * (o_rd_max - mean_rd_max)^2 )

combined$mean_rd_max = mean_rd_max
combined$sd_rd_max   = sd_rd_max

###
### p orbital radius (Angstrom) for free atoms calculated using 
###  the local density approximation (LDA)
###

al_rp_max = 1.5
ga_rp_max = 1.33
in_rp_max = 1.39
o_rp_max  = 1.406666667

mean_rp_max = (al_ratio * al_rp_max) + (ga_ratio * ga_rp_max) + (in_ratio * in_rp_max) + (o_ratio * o_rp_max)
sd_rp_max   =  sqrt( al_ratio * (al_rp_max - mean_rp_max)^2 +  ga_ratio * (ga_rp_max - mean_rp_max)^2 + 
             in_ratio * (in_rp_max - mean_rp_max)^2 + o_ratio * (o_rp_max - mean_rp_max)^2 )

combined$mean_rp_max = mean_rp_max
combined$sd_rp_max   = sd_rp_max

###
### s orbital radius (Angstrom) for free atoms calculated using 
###  the local density approximation (LDA)
###

al_rs_max = 1.5
ga_rs_max = 1.33
in_rs_max = 1.39
o_rs_max  = 1.406666667

mean_rs_max = (al_ratio * al_rs_max) + (ga_ratio * ga_rs_max) + (in_ratio * in_rs_max) + (o_ratio * o_rs_max)
sd_rs_max   =  sqrt( al_ratio * (al_rs_max - mean_rs_max)^2 +  ga_ratio * (ga_rs_max - mean_rs_max)^2 + 
             in_ratio * (in_rs_max - mean_rs_max)^2 + o_ratio * (o_rs_max - mean_rs_max)^2 )

combined$mean_rs_max = mean_rs_max
combined$sd_rs_max   = sd_rs_max

###
### atomic radius from www.ptable.com
###

al_aradius = 143
ga_aradius = 135
in_aradius = 167
o_aradius  = 66  ### used covalent radius

mean_aradius = (al_ratio * al_aradius) + (ga_ratio * ga_aradius) + (in_ratio * in_aradius) + (o_ratio * o_aradius)
sd_aradius   =  sqrt( al_ratio * (al_aradius - mean_aradius)^2 +  ga_ratio * (ga_aradius - mean_aradius)^2 + 
             in_ratio * (in_aradius - mean_aradius)^2 + o_ratio * (o_aradius - mean_aradius)^2 )

combined$mean_aradius = mean_aradius
combined$sd_aradius   = sd_aradius

###
### thermal conductivity from www.ptable.com
###

al_thermal = 237
ga_thermal = 40.6
in_thermal = 81.8
o_thermal  = 0.02658

mean_thermal = (al_ratio * al_thermal) + (ga_ratio * ga_thermal) + (in_ratio * in_thermal) + (o_ratio * o_thermal)
sd_thermal   =  sqrt( al_ratio * (al_thermal - mean_thermal)^2 +  ga_ratio * (ga_thermal - mean_thermal)^2 + 
             in_ratio * (in_thermal - mean_thermal)^2 + o_ratio * (o_thermal - mean_thermal)^2 )

combined$mean_thermal = mean_thermal
combined$sd_thermal   = sd_thermal

###
### Treat spacegroup as an ordered factor:
###

combined$spacegroup = factor(as.factor(combined$spacegroup), 
                      levels = c("12","33","167","194","206", "227"),
                      ordered = TRUE)

###
### The following are just ratios of element contents, which are in columns 3 to 5.
### 1 is added to prevent division by zero.
###

combined$al_over_ga = (combined$percent_atom_al+1)/(combined$percent_atom_ga+1)
combined$al_over_in = (combined$percent_atom_al+1)/(combined$percent_atom_in+1)
combined$ga_over_in = (combined$percent_atom_ga+1)/(combined$percent_atom_in+1)

###
### This feature is meant to capture variability in the content:
###

combined$sd_content = apply(combined[,3:5], 1, sd)

###
### I hope this captures more interactions:
###

combined$al_times_ga = (combined$percent_atom_al+1)*(combined$percent_atom_ga+1)
combined$al_times_in = (combined$percent_atom_al+1)*(combined$percent_atom_in+1)
combined$ga_times_in = (combined$percent_atom_ga+1)*(combined$percent_atom_in+1)

###
### These extract features from lattice vectors:
### lattice_vector_1_ang...lattice_vector_3_ang are columns 6:8
###

combined$size_of_lattice          = apply( combined[,6:8], 1, gmean)
combined$sd_of_lattice_size       = apply( combined[,6:8], 1, sd)

###
### These extract features from lattice angles:
### lattice_angle_alpha_degree, lattice_angle_beta_degree, and
###  lattice_angle_gamma_degree
###
### These are columns 9:11
###

combined$mean_angle = apply(combined[,9:11], 1, gmean)
combined$sd_angle   = apply(combined[,9:11], 1, sd)


###
### Here is a function that takes the geometry file and extracts features:
###

extract_features = function(tmp)
{
  start_pos = min(which(tmp == "lattice_vector"))
  tmp_new   = tmp[start_pos:length(tmp)]
  tmp_new   = tmp_new[-which(tmp_new == "atom")]
  len       = length(tmp_new)/4
  val       = rep(NA,len)
  tmp_tab   = data.frame(v1 = val, v2 = val, v3 = val, v4 = val)
  stat_ind  = 1
  end_ind   = 4
  for (i in 1:len)
  {
    tmp_tab[i,] = tmp_new[stat_ind:end_ind]
    stat_ind    = stat_ind + 4
    end_ind     = end_ind + 4
  }

  tmp_val_1 = tmp_tab[1:3,4]
  tmp_val_2 = tmp_tab[1:3,1]
  tmp_tab[1:3,4] = tmp_val_2
  tmp_tab[1:3,1] = tmp_val_1

  ###
  ### Convert to numeric
  ###

  for (i in 1:3) tmp_tab[,i] = as.numeric(tmp_tab[,i])


  ###
  ### It seems that I can ignore the lattice vectors
  ###

  ###
  ### Get average position of Al atoms:
  ###

  al_loc = which(tmp_tab$v4 == "Al")
  if ( length(al_loc) > 0 ) {

    al_table = tmp_tab[al_loc,1:3]

    al_avg_loc   = apply(al_table, 2, mean)
    al_loc_var   = apply(al_table, 2, pop_sd)
    al_f_norm    = sqrt(mean( cor(al_table)^2 ))
    if (is.na(al_f_norm)) { al_f_norm = 0 }

  } else {

    al_avg_loc  = c(0,0,0)
    al_loc_var  = c(0,0,0)
    al_f_norm   = 0
  }

  ###
  ### Get average position of Ga atoms:
  ###

  ga_loc = which(tmp_tab$v4 == "Ga")
  if ( length(ga_loc) > 0 ) {

    ga_table = tmp_tab[ga_loc,1:3]

    ga_avg_loc   = apply(ga_table, 2, mean)
    ga_loc_var   = apply(ga_table, 2, pop_sd)
    ga_f_norm    = sqrt(mean( cor(ga_table)^2 ))
    if (is.na(ga_f_norm)) { ga_f_norm = 0 }

  } else {

    ga_avg_loc  = c(0,0,0)
    ga_loc_var  = c(0,0,0)
    ga_f_norm   = 0

  }

  ###
  ### Get average position of In atoms:
  ###

  in_loc = which(tmp_tab$v4 == "In")
  if ( length(in_loc) > 0 ) {

    in_table = tmp_tab[in_loc,1:3]

    in_avg_loc = apply(in_table, 2, mean)
    in_loc_var = apply(in_table, 2, pop_sd)
    in_f_norm    = sqrt(mean( cor(in_table)^2 ))
    if (is.na(in_f_norm)) { in_f_norm = 0 }

  } else {

    in_avg_loc  = c(0,0,0)
    in_loc_var  = c(0,0,0)
    in_f_norm   = 0
  }

  ###
  ### Get average position of O2 atoms:
  ###

  o_loc = which(tmp_tab$v4 == "O")
  if ( length(o_loc) > 0 ) {

    o_table = tmp_tab[o_loc,1:3]

    o_avg_loc = apply(o_table, 2, mean)
    o_loc_var = apply(o_table, 2, pop_sd)
    o_f_norm    = sqrt(mean( cor(o_table)^2 ))
    if (is.na(o_f_norm)) { o_f_norm = 0 }

  } else {

    o_avg_loc  = c(0,0,0)
    o_loc_var  = c(0,0,0)
    o_f_norm   = 0

  }

  ###
  ### Find distances
  ###

  o_to_al    = o_avg_loc  -  al_avg_loc
  d_o_to_al  = sqrt(sum(o_to_al^2))

  o_to_ga    = o_avg_loc  -  ga_avg_loc
  d_o_to_ga  = sqrt(sum(o_to_ga^2))

  o_to_in    = o_avg_loc  -  in_avg_loc
  d_o_to_in  = sqrt(sum(o_to_in^2))

  al_to_ga   = al_avg_loc -  ga_avg_loc
  d_al_to_ga = sqrt(sum(al_to_ga^2))

  al_to_in   = al_avg_loc -  in_avg_loc
  d_al_to_in = sqrt(sum(al_to_in^2))

  ga_to_in   = ga_avg_loc -  in_avg_loc
  d_ga_to_in = sqrt(sum(ga_to_in^2))

  ###
  ### Here are the features
  ###


  out = c(al_avg_loc, al_loc_var, al_f_norm,  ga_avg_loc, ga_loc_var, ga_f_norm,
           in_avg_loc, in_loc_var, in_f_norm,  o_avg_loc, o_loc_var, 
           o_f_norm, o_to_al ,  d_o_to_al ,  o_to_ga  ,  d_o_to_ga,  o_to_in  ,  d_o_to_in ,  
            al_to_ga ,  d_al_to_ga ,  al_to_in ,  d_al_to_in,  ga_to_in ,  d_ga_to_in )
  names(out) = NULL
 
  out

}


###
### This is a loop for reading in the geometry files for the train.
### Make note of the locations. 52 comes from the 52 features extracted.
###

u_train = matrix(NA, 2400, 52)

system.time({
for (i in 1:2400)
{ 
  lead_path    = "C:/prohith/UT Austin/Spring 2018/ISL/ISL Assignment 4/train/train/" 
  the_path     = paste(lead_path, i, sep = "")
  setwd(the_path)
  tmp          = scan("geometry.xyz", what = "raw", skip = 1)
  u_train[i,]  = extract_features(tmp)
}
})

u_train = as.data.frame(u_train)


###
### Do the same for test:
###

u_test = matrix(NA, 600, 52)

system.time({
for (i in 1:600)
{ 
  lead_path    = "C:/prohith/UT Austin/Spring 2018/ISL/ISL Assignment 4/test/test/" 
  the_path     = paste(lead_path, i, sep = "")
  setwd(the_path)
  tmp          = scan("geometry.xyz", what = "raw", skip = 1)
  u_test[i,]   = extract_features(tmp)
}
})

u_test = as.data.frame(u_test)


###
### Combine the new features
###

u_combined = rbind(u_train, u_test)

###
### Now add it to the combined:
###

big_combined = cbind(combined, u_combined)

###
### First delete the originals and separate:
###

rm(train)
rm(test)
rm(i)
rm(tmp)

train = big_combined[rows_for_train, ] 
test  = big_combined[rows_for_test, ] 

###
### Add back the response:
###

train$formation_energy_ev_natom  = formation_energy_ev_natom
train$bandgap_energy_ev          = bandgap_energy_ev


###
### Reset the working directory & Save the work so far:
###

setwd("C:/prohith/UT Austin/Spring 2018/ISL/ISL Assignment 4")

###
### Try a simple multiple regression model:
###
### train[,-98] excludes bandgap_energy_ev
### train[,-97] excludes formation_energy_ev_natom 
###
### ranger is faster than randomForest.
###
#install.packages("ranger")
library(ranger)

set.seed(1)

model_for_formation = ranger(formation_energy_ev_natom ~ ., data = train[,-98])
model_for_bandgap   = ranger(bandgap_energy_ev ~ ., data = train[,-97])

formation_pred = exp(predict(model_for_formation, test)$predictions) - 1
bandgap_pred   = exp(predict(model_for_bandgap, test)$predictions) - 1

###
### Now submit:
###


submission = data.frame(id = 1:600,
                        formation_energy_ev_natom = formation_pred,
                        bandgap_energy_ev = bandgap_pred)

write.csv(submission, "Multiple linear regression submission 1.csv", row.names = F)

###
### This got me a decent score.
###

## Random Forest
model_for_formation2 = ranger(formation_energy_ev_natom ~ ., data = train[,-98], num.trees = 1000)
model_for_bandgap2   = ranger(bandgap_energy_ev ~ ., data = train[,-97], num.trees = 1000)

formation_pred2 = exp(predict(model_for_formation2, test)$predictions) - 1
bandgap_pred2   = exp(predict(model_for_bandgap2, test)$predictions) - 1

submission2 = data.frame(id = 1:600,
                        formation_energy_ev_natom = formation_pred2,
                        bandgap_energy_ev = bandgap_pred2)

write.csv(submission2, "Random Forest submission 4.csv", row.names = F)

## Gradient Boosting
## install.packages("gbm")

library(DAAG)
library(gbm)

set.seed(17)
gbm_model_formation = gbm(formation_energy_ev_natom ~ ., data = train[,-98], n.trees = 10000, cv.folds = 3, 
                          interaction.depth = 2)
gbm_model_bandgap = gbm(bandgap_energy_ev ~ ., data = train[,-97], n.trees = 10000, cv.folds = 3, interaction.depth = 2)

gbm_formation_pred = exp(predict(gbm_model_formation, test, n.trees = gbm_model_formation$n.tree)) - 1
gbm_bandgap_pred  = exp(predict(gbm_model_bandgap, test, n.trees  = gbm_model_bandgap$n.tree)) - 1

gbm_submission = data.frame(id = 1:600,
                         formation_energy_ev_natom = gbm_formation_pred,
                         bandgap_energy_ev = gbm_bandgap_pred)

write.csv(gbm_submission, "GBM submission 4.csv", row.names = F)

## Lasso and Ridge regression
##install.packages("glmnet")
library(glmnet)

x = subset(train[,-98], select = -c(formation_energy_ev_natom))
x = model.matrix(~(.)^2 , x)

y = subset(train[,-98], select = c(formation_energy_ev_natom))
y = as.matrix(y)

ridge_model = cv.glmnet(x = x, y = y, alpha = 0)  #alpha = 0 means ridge regression penalty

ridge_coef  = coef(ridge_model, s = ridge_model$lambda.min)

lasso_model = cv.glmnet(x = x, y = y, alpha = 1)  #alpha = 1 means lasso regression penalty

lasso_coef  = coef(lasso_model, s = lasso_model$lambda.min)

elasticnet_model = cv.glmnet(x = x, y = y, alpha = 0.5)

x = model.matrix(~(.)^2 , test)
y_hat = predict(lasso_model,  newx = x, s = lasso_model$lambda.min)
summary(y_hat)

##

x1 = subset(train[,-97], select = -c(bandgap_energy_ev))
x1 = model.matrix(~(.)^2 , x1)

y1 = subset(train[,-97], select = c(bandgap_energy_ev))
y1 = as.matrix(y1)

ridge_model1 = cv.glmnet(x = x1, y = y1, alpha = 0)  #alpha = 0 means ridge regression penalty

ridge_coef1  = coef(ridge_model1, s = ridge_model1$lambda.min)

lasso_model1 = cv.glmnet(x = x1, y = y1, alpha = 1)  #alpha = 1 means lasso regression penalty

lasso_coef1  = coef(lasso_model1, s = lasso_model1$lambda.min)

## elasticnet_model = cv.glmnet(x = x, y = y, alpha = 0.5)

x1 = model.matrix(~(.)^2 , test)
y_hat1 = predict(lasso_model1,  newx = x1, s = lasso_model1$lambda.min)
summary(y_hat1)

lasso_submission = data.frame(id = 1:600,
                            formation_energy_ev_natom = exp(y_hat)-1,
                            bandgap_energy_ev = exp(y_hat1)-1)

write.csv(lasso_submission, "lasso 2.csv", row.names = F)

## Tree model
install.packages("rpart")
library(rpart)
tree.model_formation = rpart(formation_energy_ev_natom ~ ., data = train[,-98])
tree.model_bandgap = rpart(bandgap_energy_ev ~ ., data = train[,-97])

tree_predict_formation = exp(predict(tree.model_formation, newdata = test ))-1
tree_predict_bandgap = exp(predict(tree.model_bandgap, newdata = test ))-1
tree_submission = data.frame(id = 1:600,
                            formation_energy_ev_natom = tree_predict_formation,
                            bandgap_energy_ev = tree_predict_bandgap)

write.csv(tree_submission, "Tree submission 1.csv", row.names = F)


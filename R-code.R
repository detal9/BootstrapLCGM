set.seed(347921) 
for (i in 1:N) {  
##############################################################################
### STEP 1: Simulation of groups: Scenario 1
##############################################################################  
  p1= 0.8 * 1.5* tso2 + rnorm(200, 0, 2) +  rnorm(200, 0, 5);
  p2 =  0.8 * ((100-tso+10*tso^2)^1/tso) + rnorm(150, 0, 2) +  rnorm(150, 0, 10); 
  p3 = 0.8* (180 * tso -10*tso^2-1.5*tso^3) + rnorm(150, 0, 2) +  rnorm(150,0, 20);
  Y= c(p1,p2, p3);
  Y <- ifelse(Y<0, 0, Y); # All negative values are replaced by 0
  dat_value = data.frame( t=1:5, Y, id= rep(1:100,each=5));
##############################################################################
### STEP 2: Selection of the number of groups
##############################################################################  
  lin = matrix(NA, nrow = 5, ncol = 2);
  models = list(); 
  output <-  matrix();
  x = matrix(NA, nrow = 1, ncol = 5);
  k=1;
  for (j in 1:5){
    mj <-  stepFlexmix(Y ~ t + I(t^2)+ I(t^3)| id, data=dat_value, k = j, nrep = 5, control=list(iter.max = 1000, minprior = 0))
    if(j == length(table(clusters(mj)))){
      lin[k,] = c(j, BIC(mj));
      x[,k] = BIC(mj);
      models[[k]] = mj;
    }
    k = k + 1;
  }
  lin2 <- as.data.frame(lin);
###############################################################################
### STEP 3: Model adequacy criteria
############################################################################### 
  n <- 500;
  K <- which.min(x);# number of groups for each replication
  m3 = models[[K]]; 
  fit <- cbind(m3@group, m3@cluster, posterior(m3));
  fit<-unique(fit); # one line per person
  #Selection of the maximum posterior probability considering the number of groups identified by replication
  if (  K == 2) { 
    group <-  cbind(fit[,1], fit[,2], pmax(fit[,3], fit[,4]));  
  } else
    if (  K == 3) {
      group <-  cbind(fit[,1],fit[,2], pmax(fit[,3],fit[,4],fit[,5]));
    } else
      if (  K == 4) { 
        group <-  cbind(fit[,1],fit[,2], pmax(fit[,3], fit[,4], fit[,5],  fit[,6]));
      } else
        if (  K == 5) {
          group <-  cbind(fit[,1],fit[,2], pmax(fit[,3], fit[,4],fit[,5], fit[,6], fit[,7]));
        } 
  #Average posterior probability
  appa <- cbind(tapply(group[,3],group[,2],mean));
  appa_min <- min(appa);
  #Entropy
  EIC <- EIC(m3);
  #Odds of correct classification
  pi <- prior(m3); #Estimated proportion of group membership
  occ <- (appa/(1 - appa))/(pi/(1 - pi));
  occ_min <- min(occ);
  output <- c(i, K, appa_min, EIC, occ_min) ;
  # Saving the output
  capture.output(output, file = "C:/. /test_rep_100.txt", append = TRUE)
}  
##############################################################################
### STEP 4: Bootstrap approach
##############################################################################  
#Transformation of long format data to wide format to resample subjects
library(reshape2)
data_large<- reshape2::dcast(dat_value, id~ t, value.var= "Y")
#Bootstrap
boot_fun <- function(data, index){
  d <- data[index,]
  d$id = 1:100; 
  #Format wide to long
  d_long <- reshape2::melt(d, id.var="id")
  names(d_long)[names(d_long)=="variable"] <- "t" 
  names(d_long)[names(d_long)=="value"] <- "Y" 
  d_long$t <- as.numeric(d_long$t)
  lin = matrix(NA, nrow = 5, ncol = 2);
  models = list(); 
  output <-  matrix();
  x = matrix(NA, nrow = 1, ncol = 5);
  k=1;
  for (j in 1:5){
    mj <-  stepFlexmix(Y ~ t + I(t^2)+ I(t^3)| id, data=dat_value, k = j, nrep = 5, control=list(iter.max = 1000, minprior = 0))
    if(j == length(table(clusters(mj)))){
      lin[k,] = c(j, BIC(mj));
      models[[k]] = mj;
      x[,k] = BIC(mj);
    }
    k = k + 1;
  }
  K <- which.min(x);
  output <- K 
  return(output)
}
res <- boot(data=data_large, statistic=boot_fun, R=100) 
results <- as.data.frame(res[["t"]], header=F)
# Saving the output
capture.output(results, file = "C:/./test_boot_100.txt", append = TRUE)

#Outlier-Dependent Comparison
#Ben Rubin
#April 2026
#This function takes a list of sample means, compares them, and returns
#which pairs are significantly difference. It uses a combination of Tukey's
#HSD and a modified Bonferroni correction depending on outlier status of
#each sample mean.

odc <- function(means, J, MSE, alpha = 0.05) {
  # means: vector of sample means
  # J: number of observations per group
  # MSE: mean square error from ANOVA
  # alpha: significance level
  
  # Returns: a data frame with columns:
  # - group1, group2: the pair being compared
  # - diff: the difference in means
  # - significant: TRUE/FALSE
  
  #vectors to store results 
  X1 <- c() #this is group 1
  X2 <- c() #this is group 2
  diff <- c()
  significant <- c()
  
#SORT SAMPLE MEANS BETWEEN OUTLIERS, EXTREME OUTLIERS, NONOUTLIERS
  
  iqr <- IQR(means) #calculate interquartile range
  med <- median(means) #calculate median
  
  #create vectors to sort sample means between outliers, extreme outliers, and non-outliers
  extreme_outliers <- c()
  outliers <- c()
  non_outliers <- c()
  
  I <- length(means) #This is the number of means
  
  for (x in means){ #check each sample mean
    if (x < med - 3 * iqr || x > med + 3 * iqr){ #if it is further than 3 times the IQR from the median
      extreme_outliers <- c(extreme_outliers, x) #add it to the vector of extreme outliers
    }
    else if (x < med - 1.5 * iqr || x > med + 1.5 * iqr){ #if it is more than 1.5 times the IQR fromt the median
      outliers <- c(outliers, x) #add it to the vector of outliers
    }
    else{
      non_outliers <- c(non_outliers, x) #otherwise add it to the vector of non-outliers
    }
  }
  
#COMPARE EXTREME OUTLIERS USING T-TEST WITH ALPHA/n^2
  
  n <- choose(I,2) #n is the total number of comparisons being made
  
  if (length(extreme_outliers) > 0){#make sure the vector is not empty
    #compare extreme outliers to other extreme outliers
    if (length(extreme_outliers) > 1){#make sure there is something else to compare to
      for (a in 1:(length(extreme_outliers)-1)){
        for (b in (a+1):length(extreme_outliers)){
          A <- extreme_outliers[a] #A and B are the
          B <- extreme_outliers[b] #sample means being compared
          X1 <- c(X1, A) #record the pair
          X2 <- c(X2, B)
          diff <- c(diff, abs(A-B)) #and difference
          t <- (A - B) / sqrt (2*MSE / J) #calculate test statistic
          p <- 2*(1-pt(t, I*(J-1))) #calculate p value
          #determine if it is significant
          if (p < alpha / n){
            significant <- c(significant, TRUE)
          }
          else{
            significant <- c(significant, FALSE)
          }
        }
      }
    }
    #compare extreme outliers to outliers
    if(length(outliers) > 0){#make sure there is something to compare to
      for (a in 1:length(extreme_outliers)){
        for (b in 1:length(outliers)){
          A <- extreme_outliers[a] #A and B are the sample
          B <- outliers[b] #means being compared
          X1 <- c(X1, A) #record the pair
          X2 <- c(X2, B)
          diff <- c(diff, abs(A-B)) #and difference
          t <- (A - B) / sqrt (2*MSE / J) #calculate test statistic
          p <- 2*(1-pt(t, I*(J-1))) #calculate p value
          #determine if it is significant
          if (p < alpha / n){
            significant <- c(significant, TRUE)
          }
          else{
            significant <- c(significant, FALSE)
          }
        }
      }
    }
    #compare extreme outliers to non_outliers
    for (a in 1:length(extreme_outliers)){
      for (b in 1:length(non_outliers)){
        A <- extreme_outliers[a] #A and B are the sample
        B <- non_outliers[b] #means being compared
        X1 <- c(X1, A) #record the pair
        X2 <- c(X2, B)
        diff <- c(diff, abs(A-B)) #and difference
        t <- (A - B) / sqrt (2*MSE / J) #calculate test statistic
        p <- 2*(1-pt(t, I*(J-1))) #calculate p value
        #determine if it is significant
        if (p < alpha / n){
          significant <- c(significant, TRUE)
        }
        else{
          significant <- c(significant, FALSE)
        }
      }
    }
  }
  
#COMPARE OUTLIERS USING T-TEST WITH ALPHA / n  
  
  if (length(outliers) > 0){#make sure the vector is not empty
    #compare with other outliers
    if(length(outliers) > 1){ #make sure there is something to compare to
      for (a in 1:(length(outliers)-1)){
        for (b in (a+1):length(outliers)){
          A <- outliers[a] #A and B are the sample
          B <- outliers[b] #means being compared
          X1 <- c(X1, A) #record the pair
          X2 <- c(X2, B)
          diff <- c(diff, abs(A-B)) #and difference
          t <- (A - B) / sqrt (2*MSE / J) #calculate test statistic
          p <- 2*(1-pt(t, I*(J-1))) #calculate p value
          #determine if it is significant
          if (p < (alpha / sqrt(n))){
            significant <- c(significant, TRUE)
          }
          else{
            significant <- c(significant, FALSE)
          }
        }
      }
    }
    #compare with non-outliers
    for (a in 1:length(outliers)){
      for (b in 1:length(non_outliers)){
        A <- outliers[a] #A and B are the sample
        B <- non_outliers[b]#means being compared
        X1 <- c(X1, A) #record the pair
        X2 <- c(X2, B)
        diff <- c(diff, abs(A-B)) #and difference
        t <- (A - B) / sqrt (2*MSE / J) #calculate test statistic
        p <- 2*(1-pt(t, 2*(J-1))) #calculate p value
        #determine if it is significant
        if (p < (alpha / sqrt(n))){
          significant <- c(significant, TRUE)
        }
        else{
          significant <- c(significant, FALSE)
        }
      }
    }
  }
  
#COMPARE NON OUTLIERS USING TUKEY
  
  if (length(non_outliers) > 0){#make sure the vector is non-empty
    Q <- qtukey(1-alpha, I, I * (J-1)) #calculate Q
    HSD <- Q * sqrt(MSE / J) #calculate the HSD
    
    for (a in 1:(length(non_outliers)-1)){
      for (b in (a+1):length(non_outliers)){
        A <- non_outliers[a] #A and B are the sample
        B <- non_outliers[b] #means being compared
        X1 <- c(X1, A) #record the pair
        X2 <- c(X2, B)        
        diff <- c(diff, abs(non_outliers[a]-non_outliers[b])) #and difference
        #determine if it is significant
        if (abs(non_outliers[a]-non_outliers[b]) > HSD){
          significant <- c(significant, TRUE)
        }
        else{
          significant <- c(significant, FALSE)
        }
      }
    }
  }

#RETURN RESULTS
  
  #create the data frame to return
  ret <- data.frame(X1, X2,diff,significant)
  return(ret) #return the data frame
}

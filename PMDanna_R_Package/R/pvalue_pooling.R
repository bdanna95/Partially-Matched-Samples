#library(roxygen2); # Read in the roxygen2 R package uncomment these two lines to modify the help files for the package
#roxygenise();      # Builds the help files
#' Liptak's weighted Z-test
#' 
#' Calculates a combined p-value by the weighted Z-test
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of the weighted Z-statistic \cr
#'    \tab \cr
#'    \code{p.value} \tab the combined p.value for the test \cr
#'    \tab \cr
#'    \code{alternative} \tab the alternative hypothesis for the test \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat") 
#' liptak.z.test(colondat$`araB 1bn`,colondat$`araB 1bT`,alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
liptak.z.test <- function(x,y,alternative = "two.sided"){
  
  if (is.null(x)) {
    stop("'x' is missing for paired test")
  }
  if (is.null(y)) {
    stop("'y' is missing for paired test")
  }
  if (length(x) < 2) {
    stop("not enough 'x' observations")
  }
  if (length(y) < 2) {
    stop("not enough 'y' observations")
  }
  
  sbset <- cbind(x,y)
  s1 <- sbset[!is.na(sbset[,1]) & !is.na(sbset[,2]),] #matched sample
  s2 <- sbset[is.na(sbset[,1]) & !is.na(sbset[,2]),2] #independent group only n2
  s3 <- sbset[!is.na(sbset[,1]) & is.na(sbset[,2]),1] #independent group only n3
  
  n1 <- nrow(s1) #matched sample size
  n2 <- length(s2) #independent group n2 sample size
  n3 <- length(s3) #independent group n3 sample size
  
  if ((n2==0)&(n3==0)){
    stop('all samples are matched pairs, use t.test with paired = True')
  } else if (n1==0){
    stop('all samples are unmatched, use t.test with paired = False')
  }
  if (n1 < 2) {
    stop("not enough paired observations")
  }
  if (n2 < 2 || n3 < 2) {
    stop("not enough unpaired observations")
  }
  if(alternative == "two.sided"){
    alt <- "greater"
  } else{
    alt <- alternative
  }
  
  w1 = sqrt(2*n1) #weight of matched sample
  w2 = sqrt(n2+n3) #weight of independent groups
  
  p1 <- t.test(s1[,1],s1[,2],paired = T,alternative = alt)$p.value
  z1 <- qnorm(1-p1)
  p2 <- t.test(s2,s3,paired = F,alternative = alt)$p.value
  z2 <- qnorm(1-p2)
  
  pc.one.sided <- 1 - pnorm(((w1*z1)+(w2*z2))/sqrt(w1^2+w2^2)) #calculate one-sided p-value
  
  if(alternative != "two.sided"){
    pc <- pc.one.sided
  } else{
    pc <- if (pc.one.sided < 0.5){ #calculate two-sided p-value
      2 *pc.one.sided
    } else {
      (2 *(1- pc.one.sided))
    }
  }
  
  method <- "Liptak's weighted Z-test"
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  tststat <- ((w1*z1)+(w2*z2))/sqrt(w1^2+w2^2)
  names(tststat) <- "weighted Z"
  
  returnlst <- list(statistic = tststat, p.value = pc, alternative = alternative, method = method, data.name = dname)
  class(returnlst) <- "htest"
  returnlst
}

#' Kim et al. modified t-statistic
#' 
#' Calculates modified t-statistic and corresponding p-value from the null distribution of t3 which is approximated with a standard Gaussian distribution.
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of the modified t-statistic \cr
#'    \tab \cr
#'    \code{p.value} \tab the p.value for the statistic \cr
#'    \tab \cr
#'    \code{alternative} \tab the alternative hypothesis for the test \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat") 
#' kim.t.test(colondat$`araB 1bn`,colondat$`araB 1bT`,alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
kim.t.test <- function(x,y,alternative = "two.sided"){
  
  if (is.null(x)) {
    stop("'x' is missing for paired test")
  }
  if (is.null(y)) {
    stop("'y' is missing for paired test")
  }
  if (length(x) < 2) {
    stop("not enough 'x' observations")
  }
  if (length(y) < 2) {
    stop("not enough 'y' observations")
  }
  
  sbset <- cbind(x,y)
  s1 <- sbset[!is.na(sbset[,1]) & !is.na(sbset[,2]),] #matched sample
  s2 <- sbset[is.na(sbset[,1]) & !is.na(sbset[,2]),2] #independent group only n2
  s3 <- sbset[!is.na(sbset[,1]) & is.na(sbset[,2]),1] #independent group only n3
  
  n1 <- nrow(s1) #matched sample size
  n2 <- length(s2) #independent group n2 sample size
  n3 <- length(s3) #independent group n3 sample size
  nh <- 2/((1/n2)+(1/n3)) #harmonic mean of n2 and n3
  
  if ((n2 < 2) || (n3 < 2)){
    stop('not enough unmatched observations')
  } else if (n1 < 2){
    stop('not enough matched observations')
  }
  
  dbar <- mean(s1[,1]-s1[,2]) #mean difference of the n1 paired samples
  tbar <- mean(s2) #mean for the n2 unmatched samples
  nbar <- mean(s3) #mean for the n3 unmatched samples
  
  sdd <- sd(s1[,1]-s1[,2]) #standard deviation of the difference of the n1 paired samples
  sdt <- sd(s2) #standard deviation of the n2 unmatched samples
  sdn <- sd(s3) #standard deviation of the n3 unmatched samples
  
  t3 <- (n1*dbar+nh*(tbar-nbar))/sqrt(n1*sdd^2+nh^2*(sdn^2/n3+sdt^2/n2)) #modified t-statistic
  
  if(alternative != "two.sided"){
    p.value <- pnorm(abs(t3),lower.tail = F) #null distribution of t3 is approximated with a standard normal distribution
  } else{
    p.value <- 2*pnorm(abs(t3),lower.tail = F)
  }
  
  method <- "Kim et al. modified t-statistic"
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  names(t3) <- "t3"
  
  returnlst <- list(statistic = t3, p.value = p.value, alternative = alternative, method = method,data.name = dname)
  class(returnlst) <- "htest"
  returnlst
}

#' Looney and Jones corrected Z-test
#' 
#' Calculates corrected Z-test and corresponding p-value based on a modified variance estimation of the standard Z-test by accounting for the correlation among the n1 matched pairs.
#' The corrected Z-test reduces to paired sample or two-sample Z-test when n2=n3=0 or n1=0, respectively.
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of the corrected Z-test \cr
#'    \tab \cr
#'    \code{p.value} \tab the p.value for the statistic \cr
#'    \tab \cr
#'    \code{alternative} \tab the alternative hypothesis for the test \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat") 
#' looneyjones.z.test(colondat$`araB 1bn`,colondat$`araB 1bT`,alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
looneyjones.z.test <- function(x,y,alternative = "two.sided"){
  
  if (is.null(x)) {
    stop("'x' is missing for paired test")
  }
  if (is.null(y)) {
    stop("'y' is missing for paired test")
  }
  if (length(x) < 2) {
    stop("not enough 'x' observations")
  }
  if (length(y) < 2) {
    stop("not enough 'y' observations")
  }
  
  sbset <- cbind(x,y)
  s1 <- sbset[!is.na(sbset[,1]) & !is.na(sbset[,2]),] #matched sample
  s2 <- sbset[is.na(sbset[,1]) & !is.na(sbset[,2]),2] #independent group only n2
  s3 <- sbset[!is.na(sbset[,1]) & is.na(sbset[,2]),1] #independent group only n3
  
  n1 <- nrow(s1) #matched sample size
  n2 <- length(s2) #independent group n2 sample size
  n3 <- length(s3) #independent group n3 sample size
  
  if ((n2==0)&(n3==0)){
    warning('Zcorr is reduced to a paired sample Z-test since there are no unmatched pairs')
  } else if (n1==0){
    warning('Zcorr is reduced to a two-sample Z-test since there are no matched pairs')
  }

  tbar.star <- mean(c(s1[,2],s2)) #mean for the n1+n2 matched and unmatched pairs combined
  nbar.star <- mean(c(s1[,1],s3)) #mean for the n1+n3 matched and unmatched pairs combined
  
  sdt.star <- sd(c(s1[,2],s2)) #standard deviation for the n1+n2 matched and unmatched pairs combined
  sdn.star <- sd(c(s1[,1],s3)) #standard deviation for the n1+n2 matched and unmatched pairs combined
  stn1 <- if (n1!=0){cov(s1[,1],s1[,2])} else {0} #sample covariance of the n1 complete pairs
  
  zcorr <- (tbar.star-nbar.star)/sqrt((((sdt.star^2)/(n1+n2))+((sdn.star^2)/(n1+n3)))-(2*n1*(stn1/((n1+n2)*(n1+n3)))))
  
  if(alternative != "two.sided"){
    p.value <- pnorm(abs(zcorr),lower.tail = F) #null distribution of t3 is approximated with a standard normal distribution
  } else{
    p.value <- 2*pnorm(abs(zcorr),lower.tail = F)
  }
  
  method <- "Looney and Jones corrected Z-test"
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  names(zcorr) <- "corrected Z"
  
  returnlst <- list(statistic = zcorr, p.value = p.value, alternative = alternative, method = method, data.name = dname)
  class(returnlst) <- "htest"
  returnlst
}

#' Lin and Stiver's MLE based test statistic under heteroscedasticity
#' 
#' Calculates MLE based test statistic and corresponding p-value based on a modified maximum likelihood estimator and simple mean difference for testing
#' the equality of two correlated means with incomplete data under the bivariate Gaussian assumption.
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of the MLE based statistic \cr
#'    \tab \cr
#'    \code{p.value} \tab the p.value for the statistic \cr
#'    \tab \cr
#'    \code{alternative} \tab the alternative hypothesis for the test \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat")  
#' linsilver.mle.test(colondat$`araB 1bn`,colondat$`araB 1bT`,alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
linsilver.mle.test <- function(x,y,alternative = "two.sided"){
  
  if (is.null(x)) {
    stop("'x' is missing for paired test")
  }
  if (is.null(y)) {
    stop("'y' is missing for paired test")
  }
  if (length(x) < 2) {
    stop("not enough 'x' observations")
  }
  if (length(y) < 2) {
    stop("not enough 'y' observations")
  }
  
  sbset <- cbind(x,y)
  s1 <- sbset[!is.na(sbset[,1]) & !is.na(sbset[,2]),] #matched sample
  s2 <- sbset[is.na(sbset[,1]) & !is.na(sbset[,2]),2] #independent group only n2
  s3 <- sbset[!is.na(sbset[,1]) & is.na(sbset[,2]),1] #independent group only n3
  
  n1 <- nrow(s1) #matched sample size
  n2 <- length(s2) #independent group n2 sample size
  n3 <- length(s3) #independent group n3 sample size

  if ((n2 < 2) || (n3 < 2)){
    stop('not enough unmatched observations')
  } else if (n1 < 2){
    stop('not enough matched observations')
  }
  
  tbar <- mean(s2) #mean for the n2 unmatched samples
  nbar <- mean(s3) #mean for the n3 unmatched samples
  tbar.1 <- mean(s1[,2]) #mean for the n1 paired samples
  nbar.1 <- mean(s1[,1]) #mean for the n2 paired samples
  
  sdt.1 <- sd(s1[,2]) #standard deviation for the n1 paired samples
  sdn.1 <- sd(s1[,1]) #standard deviation for the n2 paired samples
  stn1 <- cov(s1[,1],s1[,2]) #sample covariance of the n1 complete pairs
  
  r <- stn1/(sdt.1*sdn.1)
  f <- n1*(n1+n3+(n2*(stn1/sdt.1^2)))*(((n1+n2)*(n1+n3)-n2*n3*r^2)^-1)
  g <- n1*(n1+n2+(n3*(stn1/sdn.1^2)))*(((n1+n2)*(n1+n3)-n2*n3*r^2)^-1)
  V1 = ((((f^2)/n1+((1-f)^2)/n2)*(sdt.1^2)*(n1-1))+(((g^2)/n1+((1-g)^2)/n3)*(sdn.1^2)*(n1-1))-(2*f*g*stn1*(n1-1)/n1))/(n1-1)
  
  zls <- ((f*(tbar.1-tbar))-(g*(nbar.1-nbar))+(tbar-nbar))/sqrt(V1) #MLE based test statistic under heteroscedasticity
  
  if(alternative != "two.sided"){
    p.value <- pnorm(abs(zls),lower.tail = F) #null distribution of t3 is approximated with a standard normal distribution
  } else{
    p.value <- 2*pnorm(abs(zls),lower.tail = F)
  }
  
  method <- "Lin and Stiver's MLE based test statistic under heteroscedasticity"
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  names(zls) <- "MLE based statistic"
  
  returnlst <- list(statistic = zls, p.value = p.value, alternative = alternative, method = method, data.name = dname)
  class(returnlst) <- "htest"
  returnlst
}

#' Ekbohm's MLE based test statistic under homoscedasticity
#' 
#' Calculates MLE based test statistic and corresponding p-value based on a modified maximum likelihood estimator and simple mean difference for testing
#' the equality of two correlated means with incomplete data under the bivariate Gaussian assumption. Suggested use when the variances are equal.
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of the MLE based statistic \cr
#'    \tab \cr
#'    \code{p.value} \tab the p.value for the statistic \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{alternative} \tab the alternative hypothesis for the test \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat")  
#' ekbohm.mle.test(colondat$`araB 1bn`,colondat$`araB 1bT`,alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
ekbohm.mle.test <- function(x,y,alternative = "two.sided"){
  
  if (is.null(x)) {
    stop("'x' is missing for paired test")
  }
  if (is.null(y)) {
    stop("'y' is missing for paired test")
  }
  if (length(x) < 2) {
    stop("not enough 'x' observations")
  }
  if (length(y) < 2) {
    stop("not enough 'y' observations")
  }
  
  sbset <- cbind(x,y)
  s1 <- sbset[!is.na(sbset[,1]) & !is.na(sbset[,2]),] #matched sample
  s2 <- sbset[is.na(sbset[,1]) & !is.na(sbset[,2]),2] #independent group only n2
  s3 <- sbset[!is.na(sbset[,1]) & is.na(sbset[,2]),1] #independent group only n3
  
  n1 <- nrow(s1) #matched sample size
  n2 <- length(s2) #independent group n2 sample size
  n3 <- length(s3) #independent group n3 sample size
  
  if ((n2 < 2) || (n3 < 2)){
    stop('not enough unmatched observations')
  } else if (n1 < 2){
    stop('not enough matched observations')
  }
  
  tbar <- mean(s2) #mean for the n2 unmatched samples
  nbar <- mean(s3) #mean for the n3 unmatched samples
  tbar.1 <- mean(s1[,2]) #mean for the n1 paired samples
  nbar.1 <- mean(s1[,1]) #mean for the n2 paired samples
  
  sdt <- sd(s2) #standard deviation of the n2 unmatched samples
  sdn <- sd(s3) #standard deviation of the n3 unmatched samples
  sdt.1 <- sd(s1[,2]) #standard deviation for the n1 paired samples
  sdn.1 <- sd(s1[,1]) #standard deviation for the n2 paired samples
  stn1 <- cov(s1[,1],s1[,2]) #sample covariance of the n1 complete pairs
  
  r <- stn1/(sdt.1*sdn.1)
  f.star <- n1*(n1+n3+(n2*r))*(((n1+n2)*(n1+n3)-n2*n3*r^2)^-1)
  g.star <- n1*(n1+n2+(n3*r))*(((n1+n2)*(n1+n3)-n2*n3*r^2)^-1)
  sigma.hat <-(sdt.1^2*(n1-1)+sdn.1^2*(n1-1)+(1+r^2)*(sdt^2*(n2-1)+sdn^2*(n3-1)))/(2*(n1-1)+(1+r^2)*(n2+n3-2))
  V1.star <- sigma.hat*((2*n1*(1-r)+(n2+n3)*(1-r^2))/((n1+n2)*(n1+n3)-(n2*n3*(r^2))))
  
  ze <- ((f.star*(tbar.1-tbar))-(g.star*(nbar.1-nbar))+(tbar-nbar))/sqrt(V1.star) #MLE based test statistic under homoscedasticity

  if(alternative != "two.sided"){
    p.value <- pnorm(abs(ze),lower.tail = F) #ZE is approximately distributed as t with n1 degrees of freedom under the null hypothesis
  } else{
    p.value <- 2*pnorm(abs(ze),lower.tail = F)
  }
  
  method <- "Ekbohm's MLE based test statistic under homoscedasticity"
  
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  
  names(ze) <- "MLE based statistic"
  
  returnlst <- list(statistic = ze, p.value = p.value, alternative = alternative, method = method, data.name = dname)
  class(returnlst) <- "htest"
  returnlst
}

#' P-value Pooling Strategies for Partially Matched Samples
#' 
#' Calculates the test statistic and corresponding p-value for partially matched samples based on the user's preferred pooling method.
#' Partially matched samples are viewed as data generated from two experimental designs; n1 matched pairs, and independent groups with n2 and n3 observations
#' per group.
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param method a character string specifying which pooling method to perform, must be one of "liptak" (default), "kim", "looneyjones", "linsilver" or "ekbohm"
#' @return a list with class "htest" containing the below:\tabular{ll}{
#'    \code{statistic} \tab the value of chosen method's statistic \cr
#'    \tab \cr
#'    \code{p.value} \tab the p.value for the statistic \cr
#'    \tab \cr
#'    \code{method} \tab indicates what type of pooling approach was performed \cr
#'    \tab \cr
#'    \code{data.name} \tab the names of the data \cr
#' }
#' @examples
#' # Using colondat sample dataset
#' data("colondat") 
#' p.values.pooling(colondat$`araB 1bn`,colondat$`araB 1bT`,method="liptak",alternative="two.sided")
#' @references 
#' Kuan, P., & Huang, B. (2013). A simple and robust method for partially matched samples using the p-values pooling approach. Statistics In Medicine, 32(19), 3247-3259. doi: 10.1002/sim.5758
#' @export
p.values.pooling <- function(x,y,method = "liptak",alternative = "two.sided"){
  
  if (!(method %in% c("liptak","kim","linsilver","looneyjones","ekbohm"))){
    stop("invalid method assignment")
  }
  if (!(alternative %in% c("two.sided","greater","less"))){
    stop("invalid alternative assignment")
  }

  if (method == "liptak"){
    return(liptak.z.test(x=x,y=y,alternative=alternative))
  } else if (method == "kim"){
    kim.t.test(x=x,y=y,alternative=alternative)
  } else if (method == "linsilver"){
    linsilver.mle.test(x=x,y=y,alternative=alternative)
  } else if (method == "looneyjones"){
    looneyjones.z.test(x=x,y=y,alternative=alternative)
  } else if (method == "ekbohm"){
    ekbohm.mle.test(x=x,y=y,alternative=alternative)
  }
}

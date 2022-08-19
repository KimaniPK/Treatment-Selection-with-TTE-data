
###########################################################################################################
###########################################################################################################
#### THIS IS THE R CODE FOR PERFORMING THE COMPUTATIONS IN SECTION 11.1 IN THE SUPPLEMENTARY MATERIAL #####
###########################################################################################################
###########################################################################################################


###################################################################################################################
######################################## Loading the required R libraries #########################################
###################################################################################################################

library(survival)
library(mvtnorm)


###################################################################################################################
############################ Importing example data in the csv file "ExampleFinalData" ############################
###################################################################################################################

y <- read.csv("C:/Documents/ExampleFinalData.csv", header=F) # Replace "C:/Documents" with the folder containing the example data
yy <- cbind(y[,3], y[,1], y[,2], y[,8], y[,7]) ## Creating a copy for computing correlation of stage 1 estimates
colnames(y) <- c("EnterTime", "Stage", "Arm", "CalendarTime", "T.date", "death", "T.date1", "death1")
y <- data.frame(y)


#########################################
## Stage 1 estimates and bounds values ##
#########################################

## bounds ##

Aj <- qnorm(0.2) # This is $\Phi^{-1}(a)$ where $a=0.2$.


############################
## STAGE ONE COMPUTATIONS ##
############################

## Selecting stage 1 data ##

y1 <- y[y$Stage == 1,]
y1$SurvTime <- with(y1, Surv(T.date1, death1))


## Treatment 1 estimates ##

y.onepartition <- y1[y1$Arm != 2,] # Excluding treatment 2 patients
y.onepartition$Group <- 1 * (y.onepartition$Arm == 1)
a <- survdiff(SurvTime ~ Group,data=y.onepartition, rho = 0)
T1.stage1 <- c(((a$exp[2]-a$obs[2])>0)*-1*sqrt(a$chisq / a$var[2,2]) + ((a$exp[2]-a$obs[2])<=0)*sqrt(a$chisq / a$var[2,2]), 1/a$var[2,2]) # Treatment 1 log hazard ratio estimate $\hat{\theta}_{1,1}$ and its variance
T1.stage1LHR <- T1.stage1[1] # T1.stage1LHR is Treatment 1 log hazard ratio estimate $\hat{\theta}_{1,1}$
T1.stage1Fisher <- 1/T1.stage1[2] # Treatment 1 Fisher Information $V_{1,1}$
T1.pvalue <- pnorm(T1.stage1LHR/sqrt(T1.stage1[2])) # Stage 1 Treatment 1 "One - sided" p-value
W1 <- Aj * sqrt(T1.stage1[2]) # This is $W_1$


## Treatment 2 estimates ##

y.onepartition <- y1[y1$Arm != 1,] # Excluding treatment 1 patients
y.onepartition$Group <- 1 * (y.onepartition$Arm == 2)
a <- survdiff(SurvTime ~ Group,data=y.onepartition, rho = 0)
T2.stage1 <- c(((a$exp[2]-a$obs[2])>0)*-1*sqrt(a$chisq / a$var[2,2]) + ((a$exp[2]-a$obs[2])<=0)*sqrt(a$chisq / a$var[2,2]), 1/a$var[2,2]) # Treatment 2 log hazard ratio estimate $\hat{\theta}_{1,2}$ and its variance
T2.stage1LHR <- T2.stage1[1] # Treatment 2 log hazard ratio estimate $\hat{\theta}_{1,2}$
T2.stage1Fisher <- 1/T2.stage1[2] # Treatment 2 Fisher Information $V_{1,2}$
T2.pvalue <- pnorm(T2.stage1[1]/sqrt(T2.stage1[2])) # Stage 1 Treatment 2 "One - sided" p-value
W2 <- Aj * sqrt(T2.stage1[2]) # This is $W_2$


## Covariance between log hazards ratios corresponding to treatments 1 and 2 ##

T0 <- 1 * (yy[,1] == 0) # Indicator for control
T1 <- 1 * (yy[,1] == 1) # Indicator for Treatment 1
T2 <- 1 * (yy[,1] == 2) # Indicator for Treatment 2
z <- cbind(yy, T0, T1, T2) # New matrix with indicator columns for different treatment groups
z <- z[z[,3] == 1,] # Selecting stage 1 patients
T0n <- sum(z[,6]) # Control stage 1 sample size
T1n <- sum(z[,7]) # Treatment 1 stage 1 sample size
T2n <- sum(z[,8]) # Treatment 2 stage 1 sample size

z <- z[order(z[,5]),] # Ordering stage 1 data by survival time
z1 <- cumsum(z[,6]) 
z1 <- T0n - z1 # Number of control patients at risk
z2 <- cumsum(z[,7])
z2 <- T1n - z2 # Number of treatment 1 patients at risk
z3 <- cumsum(z[,8])
z3 <- T2n - z3 # Number of treatment 2 patients at risk
z <- cbind(z, z1, z2, z3) # Adding columns for number of patients at risk
z <- z[z[,4] == 1,] # Selecting cases where there was a relapse. That is, not censored survival time.

z <- z[z[,1] == 0,] # Selecting cases where relapse was from the control
z4 <- z[,9] + z[,10] # Number of patients at risk in control and treatment 1 groups
z5 <- z[,9] + z[,11] # Number of patients at risk in control and treatment 2 groups
z7 <- z[,10] * z[,11] # Multiplying numbers of patients at risk in treatment 1 and 2 groups

z7 <- z7/z4
z7 <- z7/z5 # This is 
covT1T2 <- sum(z7)
covT1T2 <- 	covT1T2 * T1.stage1[2] * T2.stage1[2]
covT1T2
rm(T0, T1, T2, z, T0n, T1n, T2n, z1, z2, z3, z4, z5, z7)


##################################
## NAIVE COMPUTATIONS/ESTIMATES ##
##################################

y1 <- y
y1$SurvTime <- with(y1, Surv(T.date, death))

## Treatment 1 estimates ##

y.onepartition <- y1[y1$Arm != 2,] # Excluding treatment 2 patients
y.onepartition$Group <- 1 * (y.onepartition$Arm == 1)
a <- survdiff(SurvTime ~ Group,data=y.onepartition, rho = 0)
T1.naive <- c(((a$exp[2]-a$obs[2])>0)*-1*sqrt(a$chisq / a$var[2,2]) + ((a$exp[2]-a$obs[2])<=0)*sqrt(a$chisq / a$var[2,2]), 1/a$var[2,2]) # Vector of Treatment 1 log hazard ratio naive estimate $\hat{\theta}_1$ and its variance
T1.naiveLHR <- T1.naive[1] # Treatment 1 log hazard ratio naive estimate $\hat{\theta}_1$
T1.naiveFisher <- 1/T1.naive[2] # Treatment 1 stage 1 and stage 2 data Fisher Information $V_1$


## Treatment 2 estimates ##

y.onepartition <- y1[y1$Arm != 1,] # Excluding treatment 1 patients
y.onepartition$Group <- 1 * (y.onepartition$Arm == 2)
a <- survdiff(SurvTime ~ Group,data=y.onepartition, rho = 0)
T2.naive <- c(((a$exp[2]-a$obs[2])>0)*-1*sqrt(a$chisq / a$var[2,2]) + ((a$exp[2]-a$obs[2])<=0)*sqrt(a$chisq / a$var[2,2]), 1/a$var[2,2]) # Vector of Treatment 2 log hazard ratio naive estimate $\hat{\theta}_2$ and its variance
T2.naiveLHR <- T2.naive[1] # Treatment 2 log hazard ratio naive estimate $\hat{\theta}_2$
T2.naiveFisher <- 1/T2.naive[2] # Treatment 2 stage 1 and stage 2 data Fisher Information $V_2$



########################################
## INDEPENDENT INCREMENTS COMPUTATION ##
########################################

## Treatment 1 increments ##

a <- T1.naiveFisher - T1.stage1Fisher # This is $V_1 - V_{1,1}$
T1.increment <- c(((T1.naive[1]/T1.naive[2]) - (T1.stage1[1]/T1.stage1[2]))/a, 1/a) # Vector of $\hat{\theta}_{2,1}$ and $\sigma^2_{2,1}$
T1.incrementLHR <- T1.increment[1] # This is $\hat{\theta}_{2,1}$
T1.incrementVar <- T1.increment[2] # This is $\sigma^2_{2,1}$


## Treatment 2 increments ##

a <- T2.naiveFisher - T2.stage1Fisher # This is $V_2 - V_{1,2}$
T2.increment <- c(((T2.naive[1]/T2.naive[2]) - (T2.stage1[1]/T2.stage1[2]))/a, 1/a) # Vector of $\hat{\theta}_{2,2}$ and $\sigma^2_{2,2}$
T2.incrementLHR <- T2.increment[1] # This is $\hat{\theta}_{2,2}$
T2.incrementVar <- T2.increment[2] # This is $\sigma^2_{2,2}$



############################################
## COMPUTING UNBIASED ESTIMATES (UMVCUES) ##
############################################

## Treatment 1 estimates ##

b1 <- T1.stage1[2] # This is $\sigma^2_{1,1}$
b2 <- T1.increment[2] # This is $\sigma^2_{2,1}$
T1naive <- ((b2*T1.stage1[1])+(b1*T1.increment[1]))/(b1+b2) # Naive estimate using expression (2) in the main paper. T1.stage1[1] is $\hat{\theta}_{1,1}$ and T1.increment[1] is $\hat{\theta}_{2,1}$.
var1 <- b2/sqrt(b1+b2) # This is $\sigma_{2,1}^2 / \sqrt{\sigma_{1,1}^2 + \sigma_{2,1}^2}$ in expression (12) in main paper.
var2 <- sqrt(b1+b2)/b1 # This is $\sqrt{\sigma_{1,1}^2 + \sigma_{2,1}^2}/\sigma_{1,1}^2$ in $g(W_1)$ in expression (12) in main paper.
gW1 <- var2 * (T1naive - W1) # This is $g_(W_1)$ in expression (12) in the main paper.
T1.unbiased <- T1naive + var1 * ((dnorm(gW1))/(1-pnorm(gW1))) # This is treatment 2 UMVCUE $\hat{\theta}_{1,UMV}$ using expression (12) in the main paper.


## Treatment 2 estimates ##

b1 <- T2.stage1[2] # This is $\sigma^2_{1,2}$
b2 <- T2.increment[2] # This is $\sigma^2_{2,2}$
T2naive <- ((b2*T2.stage1[1])+(b1*T2.increment[1]))/(b1+b2) # Naive estimate using expression (2) in the main paper. T2.stage1[1] is $\hat{\theta}_{1,2}$ and T2.increment[1] is $\hat{\theta}_{2,2}$.
var1 <- b2/sqrt(b1+b2) # This is $\sigma_{2,2}^2 / \sqrt{\sigma_{1,2}^2 + \sigma_{2,2}^2}$ in expression (12) in main paper.
var2 <- sqrt(b1+b2)/b1 # This is $\sqrt{\sigma_{1,2}^2 + \sigma_{2,2}^2}/\sigma_{1,2}^2$ in $g(W_2)$ in expression (12) in main paper.
gW2 <- var2 * (T2naive - W2) # This is $g_(W_2)$ in expression (12) in the main paper.
T2.unbiased <- T2naive + var1 * ((dnorm(gW2))/(1-pnorm(gW2))) # This is treatment 2 UMVCUE $\hat{\theta}_{2,UMV}$ using expression (12) in the main paper.



##########################################################
## COMPUTING SINGLE ITERATION BIAS SUBTRACTED ESTIMATES ##
##########################################################


# The expression for bias for the selection #
# rule and selection made in the example is #
# given in supplementary material Section 7.#


# Function to compute estimates: Code for expression (5) in supplementary material. #

eSIGMAc.TwoTrtsExample <- function(TransVarcov, TransMeans, FutilityThreshold, eMatrix){
	TildeSigmaOne <- TransVarcov[-1,-1] - ((1/TransVarcov[1,1]) * TransVarcov[-1,1] * TransVarcov[1,-1])
	TildeSigmaTwo <- TransVarcov[-2,-2] - ((1/TransVarcov[2,2]) * TransVarcov[-2,2] * TransVarcov[2,-2])

	deltaOneStar <- (FutilityThreshold[1] -  TransMeans[1])/TransVarcov[1,1]
	deltaTwoStar <- (FutilityThreshold[2] -  TransMeans[2])/TransVarcov[2,2]

	TildeDeltaMinusOne <- TransMeans[-1] + (TransVarcov[-1,1] * deltaOneStar)
	TildeDeltaMinusTwo <- TransMeans[-2] + (TransVarcov[-2,2] * deltaTwoStar)

	c1 <- pnorm(FutilityThreshold[2], TildeDeltaMinusOne, sqrt(TildeSigmaOne))
	c1 <- -c1
	c1 <- dnorm(FutilityThreshold[1], TransMeans[1], sqrt(TransVarcov[1,1])) * c1

	c2 <- pnorm(FutilityThreshold[1], TildeDeltaMinusTwo, sqrt(TildeSigmaTwo))
	c2 <- -c2
	c2 <- dnorm(FutilityThreshold[2], TransMeans[2], sqrt(TransVarcov[2,2])) * c2

	c <- c(c1, c2)
	c <- matrix(c, ncol=1)

	bb <- as.numeric(eMatrix %*% TransVarcov %*% c)
	return(bb)
}

# Stage 1 variance-covariance matrix #

SigmaVarcov <- matrix(c(T1.stage1[2], covT1T2, covT1T2, T2.stage1[2]), ncol=2, byrow=T)
NewSigmaVarcov <- SigmaVarcov


# Treatment 1 and Treatment 2 naive estimates as a column vector #

means <- c(T1.naive[1], T2.naive[1])
NewMeans <- matrix(means, ncol=1)
NewMeans <- as.vector(NewMeans)


# Probability of selecting Treatment 1 and Treatment 2. Using second expression in Section 7 in the supplementary material. #

a <- pmvnorm(lower=c(-Inf, -Inf), upper=c(W1,W2), mean=NewMeans, sigma=NewSigmaVarcov) # W1 is $W_1$ and W2 is $W_2$
ProbSel <- a[[1]]


### Treatment 1 single iteration estimate ###

e1 <- matrix(c(1,0), nrow=1) # This is $e_1$ as defined in expression (7) in the main paper
bb1 <- eSIGMAc.TwoTrtsExample(NewSigmaVarcov, NewMeans, c(W1, W2), e1)
t1 <- T1.naive[2]/T1.stage1[2] # This is $t_1$ as defined in expression (3) in the main paper
bias.trt1 <- t1*(bb1/ProbSel) # Estimate of bias for Treatment 1 estimate
T1.singleIter <- T1.naive[1] - bias.trt1 # This is Treatment 1 single iteration bias adjusted estimate


### Treatment 2 single iteration estimate ###

e2 <- matrix(c(0,1), nrow=1) # This is $e_2$ as defined in expression (7) in the main paper
bb2 <- eSIGMAc.TwoTrtsExample(NewSigmaVarcov, NewMeans, c(W1, W2), e2)
t2 <- T2.naive[2]/T2.stage1[2] # This is $t_2$ as defined in expression (3) in the main paper
bias.trt2 <- t2*(bb2/ProbSel) # Estimate of bias for Treatment 2 estimate
T2.singleIter <- T2.naive[1] - bias.trt2 # This is Treatment 2 single iteration bias adjusted estimate



#############################################################
## COMPUTING MULTIPLE ITERATIONS BIAS SUBTRACTED ESTIMATES ##
#############################################################

#  The function for computing multiple iterations estimates involves   #
#  adding iteration step to the single iteration function above. This  #
# includes updating probability of continuing with the two treatments. #
######  More details are in the supplementary material Section 7. ######

# Multiple iterations function #

MultipleIterationsTwoTrtsExample <- function(NaiveEsts, t, TransVarcov1, FutilityThreshold1){
	e1 <- matrix(c(1,0), nrow=1)
	e2 <- matrix(c(0,1), nrow=1)

	y <- NaiveEsts
	zrminus <- NaiveEsts # This object starts with naive estimates and is updated in each iteration step
	k1 <- 0

	for (k in 1:20){
		gg <- sum(is.finite(zrminus) == 0)
		if (gg > 0) {break()}

		zr <- rep(NA, 2)

		NewMeans <- zrminus

		a <- pmvnorm(lower=c(-Inf, -Inf), upper=c(FutilityThreshold1[1], FutilityThreshold1[2]), mean=NewMeans, sigma=TransVarcov1, algorithm = GenzBretz())
		ProbSel <- a[1] # Updated probability of selecting both treatments

		bb <- eSIGMAc.TwoTrtsExample(TransVarcov=TransVarcov1, TransMeans=NewMeans, FutilityThreshold=FutilityThreshold1, eMatrix=e1)
		Trt1Bias <- t[1] * (bb/ProbSel) # Updated bias for treatment 1 naive estimate

		bb1 <- eSIGMAc.TwoTrtsExample(TransVarcov=TransVarcov1, TransMeans=NewMeans, FutilityThreshold=FutilityThreshold1, eMatrix=e2)
		Trt2bias <- t[2] * (bb1/ProbSel) # Updated bias for treatment 2 naive estimate

		zr[1] <- y[1] - Trt1Bias
		zr[2] <- y[2] - Trt2bias

		euc.dis <- sqrt(sum((zrminus - zr)^2))
		if (euc.dis <= 0.001) {zrminus <- zr; k1 <- 1; break()}
			else {zrminus <- zr}
	}
	return(c(k1, zrminus))
}


# Using the function above to compute multiple iterations bias subtracted estimates #

T <- c(t1, t2) # Vector of $t_1$ and $t_2$
MultiEstimates <- MultipleIterationsTwoTrtsExample(NaiveEsts=c(T1.naive[1], T2.naive[1]), t=T, TransVarcov1=NewSigmaVarcov, FutilityThreshold1=c(W1, W2))
T1.multipleIter <- MultiEstimates[2] # This is Treatment 1 multiple iterations bias subtracted estimate
T2.multipleIter <- MultiEstimates[3] # This is Treatment 2 multiple iterations bias subtracted estimate



###################################
## COMPUTING SHRINKAGE ESTIMATES ##
###################################

# Function to implement iteration procedure in Section 6 in the supplementary material to estimate $\nu^2$ and then compute $\nu^2_{+} = \max \{0, \nu^2 \}$ #

EstimatingTauTwoTrts <- function(EigenValues, NaiveEstimates, Overall){
	d1 <- (EigenValues[1])*(EigenValues[1]) # $D_{1,1}^2$ in Section 6 in supplementary material, which is square of first eigenvalue. Vector of eigenvalues is first input in the function.
	d2 <- (EigenValues[2])*(EigenValues[2]) # $D_{2,2}^2$ in Section 6 in supplementary material, which is square of second eigenvalue.

	y <- c(NaiveEstimates[1], NaiveEstimates[2]) # Vector $(\hat{\theta}_{N_1}, \hat{\theta}_{N_1})^\prime$ which is the second input in the function.

	e1 <- ((y[1] - Overall)*(y[1] - Overall)) - d1 # This is $(\hat{\theta}_{1,1} - \hat{\theta}_{1,all}) - D_{1,1}^2$ in Section 6 in supplementary material
	e2 <- ((y[2] - Overall)*(y[2] - Overall)) - d2 # This is $(\hat{\theta}_{1,2} - \hat{\theta}_{1,all}) - D_{2,2}^2$ in Section 6 in supplementary material

	zrminus <- 0.05 # Initial guess for $\nu^2$

		for (k in 1:200){
			zr <- NA

			w1 <- 1/(zrminus + d1) # This is $w_1 = 1/(\nu^2+D_{1,1}^2)$
			w2 <- 1/(zrminus + d2) # This is $w_2 = 1/(\nu^2+D_{2,2}^2)$

			w <- w1 + w2 # This is $\sum_{i=1}^{2} w_i = w_1 + w_2$

			ww1 <- (w1 * e1) + (w2 * e2) # This is $w_1[(\hat{\theta}_{1,1} - \hat{\theta}_{1,all}) - D_{1,1}^2] + w_2[(\hat{\theta}_{1,2} - \hat{\theta}_{1,all}) - D_{2,2}^2]$

			zr <- ww1/w # This is updated estimated value for $\nu^2$

			euc.dis <- sqrt((zrminus - zr)*(zrminus - zr))
				if (euc.dis <= 0.0001) {zrminus <- zr; break()}
					else {zrminus <- zr}
	}
	zrminus <- max(0,zrminus)
	return(c(zrminus,k))
}


# Computing log hazard ratio using stage 1 data when treatment 1 and treatment 1 are combined into a single group #

y1 <- y[y$Stage == 1,] # Selecting stage 1 patients
y1$SurvTime <- with(y1, Surv(T.date1, death1))
y1$Group1 <- 1 * (y1$Arm > 0) # Indicator for patients randomised to treatments 1 & 2

a <- survdiff(SurvTime ~ Group1, data=y1, rho = 0)
stage1overall <- ((a$exp[2]-a$obs[2])>0)*-1*sqrt(a$chisq / a$var[2,2]) + ((a$exp[2]-a$obs[2])<=0)*sqrt(a$chisq / a$var[2,2])


# Computing the eigenvalues #

test <- svd(NewSigmaVarcov) # Single value decomposition of variance-covariance matrix for stage 1 log hazard ratios. NewSigmaVarcov is defined assigned values earlier while computing single iteration bias subtracted estimates.
EigenValues1 <- test$d # Extracting eigenvalues


# Computing $C$ in main paper Section 3.2 and then stage 1 data shrinkage estimate as $C \boldsymbol{\hat{\theta}}_1 + (\textbf{\emph{I}}_{K} - C) (\hat{\theta}_{1, all}, \hat{\theta}_{1, all})^\prime$ where $\boldsymbol{\hat{\theta}}_1 = (\hat{\theta}_{1,1}, \hat{\theta}_{1,2})^\prime$ #

nu.squared <- EstimatingTauTwoTrts(EigenValues1, c(T1.naive[1], T2.naive[1]), stage1overall)[1] # nu.squared is $\tilde{\nu}^2_{+}$ computed using function above.
new.matrix <- nu.squared * diag(2) + NewSigmaVarcov # This is $\nu^2 \textbf{\emph{I}}_2 + \Sigma_{\hat{\theta}_1}$
new.matrix <- solve(new.matrix) # This is $\left ( \nu^2 \textbf{\emph{I}}_2 + \Sigma_{\hat{\theta}_1} \right )^{-1}$
new.matrix <- diag(2) - NewSigmaVarcov %*% new.matrix # # This is $C = \textbf{\emph{I}}_2 - \Sigma_{\hat{\theta}_1} \left ( \nu^2 \textbf{\emph{I}}_2 + \Sigma_{\hat{\theta}_1} \right )^{-1}$ defined in Section 3.2 in the main paper.

ShrinkageEstimate <- new.matrix %*% matrix(c(T1.stage1[1], T2.stage1[1]), ncol=1) + ((diag(2) - new.matrix) %*% matrix(rep(stage1overall,2), ncol=1)) # This is $C \boldsymbol{\hat{\theta}}_1 + (\textbf{\emph{I}}_{2} - C) (\hat{\theta}_{1, all}, \hat{\theta}_{1, all})^\prime$
ShrinkageEstimate <- as.vector(ShrinkageEstimate) # Vector of stage 1 shrinkage estimates $\left ( \hat{\theta}^{(1)}_{1,SH}, \hat{\theta}^{(1)}_{2,SH} \right )^\prime$


# Computing two-stage shrinkage estimates #

T1.ShrinkWeight <- 0.5 # This is $\omega$ in two-stage shrinkage estimator
T1.ShrinkageEstimate <- (T1.ShrinkWeight * ShrinkageEstimate[1]) + ((1 - T1.ShrinkWeight) * T1.increment[1]) # Computing $\omega \hat{\theta}^{(1)}_{1,SH} + (1 - \omega) \hat{\theta}_{2,1}$, treatment 1 shrinkage estimate $\hat{\theta}_{1,SH}$

T2.ShrinkWeight <- 0.5 # This is $\omega$ in two-stage shrinkage estimator
T2.ShrinkageEstimate <- (T2.ShrinkWeight * ShrinkageEstimate[2]) + ((1 - T2.ShrinkWeight) * T2.increment[1]) # Computing $\omega \hat{\theta}^{(1)}_{2,SH} + (1 - \omega) \hat{\theta}_{2,2}$, treatment 2 shrinkage estimate $\hat{\theta}_{2,SH}$




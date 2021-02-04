#functions used to produce statistics reported in Baird et al.,
#"Developmental basis for climatic distribution of grass leaf size worldwide"

#required packages
library(nlme)
library(ape)
library(phytools)


#function to obtaining likelihood of pgls, used with optimize when estimating Lambda
opt.Lambda <- function(Lambda, form, tree, data){
							m <- gls(form,
												correlation = corPagel(Lambda, tree, fixed = TRUE),
												data = data,
												method = "ML")
							m$logLik
							}

#function for anova output comparing C3 and C4 mean values
#trt is a character vector giving trait names
#data is a dataframe, pre-formatted so that row.names matches tip.labels in phylo object,
#and with column "C3C4" giving C3/C4 categorisation
#tree is a phylo object matching the species list included in data
C3C4.ANOVA <- function(trt, data, tree){
	
	#create and format output dataframe
	anova.hds<-c(
		#least squares outputs
		"LS.C3.mean", "LS.C3.sem", "LS.C4.mean", "LS.C4.sem",
		"LS.F", "LS.P",
		"LS.resid.norm.stat", "LS.resid.norm.P",
		"LS.resid.homo.stat", "LS.resid.homo.P",
		
		#PGLS outputs using lambda to adjust for phylogenetic covariance
		"PGLS.C3.mean", "PGLS.C3.sem",
		"PGLS.C4.mean", "PGLS.C4.sem",
		"PGLS.F", "PGLS.P",
		"PGLS.lambda", "PGLS.lambda.P",
		"PGLS.phyloresid.norm.stat", "PGLS.phyloresid.norm.P",
		"PGLS.phyloresid.homo.stat", "PGLS.phyloresid.homo.P",
		"PGLS.phyloresid.lambda", "PGLS.phyloresid.lambda.P",
		
		#Garland et al 1993 phylogenetic ANOVA, and Blomberg's K in case a trait fails all assumptions
		"phylANOVA.F", "phylANOVA.P", "phylosig.K", "phylosig.K.P"
		)

	anova.fr <- data.frame(matrix(NA, nrow = length(trt), ncol = length(anova.hds)))
	row.names(anova.fr) <- trt
	names(anova.fr) <- anova.hds

	#create dataframe to hold residuals from LS fits
	LS.res <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(trt)))
	names(LS.res) <- trt
	row.names(LS.res) <- row.names(data)

	#same for PGLS fits
	PGLS.res <- LS.res

	#loop fitting both LS and PGLS models,
	#extracting coefficients, residuals, and testing residuals for assumptions
	for (i in 1:nrow(anova.fr)){

		ti <- row.names(anova.fr)[i]
		print(ti)#this just lets you watch the loop
		di = data[ , c("C3C4", ti)]
		
		#model formulae
		#1) intercept not estimated so that mean/se can be easily extracted from output
		f.c <- formula(paste0(ti, " ~ factor(C3C4) - 1"))		
		#2) with intercept to allow extraction of appropriate F & P
		f.1 <- formula(paste0(ti, " ~ factor(C3C4)"))

		#corStruct for lambda=0
		cor.0 = corPagel(0, tree, fixed = T)
		#use optimize to constrain an estimate of Lambda between 0 & 1
		opt.l <- optimize(f = opt.Lambda,
											interval = c(0, 1),
											maximum = TRUE,
											form = f.c,
											tree = tree,
											data = di
											)
		#create a corStruct for this lambda
		cor.l = corPagel(opt.l$maximum, tree, fixed = T)
		
		#fit coefficient extraction model with optimized lambda (PGLS)
		PGLS.coef <- gls(f.c, correlation = cor.l, data = di, method = "ML")
		#fit coefficient extraction model with lambda=0 (LS)
		LS.coef <- gls(f.c, correlation = cor.0, data = di, method = "ML")					

		#fit statistic extraction models
		PGLS.FP <- gls(f.1, correlation = cor.l, data = di, method = "ML")
		LS.FP <- gls(f.1, correlation = cor.0, data = di, method = "ML")					

		#extract coefficients
		anova.fr[i, c("LS.C3.mean", "LS.C4.mean")] <- coef(LS.coef)[c(1, 2)]
		anova.fr[i, c("LS.C3.sem", "LS.C4.sem")] <- summary(LS.coef)$tTable[c(1, 2), 2]
		anova.fr[i, c("LS.F")] <- round(anova(LS.FP)$F[2], 2)
		anova.fr[i, c("LS.P")] <- round(anova(LS.FP)$p[2], 4)

		anova.fr[i, c("PGLS.C3.mean", "PGLS.C4.mean")] <- coef(PGLS.coef)[c(1, 2)]
		anova.fr[i, c("PGLS.C3.sem", "PGLS.C4.sem")] <- summary(PGLS.coef)$tTable[c(1, 2), 2]
		anova.fr[i, c("PGLS.F")] <- round(anova(PGLS.FP)$F[2], 2)
		anova.fr[i, c("PGLS.P")] <- round(anova(PGLS.FP)$p[2], 4)
		anova.fr[i, c("PGLS.lambda")] <- round(opt.l$maximum, 3)
		#approx p-value for lambda == 0
		pL0 <- 1 - pchisq(2 * (PGLS.coef$logLik - LS.coef$logLik), 1)
		anova.fr[i, c("PGLS.lambda.P")] <- round(pL0, 4)
		
		#extract residuals
		LS.res[ , ti] <- resid(LS.coef)
		#the below are corrected for effect of phylogeny
		PGLS.res[ , ti] <- resid(PGLS.coef, type = "normalized")
 		#normality testing for LS.res using Anderson Darling
		LS.AD <- ad.test(scale(LS.res[ , ti]))
		anova.fr[i, "LS.resid.norm.stat"] <- round(LS.AD$statistic, 2)
		anova.fr[i, "LS.resid.norm.P"] <- round(LS.AD$p.value, 4)

		#homogeneity testing for LS.res using Bartlett
		LS.B <- bartlett.test(scale(LS.res[ , ti]) ~ factor(di[ , "C3C4"]))
		anova.fr[i, "LS.resid.homo.stat"] <- round(LS.B$statistic, 2)
		anova.fr[i, "LS.resid.homo.P"] <- round(LS.B$p.value, 4)

		#normality testing for PGLS.res using Anderson Darling
		PGLS.AD <- ad.test(scale(PGLS.res[ , ti]))
		anova.fr[i, "PGLS.phyloresid.norm.stat"] <- round(PGLS.AD$statistic, 2)
		anova.fr[i, "PGLS.phyloresid.norm.P"] <- round(PGLS.AD$p.value, 4)

		#homogeneity testing for PGLS.res using Bartlett
		PGLS.B <- bartlett.test(scale(PGLS.res[ , ti]) ~ factor(di[ , "C3C4"]))
		anova.fr[i, "PGLS.phyloresid.homo.stat"] <- round(PGLS.B$statistic, 2)
		anova.fr[i, "PGLS.phyloresid.homo.P"] <- round(PGLS.B$p.value, 4)
		
		#is lambda 0 for the PGLS residuals?
		f.r1 <- formula(paste0(ti, " ~ 1"))
		#estimate lambda		
		opt.l.res <- optimize(f = opt.Lambda,
													interval = c(0, 1),
													maximum = TRUE,
													form = f.r1,
													tree = tree,
													data = PGLS.res
													)
		#create corStruct		
		cor.l.res = corPagel(opt.l.res$maximum, tree, fixed = TRUE)
		#fit PGLS (with lambda)		
		PGLS.r1 <- gls(f.r1,
										correlation = cor.l.res,
										data = PGLS.res,
										method = "ML"
										)	
		#fit LS (lambda = 0)
		LS.r1 <- gls(f.r1,
									correlation = cor.0,
									data = PGLS.res,
									method = "ML"
									)	

		anova.fr[i, "PGLS.phyloresid.lambda"] <- round(opt.l.res$maximum, 3)
		#approx p-value for lambda == 0
		pL0.resids <- 1 - pchisq(2 * (PGLS.r1$logLik - LS.r1$logLik), 1)
		anova.fr[i, "PGLS.phyloresid.lambda.P"] <- round(pL0.resids, 4 )
		
		#non parametric tests, just in case
		#formatting the data so these tests work
		xf <- di[ , "C3C4"]
		names(xf) <- row.names(di)
		yv <- di[ , ti]
		names(yv) <- row.names(di)
		
		#non-parametric, i.e., permutation/simulation based phylogenetic ANOVA: phylANOVA (phytools)		
		pA <- phylANOVA(tree, xf, yv, posthoc = FALSE)
		anova.fr[i, "phylANOVA.F"] <- pA$F
		anova.fr[i, "phylANOVA.P"] <- pA$P
		
		#non-parametric, i.e., permutation/simulation based phylogenetic signal
		#via Blomberg's K: phylosig (phytools)		
		BK <- phylosig(tree, yv, method = "K", test = TRUE)
		anova.fr[i, "phylosig.K"] <- BK$K
		anova.fr[i, "phylosig.K.P"] <- BK$P		

		}
	
	return(
					list(summary = anova.fr,
								LS.resids = LS.res,
								PGLS.resids = PGLS.res,
								data = data
								)
					)

}

#functions called by phylo.REGRESSION
	
#correlation coefficient, based on the covariance matrix of phyl.RMA
#rma is a fitted phyl.RMA object
r.RMA <- function(rma){

	r <- rma$V[1, 2] / (sqrt(rma$V[1, 1]) * sqrt(rma$V[2, 2]))
	lt = sign(coef(rma)[2]) == -1
	r.q = r * sqrt(length(rma$resid) - 2) / sqrt(1 - r^2)
	r.df = length(rma$resid) - 2
	p.r <- 2 * pt(q = r.q, df = r.df, lower.tail = lt)

	return(list(r = r,	t.stat = r.q,	df = r.df, P = p.r))

} 
		
#95% confidence intervals for slope of RMA (following Warton et al 2006)
#relies on r.RMA
cislope.RMA <- function(rma){

	r <- r.RMA(rma)		
	B <- (1 - r$r^2) / r$df * qf(1 - 0.05, 1, r$df)
	ci.u <- rma$RMA.beta[2] * (sqrt(B + 1) + sqrt(B))
	ci.l <- rma$RMA.beta[2] * (sqrt(B + 1) - sqrt(B))
	
	return(c(ci.l, ci.u))

	}

#95% confidence intervals for intercept (int) of RMA fit (following Warton et al 2006)
#requires:
#r.RMA
#X is x variable data from the original fit
#output not used in paper - see note below
ciint.RMA <- function(rma, X){

	r <- r.RMA(rma)
	V.s <- (1 / r$df) * (rma$V[2, 2] / rma$V[1, 1]) * (1 - r$r^2)
	#using raw residuals, which may be strictly incorrect when lambda is used
	V.i <- (var(rma$resid) / length(X)) + mean(X)^2 * V.s
	s.i <- sqrt(V.i)
	ci.u <- rma$RMA.beta[1] - s.i * qt(p = 1 - (0.05 / 2), df = r$df)
	ci.l <- rma$RMA.beta[1] + s.i * qt(p = 1 - (0.05 / 2), df = r$df )
	
	return(c(ci.l, ci.u))

	}
		
#transform residuals from phyl.RMA to phylogenetically corrected residuals
#these are equivalent to resid.gls(... , type = "n")
#adapted from code originally obtained from R.P. Freckleton
p.resids<-function(resids, lambda, tree){
				
	V <- vcv(tree)
	V1 <- V
	diag(V1) <- 0
	V2 <- diag(diag(V), ncol = length(V[1, ]), nrow = length(V[, 1]))
	Vmat <- V1 * lambda + V2
	phylo.resids <- t(solve(chol(Vmat))) %*% resids
	
	return(phylo.resids)
	
	}

#function that returns rounded regression output using both LS/PGLS and RMA/phylRMA
#y.trt is a character vector, naming traits to treat as dependent
#x.trt is a character vector, naming traits to treat as independent
#data is a data matrix containing x and y variates
#tree is a phylo object with tip.labels matching the row.names of data
phylo.REGRESSION<-function(y.trt, x.trt, data, tree){

	R.y <- data[ , y.trt]
	R.x <- data[ , x.trt]

	R.hds <- c(
		#variables under test
		"x", "y",
		#phylogenetic correlation coefficient after Revell, and t-test
		"r", "r.t", "r.df", "r.P",
		"phylo.r", "phylo.r.t", "phylo.r.df", "phylo.r.P",
		#LS regression
		"LS.intercept", "LS.intercept.ci95.lo", "LS.intercept.ci95.hi",
		"LS.intercept.F", "LS.intercept.P",
		"LS.slope", "LS.slope.ci95.lo", "LS.slope.ci95.hi",
		"LS.slope.F", "LS.slope.P",
		"LS.resid.norm.stat", "LS.resid.norm.P",
		#PGLS regression
		"PGLS.intercept", "PGLS.intercept.ci95.lo", "PGLS.intercept.ci95.hi",
		"PGLS.intercept.F", "PGLS.intercept.P",
		"PGLS.slope", "PGLS.slope.ci95.lo", "PGLS.slope.ci95.hi",
		"PGLS.slope.F", "PGLS.slope.P",
		"PGLS.lambda",
		"PGLS.lambda.L", "PGLS.lambda.P",
		"PGLS.phyloresid.norm.stat", "PGLS.phyloresid.norm.P",
		"PGLS.phyloresid.lambda",
		"PGLS.phyloresid.lambda.L", "PGLS.phyloresid.lambda.P",
		#RMA
		"RMA.intercept", "RMA.intercept.ci95.lo", "RMA.intercept.ci95.hi",
		#"RMA.intercept.F","RMA.intercept.P",
		"RMA.slope", "RMA.slope.ci95.lo", "RMA.slope.ci95.hi",
		#"RMA.slope.F","RMA.slope.P",
		"RMA.resid.norm.stat", "RMA.resid.norm.P",
		#phyloRMA
		"PRMA.intercept", "PRMA.intercept.ci95.lo", "PRMA.intercept.ci95.hi",
		"PRMA.slope", "PRMA.slope.ci95.lo", "PRMA.slope.ci95.hi",
		"PRMA.lambda",
		"PRMA.lambda.L", "PRMA.lambda.P",
		"PRMA.phyloresid.norm.stat", "PRMA.phyloresid.norm.P",
		"PRMA.phyloresid.lambda",
		"PRMA.phyloresid.lambda.L", "PRMA.phyloresid.lambda.P"
		)

	#create dataframe for output
	R.fr <- data.frame(matrix(NA,	
														nrow = length(y.trt) * length(x.trt),
														ncol = length(R.hds)
														)
											)
	names(R.fr) <- R.hds
	R.fr$y <- rep(y.trt, length(x.trt))
	R.fr$x <- rep(x.trt, each = length(y.trt))
	#remove the self-self relationships
	R.fr <- R.fr[!(R.fr$x == R.fr$y), ]

	#create dataframe to hold residuals from LS fits
	LS.res <- data.frame(matrix(NA,
															nrow = nrow(data),
															ncol = nrow(R.fr)
															)
												)
	names(LS.res) <- paste(R.fr$y, R.fr$x, sep = "_")
	row.names(LS.res) <- row.names(data)
		
	#repeat for PGLS, RMA, phylRMA
	PGLS.res <- LS.res
	RMA.res <- LS.res
	PRMA.res <- LS.res

	#loop over the x-y combinations and fit all the alternative models
	for (i in 1:nrow(R.fr)) {
	
		#extract x and y for use in all fits
	
		if (!is.null(ncol(R.y))) { y = R.y[ , R.fr$y[i]] } else { y = R.y }
		names(y) <- row.names(data)

		if (!is.null(ncol(R.x))) { x = R.x[ , R.fr$x[i]] } else { x = R.x }
		names(x)<-row.names(data)

		#re-format for LS/PGLS
		R.d <- data.frame(y, x)
		names(R.d) <- c(R.fr$y[i], R.fr$x[i])
		row.names(R.d) <- row.names(data)

		#formula for y ~ x
		f.xy = formula(paste(R.fr$y[i], "~", R.fr$x[i], sep = " "))

		#produce relevant corStruct objects, including estimation of lambda
		#corStruct for lambda == 0
		cor.0 = corPagel(0, tree, fixed = TRUE)
		#use optimize to constrain an estimate of Lambda between 0 & 1
		opt.l <- optimize(f = opt.Lambda,
											interval = c(0, 1),
											maximum = TRUE,
											form = f.xy,
											tree = tree,
											data = R.d
											)
		#create a corStruct for this lambda
		cor.l = corPagel(opt.l$maximum, tree, fixed = TRUE)

		#fit model with optimized lambda (PGLS)
		PGLS.xy <- gls(f.xy, correlation = cor.l, data = R.d, method = "ML")
		#fit model with lambda == 0 (LS)						
		LS.xy <- gls(f.xy, correlation = cor.0, data = R.d, method = "ML")					

		#what are predictions, ci, significance?
		#LS
		R.fr[i, c("LS.intercept", "LS.slope")] <- summary(LS.xy)$tTable[ , "Value"]
		R.fr[i, c("LS.intercept.F", "LS.slope.F")] <- round(anova(LS.xy)[ , "F-value"], 2)
		R.fr[i, c("LS.intercept.P", "LS.slope.P")] <- round(anova(LS.xy)[ , "p-value"], 5)
		#use confint here
		LSci.nms <- c("LS.intercept.ci95.lo",
									"LS.slope.ci95.lo",
									"LS.intercept.ci95.hi",
									"LS.slope.ci95.hi"
									)
		R.fr[i, LSci.nms] <- c(confint(LS.xy))
		#PGLS
		R.fr[i, "PGLS.lambda"] <- round(opt.l$maximum, 3)
		R.fr[i, "PGLS.lambda.L"] <- round(2 * (PGLS.xy$logLik - LS.xy$logLik), 2)
		R.fr[i, "PGLS.lambda.P"] <- round(1 - pchisq(2 * (PGLS.xy$logLik - LS.xy$logLik), 1), 5)
		R.fr[i, c("PGLS.intercept", "PGLS.slope")] <- summary(PGLS.xy)$tTable[ , "Value"]
		R.fr[i, c("PGLS.intercept.F", "PGLS.slope.F")] <- round(anova(PGLS.xy)[ , "F-value"], 2)
		R.fr[i, c("PGLS.intercept.P", "PGLS.slope.P")] <- round(anova(PGLS.xy)[ , "p-value"], 5)
		PGLSci.nms <- c("PGLS.intercept.ci95.lo",
										"PGLS.slope.ci95.lo",
										"PGLS.intercept.ci95.hi",
										"PGLS.slope.ci95.hi"
										)
		R.fr[i, PGLSci.nms] <- c(confint(PGLS.xy))

		#output residuals and test normality		
		#extract residuals
		LS.res[,i]<-resid(LS.xy)
		PGLS.res[,i]<-resid(PGLS.xy,type="normalized")#these are corrected for effect of phylogeny
		
 		#normality testing for LS.res using Anderson Darling
		LS.AD<-ad.test( scale( LS.res[,i] ) )
		R.fr[ i, "LS.resid.norm.stat" ]<-round( LS.AD$statistic, 2 )
		R.fr[ i, "LS.resid.norm.P" ]<-round( LS.AD$p.value, 5 )
		
		#normality testing for PGLS.res using Anderson Darling
		PGLS.AD <- ad.test(scale(PGLS.res[ , i]))
		R.fr[i, "PGLS.phyloresid.norm.stat"] <- round(PGLS.AD$statistic, 2)
		R.fr[i, "PGLS.phyloresid.norm.P"] <- round(PGLS.AD$p.value, 5)

		#is lambda 0 for the PGLS residuals?
		#null formula for testing lambda of resids	
		f.r1 <- formula(paste(names(PGLS.res)[i], " ~ 1", sep = ""))
		#estimate lambda		
		PGLS.opt.l.res <- optimize(f = opt.Lambda,
																interval = c(0, 1),
																maximum = TRUE,
																form = f.r1,
																tree = tree,
																data = PGLS.res
																)
		#create corStruct		
		PGLS.cor.l.res = corPagel(PGLS.opt.l.res$maximum, tree, fixed = TRUE)
		#fit PGLS (lambda from residuals) 		
		PGLS.r1 <- gls(f.r1,
										correlation = PGLS.cor.l.res,
										data = PGLS.res,
										method = "ML"
										)	
		#fit LS (lambda = 0)
		LS.r1 <- gls(f.r1,
									correlation = cor.0,
									data = PGLS.res,
									method = "ML"
									)
		#lambda for residuals
		R.fr[i, "PGLS.phyloresid.lambda"] <- round(PGLS.opt.l.res$maximum, 3)
		#likelihood ratio
		R.fr[i, "PGLS.phyloresid.lambda.L"] <- round(2 * (PGLS.r1$logLik - LS.r1$logLik), 2)
		#approx p-value for lambda==0
		p.PGLS.resid <- 1 - pchisq(2 * (PGLS.r1$logLik - LS.r1$logLik), 1)
		R.fr[i, "PGLS.phyloresid.lambda.P"] <- round(p.PGLS.resid, 5)
		
		#there will be a series of warnings here... linked with P == 1 test for slope
		#output from that test is not relied upon
		#fit RMA with optimised lambda
		PRMA.xy <- phyl.RMA(x,
												y,
												tree,
												method = "lambda",
												lambda = NULL,
												fixed=FALSE,
												h0 = 1
												)
		#RMA with lambda == 0
		RMA.xy <- phyl.RMA(x,
											 y,
											 tree,
											 method = "lambda",
											 lambda = 0,
											 fixed = TRUE,
											 h0 = 1
											 )

		#what are predictions, ci, significance?
		#RMA
		R.fr[i, c("RMA.intercept", "RMA.slope")] <- RMA.xy$RMA.beta
		R.fr[i, c("RMA.intercept.ci95.lo", "RMA.intercept.ci95.hi")] <- ciint.RMA(RMA.xy, x)
		R.fr[i, c("RMA.slope.ci95.lo", "RMA.slope.ci95.hi")] <- cislope.RMA(RMA.xy)
		#PRMA
		R.fr[i, "PRMA.lambda"] <- round(opt.l$maximum, 3)
		R.fr[i, "PRMA.lambda.L"] <- round(2 * (PRMA.xy$logL - RMA.xy$logL), 2)
		p.PRMA.lambda <- 1 - pchisq(2 * (PRMA.xy$logL - RMA.xy$logL), 1)
		R.fr[i, "PRMA.lambda.P"] <- round(p.PRMA.lambda, 5)
		R.fr[i, c("PRMA.intercept", "PRMA.slope")] <- PRMA.xy$RMA.beta
		R.fr[i, c("PRMA.intercept.ci95.lo", "PRMA.intercept.ci95.hi")] <- ciint.RMA(PRMA.xy, x)
		R.fr[i, c("PRMA.slope.ci95.lo", "PRMA.slope.ci95.hi")] <- cislope.RMA(PRMA.xy)

		#RMA residuals and tests
		RMA.res[ , i] <- RMA.xy$resid
		PRMA.res[ , i] <- p.resids(PRMA.xy$resid, PRMA.xy$lambda, tree)

 		#normality testing for RMA.res using Anderson Darling
		RMA.AD <- ad.test(scale(RMA.res[ , i]))
		R.fr[i, "RMA.resid.norm.stat"] <- round(RMA.AD$statistic, 2)
		R.fr[i, "RMA.resid.norm.P"] <- round(RMA.AD$p.value, 5)

		#normality testing for PRMA.res using Anderson Darling
		PRMA.AD <- ad.test(scale(PRMA.res[ , i]))
		R.fr[i, "PRMA.phyloresid.norm.stat"] <- round(PRMA.AD$statistic, 2)
		R.fr[i, "PRMA.phyloresid.norm.P"] <- round(PRMA.AD$p.value, 5)
		
		#is lambda 0 for the PRMA residuals?
		#null formula for testing lambda of resids	is as above	
		
		#estimate lambda		
		PRMA.opt.l.res <- optimize(f = opt.Lambda,
																interval = c(0, 1),
																maximum = TRUE,
																form = f.r1,
																tree = tree,
																data = PRMA.res
																)
		#create corStruct		
		PRMA.cor.l.res = corPagel(PRMA.opt.l.res$maximum, tree, fixed = TRUE)
		#fit PGLS (lambda from residuals) 		
		PGLS.PRMAr1 <- gls(f.r1,
												correlation = PRMA.cor.l.res,
												data = PRMA.res,
												method = "ML"
												)	
		#fit LS (lambda == 0)
		LS.PRMAr1 <- gls(f.r1,
											correlation = cor.0,
											data = PRMA.res,
											method = "ML"
											)

		R.fr[i, "PRMA.phyloresid.lambda"] <- round(PRMA.opt.l.res$maximum, 3)
		#likelihood ratio
		PRMA.LR <- 2 * (PGLS.PRMAr1$logLik - LS.PRMAr1$logLik)
		R.fr[i, "PRMA.phyloresid.lambda.L"] <- round(PRMA.LR, 2)
		#approx p-value for lambda == 0
		p.PRMA.lambda <- 1 - pchisq(2 * (PGLS.PRMAr1$logLik - LS.PRMAr1$logLik), 1)
		R.fr[i, "PRMA.phyloresid.lambda.P"] <- round(p.PRMA.lambda, 5)
	
		r.nRMA <- r.RMA(RMA.xy)
		R.fr[i, "r"] <- round(r.nRMA$r, 3)
		R.fr[i, "r.t"] <- round(r.nRMA$t.stat, 2)
		R.fr[i, "r.df"] <- round(r.nRMA$df, 0)
		R.fr[i, "r.P"] <- round(r.nRMA$P, 5)

		r.PRMA <- r.RMA(PRMA.xy)
		R.fr[i, "phylo.r"] <- round(r.PRMA$r, 3)
		R.fr[i, "phylo.r.t"] <- round(r.PRMA$t.stat, 2)
		R.fr[i, "phylo.r.df"] <- round(r.PRMA$df, 0)
		R.fr[i, "phylo.r.P"] <- round(r.PRMA$P, 5)
		
		}

	#generate output	
	return(
					list(summary = R.fr,
								LS.resids = LS.res,
								PGLS.resids = PGLS.res,
								RMA.resids = RMA.res,
								PRMA.resids = PRMA.res,
								data = data
								)
					)
	
	}

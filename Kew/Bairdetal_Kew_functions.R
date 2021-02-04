##Functions for analysis of Baird et al Kew leaf trait dataset,  climate traits and phylogeny

#function for obtaining optimal value of lambda
opt.Lambda <- function(Lambda, form, tree, dat){
							m <- gls(form,
												correlation = corPagel(Lambda, tree, fixed = TRUE),
												data = dat,
												method = "ML")
							m$logLik
							}

#function that fits PGLS and OLS models
#and extracts a couple of key stats for convenience
#ouputs a named list of models and stats
#takes a single formula as input
pullmr <- function(ffi, tree, data){
	
	print(ffi)
	
	#dependent variable
	trt <- all.vars(ffi)[1]
	
	#subset data and tree
	dat <- data[, all.vars(ffi)]
	dat <- na.omit(dat)
	
	drt <- tree$tip.label[!tree$tip.label %in% rownames(dat)]
	tree.dat <- drop.tip(tree, drt)
	
	#using formulae with gls/lme doesn't work with update
	#hence, below setting out alternative models more explicitly
	
	#find lambda
	L.opt <- optimize(f = opt.Lambda,
										interval = c(0, 1),
										maximum = TRUE,
										form = ffi,
										tree = tree.dat,
										dat = dat
										)

	#model with lambda
	opt <- gls(ffi,
							correlation = corPagel(L.opt$maximum, tree.dat, fixed = TRUE),
							data = dat,
							method = "ML"
							)
							
	#model with lambda = 0
	opt.L0 <- gls(ffi,
								correlation = corPagel(0, tree.dat, fixed = TRUE),
								data = dat,
								method = "ML"
								)
	
	#Likelihood Ratio Test for lambda
	L.opt.LRT <- 1 - pchisq(2 * (opt$logLik - opt.L0$logLik), 1)
	
	#test if lambda is > zero and the distribution is normal for phylogenetic residuals
	
	#get phylogenetic residuals
	phy.resid <- data.frame(resid(opt, "normalized"))
	names(phy.resid) = trt

	#find their optimal lambda
	L.opt.phy.resid <- optimize(f = opt.Lambda,
															interval = c(0, 1),
															maximum = TRUE,
															form = formula(terms(ffi)[0]),
															tree = tree.dat,
															dat = phy.resid
															)

	#model with lambda
	phy.resid.opt <- gls(formula(terms(ffi)[0]),
												correlation = corPagel(L.opt.phy.resid$maximum, tree.dat, fixed = TRUE),
												data = phy.resid,
												method = "ML"
												)	
	#model with lambda = 0
	phy.resid.0 <- gls(formula(terms(ffi)[0]),
											correlation = corPagel(0, tree.dat, fixed = TRUE),
											data = phy.resid,
											method = "ML"
											)
	
	#Likelihood Ratio Test
	L.opt.phy.resid.LRT <- 1 - pchisq(2*(phy.resid.opt$logLik - phy.resid.0$logLik), 1)
	
	#normality test
	phy.resid.norm <- ad.test(scale(phy.resid))

	#get OLS residuals and test if they are normal 
	OLS.resid <- resid(opt.L0, "normalized")
	resid.norm <- ad.test(scale(OLS.resid))
	
	list(	opt = opt,
				opt.L0 = opt.L0,
				L.opt = L.opt,
				L.opt.LRT = L.opt.LRT,
				L.opt.phy.resid = L.opt.phy.resid,
				L.opt.phy.resid.LRT = L.opt.phy.resid.LRT,
				phy.resid.norm = phy.resid.norm,
				resid.norm = resid.norm
				)
	}

#function used in pullmr.stats to grab tTable output from pullmr output
#tTname = term name from a model
#pmr = pullmr output
tTab.extract <- function(tTname, pmr){
	
	ols.jnms <- paste0("OLS.", c("coef", "se", "t", "t.P", "F", "F.P"), ".", tTname)
	pgls.jnms <- paste0("PGLS.", c("coef", "se", "t", "t.P", "F", "F.P"), ".", tTname)
	
	scrap <- rep(NA, length(ols.jnms) + length(pgls.jnms))
	names(scrap) <- c(ols.jnms, pgls.jnms)
	
	scrap[ols.jnms] <- c(summary(pmr$opt.L0)$tTable[tTname, ],
													as.data.frame(anova(pmr$opt.L0))[tTname, c("F-value", "p-value")]
													)
	
	scrap[pgls.jnms] <- c(summary(pmr$opt)$tTable[tTname, ],
													as.data.frame(anova(pmr$opt))[tTname, c("F-value", "p-value")]
													)
	
	scrap

}

#function to format select outputs from a single pullmr item
#will produce a single row dataframe
#pmr = pullmr output
#tt = complete list of terms for the maximal model
#note that terms in tt may not all be included in the model being processed
#so this function only works efficiently within phylo.MR
pullmr.stats <- function(pmr, tt){	
	
	#names for ouput vector
	out.nms <- c(
								#OLS
								#AIC
								"AIC.OLS", # = NA, #AIC(opt.L0), 
								#coefs
								paste0("OLS.", c("coef", "se", "t", "t.P"), ".", "(Intercept)"), 
								paste0("OLS.", c("coef", "se", "t", "t.P"), ".", tt[1]), 
								paste0("OLS.", c("coef", "se", "t", "t.P"), ".", tt[2]), 
								paste0("OLS.", c("coef", "se", "t", "t.P"), ".", tt[3]), 
								#F & P
								paste0("OLS.", c("F", "F.P"), ".", "(Intercept)"), 
								paste0("OLS.", c("F", "F.P"), ".", tt[1]), 
								paste0("OLS.", c("F", "F.P"), ".", tt[2]), 
								paste0("OLS.", c("F", "F.P"), ".", tt[3]), 
								"OLS.resid.norm.P", 
								#PGLS
								#AIC
								"AIC.PGLS", #=NA, #AIC(opt), 
								#what is lambda and P lambda?
								"lambda", 
								"lambda.P", 
								#coefs
								paste0("PGLS.", c("coef", "se", "t", "t.P"), ".", "(Intercept)"), 
								paste0("PGLS.", c("coef", "se", "t", "t.P"), ".", tt[1]), 
								paste0("PGLS.", c("coef", "se", "t", "t.P"), ".", tt[2]), 
								paste0("PGLS.", c("coef", "se", "t", "t.P"), ".", tt[3]), 
								#F & P? 
								paste0("PGLS.", c("F", "F.P"), ".", "(Intercept)"), 
								paste0("PGLS.", c("F", "F.P"), ".", tt[1]), 
								paste0("PGLS.", c("F", "F.P"), ".", tt[2]), 
								paste0("PGLS.", c("F", "F.P"), ".", tt[3]), 
								#validation
								"PGLS.phy.resid.lambda", 
								"PGLS.phy.resid.lambda.P", 
								"PGLS.phy.resid.norm.P"
								)
								
		#establish the output vector
		ff.outsi <- rep(NA, length(out.nms))
		names(ff.outsi) <- out.nms		
	
		#add the various outputs

		#things that are always there
		ff.outsi["AIC.OLS"] <- AIC(pmr$opt.L0)
		ff.outsi["OLS.resid.norm.P"] <- pmr$resid.norm$p.value
		
		ff.outsi["AIC.PGLS"] <- AIC(pmr$opt)		
		ff.outsi[c("lambda", "lambda.P")] <- c(pmr$L.opt$maximum, pmr$L.opt.LRT)		
		pgls.nms <- paste0("PGLS.phy.resid.", c("lambda", "lambda.P", "norm.P"))
		ff.outsi[pgls.nms] <- c(pmr$L.opt.phy.resid$maximum,
														pmr$L.opt.phy.resid.LRT,
														pmr$phy.resid.norm$p.value
														)

		#things that depend on the terms included - see tTab.extract function above
		tTable.nms <- row.names(summary(pmr$opt.L0)$tTable)
		tTab.outs <- unlist(lapply(tTable.nms, tTab.extract, pmr = pmr))
		ff.outsi[names(tTab.outs)] <- tTab.outs
		
		data.frame(ff.outsi)
	}

#function to fit mutiple regressions
#providing PGLS (with lambda) and OLS (no lambda) output via pullmr
#form = a maximal formula: trait~independent variables
#tree = the phylogeny
#data = dataframe
#needs:
#tip.labels for tree == rownames for data
#pullmr()
#will evaluate all simplified submodels from 'form'
phylo.MR <- function(form, tree, data){

	#obtain names of independent variables from formula
	#and establish the sub-models to be tested
 	f1 = formula(form)
	#list of sub-models
	av <- all.vars(f1)[1]
	at <- all.vars(f1)[2:3]
	#term names
	tn <- paste(av,
							c(paste(at, collapse = "_by_"), paste(at, collapse = "_plus_"), at[2], at[1]),
							sep = "_")

	ff <- list(
							f1, 
							formula(terms(f1)[1:2]), 
							formula(terms(f1)[2]), 
							formula(terms(f1)[1])
							)
	names(ff) <- tn
	#terms from maximal model
	tt <- attr(terms(f1), "term.labels")
	
	#list for outputs with components set to NA
	ff.outs <- lapply(ff, function(.){ . <- NA })
	#ff.outs <- ff
	#for (h in 1:length(ff.outs)){
	#	ff.outs[[h]] <- NA
	#}
	
	#extract stats one model at a time
	ls.pmr <- lapply(ff, pullmr, tree = tree, data = data)
	#create a dataframe output from these
	df.pmr <- lapply(ls.pmr, pullmr.stats, tt = tt)
	
	list(terms = tt, models = ls.pmr, stats = df.pmr)
	
	}

#a function that processes a series of plant traits (dependent vars)
#in combination with a set of climate traits (indepedent vars)
#a phylogeny = tree
#and a dataframe in which the traits are held
#produces a large list object with one complete set of models per trait
#and a summary table for the full set of models
allmods <- function(pl.trts, cl.trts, tree, data){
	behemoth <- as.list(rep(NA,length(pl.trts)))
	names(behemoth) <- pl.trts
	
	for (i in 1:length(pl.trts)){
		#formula for trait i
		tmsf <- paste0(pl.trts[i], "~", paste(cl.trts, collapse = "*"))
		models <- phylo.MR(form = tmsf, tree = tree, data = data)
		behemoth[[pl.trts[i]]] <- models
		
		if (i == 1){
			tms <- do.call(cbind, models$stats)
			names(tms) <- names(models$stats)
			} else {
				tms2 <- do.call(cbind, models$stats)
				names(tms2) <- names(models$stats)
				tms <- data.frame(cbind(tms, tms2))
				}
	}

	list(models = behemoth, stats_table = tms)

}

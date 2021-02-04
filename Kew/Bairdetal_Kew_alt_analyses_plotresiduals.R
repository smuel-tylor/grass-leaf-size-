#plotting residuals for leaf area to demonstrate how our analysis accounts for assumptions
setwd(choose.dir())

load("Bairdetal_Kew_alt_analyses.Rdata")

#a function for plotting residuals from the R image 'Bairdetal_Kew_alt_analyses.Rdata'
residuals_plot <- function(object, data, trait, dep1, dep2 = NA, iact = "by", PGLS = FALSE){
	if (!is.na(dep2)) { model <- paste(trait, dep1, iact, dep2, sep = "_") }
		else { model <- paste(trait, dep1, sep = "_") }
	is.PGLS <- ifelse(PGLS, "opt", "opt.L0")
	res <- residuals(object$models[[trait]]$models[[model]][[is.PGLS]],
										"normalized"
										)
	fit <- fitted(object$models[[trait]]$models[[model]][[is.PGLS]])

	idep.str <- function(idep){
		idep <- gsub("over" , " / ", idep)
		idep <- gsub("plus" , " + ", idep)
		idep
	}
	
	#probably tidier to add units to legends
	#dep.unit <- function(dep){
	#	if (grepl("P", dep)) { unt <- "mm" }
	#	if (grepl("T", dep)) { unt <- "\u00b0 C" }
	#	unt
	#}
	
	if (!is.na(dep2)) {
		par(mfrow = c(2, 2),
				mgp = c(3, 1, 0), 
				mar = c(5, 4.5, 1.5, 2),
				cex = 1,
				cex.main = 1,
				las = 1
				)
		plot(density(res),
					main = bquote(normalized~~residuals~~log[10](.(trait)))
					)
		par(mar = c(5, 5, 1, 1),
				mgp = c(2, 0.6, 0),
				tcl = -0.3 				
				)
		plot(res ~ fit,
					ylab = bquote(atop(normalized~~residuals,log[10](.(trait)))),
					xlab = bquote(fitted~~values~~log[10](.(trait))),
					col = rgb(0, 0, 0, alpha = 0.2),
					pch = 19
					)
		plot(res ~ data[, dep1],
					ylab = bquote(atop(normalized~~residuals,log[10](.(trait)))),
					xlab = bquote(centred~~log[10](.(idep.str(dep1)))),
					col = rgb(0, 0, 0, alpha = 0.2),
					pch = 19
					)
		plot(res ~ data[, dep2],
					ylab = bquote(atop(normalized~~residuals,log[10](.(trait)))),
					xlab = bquote(centred~~log[10](.(idep.str(dep2)))),
					col = rgb(0, 0, 0, alpha = 0.2),
					pch = 19
					)

	}
		else {
		par(mfrow = c(2, 2),
				mgp = c(3, 1, 0), 
				mar = c(5, 5, 1.5, 1.5),
				cex = 1,
				cex.main = 1,
				las = 1
				)
		plot(density(res),
					main = bquote(normalized~~residuals~~log[10](.(trait)))
					)
		par(mar = c(5, 5, 1, 1),
				mgp = c(2, 0.6, 0),
				tcl = -0.3 				
				)
		plot(res ~ fit,
					ylab = bquote(atop(normalized~~residuals,log[10](.(trait)))),
					xlab = bquote(fitted~~values~~log[10](.(trait))),
					col = rgb(0, 0, 0, alpha = 0.2),
					pch = 19
					)
		plot(res ~ data[, dep1],
					ylab = bquote(atop(normalized~~residuals,log[10](.(trait)))),
					xlab = bquote(centred~~log[10](.(idep.str(dep1)))),
					col = rgb(0, 0, 0, alpha = 0.2),
					pch = 19
					)
		plot(1, 1, type = "n", bty = "n", axes = FALSE, ylab = "", xlab = "")
		}	
}

#function to insert blank page in pdf output containing a title (main is an expression or character vector)
header.page <- function(main){
	par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
	plot(1, 1, type = "n",
				bty = "n", axes = FALSE,
				xlab = "", ylab = "",
				main = main
				)
}

fig.leg <- function(number, text){
	mtext(bquote(bold(Figure~~.(number))),
				side = 1,
				line = 0,
				outer = TRUE,
				adj = 0,
				padj = 1,
				family = "ArialMT"
				)
	txt <- strwrap(text, width = 75) #par()$din * nchar(text) / strwidth(text, units = "inches"))
	txt <- paste(txt, collapse = "\n")
	mtext(txt,
				side = 1,
				line = 1,
				outer = TRUE,
				adj = 0,
				padj = 1,
				family = "ArialMT"
				)
}

png("Bairdetal_Kew_alt_analyses_plotresiduals_MA1752.png",
		width = 6, height = 8, units = "in", res = 300, type = "cairo")

#MAP MAT, 1752		
par(oma = c(10, 0, 0, 0))
residuals_plot(MA.1752, Kew.Spriggs, "LA", MA[1], MA[2], PGLS = TRUE)
fig.leg(1,
				"Density plot showing distribution of residuals (normalized, i.e. phylogenetic residuals) for a multiple regression using Pagel's \u03bb,
				and plots of normalized residuals against: fitted values for log-log scaling of leaf area (LA, cm\u00b2), and primary climate axes
				Mean annual precipitation (MAP, mm) and mean annual temperature (MAT, \u00b0C). In this 1752 taxa dataset,
				notable heteroskedasticity results from left skewed log(MAT + 20), linked with smaller residuals at low MAT"
				)
dev.off()

#MAP MAT > 0		
png("Bairdetal_Kew_alt_analyses_plotresiduals_MAMAT0.png",
		width = 6, height = 8, units = "in", res = 300, type = "cairo")
par(oma = c(10, 0, 0, 0))
residuals_plot(MA.MAT0, Kew.Spriggs.MAT0, "LA", MA[1], MA[2], PGLS = TRUE)
fig.leg(2,
				"As Figure 1, but taking only those taxa with MAT \u2265 0 (MAT + 20 \u2265 20) prior to log transformation (1723 taxa).
				Heteroskedasticity for log(MAT + 20) has been substantially reduced"
				)
dev.off()
				
#GSP GST, 1752		
png("Bairdetal_Kew_alt_analyses_plotresiduals_GS1752.png",
		width = 6, height = 8, units = "in", res = 300, type = "cairo")
par(oma = c(10, 0, 0, 0))
residuals_plot(GS.1752, Kew.Spriggs, "LA", GS[1], GS[2], PGLS = TRUE)
fig.leg(3,
				"Density plot showing distribution of residuals (normalized, i.e. phylogenetic residuals) for a multiple regression using Pagel's \u03bb,
				and plots of normalized residuals against: fitted values for log-log scaling of leaf area (LA, cm\u00b2) and growing season climate as
				mean growing season precipitation (GSP, mm) and mean growing season temperature (GST, \u00b0C).
				In this 1752 taxa dataset there is less heteroskedasticity than is observed when using MAP and MAT"
				)
dev.off()

#GSP GST, MAT > 0		
png("Bairdetal_Kew_alt_analyses_plotresiduals_GS.MAT0.png",
		width = 6, height = 8, units = "in", res = 300, type = "cairo")
par(oma = c(10, 0, 0, 0))
residuals_plot(GS.MAT0, Kew.Spriggs.MAT0, "LA", GS[1], GS[2], PGLS = TRUE)
fig.leg(4,
				"As Figure 3, but taking only those taxa with MAT \u2265 0 (MAT + 20 \u2265 20) prior to log transformation (1723 taxa).
				Removal of low MAT points has limited effects on the spread of residuals against GST"
				)
dev.off()

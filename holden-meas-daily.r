#------------------------------------------------------------------------------------------
# Holden low-frequency liquity measures model script
# For a measure of comparison, the Roll effective spread proxy is first estimated
# for comparision to the Holden and Holden2 effective spread proxies
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
# Initialization and setup
require("DEoptim")
rm(list = ls(all =T ))
set.seed(12345)

DEoptim.cntrl <- list(itermax = 150, F = 0.7, NP = 70, trace = FALSE)	

# Fractional grid cluster probabilities, from hand calculations, see excel spreadsheet 'Spread Coefficients.xls'
fpp.mat <- matrix(c(2*c(0.5,0,0,0,0,0.25,0.5,0,0,0,0.125,0.25,0.5,0,0,0.0625,0.125,0.25,0.5,0,0.0625,0.125,0.25,0.5,1),
		  1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),10,5, byrow = TRUE)

# Decimal grib cluster probabilities, from hand calculations, see excel spreadsheet 'Spread Coefficients.xls'
dpp.mat <- matrix(c(2*c(0.8,0,0,0,0,0.08,0.4,0,0,0,0.08,0.4,0.8,0,0,0.03,0.15,0.1,0.75,0,0.01,0.05,0.1,0.25,1),
		  0.8,0,0,0,0,0.16,0.8,0,0,0,0,0,1,0,0,0.04,0.2,0,1,0,0,0,0,0,1),10,5, byrow = TRUE)
	  
# Log-likelihood function for determining Holden spread estimations via MLE with DEoptim. 
clus.lik <- function(parms, llrow, hsprow)
{
	clus1 <- clus2 <- clus3 <- vector("numeric", length = nrow(llrow))
	g <- parms[1:5]
	mu <- parms[6]
	lam <- parms[7]

	if(sum(g) < 0 || sum(g) > 1) return(Inf)

	for(ii in 3:nrow(llrow))
	{
		cofs <- llrow[ii,1:15]
		norms1 <- llrow[ii,17:18]
		norms2 <- hsprow[ii,]
		
		clus1[ii] <- sum(mu*cofs[1:5]*g)
		clus2[ii] <- sum(mu*cofs[6:10]*g*dnorm( norms1[1]*norms2[1:5] - 2*(1-lam)*norms2[6:10] ))
		clus3[ii] <- sum(mu*cofs[11:15]*g*dnorm( norms1[2]*norms2[6:10] - 2*(1-lam)*norms2[11:15] ))
	}

	return( -(sum(log((clus3*clus2*clus1)[-(1:2)]))) )
}
  
#------------------------------------------------------------------------------------------
# User selects raw CRSP data file to load, processes data for calculations
sel.lab <- "- select CRSP data file -"
selected <- select.list(list.files(pattern = ".txt"), preselect = NULL, multiple = FALSE, title = sel.lab)
resultsfilename <- gsub('(.txt)|(.dat)|(.csv)', "-HOLDEN.csv", selected)

# expected format, order is ambiguous
# DATE		PERMNO	BIDLO		ASKHI		PRC			VOL		RET
# 19890103	10107	52.75000	53.75000	53.62500	408767	0.007042
CRSP <- read.table(selected, sep = '\t', header = TRUE)
colnames(CRSP) <- toupper(colnames(CRSP))
#PERMNOsb4 <- unique(CRSP$PERMNO)
cleaned <- is.na(CRSP$PRC)
CRSP <- CRSP[!cleaned, ]
PERMNOs <- unique(CRSP$PERMNO)
#stop('.......................')

cat(sprintf(" Processing file: %s \n", selected))
cat(sprintf(" Removed %d items due to missing data \n", length(cleaned)))
cat(sprintf(" Beginning Holden estimate calculations for %d PERMNO...\n", length(PERMNOs)))

#------------------------------------------------------------------------------------------
# Loop thru PERMNOs and obtain daily spread estimates 
#HOLD <- NULL 
cnt <- 1

for(p in PERMNOs)
{
	crsp <- CRSP[CRSP$PERMNO == p,]
	
	if(nrow(crsp) < 5)
	{
		cat(sprintf("Excluded PERMNO:%d not enough data to calculate Holden measure", p))
		next
	}
		
	crsp$PRC <- ifelse(crsp$PRC < 0, -crsp$PRC, crsp$PRC)
 
	cvl <- length(which(crsp$VOL == 0))
	td.spd <- crsp$PRC
	
	# prob of TD: 1 - NTD/total days
	prob.mu <- ifelse(crsp$VOL == 0, 1, 1 - cvl/nrow(crsp))
	# for the empirical spread
	#ntd.spd <- ifelse(crsp$VOL == 0, 0, sum(crsp$ASKHI-crsp$BIDLO, na.rm = TRUE)/cvl) 
	# indices of TDs
	td.ind <- ifelse(crsp$VOL == 0, which(crsp$VOL != 0), 1:length(crsp$PRC))
	
	# Define the price spread grids
	td.spd.mod <- round(td.spd%%1, 3)
	td.spd.mprc <- td.spd.mod[td.ind]
	td.spd.mmid <- td.spd.mod[-td.ind]
	
	# 1/16 Fractional price spread grid
	frac1 <- seq(1,16,2)/16
	frac2 <- seq(1,8,2)/8
	frac3 <- c(1,3)/4
	frac4 <- c(1,2)/2
	frac5 <- c(0,1)  
	spf1 <- td.spd.mprc %in% frac1
	spf2 <- td.spd.mprc %in% frac2
	spf3 <- td.spd.mprc %in% frac3
	spf4 <- td.spd.mprc %in% frac4
	spf5 <- td.spd.mprc %in% frac5
	
	# Get the midpoint probabilities
	mpf1 <- td.spd.mmid %in% (frac1*0.5)
	mpf2 <- td.spd.mmid %in% (frac2*0.5)
	mpf3 <- td.spd.mmid %in% (frac3*0.5)
	mpf4 <- td.spd.mmid %in% (frac4*0.5)
	mpf5 <- td.spd.mmid %in% (frac5*0.5)

	# Unconstrained probabilities
	fprob.vec <- c(length(which(spf1)),length(which(spf2)),length(which(spf3)),length(which(spf4)),length(which(spf5)),
                     length(which(mpf1)),length(which(mpf2)),length(which(mpf3)),length(which(mpf4)),length(which(mpf5)))
	fprob.vec <- fprob.vec/length(td.spd)
	Uf1 <- 2*fprob.vec[1] - fprob.vec[6]
	Uf2 <- 2*fprob.vec[2] - fprob.vec[1] + fprob.vec[7]
	Uf3 <- 2*fprob.vec[3] - fprob.vec[2] + fprob.vec[8]
	Uf4 <- 2*fprob.vec[4] - fprob.vec[3] + fprob.vec[9]
	Uf5 <- fprob.vec[5] - fprob.vec[4] - 0*fprob.vec[10]
	
	# Constrained probabilities
	Gf1 <- min(max(Uf1,0),1)
	Gf2 <- min(max(Uf2,0),(1-Gf1))
	Gf3 <- min(max(Uf3,0),(1-Gf1-Gf2))
	Gf4 <- min(max(Uf4,0),(1-Gf1-Gf2-Gf3))
	Gf5 <- min(max(Uf5,0),(1-Gf1-Gf2-Gf3-Gf4))

	# Decimal price spread grid
	deci1 <- seq(1,100,1)[-seq(5,100,5)]/100
	deci2 <- seq(5,100,5)/100
	deci3 <- seq(10,100,10)/100
	deci4 <- seq(25,100,25)/100
	deci5 <- c(0,1)
	spd1 <- td.spd.mprc %in% deci1
	spd2 <- td.spd.mprc %in% deci2
	spd3 <- td.spd.mprc %in% deci3
	spd4 <- td.spd.mprc %in% deci4
	spd5 <- td.spd.mprc %in% deci5
	
	# Get the midpoint probabilities
	mpd1 <- td.spd.mmid %in% (deci1*0.5)
	mpd2 <- td.spd.mmid %in% (deci2*0.5)
	mpd3 <- td.spd.mmid %in% (deci3*0.5)           
	mpd4 <- td.spd.mmid %in% (deci4*0.5)
	mpd5 <- td.spd.mmid %in% (deci5*0.5)

	# Unconstrained probabilities
	dprob.vec <- c(length(which(spd1)),length(which(spd2)),length(which(spd3)),length(which(spd4)),length(which(spd5)),
                    length(which(mpd1)),length(which(mpd2)),length(which(mpd3)),length(which(mpd4)),length(which(mpd5)))
	dprob.vec <- dprob.vec/length(td.spd)
	
	Ud1 <- 1.25*dprob.vec[1] + 1.25*dprob.vec[6]
	Ud2 <- 2.5*dprob.vec[2] - 0.25*dprob.vec[1] + 1.25*dprob.vec[7] - 0.25*dprob.vec[6]
	Ud3 <- 1.25*dprob.vec[3] - 1.25*dprob.vec[2] + dprob.vec[8]
	Ud4 <- (4/3)*dprob.vec[4] - 0.25*dprob.vec[2] - 0.25*dprob.vec[3] + dprob.vec[9] - 0.25*dprob.vec[8]
	Ud5 <- dprob.vec[5] - (1/3)*dprob.vec[4] + 0*dprob.vec[10]

	# Constrained probabilities
	Gd1 <- min(max(Ud1,0),1)
	Gd2 <- min(max(Ud2,0),(1-Gd1))
	Gd3 <- min(max(Ud3,0),(1-Gd1-Gd2))
	Gd4 <- min(max(Ud4,0),(1-Gd1-Gd2-Gd3))
	Gd5 <- min(max(Ud5,0),(1-Gd1-Gd2-Gd3-Gd4))
	
	clustID <- vector("numeric",length=length(td.spd))
	clustID[which(spd1)] <- 1 ; clustID[which(spd2)] <- 2 ; clustID[which(spd3)] <- 3 
	clustID[which(spd4)] <- 4 ; clustID[which(spd5)] <- 5 ; clustID[which(mpd1)] <- 6
	clustID[which(mpd2)] <- 7 ; clustID[which(mpd3)] <- 8 ; clustID[which(mpd4)] <- 9
	clustID[which(mpd5)] <- 10
	clustID[which(clustID == 0)] <- sample(seq(1,5), length(which(clustID == 0)), replace = TRUE)

	# -------------------------------------------------------
	# For fractional price grid    
    if(dprob.vec[1] == 0)
	{
		hf5 <- c(1,0,0,0,0, 2,1,0,0,0, 4,2,1,0,0, 8,4,2,1,0, 16,8,4,2,1, rep(0,25))/32
		hfs <- matrix(hf5, ncol = 5, nrow = 10, byrow = TRUE)
		plf <- matrix(0, nrow = length(td.spd), ncol = 6)
		
		for(iii in 3:length(td.spd))
		{
			plf[iii,] <- c(clustID[iii-2], clustID[iii-1], clustID[iii], 0, 
				crsp$PRC[iii-1] - crsp$PRC[iii-2],
				crsp$PRC[iii] - crsp$PRC[iii-1])
			
			if(sum(is.na(plf[iii,])) > 0)
				plf[iii, which(is.na(plf[iii,]))] <- 5
		}
		
		llrow <- matrix(0, nrow = length(td.spd), ncol = (3*5 + 3))
		hsprow <- matrix(0, nrow = length(td.spd), ncol = 3*5)
		for(jj in 3:length(td.spd))
		{ 
			llrow[jj,] <- c(fpp.mat[plf[jj,1],], fpp.mat[plf[jj,2],], fpp.mat[plf[jj,3],], plf[jj,4], plf[jj,5], plf[jj,6]) 
			hsprow[jj,] <- c(-hfs[plf[jj,1],], plf[jj,5] - hfs[plf[jj,2],], plf[jj,6] - hfs[plf[jj,3],])
		}   
		
		# nonlinear optimizer for the MLE estimates of the gammas, mu, and lambda
		lower <- rep(0,7) ; upper <- rep(1,7) 
		optim.out <- DEoptim(clus.lik, llrow = llrow, hsprow = hsprow, lower, upper, control = DEoptim.cntrl) 
				
		#EFFTIC[k] <- 100*sum(0.0625*Gf1, 0.125*Gf2, 0.25*Gf3, 0.5*Gf4, Gf5)/mean(td.spd)
		#HOLD[k] <- 100*sum(optim.out$optim$bestmem[1:5]*c(0.0625,0.125,0.25,0.5,1))/mean(td.spd)
		#EFFTIC.D <- 100*sum(0.0625*Gf1, 0.125*Gf2, 0.25*Gf3, 0.5*Gf4, Gf5)/td.spd
		hold <- 100*sum(optim.out$optim$bestmem[1:5]*c(0.0625,0.125,0.25,0.5,1))/td.spd
		
	}
	# end of fractional price grid
	# -------------------------------------------------------
		
	# -------------------------------------------------------
	# For Decimal Price grid        
    if(dprob.vec[1] != 0)
	{
		hd5 <- c(1,0,0,0,0, 5,1,0,0,0, 10,5,1,0,0, 25,10,5,1,0, 100,25,10,5,1, 5,1,0,0,0, 1,0,0,0,0, rep(0,15))/200
		hds <- matrix(hd5, ncol = 5, nrow = 10, byrow = TRUE)
		pld <- matrix(0, nrow = length(td.spd), ncol = 6)
      
		for(iii in 3:length(td.spd))
		{
			pld[iii,] <- c(clustID[iii-2], clustID[iii-1], clustID[iii], 0,
				crsp$PRC[iii-1] - crsp$PRC[iii-2],
				crsp$PRC[iii] - crsp$PRC[iii-1])
			
			if(sum(is.na(pld[iii,])) > 0)
				pld[iii, which(is.na(pld[iii,]))] <- 5
		}
	  
		llrow <- matrix(0,nrow = length(td.spd), ncol = (3*5 + 3))
		hsprow <- matrix(0,nrow = length(td.spd), ncol = 3*5)
		for(jj in 3:length(td.spd))
		{ 
			llrow[jj,] <- c(dpp.mat[pld[jj,1],], dpp.mat[pld[jj,2],], dpp.mat[pld[jj,3],], pld[jj,4], pld[jj,5], pld[jj,6])
			hsprow[jj,] <- c(-hds[pld[jj,1],] ,pld[jj,5] - hds[pld[jj,2],], pld[jj,6] - hds[pld[jj,3],]) 
		}
		
        # nonlinear optimizer for the MLE estimates of the gammas, mu, and lambda
		lower <- rep(0,7) ; upper <- rep(1,7) 
		optim.out <- DEoptim(clus.lik, llrow = llrow, hsprow = hsprow, lower, upper, control = DEoptim.cntrl) 
	  
		#EFFTIC[k] <- 100*sum(0.01*Gd1, 0.05*Gd2, 0.1*Gd3, 0.25*Gd4, Gd5)/mean(td.spd)
		#HOLD[k] <- 100*sum(optim.out$optim$bestmem[1:5]*c(0.01,0.05,0.1,0.25,1))/mean(td.spd)
		#EFFTIC.D <- 100*sum(0.01*Gd1, 0.05*Gd2, 0.1*Gd3, 0.25*Gd4, Gd5)/td.spd
		hold <- 100*sum(optim.out$optim$bestmem[1:5]*c(0.01,0.05,0.1,0.25,1))/td.spd  
	  
    }
	# end of decimal price grid
	# -------------------------------------------------------

	# Now the Holden2 and Effective-Tick2 estimates
    #EFFTIC2[k] <- prob.mu*EFFTIC[k] + (1-prob.mu)*ntd.spd/mean(td.spd)    
    #HOLD2[k] <- optim.out$optim$bestmem[6]*HOLD[k] + (1-optim.out$optim$bestmem[6])*ntd.spd/mean(td.spd)
    #EFFTIC2.D <- prob.mu*EFFTIC.D + (1-prob.mu)*ntd.spd/td.spd    
    #HOLD2.D <- optim.out$optim$bestmem[6]*HOLD.D + (1-optim.out$optim$bestmem[6])*ntd.spd/td.spd
	
	#avgsprd <- mean(abs(crsp$ASKHI-crsp$BIDLO))
	#avgholden <- mean(hold)
	#cat(sprintf(" Progress [prcnt]: %2.2f| avg sprd:%3.4f avg holden:%3.4f\n", 100*cnt/length(PERMNOs), avgsprd, avgholden)) 
	
	write.table(as.data.frame(cbind(DATE = crsp$DATE, PERMNO = p, HOLDEN = hold)), file = resultsfilename, 
		append = TRUE, quote = FALSE, sep = ',', row.names = FALSE)
	cat(sprintf(" Progress [prcnt]: %2.2f, results incrementially saved to %s\n", 100*cnt/length(PERMNOs), resultsfilename))
	cnt <- cnt + 1
	#HOLD <- rbind(HOLD, cbind(DATE = crsp$DATE, PERMNO = p, HOLDEN = hold))
}

cat(sprintf(" Finished, saved results to: %s \n", sname))
#' optM function
#'
#' Load a folder of .llik files from the program Treemix and determine the optimal number of migration edges to include
#' @param folder A character string of the path to a directory containing .llik, .cov.gz and .modelcov.gz files produced by Treemix
#' @param tsv a string defining the name of the tab-delimited output file.
#' If NULL (default), then no data file is produced.
#' @param method a string containing the method to use, either "Evanno", "linear", or "SiZer".  Default is "Evanno".
#' @param thresh a numeric value between 0 and 1 for the threshold to use for the proportion of increase
#' in likelihood that defines when a plateau is reached.  Default is 0.05 (5\%), only applicable for method = "linear".
#' @param ... other options sent to the function "SiZer" - see the R package 'SiZer'
#' @keywords optM
#' @return If method = "Evanno": A data frame with 17 columns summarizing the results for each migration edge (rows).
#' @return The columns are: "m" - number of migration edges from the model; "runs" = number of iterations for "m";
#' "mean(Lm)" - mean log likelihood across runs; "sd(Lm)" - standard deviation of log likelihood across runs;
#' "min(Lm)" - minimum log likelihood across runs; "max(Lm)" - maximum log likelihood across runs;
#' "L'(m)" - first-order rate of change in log likelihood; "sdL'(m)" - standard deviation of first-order rate of change in log likelihood;
#' "minL'(m)" - minimum first-order rate of change in log likelihood; "maxL'(m)" - maximum first-order rate of change in log likelihood;
#' "L''(m)" - second-order rate of change in log likelihood; "sdL''(m)" - standard deviation of the second-order rate of change in log likelihood;
#' "minL''(m)" - minimum second-order rate of change in log likelihood; "maxL''(m)" - maximum second-order rate of change in log likelihood;
#' "Deltam" - the ad hoc deltaM statistic (secord order rate of change in log likelihood);
#' "mean(f)" - mean proportion of variation explained by the models; "sd(f)" - standard deviation of the proportion of variation explained by the models
#' @return If method = "linear": A list containing 5 elements:
#' @return $out - a data frame with the name of each model, the degrees of freedom (df), the Akaike information criterion (AIC), the deltaAIC, and the optimal estimate for m based on the model.
#' @return $PiecewiseLinear - the piecewise linear model object
#' @return $BentCable - the bent cable model object
#' @return $SimpleExponential - the simple exponential model object
#' @return $NonLinearLeastSquares - the NLS model object
#' @return If method = "SiZer": an object of class "SiZer" (see the R package 'SiZer' for more information)
#' @import SiZer
#' @export
#' @examples
#' # Load a folder of simulated test data for m = 3
#' folder <- system.file("extdata", package = "OptM")
#' test.optM = optM(folder)
#'
#' # To view the various linear modeling estimates:
#'    # test.linear = optM(folder, method = "linear")
#'
#' # To view the results from the SiZer package:
#'    # test.sizer = optM(folder, method = "SiZer")

optM <- function(folder, tsv = NULL, method = "Evanno", thresh = 0.05, ...){
 
    # ... are option to pass to the function 'SiZer'

	# Load .llik input files
	tbl = read.treemix(folder)
	#Sys.sleep(2)
	
	# Setup output file name	
		if (is.null(tsv)){
		message("No output file will be saved. To save an output file, run with 'tsv = \"file.tsv\"'\n")
	} else {
		if (length(tsv) != 1 | !is.character(tsv)) stop("Output tsv file incorrectly specified.\n")
		tsv = tsv
	}
	
	# Check for correct 'method' formatting
	if (missing(method)) 
        method = "Evanno"
    else method = method
    if (!is.character(method) | length(method) > 1) stop("The 'method' argument was not set correctly\n")
    methods = c("SiZer", "linear", "Evanno") 
    if(!(method %in% methods)) stop("Could not find the selected 'method'.  Please check.\n")
    
    # Check for correct setting of the 'thresh' parameter
       # 'Percent increase in log likelihood less than thresh is the optimal number of M for exponential models
    if(missing(thresh) & method == "linear"){
    	thresh = 0.05
    	message("'thresh' has been set to the default value of 0.05\n")
    } else if(!is.numeric(thresh) | length(thresh) > 1 | thresh > 1 | thresh < 0) stop("The 'thresh' parameter has not been set correctly to a number between 0 and 1.  Please check.\n")
    else thresh = thresh
   
	# Reduce table to two columns: m and lnpd
	M <- NULL  # avoids variable scope error with 'subset' and CRAN submission
	LnPD <- NULL  # avoids variable scope error with 'subset' and CRAN submission
	tbl = tbl[,c(4,7)]
	if (any(is.na(tbl))) stop("Error: One or more likelihoods are \"NaN\", please check datafiles and/or repeat the treemix run!\n")
	colnames(tbl) <- c("M", "LnPD")
	m = max(tbl[,1], na.rm = T)
	low = min(tbl[,1], na.rm = T)
	
	# Get number of iterations per value of m
	runs = vector()
	for (i in 1:m){
		runs = c(runs,nrow(subset(tbl, M == i)))
	}
		mean.runs = mean(runs, na.rm = T)
		if(ceiling(mean.runs) != mean.runs) {
		warning(paste0("The mean number of runs detected is ", mean.runs, ", but this must be a whole number.  Rounding up and continuing anyway.\n"))
		# Sys.sleep(3)
		mean.runs = round(mean.runs)
	}
	
	message(paste("m ranges between ", low, " and ", m, ".\n", sep = ""))
	#Sys.sleep(2)
	message(paste("The average number of iterations per m (not including m = 0) was ", mean(runs, na.rm = T), ".\n", sep = ""))
	#Sys.sleep(2)
	message("Make sure these values are correct...\n")
	#Sys.sleep(2)
	
	# Check the cov.gz and modelcov.gz files
	cov.files = list.files(path = folder, pattern = "\\.cov.gz", full.names = T)
	modelcov.files = list.files(path = folder, pattern = "\\.modelcov.gz", full.names = T)
	if(length(cov.files) != length(modelcov.files)) stop("Error: Could not find the same number of .cov.gz and .modelcov.gz files.\n")
	cov.files = sort(unlist(lapply(cov.files, find.cov)))
	modelcov.files =  sort(unlist(lapply(modelcov.files, find.modelcov)))
	if(!all(cov.files == modelcov.files)) stop("Error:  The .cov.gz and .modelcov.gz files do not all match up properly")
	stem = unique(sub("\\..*\\..*", "", c(cov.files, modelcov.files)))
	if (is.null(stem) | length(stem) != 1 | !is.character(stem)) stop("Error: The file stem name was not correctly identified. Check naming convention\n")

    # Get the variation explained by each model
	var.expl = data.frame()
	for (i in 1:m){
		for (n in 1:runs[m]){
		stem1 = paste(folder, "/", stem, ".",n, ".",i, sep = "")
		if (!file.exists(paste0(stem1, ".cov.gz"))) next
		if (!file.exists(paste0(stem1, ".modelcov.gz"))) next
		var.expl = rbind(var.expl,c(n, i, get_f(stem1)))
		}
	}
	colnames(var.expl) = c("run", "m", "f")
	
	# Get the mean and standard deviation of % variation explained across runs per m
	f = c(NA, stats::aggregate(var.expl, by = list(var.expl$m), mean, na.rm = T)$f)
	sdf = c(NA, stats::aggregate(var.expl, by = list(var.expl$m), stats::sd, na.rm = T)$f)

	
	######################  Now separate based on method chosen  #############################################
	
	if (method == "SiZer"){
	   message(paste0("Analyzing the treemix results using the ", method, " method.\n"))
	   out = SiZer::SiZer(tbl$M, tbl$LnPD, ...)
	}

	else if (method == "linear"){
	   message("Analyzing the treemix results using various linear models.\n")
	   x = tbl$M
       y = tbl$LnPD
       data = data.frame(x = x, y = y)
       x.grid <- seq(min(x), max(x), length = 200)
       
	   # Fit piecewise linear
	   	  message("Fitting piecewise linear model...\n")
	      pl = SiZer::piecewise.linear(x, y, middle = 1, CI = TRUE, bootstrap.samples = 1000, sig.level = 0.05)
	      
	   # Fit bent cable model
	   	  message("Fitting bent cable model...\n")
	      bc = SiZer::bent.cable(x, y, grid.size = 100)
	      
	   # Fit simple exponential
	   	  message("Fitting a simple exponential model...\n")
	   	  z = max(y) + 1
          y2 = log(-y + z)
          sim.exp = stats::lm(y2 ~ x)
          lnPD.exp = z - exp(stats::predict(sim.exp, newdata = data.frame(x = c(low:m))))
          lnPD.exp.diff = c(0, diff(lnPD.exp)/lnPD.exp[-length(lnPD.exp)])
          opt.sim.exp = (low:m)[which(lnPD.exp.diff[-1] < thresh)[1] + 1]
          
          ####### Old code to find radius of curvature  #######
          # Get first and second derivative
          # m = coef(sim.exp)[2]
          # b = coef(sim.exp)[1]
          # dv = D(expression(z-exp(m*x.grid+b)), 'x.grid')
          # slopes = eval(dv)
          # dv2 = D(D(expression(z-exp(m*x.grid+b)), 'x.grid'), 'x.grid')
          # slopes2 = eval(dv2)
          # Make a function for the derivatives
          # fdv1 = function(x, m, b) return(-(exp(m * x + b) * m))
          # fdv2 = function(x, m, b) return(-(exp(m * x + b) * m * m))
          # Radius of curvature K is defined as:
          # R = abs(((1 + (fdv1(x.grid, m, b)^2))^(3/2)) / fdv2(x.grid, m, b))
          
       # Fit non-linear least squares (NLS) model
	   	  message("Fitting a NLS model...\n")
          nl_ls = stats::nls(y2 ~ I(exp(1)^(a + b * x)), data = data.frame(x = x, y2 = y2), start = list(a = 0, b = 1), trace = T)
          message("Done fitting a NLS model...\n")
          lnPD.nls = z - exp(stats::predict(nl_ls, newdata = data.frame(x = c(low:m))))
          lnPD.nls.diff = c(0, diff(lnPD.nls)/lnPD.nls[-length(lnPD.nls)])
          opt.nls = (low:m)[which(lnPD.nls.diff[-1] < thresh)[1] + 1]
          
       # Fit Michaelis-Menten equation
	   	  # message("Fitting the Michaelis-Menten model...\n")
          # mme = nls(y ~ (Vm * x) / (K + x), data = data, start = list(K = max(data$y)/2, Vm = max(data$y)), trace = T, control = nls.control(warnOnly = TRUE, minFactor = 1/100000))
          # message("Done fitting the Michaelis-Menten model...\n")
       
       # Get AIC and deltaAIC
       aic = stats::AIC(pl, bc, sim.exp, nl_ls)
       delta.aic = aic$AIC - min(aic$AIC)
       out = cbind(aic, delta.aic, ChangePoints = c(pl$change.point, bc$alpha, opt.sim.exp, opt.nls))
       rownames(out) <- c("PiecewiseLinear", "BentCable", "SimpleExponential", "NonLinearLeastSquares")
       out = out[order(out$delta.aic),]
       
       # Combine output estimates and models
       out = list(out = out, PiecewiseLinear = pl, BentCable = bc, SimpleExponential = sim.exp, NonLinearLeastSquares = nl_ls)
	}

	else if (method == "Evanno"){
	   message(paste0("Analyzing the treemix results using the ", method, " method.\n"))
	   
	   # Get mean and sd across runs for each m
	   data.summary = data.frame()
	   for (i in low:m){
		   probs = subset(tbl, M == i, select = LnPD)
		   probs = as.vector(probs[,1])
		   tmp = c(i, length(probs), mean(probs, na.rm = T), stats::sd(probs, na.rm = T), min(probs, na.rm = T), max(probs, na.rm = T))
		   data.summary = rbind(data.summary, tmp)
	   }
	   colnames(data.summary)<-c("m", "runs", "meanLm", "sdLm", "minLm", "maxLm")
	   data.summary <- as.data.frame(data.summary, stringsAsFactors = F)
	
	   # Check data for errors
	   if (length(data.summary$m) < 3) stop("Error: The number of migration edges, m, must be > 3.\n")
	   diffs = abs(diff(data.summary$m))
	   if (all(diffs != 1) == TRUE) stop("Error: The method requires sequential values of m.\n")
	   if (any(data.summary$runs < 2)) stop("Error: It is recommended to run more than 2 iterations for each m.\n")
       if (all(data.summary$sdLm == 0)) {
    	   #Sys.sleep(3)
    	   stop("SD is zero for all runs!  Check your analysis.\n", immediate.=T)
       }
       if (any(data.summary$sdLm == 0) && !all(data.summary$sdLm == 0)) {
          #Sys.sleep(3)
          stop("SD is zero for 1 or more runs!  Check your analysis.\n", immediate.=T)
          #data.summary$sdLm[which(data.summary$sdLm == 0)] <- min(data.summary$sdLm[which(data.summary$sdLm != 0)], na.rm = T)
       }
	
	   # Convert data to a list
	   data.list = as.list(data.summary)
	
	   # Get l'(m)
	   l1m = vector()
	   l1msd = vector()
	   n = length(data.list$meanLm) - 1
	   for (i in 1:n) {
          l1m[i] = data.list$meanLm[i+1] - data.list$meanLm[i]
          l1msd[i] = abs(data.list$sdLm[i+1] - data.list$sdLm[i])
	   }
	   l1m = as.numeric(l1m)
	   l1m = c(NA, l1m)
	   l1msd = as.numeric(l1msd)
	   l1msd = c(NA, l1msd)
	
	   # Get l''(m)
	   l2m = vector()
	   l2msd = vector()
	   n = length(l1m)-1
	   for (i in 2:n){
	      l2m[i-1] = abs(l1m[i+1] - l1m[i])
          l2msd[i-1] = abs(l1msd[i+1] - l1msd[i])
	   }
	   l2m = as.numeric(l2m)
	   l2m = c(NA, l2m, NA)
	   l2msd = as.numeric(l2msd)
	   l2msd = c(NA, l2msd, NA)
	
	   # Get the min/max for each derivative
	   l1m.max = l1m + l1msd
	   l1m.min = l1m - l1msd
	   l2m.max = l2m + l2msd
	   l2m.min = l2m - l2msd
	
	   # Get delta(m)
	   delta.m = abs(l2m)/data.list$sdLm
	
	   message("Finished calculating delta m.\n")
	   message(paste("The maximum value for delta m was ", 
		   round(max(delta.m, na.rm  = T), 4)," at m = ", which.max(delta.m) - 1, 
		   " edges.\n", sep = ""))
	   #Sys.sleep(2)
	
	   # Create output table
	   out = cbind(data.summary, l1m, l1msd, l1m.min, l1m.max, l2m, l2msd, l2m.min, l2m.max, delta.m, f, sdf)
	   colnames(out) = c("m", "runs", "mean(Lm)", "sd(Lm)", "min(Lm)", "max(Lm)", "L'(m)", "sdL'(m)", "minL'(m)", "maxL'(m)", "L''(m)", "sdL''(m)", "minL''(m)", "maxL''(m)", "Deltam", "mean(f)", "sd(f)")
	   out = as.data.frame(out)
	   if (!is.null(tsv)){
	   utils::write.table(out, file = tsv, quote = F, row.names = F, sep = "\t", dec = ".")
	   message(paste("Output table ", tsv, " was written to the current directory.\n", sep = ""))
	   }

	   message("Finished all calculations for the Evanno method.\n\n\n")
	   #Sys.sleep(2)
   }
	   return(out)
}

###############################
####### Other Functions #######
###############################

# Read in treemix .llik files
read.treemix <- function(folder){
	# Check input parameters
	if (is.null(folder) | length(folder) != 1 | !is.character(folder)) stop("No input folder correctly provided.\n")	

	# Gather output files into table
	fileList = list.files(path = folder, pattern = ".llik", full.names = T)
	tbl = data.frame()

	for (i in 1:length(fileList)){
		info = file.info(fileList[i])
		if (info$size == 0) stop("At least one of the .llik files is empty. Check results\n")
		tbl = rbind(tbl, utils::read.table(fileList[i], header = F, sep = " "))
	}
	
	if (length(tbl) == 0){
		stop("Could not properly load .llik files\n")
		} else {
			message("Finished reading .llik (likelihood) files.\n")
		}
	
	return(tbl)
}

# Read treemix cov.gz and modelcov.gz files
find.cov = function(x){
	z = sub("\\.cov.gz", "", utils::tail(unlist(strsplit(x, "/")), n = 1))
	return(z)
}
find.modelcov = function(x){
	z = sub("\\.modelcov.gz", "", utils::tail(unlist(strsplit(x, "/")), n = 1))
	return(z)
}


# get_f function
# Calculate proportion of explained variance from Treemix.
# This function is slightly modified from that available currently as part of
# the "plotting_funcs.R" available in the Treemix package to calculate the
# % variation explained by the model with m migration edges to a model without
# migration.  This code can be orignally attributed to J. Pickrell.
# Input: stem - a character string of the path to a directory containing stem.cov.gz and stem.modelcov.gz files produced by Treemix
# Returns: the proportion of explained variance
get_f = function(stem){
	d = paste(stem, ".cov.gz", sep = "")
	if (file.exists(d) == FALSE) stop("Cannot find the .cov.gz file. Check naming convention or path.\n")
	d2 = paste(stem, ".modelcov.gz", sep = "")
	if (file.exists(d2) == FALSE) stop("Cannot find the .modelcov.gz file. Check naming convention or path.\n")
	d = utils::read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
	d2 = utils::read.table(gzfile(d2), as.is = T, comment.char = "", quote = "")
	d = d[order(names(d)), order(names(d))]
	d2 = d2[order(names(d2)), order(names(d2))]
	tmpcf = vector()
        tmpmcf = vector()
        for (j in 1:nrow(d)){
                for (k in (j+1):nrow(d)){
                        tmpcf = append(tmpcf, d[j,k])
                        tmpmcf = append(tmpmcf, d[j,k] - d2[j,k])
                }
        }
        tmpv = stats::var(tmpmcf)/stats::var(tmpcf)
	return(1-tmpv)
}
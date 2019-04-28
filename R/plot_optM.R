#' plot_optM function
#'
#' Plotting the optM results.
#' This function visualizes the output of optM, including the amount of total variation
#' explained across each value of the migration rate
#' @param input an object produced by the fucntion 'optM'
#' @param method a string containing the method to use, either "Evanno", "linear", or "SiZer".  Default is "Evanno", but needs to match that used in 'optM'
#' @param plot logical of whether or not to display the plot
#' @param pdf string of the file name to save the resulting pdf plot.  If NULL, no file is saved.  Default is NULL
#' @keywords plot_optM
#' @return a plot or pdf of a plot
#' @export
#' @examples
#' # Load a folder of simulated test data for m = 3
#' folder <- system.file("extdata", package = "OptM")
#' # Run the Evanno method and plot the results
#' test.optM = optM(folder)
#' plot_optM(test.optM, method = "Evanno")
#'
#' # To view the various linear modeling estimates and plot:
#'    # test.linear = optM(folder, method = "linear")
#'    # plot_optM(test.linear, method = "linear")
#'
#' # To view the results from the SiZer package:
#'    # test.sizer = optM(folder, method = "SiZer")
#'    # plot_optM(test.sizer, method = "SiZer")

plot_optM <- function(input, method = "Evanno", plot = TRUE, pdf = NULL ){
	
	# Check for correct 'method' formatting
	if (missing(method)) {
		method = "Evanno"
	} else method = method
    if (!is.character(method) | length(method) > 1) stop("The 'method' argument was not set correctly\n")
    methods = c("SiZer", "linear", "Evanno") 
    if(!(method %in% methods)) stop("Could not find the selected 'method'.  Please check.\n")

	# Setup output file name	
		if (is.null(pdf)){
		message("No output file will be saved. To save an output file, run with 'pdf = \"file.pdf\"'\n")
	} else {
		if (length(pdf) != 1 | !is.character(pdf)) stop("Output pdf file is incorrectly specified.\n")
		pdf = pdf
	}

	# Do you want a plot opened up?
	if(!is.logical(plot)) stop("Please set 'plot' as either TRUE or FALSE.\n")
	if(!plot & is.null(pdf)) stop(" You want neither a plot opened or an input file saved.  No reason to continue. Exiting.\n")

	# Now separate by method
	if (method == "Evanno"){
	   message(paste0("Plotting the treemix results using the ", method, " method.\n"))
	   if (is.null(input) | !is.data.frame(input)) stop("Proper input data frame was not detected.\n")	   
	   c = ncol(input)
	   if (c != 17) warning(paste0("Warning: This function expects 17 columns in the table but detected ", c, " columns.  Proceeding anyways even if this is incorrect!\n"))
	   # Sys.sleep(2)
	   
	   # Get the number of migration edges and iterations
	   m = max(input$m, na.rm = T)
	   low = min(input$m, na.rm = T)
	   runs = mean(input$runs[2:length(input$runs)])
	   if(ceiling(runs) != runs) {
		warning(paste0("The mean number of runs detected is ", runs, ", but this must be a whole number.  Rounding up and continuing anyway.\n"))
		# Sys.sleep(2)
		runs = round(runs)
	   }
	   
	   # Plotting
	   if (!plot){
	      pdf(pdf, width = 7, height = 7)
	      graphics::par(mfrow = c(2, 1),              # 2x1 layout
   		     mar = c(4.1, 4.1, 1.1, 5.1),   # space for one row of text at ticks and to separate plots
   		     mgp = c(3, 1, 0))              # axis label at 2 rows distance, tick labels at 1 row
   	
		  # Draw plot of likelihoods with SD bars
		  plot(input$m, input$'mean(Lm)', pch = 1, axes = F, ann = F)
		  graphics::axis(2, las = 1)
		  graphics::axis(1)
		  graphics::box()
		  graphics::segments(input$m, input$'mean(Lm)' - input$'sd(Lm)', input$m,
		     input$'mean(Lm)' + input$'sd(Lm)')
		  epsilon = 0.1
		  graphics::segments(input$m - epsilon, input$'mean(Lm)' - input$'sd(Lm)',
		     input$m + epsilon, input$'mean(Lm)' - input$'sd(Lm)')
		  graphics::segments(input$m - epsilon, input$'mean(Lm)' + input$'sd(Lm)',
		     input$m + epsilon, input$'mean(Lm)' + input$'sd(Lm)')
		  graphics::title(ylab = "Mean L(m) +/- SD")
	
		  # Calculate % variation explained means/SD
		  f.means = input$'mean(f)'
	 	  f.sd = input$'sd(f)'
	
		  # Plot variation explained on second Y axis
		  graphics::par(new = T)
		  # plot(input$m, f.means, pch = 1, col = "red", axes = F, ann = F) # v0.1.1
		  plot(input$m, f.means, pch = 19, col=grDevices::rgb(255/255,0,0,89.25/255), axes = F, ann = F)
		  graphics::axis(4, las = 1)
		  graphics::mtext("Variance Explained", side = 4, line = 3.5)
	
    	  # Get Y scale to see if the horizontal cutoff line will fit
    	  y.limits = graphics::par("usr")[3:4]
    	  if ((y.limits[1]) < 0.998 && (y.limits[2] > 0.998)){
		     graphics::abline(h = 0.998, col = "black", lty = "dotted")
    	  } else { 
    	     warning("Horizontal line at 99.8% variation cutoff is out of bounds. This is not a big deal and the program is continuing anyway without plotting the line.\n", immediate.=T)
    		 # Sys.sleep(2)
    	  }
		  graphics::segments(input$m, f.means - f.sd, input$m,
		     f.means + f.sd, col = "red")
		  epsilon = 0.1
		  graphics::segments(input$m - epsilon, f.means - f.sd,
		     input$m + epsilon, f.means - f.sd, col = "red")
		  graphics::segments(input$m - epsilon, f.means + f.sd,
		     input$m + epsilon, f.means + f.sd, col = "red")
		  graphics::legend("bottomright", legend = c("likelihoods +/- SD", "% variance", "99.8% threshold"), col = c("black", grDevices::rgb(255/255,0,0,89.25/255), "black"), bty = "n", pch = c(1, 19, NA), lty = c(NA, NA, "dotted"))
		
		  # Draw second plot of delta m
		  # plot(input$m, input$Deltam, col = "blue", pch = 19, xlab = "m", ylab = "Delta m") #v0.1.1
		  plot(input$m, input$Deltam, col = "blue", pch = 19, xlab = "m", ylab = expression(italic(paste(symbol(Delta),"m"))))
		  graphics::points(input$m, input$Deltam, col = "blue", type = "l")
		  grDevices::dev.off()
	   } else {
	  grDevices::dev.new(width = 7, height = 7)
	  graphics::par(mfrow = c(2, 1),              # 2x1 layout
   	     mar = c(4.1, 4.1, 1.1, 5.1),   # space for one row of text at ticks and to separate plots
  	     mgp = c(3, 1, 0))              # axis label at 2 rows distance, tick labels at 1 row
   	
	  # Draw plot of likelihoods with SD bars
	  plot(input$m, input$'mean(Lm)', pch = 1, axes = F, ann = F)
	  graphics::axis(2, las = 1)
	  graphics::axis(1)
	  graphics::box()
	  graphics::segments(input$m, input$'mean(Lm)' - input$'sd(Lm)', input$m,
		input$'mean(Lm)' + input$'sd(Lm)')
	  epsilon = 0.1
	  graphics::segments(input$m - epsilon, input$'mean(Lm)' - input$'sd(Lm)',
		input$m + epsilon, input$'mean(Lm)' - input$'sd(Lm)')
	  graphics::segments(input$m - epsilon, input$'mean(Lm)' + input$'sd(Lm)',
		input$m + epsilon, input$'mean(Lm)' + input$'sd(Lm)')
	  graphics::title(ylab = "Mean L(m) +/- SD")
	
	  # Calculate % variation explained means/SD
	  f.means = input$'mean(f)'
	  f.sd = input$'sd(f)'
	
	  # Plot variation explained on second Y axis
	  graphics::par(new = T)
	  # plot(input$m, f.means, pch = 1, col = "red", axes = F, ann = F) #v0.1.1
	  plot(input$m, f.means, pch = 19, col = grDevices::rgb(255/255,0,0,89.25/255), axes = F, ann = F)
	  graphics::axis(4, las = 1)
	  graphics::mtext("Variance Explained", side = 4, line = 3.5)
	
      # Get Y scale to see if the horizontal cutoff line will fit
      y.limits = graphics::par("usr")[3:4]
      if ((y.limits[1]) < 0.998 && (y.limits[2] > 0.998)){
	  graphics::abline(h = 0.998, col = "black", lty = "dotted")
      } else {
    	warning("Horizontal line at 99.8% variation cutoff is out of bounds. This is not a big deal and the program is continuing anyway without plotting the line.\n", immediate.=T)
    	# Sys.sleep(2)
      }
	  graphics::segments(input$m, f.means - f.sd, input$m,
		f.means + f.sd, col = "red")
	  epsilon = 0.1
	  graphics::segments(input$m - epsilon, f.means - f.sd,
		input$m + epsilon, f.means - f.sd, col = "red")
	  graphics::segments(input$m - epsilon, f.means + f.sd,
		input$m + epsilon, f.means + f.sd, col = "red")
	  graphics::legend("bottomright", legend = c("likelihoods +/- SD", "% variance", "99.8% threshold"), col = c("black", grDevices::rgb(255/255,0,0,89.25/255), "black"), bty = "n", pch = c(1, 19, NA), lty = c(NA, NA, "dotted"))
		
	  # Draw second plot of delta m
	  plot(input$m, input$Deltam, col = "blue", pch = 19, xlab = "m (migration edges)", ylab = expression(italic(paste(symbol(Delta),"m"))))
	  graphics::points(input$m, input$Deltam, col = "blue", type = "l")
	  
	  if(!is.null(pdf)){
	     grDevices::dev.copy(grDevices::pdf, file = pdf)
	     grDevices::dev.off()
	     message(paste0("Plot saved to file ", pdf, ".\n"))
	  } else if(is.null(pdf)) message("No plot file has been saved.\n")
  	}

   } else if (method == "linear"){
	   message("Plotting the treemix results using various linear models.\n")
	   if (is.null(input) | !is.list(input) | length(input) != 5) stop("Proper input list was not detected.\n")	   
	   # Sys.sleep(2)
	   if(plot){
	      x = input$PiecewiseLinear$model$model[,2]
	      y = input$PiecewiseLinear$model$model[,1]
	      pl = input$PiecewiseLinear
	      bc = input$BentCable
	      sim.exp = input$SimpleExponential
	      nl_ls = input$NonLinearLeastSquares
	      plot(x,y, ylab = "Log Likelihood", xlab = "m (migration edges)", axes = F)
	      graphics::box()
	      graphics::axis(1)
	      graphics::axis(2, las = 1)
	      x.grid <- seq(min(x), max(x), length = 200)
	      graphics::lines(x.grid, stats::predict(pl, x.grid), col='darkgreen', lwd = 2)
	      graphics::lines(x.grid, stats::predict(bc, x.grid), col='orange', lwd = 2)
	      z = max(y) + 1
	      y2 = log(-y + z)
	      graphics::lines(x.grid, z - exp(stats::predict(sim.exp, newdata = data.frame(x = x.grid))), col = 'red', lwd = 2)
	      graphics::lines(x.grid, z - exp(stats::predict(nl_ls, newdata = data.frame(x = x.grid))), col='blue', lwd = 2)
		  graphics::legend("bottomright",
		     legend = c("Observed data", "Piecewise Linear", "Bent Cable", "Simple Exponential", "Non-linear Least Squares", "change points"),
		     col = c("black", "darkgreen", "orange", "red", "blue", "black"),
		     bty = "y",
		     pch = c(1, NA, NA, NA, NA, 8),
		     lty = c(NA, 1, 1, 1, 1, NA))
	      # Plot change points
	      cp.pl = input$out[which(rownames(input$out) == "PiecewiseLinear"),4]
	      lnPD.pl = stats::predict(pl, cp.pl)
	      cp.bc = input$out[which(rownames(input$out) == "BentCable"),4]
	      lnPD.bc = stats::predict(bc, cp.bc)
	      cp.simexp = input$out[which(rownames(input$out) == "SimpleExponential"),4]
	      lnPD.simexp = z - exp(stats::predict(sim.exp, newdata = data.frame(x = cp.simexp)))
	      cp.nlls = input$out[which(rownames(input$out) == "NonLinearLeastSquares"),4]
	      lnPD.nlls = z - exp(stats::predict(nl_ls, newdata = data.frame(x = cp.nlls)))
	      graphics::points(x = c(cp.pl, cp.bc, cp.simexp, cp.nlls), y = c(lnPD.pl, lnPD.bc, lnPD.simexp, lnPD.nlls), pch = 8, col = c("darkgreen", "orange", "red", "blue"), cex = 1.5)
	      if(!is.null(pdf)){
	      	grDevices::dev.copy(grDevices::pdf, file = pdf)
	      	grDevices::dev.off()
	      	message(paste0("Plot saved to file ", pdf, ".\n"))
	      } else if(is.null(pdf)) message("No plot file has been saved.\n")
	      
	   } else if (!plot){
	   	  x = input$PiecewiseLinear$model$model[,2]
	      y = input$PiecewiseLinear$model$model[,1]
	      pl = input$PiecewiseLinear
	      bc = input$BentCable
	      sim.exp = input$SimpleExponential
	      nl_ls = input$NonLinearLeastSquares
	      # Build plot and save
	      pdf(pdf, width = 7, height = 7)
	      plot(x,y, ylab = "Log Likelihood", xlab = "m (migration edges)", axes = F)
	      graphics::box()
	      graphics::axis(1)
	      graphics::axis(2, las = 1)
	      x.grid <- seq(min(x), max(x), length = 200)
	      graphics::lines(x.grid, stats::predict(pl, x.grid), col='darkgreen', lwd = 2)
	      graphics::lines(x.grid, stats::predict(bc, x.grid), col='orange', lwd = 2)
	      z = max(y) + 1
	      y2 = log(-y + z)
	      graphics::lines(x.grid, z - exp(stats::predict(sim.exp, newdata = data.frame(x = x.grid))), col = 'red', lwd = 2)
	      graphics::lines(x.grid, z - exp(stats::predict(nl_ls, newdata = data.frame(x = x.grid))), col='blue', lwd = 2)
		  graphics::legend("bottomright",
		     legend = c("Observed data", "Piecewise Linear", "Bent Cable", "Simple Exponential", "Non-linear Least Squares", "change points"),
		     col = c("black", "darkgreen", "orange", "red", "blue", "black"),
		     bty = "y",
		     pch = c(1, NA, NA, NA, NA, "X"),
		     lty = c(NA, 1, 1, 1, 1, NA))
	      # Plot change points
	      cp.pl = input$out[which(rownames(input$out) == "PiecewiseLinear"),4]
	      lnPD.pl = stats::predict(pl, cp.pl)
	      cp.bc = input$out[which(rownames(input$out) == "BentCable"),4]
	      lnPD.bc = stats::predict(bc, cp.bc)
	      cp.simexp = input$out[which(rownames(input$out) == "SimpleExponential"),4]
	      lnPD.simexp = z - exp(stats::predict(sim.exp, newdata = data.frame(x = cp.simexp)))
	      cp.nlls = input$out[which(rownames(input$out) == "NonLinearLeastSquares"),4]
	      lnPD.nlls = z - exp(stats::predict(nl_ls, newdata = data.frame(x = cp.nlls)))
	      graphics::points(x = c(cp.pl, cp.bc, cp.simexp, cp.nlls), y = c(lnPD.pl, lnPD.bc, lnPD.simexp, lnPD.nlls), pch = 8, col = c("darkgreen", "orange", "red", "blue"), cex = 1.5)
	      grDevices::dev.off()
	   }
	   
	} else if (method == "SiZer"){
	   message("Plotting the treemix results using SiZer.\n")
	   if(class(input) != "SiZer") stop("Input object is not of class SiZer.\n")
	   if(plot){
	   	plot(input, xlab = "m (migration edges)")
	    if(!is.null(pdf)){
	      grDevices::dev.copy(grDevices::pdf, file = pdf)
	      grDevices::dev.off()
	      message(paste0("Plot saved to file ", pdf, ".\n"))
	      } else if(is.null(pdf)) message("No plot file has been saved.\n")
	   }  else if (!plot){
	   	pdf(pdf, width = 7, height = 7)
	   	plot(input, xlab = "m (migration edges)")
	   	grDevices::dev.off()
	   }
	}
	
	message("Finished plotting.  All results are saved to the current directory as requested.\n")
	# Sys.sleep(2)
}
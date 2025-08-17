rm(list = ls())
Sys.setenv(LANG = "en")
setwd(this.path::here())
source("ipmr_lib.r") # NECESSARY LIBRARY! Insert the path to that file!
library(data.table)

################# INPUT #######################
# All values are expressed in SI units.
# heights_of_profile -> Vector of layer thicknesses (in meters) for evenly divided sections of the soil sample.

  height_of_profile <- c(as.double(0.01666666667))
  tdr_filename <- "input_data.tdr"
  uhc_filename <- "uhc.csv"
  plot_filename <- "plot_uhc.png"
  working_directory <- this.path::here()
############## END OF INPUT ###################
  
          # Paths of data
          input_fname <- sprintf("%s/%s", working_directory, tdr_filename)
          output_fname <- sprintf("%s/%s", working_directory, uhc_filename)
          plot_fname <- sprintf("%s/%s", working_directory, plot_filename )
          
############# UHC Calculation #################
          tryCatch(
            {# UHC results
              result_final <- ipm_calculate(input_fname, height_of_profile, out_fname = output_fname)
              
              # Create PNG plot
              png(filename = plot_fname, width = 2000, height = 1000, res = 150)
              par(mfrow = c(1, 1), mar = c(5, 6, 3, 1), oma = c(1, 1, 3, 1), cex.axis = 2.2)
              
              # Upper layer plot
              plot(log10(-result_final$h1_average), log10(result_final$k1_ipm), type = "l", lwd = 2,
                   main = "Upper layer", xlab = "log10(|h|) [m]", ylab = "log10(UHC) [m/s]", cex.lab = 1.5)
              
              # Lower layer plot
              lines(log10(-result_final$h2_average), log10(result_final$k2_ipm), type = "l", lwd = 2, col="red",
                   main = "Lower layer", xlab = "log10(|h|) [m]", ylab = "log10(UHC) [m/s]", cex.lab = 1.5)
              # Legend
              legend("topright", legend = c("Upper Layer", "Lower Layer"),
                     col = c("black", "red"), lwd = 2, bty = "n", cex = 1.3)
              
              # Global title
              mtext("Unsaturated Hydraulic Conductivity calculated by IPM",
                    outer = TRUE, cex = 1.5, font = 1)
              dev.off()
            },
            error = function(e) {
              msg <- e$message
              cat(msg)
              write(msg, file = sprintf("%s/error_log.txt", working_directory), append = TRUE)
            }
          )
#################### END #####################
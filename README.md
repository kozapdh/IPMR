# IPMR
This R script calculates Unsaturated Hydraulic Conductivity (UHC) using the Instantaneous Profile Method (IPM).

Authors:
*Maciej Kozyra, e-mail: xkozyra@gmail.com, m.kozyra@ipan.lublin.pl
*Krzysztof Lamorski, e-mail: k.lamorski@ipan.lublin.pl

Essential files: ipmr.R; ipmr_lib.R (library contains IPM procedure)
Exemplary file: input_data.tdr

Input for impr.R: 
*TDR file (CSV format), 
*Height of each soil layer (in meters), 
*Output filename to be created.

Output:
*CSV file containg: 
	**Time in seconds, 
	**UHC results for upper and lower layer in dependence of log10(|h|),
	**Averaged Soil Water Potential(SWP) for upper and lower layer.
*PNG file with plots of UHC vs. log10(|h|) for upper and lower layer.
*TXT file containing potential errors messages from the calculations.

An exemplary TDR file is included in the repository to presents the form of input file and to demonstrate exemplary UHC results in uhc.csv.

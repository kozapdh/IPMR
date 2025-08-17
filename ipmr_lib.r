rnz <- function(h, abs_limit=0.000001) {
    if (h == 0) {
        return(abs_limit)
    } else if (abs(h) < abs_limit) {
        return(sign(h) * abs_limit)
    } else {
        return(h)
    }
}
unsaturated_conductivity <- function(filtered_data, i, height_of_profile) {   
    gravitation <- 1
    # Extract boundary condition
    bonduary_condition <- as.double(0)
        # Central difference for the rest
        weighted_moist_0_fl <- filtered_data[i - 1, "M2"]
        weighted_moist_k_fl <- filtered_data[i + 1, "M2"]

        t_0_fl <- filtered_data[i - 1, "seconds"]
        t_k_fl <- filtered_data[i + 1, "seconds"]

        h0_k_fl <- filtered_data[i - 1, "X21"]
        h0_0_fl <- filtered_data[i - 1, "X41"]
        h1_k_fl <- filtered_data[i, "X21"]
        h1_0_fl <- filtered_data[i, "X41"]
        h2_k_fl <- filtered_data[i + 1, "X21"]
        h2_0_fl <- filtered_data[i + 1, "X41"]

        waga_i <- (1 / 2)
        waga_i_minus_1 <- (1 / 2) * ((filtered_data[i, "seconds"] - filtered_data[i - 1, "seconds"]) / (filtered_data[i + 1, "seconds"] - filtered_data[i - 1, "seconds"]))
        waga_i_plus_1 <- (1 / 2) * ((filtered_data[i + 1, "seconds"] - filtered_data[i, "seconds"]) / (filtered_data[i + 1, "seconds"] - filtered_data[i - 1, "seconds"]))

        weighted_diff_h_1 <- waga_i_minus_1 * ((h0_k_fl - h0_0_fl) / height_of_profile) +
            waga_i * ((h1_k_fl - h1_0_fl) / height_of_profile) +
            waga_i_plus_1 * (h2_k_fl - h2_0_fl) / height_of_profile

        # Central difference
        k1 <- -((weighted_moist_k_fl - weighted_moist_0_fl) / (t_k_fl - t_0_fl)) * height_of_profile / rnz(weighted_diff_h_1 + gravitation) + bonduary_condition
        
    time <- filtered_data[i, "seconds"]

    # Update variables for the second layer (suffix _sl)
    weighted_moist_0_sl <- filtered_data[i - 1, "M4"]
    weighted_moist_k_sl <- filtered_data[i + 1, "M4"]

    h0_k_sl <- filtered_data[i - 1, "X41"]
    h0_0_sl <- filtered_data[i - 1, "X61"]
    h1_k_sl <- filtered_data[i, "X41"]
    h1_0_sl <- filtered_data[i, "X61"]
    h2_k_sl <- filtered_data[i + 1, "X41"]
    h2_0_sl <- filtered_data[i + 1, "X61"]

    h_average_fl <- filtered_data[i, "X31"]
    h_average_sf <- filtered_data[i, "X51"]

    waga_i <- (1 / 2)
    waga_i_minus_1 <- (1 / 2) * ((filtered_data[i, "seconds"] - filtered_data[i - 1, "seconds"]) / (filtered_data[i + 1, "seconds"] - filtered_data[i - 1, "seconds"]))
    waga_i_plus_1 <- (1 / 2) * ((filtered_data[i + 1, "seconds"] - filtered_data[i, "seconds"]) / (filtered_data[i + 1, "seconds"] - filtered_data[i - 1, "seconds"]))

    weighted_diff_h_2 <- waga_i_minus_1 * ((h0_k_sl - h0_0_sl) / height_of_profile) +
        waga_i * (h1_k_sl - h1_0_sl) / height_of_profile +
        waga_i_plus_1 * (h2_k_sl - h2_0_sl) / height_of_profile

    # Calculate conductivity for the second layer (k2)
    k2 <- ((-(weighted_moist_k_sl - weighted_moist_0_sl) / (t_k_fl - t_0_fl)) * height_of_profile + k1 * rnz(weighted_diff_h_1 + gravitation)) / rnz(weighted_diff_h_2 + gravitation)


    return(setNames(
        c(time, k1, k2, h_average_fl, h_average_sf, weighted_diff_h_1, weighted_diff_h_2, weighted_moist_k_fl - weighted_moist_0_fl, weighted_moist_k_sl - weighted_moist_0_sl, weighted_moist_k_fl, weighted_moist_0_fl, t_k_fl - t_0_fl),
        c("time", "k1", "k2", "h_average_fl", "h_average_sf", "weighted_diff_h_1", "weighted_diff_h_2", "diff_weighted_moist_1", "diff_weighted_moist_2", "weighted_moist_k_fl", "weighted_moist_0_fl", "dt")
    ))
} # Procedure of calculating UHC using IPM method.
#in_fname -> Name of TDR file used to calculation of UHC using both of methods.
#height_of_profile -> Height of soil layer. [m]
#out_fname -> Name of file containing: time, UHC by IPM for first layer, UHC by IPM for second layer, UHC by MvG for first layer,UHC by MvG for second layer, Average SWP for first layer Average SWP for second layer.
ipm_calculate <- function(in_fname, height_of_profile, out_fname = NULL) {
    in_data <- fread(in_fname, header = FALSE, skip = 1)
    #set the names of the columns
    colnames(in_data) <- c("seconds", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "X11", "X21", "X31", "X41", "X51", "X61", "X71")
    # Moisture columns
    moist_column_names <- grep("^M", names(in_data), value = TRUE)
    wp_column_names <- grep("^X", names(in_data), value = TRUE)

    for (moist_column_name in moist_column_names) {
        # Check if all elements in the column are greater than 1
        if (all(in_data[[moist_column_name]] > 1)) {
            # If true, divide the values in the column by 100
            in_data[[moist_column_name]] <- in_data[[moist_column_name]] / 100.0
        } else {
            # Else, leave the values unchanged
            in_data[[moist_column_name]] <- in_data[[moist_column_name]]
        }
    }

    start_time <- Sys.time()
    results_list <- list()
    result_row <- list()
    for (i in 2:(nrow(in_data) - 1)) { # Fixed loop range
        tryCatch(
            {
                result <- unsaturated_conductivity(in_data, i, height_of_profile)
                time <- result$time
                k1_mualem <- result$k1
                k2_mualem <- result$k2
                h_average_fl <- result$h_average_fl
                h_average_sf <- result$h_average_sf

                result_row <- setNames(
                    c(
                        time, k1_mualem, k2_mualem, h_average_fl, h_average_sf
                    ),
                    c(
                        "time", "k1_ipm", "k2_ipm", "h1_average", "h2_average"
                    )
                )
                if ((k1_mualem <= 0 || k2_mualem <= 0) && (h_average_sf < -10 || h_average_fl < -10 )) {
                    break
                } else {
                    results_list[[i - 1]] <- result_row
                }
            },
            error = function(e) {
                cat("Error: ", conditionMessage(e), "\n")
            }
        )
    }
    # Convert list of data to array
    result_final <- as.data.frame(do.call(rbind, results_list))
    if (is.null(out_fname)) {
        out_fname <- paste0(tools::file_path_sans_ext(in_fname), "_results.csv")
    }
    write.table(result_final, file = out_fname, sep = " ", row.names = FALSE)
    return(result_final)
} # Procedure conversing TDR file to data.frame, using unsaturated_conductivity procedure, writng UHC results to csv file.
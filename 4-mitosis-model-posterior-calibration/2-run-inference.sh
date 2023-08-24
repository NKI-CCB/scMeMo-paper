# This will use 8*4 threads; adjust the j and k parameters for something suitable for the machine this is run on. 
$BCM3_ROOT/bin/bcminf --prior=prior_simulated_data.xml --likelihood=likelihood_simulated_data.xml -j 8 -k 4 --output.folder=output_simdata_correctspecified_t24_pt100
$BCM3_ROOT/bin/bcminf --prior=prior_simulated_data_misspecified.xml --likelihood=likelihood_simulated_data_misspecified.xml -j 8 -k 4 --output.folder=output_simdata_misspecified_t24_pt100
$BCM3_ROOT/bin/bcminf --prior=prior_simulated_data_misspecified.xml --likelihood=likelihood_simulated_data_misspecified_with_mitosis_timings.xml -j 8 -k 4 --output.folder=output_simdata_misspecified_with_mitosis_timings_t24_n100
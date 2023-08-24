#cp ../1-gather-data/Akopyan2014/Akopyan2014.nc data

# This will use 8*4 threads; adjust the j and k parameters for something suitable for the machine this is run on. 
$BCM3_ROOT/bin/bcminf --prior=prior_akopyan.xml --likelihood=likelihood_akopyan.xml -j 8 -k 4 --output.folder=output_akopyan
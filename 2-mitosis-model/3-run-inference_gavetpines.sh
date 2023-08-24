#cp ../1-gather-data/GavetPines2010/GavetPines2010.nc data

# This will use 8*4 threads; adjust the j and k parameters for something suitable for the machine this is run on. 
$BCM3_ROOT/bin/bcminf --prior=prior_gavetpines.xml --likelihood=likelihood_gavetpines.xml -j 8 -k 4 --output.folder=output_gavetpines
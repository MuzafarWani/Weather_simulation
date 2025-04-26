# Weather_simulation
The weather phenomenon simulated in this code is a simplified model of atmospheric advection and diffusion. The code is a simple 2D weather simulation that models the movement of air parcels (advection) and the smoothing of weather properties (diffusion) over a grid. This model doesn't represent the full complexity of atmospheric processes but serves as a basic illustration of how numerical simulations can be used to study weather phenomena.

This code writes the output the field data to a NetCDF (https://www.unidata.ucar.edu/software/netcdf/) file named "output.nc." NetCDF is a file format often used for storing scientific data with metadata. It defines dimensions (time, x, y) and a variable ("field") to store the simulation data. The data is written to the NetCDF file in a loop over time steps.

NetCDF, which stands for Network Common Data Form, is a set of software libraries and data formats that provide a way to store and exchange scientific data in a self-describing and platform-independent manner. NetCDF was developed to facilitate the storage, access, and sharing of large datasets, particularly in the fields of Earth sciences, atmospheric sciences, oceanography, and climate modeling. To use it, you can simply type ncview output.nc (http://meteora.ucsd.edu/~pierce/ncview_home_page.html).

To compile the code you need netcdf and mpich(https://www.mpich.org/downloads/versions/) installed. Then use the below command 

mpicc -o weather_simulation weather_sim.c -I/path/to/netcdf/include -L/path/to//netcdf/lib -lnetcdf

To run the code

mpirun -np 4 ./weather_simulation n (Where n is the number of time steps)

To increase the output file size you can increase the number of time steps n.

To change the size of the grid you can change the size of NX and NY in the code.

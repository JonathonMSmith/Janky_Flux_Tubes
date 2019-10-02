Plots flux tubes that are "connected" by the air-glow layer above all sky observatories and a spherical orbit near an ISS altitude

Uses pysatMagVect to work out the field lines.

Uses mayavi to plot all of this in 3d. Resolution of the flux tubes is currently very low.

![Sample Output](/snapshot_pysat_magvect.png)
Using pysatMagVect

This packages also calculates the maximum delta longitude between a satellite and a connected airglow layer flux tube.

This is simply done by finding delta lon between the air-glow layer at the observatory lat/lon and the conjugate point along the field line at the orbit altitude for each observatory. 

Some print statements are currently left in to output to the terminal the latitudes of the conjugate points at the orbit altitude.

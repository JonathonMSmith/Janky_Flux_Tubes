Plots flux tubes that are "connected" by the air-glow layer above all sky observatories and a spherical orbit near an ISS altitude

Uses mayavi to plot all of this in 3d. Resolution of the flux tubes is currently very low.

![Sample Output](/snapshot.png)


This packages also calculates the maximum delta longitude between a satellite and a connected airglow layer flux tube.

This is simply done by finding all of the flux tubes within 10 degrees latitude of the ground station and finding delta lon between the air-glow layer altitude and the orbit altitudealong each of those flux tubes. The maximum difference is the max delta lon for that observatory. 

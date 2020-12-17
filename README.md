# TTK4550_project

This code implements a model of the Butterfly Robot.

This code needs two dependancies. YALMIP and SDPT3. These must be in Matlabs path to be able to calculate the PRDE.

To run:
  -  Chose either butterfly_robot_phi_varphi or butterfly_robot_center. The first simulates perpetual motion while the second creates a center in varphi = pi/2.
  -  Both the constructors takes a boolean. If this is true the solutuion to the PRDE is calculated, this is needed to run the stabilizing controller. 

The simulink files expect there to exist an object of type eiteher butterfly_robot_phi_varphi or butterfly_robot_center with the name bf.
Create this by calling bf = butterfly_robot_phi_varphi(true) %(or false)
or similarly bf = butterfly_robot_center(true) %(or false)

The output from simulink is used in plotting_from_timeseries and bf.make_movie_butterfly.

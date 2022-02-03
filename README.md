# Robotics II Project 1

A simulation of a mobile robot in a room with furniture was created, where this robot has multiple sensors including LIDAR and UWB. The robot position was estimated using linear and non-linear single-timestep multilateration, using a localization Kalman filter, and with SLAM. Another Kalman filter mapped the UWB anchor positions given the robot position, and LIDAR mapping of the room was also performed. A room covering algorithm was developed which can correctly navigate around the room even without perfect robot state knowledge.

To run the simulation, use `gui_drive.mlx`

Example paths are `manual_path.mat` and `auto_path.mat`

To generate graphs, run the following files:

`single_timestep_UWB_solution.mlx`

`run_localization_kalman.mlx`

`run_mapping_kalman.mlx`

`run_slam_kalman.mlx`

`plot_room_covering.mlx`

`lidar_mapping.mlx`

# Generic control trajectory message
std_msgs/Header header

# control vector
Control[{{ solver_options.N_horizon }}] control_sequence

# Status, 0 = OK
uint8 status
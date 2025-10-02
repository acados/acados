# Generic state message
std_msgs/Header header

# state vector
float64[{{ dims.nx }}] x

# status, 0 = OK
uint8 status

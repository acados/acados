# Generic control message
std_msgs/Header header

# control vector
float64[{{ dims.nu }}] u

# Status, 0 = OK
uint8 status

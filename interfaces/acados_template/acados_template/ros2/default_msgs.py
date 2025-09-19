from .mapping_node import RosField, RosTopicMsg



# --- HEADER MESSAGE ---
STD_MSGS_HEADER_FIELD_TREE = [
    RosField(name="stamp", ftype="builtin_interfaces/Time", children=[
        RosField(name="sec", ftype="int32"),
        RosField(name="nanosec", ftype="uint32")]),
    RosField(name="frame_id", ftype="string")]

STD_MSGS_HEADER = RosTopicMsg()
STD_MSGS_HEADER.msg_type = "std_msgs/Header"
STD_MSGS_HEADER.field_tree = STD_MSGS_HEADER_FIELD_TREE


# --- Vector3 MESSAGE ---
GEOMETRY_MSGS_VECTOR3_FIELD_TREE = [
    RosField(name="x", ftype="float64"),
    RosField(name="y", ftype="float64"),
    RosField(name="z", ftype="float64")]

GEOMETRY_MSGS_VECTOR3 = RosTopicMsg()
GEOMETRY_MSGS_VECTOR3.msg_type = "geometry_msgs/Vector3"
GEOMETRY_MSGS_VECTOR3.field_tree = GEOMETRY_MSGS_VECTOR3_FIELD_TREE


# --- POINT MESSAGE ---
GEOMETRY_MSGS_POINT_FIELD_TREE = [
    RosField(name="x", ftype="float64"),
    RosField(name="y", ftype="float64"),
    RosField(name="z", ftype="float64")]

GEOMETRY_MSGS_POINT = RosTopicMsg()
GEOMETRY_MSGS_POINT.msg_type = "geometry_msgs/Point"
GEOMETRY_MSGS_POINT.field_tree = GEOMETRY_MSGS_POINT_FIELD_TREE


# --- QUATERNION MESSAGE ---
GEOMETRY_MSGS_QUATERNION_FIELD_TREE = [
    RosField(name="x", ftype="float64"),
    RosField(name="y", ftype="float64"),
    RosField(name="z", ftype="float64"),
    RosField(name="w", ftype="float64")]

GEOMETRY_MSGS_QUATERNION = RosTopicMsg()
GEOMETRY_MSGS_QUATERNION.msg_type = "geometry_msgs/Quaternion"
GEOMETRY_MSGS_QUATERNION.field_tree = GEOMETRY_MSGS_QUATERNION_FIELD_TREE


# --- INERTIA MESSAGE ---
GEOMETRY_MSGS_INERTIA_FIELD_TREE = [
    RosField(name="m", ftype="float64"),
    RosField(name="com", ftype="geometry_msgs/Vector3", children=GEOMETRY_MSGS_VECTOR3_FIELD_TREE),
    RosField(name="ixx", ftype="float64"),
    RosField(name="ixy", ftype="float64"),
    RosField(name="ixz", ftype="float64"),
    RosField(name="iyy", ftype="float64"),
    RosField(name="iyz", ftype="float64"),
    RosField(name="izz", ftype="float64")]

GEOMETRY_MSGS_INERTIA = RosTopicMsg()
GEOMETRY_MSGS_INERTIA.msg_type = "geometry_msgs/Inertia"
GEOMETRY_MSGS_INERTIA.field_tree = GEOMETRY_MSGS_INERTIA_FIELD_TREE


# --- TWIST MESSAGE ---
GEOMETRY_MSGS_TWIST_FIELD_TREE = [
    RosField(name="linear", ftype="geometry_msgs/Vector3", children=GEOMETRY_MSGS_VECTOR3_FIELD_TREE),
    RosField(name="angular", ftype="geometry_msgs/Vector3", children=GEOMETRY_MSGS_VECTOR3_FIELD_TREE)]


GEOMETRY_MSGS_TWIST = RosTopicMsg()
GEOMETRY_MSGS_TWIST.msg_type = "geometry_msgs/Twist"
GEOMETRY_MSGS_TWIST.field_tree = GEOMETRY_MSGS_TWIST_FIELD_TREE

# --- TWIST STAMPED MESSAGE ---
GEOMETRY_MSGS_TWIST_STAMPED_FIELD_TREE = [
    RosField(name="header", ftype="std_msgs/Header", children=STD_MSGS_HEADER_FIELD_TREE),
    RosField(name="twist", ftype="geometry_msgs/Twist", children=GEOMETRY_MSGS_TWIST_FIELD_TREE)
]


# --- POSE MESSAGE ---
GEOMETRY_MSGS_POSE_FIELD_TREE = [
    RosField(name="position", ftype="geometry_msgs/Point", children=GEOMETRY_MSGS_POINT_FIELD_TREE),
    RosField(name="orientation", ftype="geometry_msgs/Quaternion", children=GEOMETRY_MSGS_QUATERNION_FIELD_TREE)]

GEOMETRY_MSGS_POSE = RosTopicMsg()
GEOMETRY_MSGS_POSE.msg_type = "geometry_msgs/Pose"
GEOMETRY_MSGS_POSE.field_tree = GEOMETRY_MSGS_POSE_FIELD_TREE

# --- POSE STAMPED MESSAGE ---
GEOMETRY_MSGS_POSE_STAMPED_FIELD_TREE = [
    RosField(name="header", ftype="std_msgs/Header", children=STD_MSGS_HEADER_FIELD_TREE),
    RosField(name="pose", ftype="geometry_msgs/Pose", children=GEOMETRY_MSGS_POSE_FIELD_TREE)
]
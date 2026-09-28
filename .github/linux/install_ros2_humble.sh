#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

set -e

# 1. install dependencies for ROS 2 Humble
echo "Setting up system locales and ROS 2 package sources..."
sudo apt-get update
sudo apt-get install -y locales software-properties-common curl
sudo locale-gen en_US en_US.UTF-8
sudo update-locale LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
export LANG=en_US.UTF-8
sudo add-apt-repository universe -y

# 2. add the ROS 2 GPG key and repository
sudo curl -sSL https://raw.githubusercontent.com/ros/rosdistro/master/ros.key -o /usr/share/keyrings/ros-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/ros-archive-keyring.gpg] http://packages.ros.org/ros2/ubuntu $(. /etc/os-release && echo $UBUNTU_CODENAME) main" | sudo tee /etc/apt/sources.list.d/ros2.list > /dev/null

# 3. update package index
sudo apt-get update

# 4. check if ROS 2 is already installed (cached), if not install it
echo "Installing ros-humble-ros-base..."
sudo apt-get install -y ros-humble-ros-base

# 5. rosdep init and update
echo "Initializing and updating rosdep..."
sudo apt-get install -y python3-rosdep ros-dev-tools
sudo rosdep init || echo "rosdep already initialized."
rosdep update
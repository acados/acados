# Pull base image
FROM acados/rpi-debian-casadi:latest

# Set PATH environment to let casadi and acados find precompiled swig (duplicate)
ENV PATH="/home/pi/acados/external/swig:${PATH}"

# Set PYTHONPATH environment for casadi and acados lib, note that user is root
ENV PYTHONPATH="${PYTHONPATH}:/usr/local/lib:/usr/local/python:${HOME}/local/lib"

# Add current source to docker
ADD . /home/pi/acados/

# Compile acados with Python interface
RUN cd /home/pi/acados && \
    rm -rf build && mkdir build && \
    cd build && \
    cmake -D SWIG_MATLAB=0 -D SWIG_PYTHON=1 .. && \
    make install

# Test
RUN python3 -c "import acados"

# Define default command
CMD ["bash"]

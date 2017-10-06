# Pull base image
FROM doanminhdang/rpi-raspbian-casadi:jessie

# Set PATH environment to let casadi and acados find swig
ENV PATH="/home/pi/acados/external/swig:${PATH}"

ENV PYTHONPATH="${PYTHONPATH}:/usr/local/lib:/usr/local/python:~/local/lib"

# Compile acados
RUN cd /home/pi/acados && \
    mkdir build && \
    cd build && \
    cmake -D SWIG_MATLAB=0 -D SWIG_PYTHON=1 .. && \
    make install

# Make port 22 available to the world outside this container
EXPOSE 22

# Define default command
CMD ["bash"]

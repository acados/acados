which python3;
# sudo apt-get update -yqq;
# sudo apt-get --allow-unauthenticated install -yqq $CXX $CC $COVERAGE python3.8 python3.8-tk python3.8-dev;

# install virtualenv
sudo pip3 install virtualenv;
# create virtualenv
virtualenv --python=python3 acadosenv;
# source virtualenv
source acadosenv/bin/activate;
which python;

pip install interfaces/acados_template
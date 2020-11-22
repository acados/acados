which python3;
python3 --version;
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

# install tera
TERA_RENDERER_VERSION='0.0.34';
_TERA_RENDERER_GITHUB_RELEASES="https://github.com/acados/tera_renderer/releases/download/v${TERA_RENDERER_VERSION}/";
TERA_RENDERER_URL="${_TERA_RENDERER_GITHUB_RELEASES}/t_renderer-v${TERA_RENDERER_VERSION}-linux";

mkdir -p bin;
pushd bin;
	wget -O t_renderer "${TERA_RENDERER_URL}";
	chmod +x t_renderer
popd;
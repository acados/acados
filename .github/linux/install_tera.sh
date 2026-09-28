#!/bin/bash
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# install tera
TERA_RENDERER_VERSION='0.2.1';
_TERA_RENDERER_GITHUB_RELEASES="https://github.com/acados/tera_renderer/releases/download/v${TERA_RENDERER_VERSION}/";
TERA_RENDERER_URL="${_TERA_RENDERER_GITHUB_RELEASES}/t_renderer-v${TERA_RENDERER_VERSION}-linux-amd64";

mkdir -p bin;
pushd bin;
	wget -O t_renderer "${TERA_RENDERER_URL}";
	chmod +x t_renderer
popd;

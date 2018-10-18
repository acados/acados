#!/bin/bash

# Taken from: https://docs.travis-ci.com/user/caching/#ccache-cache
export HOMEBREW_NO_AUTO_UPDATE=1; # do not let homebrew update
brew install ccache;
export PATH="/usr/local/opt/ccache/libexec:$PATH";
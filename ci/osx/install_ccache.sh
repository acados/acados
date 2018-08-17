#!/bin/bash -e

# Taken from: https://docs.travis-ci.com/user/caching/#ccache-cache
brew install ccache;
export PATH="/usr/local/opt/ccache/libexec:$PATH";
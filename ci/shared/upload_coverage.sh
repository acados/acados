#!/bin/bash

pushd build;
	if [ "$COVERAGE" == 'lcov' ]; then
		bash <(curl -s https://codecov.io/bash);
	else
		echo "Codecov did not collect coverage reports";
	fi
popd

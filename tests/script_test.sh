#!/usr/bin/env bash

set -eu

FILE="../README.md"
if [ -f $FILE ]; then
    echo "File found! Test passed!"
else
    echo "File not found! Test failed!"
    exit 1
fi

exit 0

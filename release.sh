#!/usr/bin/env bash

rm -r dist
python setup.py sdist

twine upload dist/*

rm -r dist
rm -r src/rnaseq_lib.egg-info
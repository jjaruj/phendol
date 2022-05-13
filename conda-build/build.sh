#!/bin/sh
set -e

mkdir -p                		"${PREFIX}/bin"
mv ./blast2out.py   			"${PREFIX}/bin/"
mv ./Snakefile      			"${PREFIX}/bin/"
mv ./phendol        			"${PREFIX}/bin/"
mv ./lysin_database			"${PREFIX}/bin/"
mv ./phanotate.py			"${PREFIX}/bin/"

#!/bin/bash

if [[ "x${1}" == "x--help" || "x${1}" == "x-h" ]] ; then
  echo ""
  echo "This script is a wrapper for mod_dependency.py,"
  echo "Purpose is to generate figures of the module dependency"
  echo "of the gkw modules."
  echo "The image creation is done with dot, which therefore must to be found"
  echo "on your machine."
  echo "If you do not have it installed yet (check with 'which dot', most"
  echo "reasonable source is the graphviz package for your distribution."
  echo ""
  echo "Usage of this script:"
  echo "call it from the gkw main direcetory, e.g. with"
  echo ""
  echo "./scripts/mod_dependency.sh"
  echo ""
  echo "(If the scripts path is in your PATH variable, you can"
  echo "just use mod_dependency.sh.)"
  echo ""
  echo "This will run the python script and run dot on the result."
  echo "A figure in svg, eps and png format is produced, named"
  echo "'mod_dep' + the corresponding extension and put into doc/fig."
  echo ""
  exit 0
fi

which dot 2> /dev/null 1> /dev/null
if [ $? -eq 1 ] ; then
  echo ""
  echo "dot is needed to create the dependency graph but was not found."
  echo "Make sure dot is installed (can be found in the graphviz package) and"
  echo "that the location is in the PATH variable."
  echo ""
  exit 1;
fi

./python/mod_dependency.py `ls src/*90` > mod_dep.dot

IMAGE_DIR=doc/fig

dot -Tsvg -o${IMAGE_DIR}/mod_dep.svg mod_dep.dot
dot -Teps -o${IMAGE_DIR}/mod_dep.eps mod_dep.dot
dot -Tpng -o${IMAGE_DIR}/mod_dep.png mod_dep.dot

rm mod_dep.dot

#!/bin/bash

set -e

if [ ! -f "include/mike_namespace.h" ]; then
  echo "Please run script from the mike root directory"
  exit 1
fi

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' 's|//#define DISABLE_NAMESPACING|#define DISABLE_NAMESPACING|' ./include/mike_namespace.h
else
    sed -i 's|//#define DISABLE_NAMESPACING|#define DISABLE_NAMESPACING|' ./include/mike_namespace.h
fi

mkdir -p build_broadwell && cd build_broadwell && cmake -Dmike_BUILD_TYPE=broadwell .. && make -j8 && cd ..
mkdir -p build && cd build && cmake .. && make -j8
find . ../build_broadwell -name '*.a' -exec nm {} \; | grep '.c.o:\|T ' | scala -nc ../scripts/Namespace.scala > mike_namespace.h

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' 's|#define DISABLE_NAMESPACING|//#define DISABLE_NAMESPACING|' ../include/mike_namespace.h
else
    sed -i 's|#define DISABLE_NAMESPACING|//#define DISABLE_NAMESPACING|' ../include/mike_namespace.h
fi

diff mike_namespace.h ../include/mike_namespace.h


# Check the exit code of diff
if [ $? -eq 0 ]; then
  echo "No change in namespace."
  exit 0
else
  echo "Namespace changed, please update."
  exit 1
fi

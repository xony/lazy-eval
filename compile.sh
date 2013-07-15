#!/bin/sh
../sequoia-core/configure
make -C ../sequoia-core -j 8
aclocal
autoreconf --install --force
automake
./configure
make -j 8

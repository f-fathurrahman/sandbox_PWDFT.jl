#!/bin/bash

LIBXC_DIR=/home/efefer/mysoftwares/libxc-5.1.5/
#LIBXC_DIR=/home/efefer/mysoftwares/libxc-4.3.4/

#LINK1="-Wl,-rpath,$LIBXC_DIR/lib: $LIBXC_DIR/lib/libxc.so"
LINK1="$LIBXC_DIR/lib/libxc.a"

if [ "$#" -eq 1 ]; then

  basnam=`basename $1 .c`
  gcc $1 -I${LIBXC_DIR}/include -o $basnam.x $LINK1 -lm
  ./$basnam.x

else

  echo "Wrong number of arguments: $#"

fi


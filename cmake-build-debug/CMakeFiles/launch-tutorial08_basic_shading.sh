#!/bin/sh
bindir=$(pwd)
cd /home/bq2139/Documents/ogl/tutorial08_basic_shading/
export 

if test "x$1" = "x--debugger"; then
	shift
	if test "xYES" = "xYES"; then
		echo "r  " > $bindir/gdbscript
		echo "bt" >> $bindir/gdbscript
		/usr/bin/gdb -batch -command=$bindir/gdbscript --return-child-result /home/bq2139/Documents/ogl/cmake-build-debug/tutorial08_basic_shading 
	else
		"/home/bq2139/Documents/ogl/cmake-build-debug/tutorial08_basic_shading"  
	fi
else
	"/home/bq2139/Documents/ogl/cmake-build-debug/tutorial08_basic_shading"  
fi


# just for legacy.
# do not use: use makefile instead

qcc -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -python -c -fpic -I/usr/include/python3.8 qg.c
#qcc -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -c -fpic -I/usr/include/python3.8 qg.c
swig -I/home/bderembl/work/basilisk/basilisk/src -python -py3 qg.i

cc -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -c -fpic -I/usr/include/python3.8 -I/usr/lib/python3.8/site-packages/numpy/core/include/ qg_wrap.c

cc -shared qg.o qg_wrap.o -o _qg.so

qcc -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -python -c -fpic -I/usr/include/python3.7m pypg.c
#cp -f pg.i.bak pg.i
swig -I/home/bderembl/work/basilisk/basilisk/src -python -py3 pypg.i


cc -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -c -fpic -I/usr/include/python3.7m -I/usr/lib/python3.7/site-packages/numpy/core/include/ pypg_wrap.c
cc -shared pypg.o pypg_wrap.o -o _pypg.so

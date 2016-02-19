NVCC_OPTS = -arch=sm_35 -Xptxas -v -maxrregcount=256 

all: a.out b.out

# link objects(binaries) together
a.out: euler3d.cu define.c ./AVI/avilib.o ./AVI/avimod.o
	nvcc -arch=sm_35 euler3d.cu -I/usr/local/cuda-5.5/samples/common/inc  -L/usr/local/cuda-5.5/samples/lib -IAVI AVI/avilib.o AVI/avimod.o

b.out: hyperQ_euler3d.cu define.c ./AVI/avilib.o ./AVI/avimod.o
	nvcc -arch=sm_35 hyperQ_euler3d.cu -I/usr/local/cuda-5.5/samples/common/inc  -L/usr/local/cuda-5.5/samples/lib -IAVI AVI/avilib.o AVI/avimod.o -o b.out

./AVI/avilib.o ./AVI/avimod.o:
	cd AVI; make;


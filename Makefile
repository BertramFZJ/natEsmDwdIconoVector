FC = nvfortran
# COMPILER_OPTIONS = -mp -acc=gpu -gpu=ptxinfo $(CUDAVER) $(SMARCH) $(LTOFLAG) -Minfo=accel,vec,inline -cpp -O2 -Mnoautoinline
COMPILER_OPTIONS = -mp -Minfo=accel,vec,inline -cpp -O2
# -Mnovect
LINKER_OPTIONS = -mp
TARGET = vectorBasicExp

SOURCES = main.F90

OBJECTS = $(patsubst %.F90, %.o, $(SOURCES))

%.o : %.F90
	@echo "compiling $^ ..."
	@$(FC) $(COMPILER_OPTIONS) -c $^ -o $@

all : $(OBJECTS)
	@echo "building $(TARGET) ..."
	@$(FC) $(LINKER_OPTIONS) $(OBJECTS) -o $(TARGET).exe

clean :
	@rm -f $(TARGET).exe $(OBJECTS) *.mod

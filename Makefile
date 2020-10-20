FC = gfortran
EXE = run

OP_COMP = -O3 -fopenmp
OP_LINK = -O3 -fopenmp

SRC = main.f90 \
constant.f90 \
declaration.f90 \
grid.f90 \
misc.f90 \
boundary.f90 \
weno.f90 \
flux.f90 \
odesolver.f90 

OBJ = main.o \
constant.o \
declaration.o \
grid.o \
misc.o  \
boundary.o \
weno.o \
flux.o \
odesolver.o

all: $(EXE)

compile: $(OBJ)

$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(OBJ) $(OP_LINK)

# $(EXE): $(OBJ)
# 	$(FC) $(OP_LINK) -o $(EXE) $(OBJ)

$(OBJ):
	$(FC) $(OP_COMP) -c $(SRC_DIR)$(@:.o=.f90) -o $@

# %.o %.mod: %.f90
# 	$(FC) -c $(OP_COMP) $<

.PHONY = clean

clean:
	rm $(OBJ) *.mod
	
cleaner:
	$(clean)
	rm $(EXE)
	rm *.dat

constant.o: constant.f90

declaration.o: declaration.f90 constant.o

grid.o: grid.f90 constant.o declaration.o

misc.o: misc.f90 constant.o declaration.o

boundary.o: boundary.f90 constant.o declaration.o

weno.o: weno.f90 constant.o

flux.o: flux.f90 constant.o declaration.o misc.o weno.o

odesolver.o: odesolver.f90 constant.o declaration.o
# odesolver.o: odesolver.f90 constant.o declaration.o

# # viscous.o: viscous.f90 constant.o declaration.o

# # gradient.o: gradient.f90 constant.o declaration.o

# # flux.o: flux.f90 constant.o declaration.o

# # residual.o: residual.f90 constant.o declaration.o flux.o 

# # wenonew.o: wenonew.f90 constant.o declaration.o flux.o

# # odesolver.o: odesolver.f90 constant.o declaration.o

# # main.o: main.f90 constant.o declaration.o grid.o misc.o \
# # 		viscous.o gradient.o residual.o wenonew.o odesolver.o

main.o: main.f90 constant.o declaration.o grid.o misc.o \
		boundary.o flux.o odesolver.o
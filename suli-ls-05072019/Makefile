#======================================================================
# Makefile for Imcompact3D
#=======================================================================

# GNU
FC = gfortran
OPTFC = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium # -fcheck=all  #  #large #-floop-optimize -funroll-loops -ftree-fre -finline-functions -fpeel-loops -funsafe-math-optimizations -fprefetch-loop-arrays 
CC = gcc

# Intel
#FC = ifort
#OPTFC =  -fpp -O3 -xHost -IPO -autodouble
#CC = icc

FLAGS = -O3
 
# List of source files
SRC =   modules.f90 principal.f90 tempo.f90 convdiff.f90 advectivo.f90 pressao.f90 pos.f90 iniciais.f90 contorno.f90 viscosidade.f90 derivadas.f90 plot.f90 level_set.f90 paraview_3D.f90 mms.f90 restart.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJ =	$(SRC:.f90=.o)

rodar_3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBS)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod 

allclean:
	rm -f *.o *.mod rodar_3d
	rm -f arquivos/*
	rm -f dados/*


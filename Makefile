CPP = g++ 
FLAGS = -O5 -Wall -W
RM  = rm -f
EXEC = test

OBJECTS = \
	 Main.o \
	 matrix.o \
	 matrixsol.o \
	 function_evaluation.o \

all: $(OBJECTS) compile touch 

Main.o : Main.cpp
					 $(CPP) $(FLAGS) -c Main.cpp
matrix.o : matrix.cpp
			   		 $(CPP) $(FLAGS) -c matrix.cpp
matrixsol.o : matrixsol.cpp
			   		 $(CPP) $(FLAGS) -c matrixsol.cpp
function_evaluation.o : function_evaluation.cpp
			   		 $(CPP) $(FLAGS) -c function_evaluation.cpp
clean:  
					 $(RM) $(OBJECTS) $(EXEC) 

compile: 
					 $(CPP) $(FLAGS) $(OBJECTS) -o $(EXEC) 

touch:
					 @echo " "
					 @echo "Compilation done successfully..................."
					 @echo " "
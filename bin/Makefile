EXE_DIR = .
SRC_DIR = ../src
INCLUDE_DIR = ../src
OUT_DIR = ../output

CC = mpicc
EXEC = Swarm

ifdef serial
  CC = icc
endif

ifdef hdf5
  CC = h5pcc
endif

ifdef tau
  CC = h5pcc
endif

CFLAGS = -O1 -c -std=c99 -g -traceback 
L_FLAGS =
C_CMD = $(CC) $(CFLAGS)
C_CMD += $(addprefix -I, $(INCLUDE_DIR))

LINKER_C = $(CC) $(LFLAGS)

MAIN = main.o

C_OBJS = constants.o \
	 distributions.o \
	 random.o \
	 mesh.o \
	 tree.o \
	 refinment.o \
	 visit_writer.o \
	 print.o \
	 $(MAIN)

Swarm: Makefile $(C_OBJS) 
		$(LINKER_C) *.o -o $(EXEC) $(LFLAGS)

constants.o: $(SRC_DIR)/constants.c
	$(C_CMD) $(SRC_DIR)/constants.c
distributions.o: $(SRC_DIR)/distributions.c
	$(C_CMD) $(SRC_DIR)/distributions.c
random.o: $(SRC_DIR)/random.c
	$(C_CMD) $(SRC_DIR)/random.c
mesh.o: $(SRC_DIR)/mesh.c
	$(C_CMD) $(SRC_DIR)/mesh.c
tree.o: $(SRC_DIR)/tree.c
	$(C_CMD) $(SRC_DIR)/tree.c
refinment.o: $(SRC_DIR)/refinment.c
	$(C_CMD) $(SRC_DIR)/refinment.c
visit_writer.o: $(SRC_DIR)/visit_writer.c
	$(C_CMD) $(SRC_DIR)/visit_writer.c
print.o: $(SRC_DIR)/print.c
	$(C_CMD) $(SRC_DIR)/print.c
main.o: $(SRC_DIR)/main.c
	$(C_CMD) $(SRC_DIR)/main.c
	
all: $(C_CMD)
	
clean:
	rm -f Makefile~
	rm -f *.mod
	rm -f *.o
	rm -f *.h~
	rm -f *.c~
	rm -f *.vtk~	
	rm -f $(EXE_DIR)/$(EXEC)
	rm -f $(OUT_DIR)/*.*
	
clean_results:
	rm -f $(OUT_DIR)/*.*

run:
	$(EXE_DIR)/$(EXEC)


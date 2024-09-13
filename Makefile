# Makefile,
CC := g++ -std=c++17

# path macros
BIN_PATH := bin
OBJ_PATH := obj
DBG_PATH := debug
SRC_PATH := src

RADIUS := 1
SMIN := 2
SMAX := 3
BMIN := 3
BMAX := 3

#RADIUS := 5
#SMIN := 40
#SMAX := 80
#BMIN := 41
#BMAX := 81
#RADIUS := 4
#SMIN := 40
#SMAX := 80
#BMIN := 41
#BMAX := 81
#RADIUS := 3
#SMIN := 40
#SMAX := 80
#BMIN := 41
#BMAX := 81

TARGET_NAME := prog

TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

MEASURE_POWER := NO
CCOBJFLAGS=-O3 -fopenmp -march=native


GLOBALDEFINES := -Iinclude/ -Isrc/ -w -DRADIUS=${RADIUS} -DSMIN=${SMIN} -DSMAX=${SMAX} -DBMIN=${BMIN} -DBMAX=${BMAX} -D${MEASURE_POWER}
CCDEFINES :=
DBGDEFINES := -DDEBUG -DVERIFY


CPP_SRC := $(shell find $(SRC_PATH) -name "*.cpp")

# Generate the list of .o files with subdirectory structure
OBJ := $(patsubst $(SRC_PATH)/%.cpp, $(OBJ_PATH)/%.o, $(CPP_SRC))

DBG_OBJ := $(patsubst $(SRC_PATH)/%.cpp, $(DBG_PATH)/%.o, $(CPP_SRC))

default: makedir all

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CCDEFINES) $(CCOBJFLAGS) $(GLOBALDEFINES) -c -o $@ $<


$(DBG_PATH)/%.o: $(SRC_PATH)/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CCDEFINES) $(CCOBJFLAGS) $(GLOBALDEFINES) $(DBGDEFINES) -c -o $@ $<

$(TARGET): $(OBJ) $(CUDA_OBJ)
	$(CC) -o $@ $(OBJ) $(CUDA_OBJ)

$(TARGET_DEBUG) : $(DBG_OBJ) $(DBG_CUDA_OBJ)
	$(CC) -o $@ $(DBG_OBJ) $(DBG_CUDA_OBJ)


# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(BIN_PATH) $(OBJ_PATH) $(DBG_PATH)

.PHONY: debug
debug: $(TARGET_DEBUG)

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	-@rm -r $(DBG_PATH)/*
	-@rm -r $(OBJ_PATH)/*
	-@rm $(TARGET)
	-@rm $(TARGET_DEBUG)

-include $(OBJ_PATH)/*.d
-include $(DBG_PATH)/*.d


CPPFLAGS = --std=c++11 -O3 -I. -pthread

OBJ = trimZeroOne trimZeroOnePercentZerosAllowed trimIntegerMean trimZeroOneZerosAllowed

all: $(OBJ)

.PHONY: clean
clean:
	rm -rf $(OBJ)
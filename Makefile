CPPFLAGS = --std=c++11 -O3 -I.

OBJ = trimZeroOne trimZeroOneZerosAllowed trimZeroOnePercentZerosAllowed trimIntegerMean

all: $(OBJ)

.PHONY: clean
clean:
	rm -rf $(OBJ)
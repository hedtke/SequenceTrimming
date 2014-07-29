CPPFLAGS = --std=c++11 -O3 -I.

OBJ = trimZeroOne trimZeroOneZerosAllowed trimZeroOnePercentZerosAllowed trimIntegerMean trimZeroOneZerosAllowedPar trimZeroOnePercentZerosAllowedPar trimIntegerMeanPar

all: $(OBJ)

.PHONY: clean
clean:
	rm -rf $(OBJ)
CPPFLAGS = --std=c++11 -O3 -I.

OBJ = trimZeroOne trimZeroOnePercentZerosAllowed trimIntegerMean trimZeroOneZerosAllowed trimZeroOnePercentZerosAllowedPar trimIntegerMeanPar

all: $(OBJ)

.PHONY: clean
clean:
	rm -rf $(OBJ)
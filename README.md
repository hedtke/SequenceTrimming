# Sequence Trimming Algorithms
Ivo Hedtke, Ioana Lemnian, Matthias Mueller-Hannemann, Ivo Grosse.  
*On Optimal Read Trimming in Next Generation Sequencing and Its Complexity*  
**29 Jul 2014**.

## VERSION HISTORY
**Version 1.2:** trimZeroOneZerosAllowed, trimZeroOnePercentZerosAllowed and trimIntegerMean are now parallized. The number of worker threads that run in parallel to the "file reading thread" can be specified.

**Version 1.1:** Speed-Ups for trimZeroOnePercentZerosAllowed and trimIntegerMean (expected runtime is now linear). Worst-case runtime remains quadratic.

## FILES
| file(s)                            | description                                  |
| ---------------------------------- | -------------------------------------------- |
| ComputeMatrices.h                  | Algorithms that are called by *.cpp          |
| ComputeMatricesParallel.h          | Parallel algorithms that are called by *.cpp |
| Results.h                          | Export output file                           |
| ConcurrentQueue.h                  | Thread-safe queue for parallel algorithms    |
| tclap/\*                           | Parsing command line arguments               |
| trimZeroOne.cpp                    | Problem 0-zeros                              |
| trimZeroOneZerosAllowed.cpp        | Problem *z*-zeros                            |
| trimZeroOnePercentZerosAllowed.cpp | Problem *p*-percent                          |
| trimIntegerMean.cpp                | Problem *m*-mean                             |

## COMPILE
`make` or `make CXX=g++-4.8`

## INPUT FORMAT
The input is a FASTQ file with a shift for
the ASCII-Char -> Integer transformation. A threshold is used to say what qualities
are "good" and "bad": Let *t* be the threshold, *s* be the shift, *c* a char ASCII
score and *I(c)* the ASCII index of *c*. We say that *c* is "bad" (a "0") if *I(c)-s<t*.
Otherwise it is good ("1").

## OUTPUT FORMAT
If `--outfile` is used, the output file is an CSV format. It consists of the
columns "left", "right" and "reads". For each pair *left* <= *right* the number of
reads *g* in the input file is given such that *g*[*left*..*right*] fulfills the desired row constraint of the problem (0-zeros, *z*-zeros, *p*-percent, *m*-mean).

If `--outfile` is not used, the output on the screen lists the left border, the
right border, the width, the number of selected reads and the number of selected nucleotides.

## USAGE
### trimZeroOne
| parameter     | short | type   | required | description                                                                                |
| ------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`    | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`   | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`     | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`    | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--threshold` | `-t`  | int    | yes      | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`     | `-s`  | int    | yes      | which ASCII index represents the "0" quality?                                              |

### trimZeroOneZerosAllowed
| parameter       | short | type   | required | description                                                                                |
| --------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`      | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`     | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`       | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`      | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--zeros`       | `-z`  | int    | yes      | number of allowed zeros per read                                                           |
| `--threshold`   | `-t`  | int    | yes      | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`       | `-s`  | int    | yes      | which ASCII index represents the "0" quality?                                              |
| `--workthreads` | `-w`  | int    | no       | number of parallel worker threads (if omitted the sequential algorithm is used)            |

### trimZeroOnePercentZerosAllowed
| parameter       | short | type   | required | description                                                                                |
| --------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`      | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`     | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`       | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`      | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--percent`     | `-p`  | double | yes      | percent of allowed zeros per read: value between 0.0 and 1.0                               |
| `--threshold`   | `-t`  | int    | yes      | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`       | `-s`  | int    | yes      | which ASCII index represents the "0" quality?                                              |
| `--workthreads` | `-w`  | int    | no       | number of parallel worker threads (if omitted the sequential algorithm is used)            |

### trimIntegerMean
| parameter       | short | type   | required | description                                                                                |
| --------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`      | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`     | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`       | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`      | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--mean`        | `-m`  | double | yes      | min. mean per selected read                                                                |
| `--shift`       | `-s`  | int    | yes      | which ASCII index represents the "0" quality?                                              |
| `--workthreads` | `-w`  | int    | no       | number of parallel worker threads (if omitted the sequential algorithm is used)            |
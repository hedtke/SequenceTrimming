# Sequence Trimming Algorithms
Ivo Grosse, Ivo Hedtke, Ioana Lemnian, Matthias Mueller-Hannemann.  
*19 Dec 2013*.

## FILES
| file(s)                            | description                         |
| ---------------------------------- | ----------------------------------- |
| ComputeMatrices.h                  | Algorithms that are called by *.cpp |
| Results.h                          | Export output file                  |
| tclap/\*                           | Parsing command line arguments      |
| trimZeroOne.cpp                    | Problem 0-zeros                     |
| trimZeroOneZerosAllowed.cpp        | Problem *z*-zeros                   |
| trimZeroOnePercentZerosAllowed.cpp | Problem *p*-percent                 |
| trimIntegerMean.cpp                | Problem *m*-mean                    |

## COMPILE
```
c++ -I. -O3 trimZeroOne.cpp -o trimZeroOne
c++ -I. -O3 trimZeroOneZerosAllowed.cpp -o trimZeroOneZerosAllowed
c++ -I. -O3 trimZeroOnePercentZerosAllowed.cpp -o trimZeroOnePercentZerosAllowed
c++ -I. -O3 trimIntegerMean.cpp -o trimIntegerMean
```

## INPUT FORMAT
If the input file is a FASTQ file, use the switch `--fastq` and give a shift for
the ASCII-Char -> Integer tranformation. Use the threshold to say what qualities
are "good" and "bad": Let *t* be the threshold, *s* be the shift, *c* a char ASCII
score and *I(c)* the ASCII index of *c*. We say that *c* is "bad" (a "0") if *I(c)-s<t*.
Otherwise it is good ("1").

If the input file consists only of the quality score lines of a FASTQ files do
the same as above but do not use the `--fastq' switch.

If neither of `--fastq`, `--shift` and `--threshold` is used, the input file is a 0/1
ASCII grid.

## OUTPUT FORMAT
If `--outfile` is used, the output file is an CSV format. It consists of the
columns "left", "right" and "reads". For each pair *left* <= *right* the number of
reads *g* in the grid is given such that *g*[*left*..*right*] doesn't contain a zero.

If `--outfile` is not used, the output on the screen lists the left border, the
right border, the width, the number of selected reads and the area of the biggest
1-block in the grid.

## USAGE
### trimZeroOne
| parameter     | short | type   | required | description                                                                                |
| ------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`    | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`   | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`     | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`    | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--fastq`     | `-f`  |        | no       | is input file in FASTQ format?                                                             |
| `--threshold` | `-t`  | int    | no\*     | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`     | `-s`  | int    | no\*     | which ASCII index represents the "0" quality?                                              |

\*: You can only use `--shift` and `--threshold` together. It is only allowed to use
non of them or both. If `--fastq` is used, you have to use `--shift` and `–threshold`.

### trimZeroOneZerosAllowed
| parameter     | short | type   | required | description                                                                                |
| ------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`    | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`   | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`     | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`    | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--zeros`     | `-z`  | int    | yes      | number of allowed zeros per read                                                           |
| `--fastq`     | `-f`  |        | no       | is input file in FASTQ format?                                                             |
| `--threshold` | `-t`  | int    | no\*     | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`     | `-s`  | int    | no\*     | which ASCII index represents the "0" quality?                                              |

\*: You can only use `--shift` and `--threshold` together. It is only allowed to use
non of them or both. If `--fastq` is used, you have to use `--shift` and `–threshold`.

### trimZeroOnePercentZerosAllowed
| parameter     | short | type   | required | description                                                                                |
| ------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`    | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`   | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`     | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`    | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--percent`   | `-p`  | double | yes      | percent of allowed zeros per read: value between 0.0 and 1.0                               |
| `--fastq`     | `-f`  |        | no       | is input file in FASTQ format?                                                             |
| `--threshold` | `-t`  | int    | no\*     | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`     | `-s`  | int    | no\*     | which ASCII index represents the "0" quality?                                              |

\*: You can only use `--shift` and `--threshold` together. It is only allowed to use
non of them or both. If `--fastq` is used, you have to use `--shift` and `–threshold`.

### trimIntegerMean
| parameter     | short | type   | required | description                                                                                |
| ------------- | ----- | ------ | -------- | ------------------------------------------------------------------------------------------ |
| `--infile`    | `-i`  | string | yes      | file name of input file                                                                    |
| `--outfile`   | `-o`  | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--reads`     | `-r`  | int    | yes      | number of reads in the input file                                                          |
| `--length`    | `-l`  | int    | yes      | length of each read in the input file                                                      |
| `--mean`      | `-m`  | double | yes      | min. mean per selected read                                                                |
| `--threshold` | `-t`  | int    | yes      | quality scores less than the threshold are "bad", others are "good"                        |
| `--shift`     | `-s`  | int    | yes      | which ASCII index represents the "0" quality?                                              |
# Sequence Trimming Algorithms
Ivo Grosse, Ivo Hedtke, Ioana Lemnian, Matthias Mueller-Hannemann.  
*17 Dec 2013*.

## FILES

| file(s) | description |
| ------- | ----------- |
| ComputeMatrices.h | Algorithms that are called by *.cpp |
| Results.h | Export output file |
| tclap/\* | Parsing command line arguments |
| trimZeroOne.cpp | Problem 0-zeros |
| trimZeroOneZerosAllowed.cpp | Problem *z*-zeros |

## COMPILE
```
c++ -I. -O3 trimZeroOne.cpp -o trimZeroOne
c++ -I. -O3 trimZeroOneZerosAllowed.cpp -o trimZeroOneZerosAllowed
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
columns "left", "right" and "rows". For each pair left <= right the number of
rows g in the grid is given such that g[left..right] doesn't contain a zero.

If `--outfile` is not used, the output on the screen lists the left border, the
right border, the width, the number of selected rows and the area of the biggest
1-block in the grid.

## USAGE
### trimZeroOne
| parameter     | type   | required | description |
| ------------- | ------ | -------- | ----------- |
| `--infile`    | string | yes      | file name of input file |
| `--outfile`   | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--rows`      | int    | yes      | number of rows in the input file |
| `--length`    | int    | yes      | length of each row in the input file |
| `--fastq`     |        | no       | is input file in FASTQ format? |
| `--threshold` | int    | no\*   | quality scores less than the threshold are "bad", others are "good" |
| `--shift`     | int    | no\*   | which ASCII index represents the "0" quality? |

\*: You can only use `--shift` and `--threshold` together. It is only allowed to use
non of them or both. If `--fastq` is used, you have to use `--shift` and `–threshold`.

### trimZeroOneZerosAllowed
| parameter     | type   | required | description |
| ------------- | ------ | -------- | ----------- |
| `--infile`    | string | yes      | file name of input file |
| `--outfile`   | string | no       | file name of output file (CSV format), if skipped, only a short summary on screen is given |
| `--rows`      | int    | yes      | number of rows in the input file |
| `--length`    | int    | yes      | length of each row in the input file |
| `--zeros`     | int    | yes      | number of allowed zeros per row |
| `--fastq`     |        | no       | is input file in FASTQ format? |
| `--threshold` | int    | no\*   | quality scores less than the threshold are "bad", others are "good" |
| `--shift`     | int    | no\*   | which ASCII index represents the "0" quality? |

\*: You can only use `--shift` and `--threshold` together. It is only allowed to use
non of them or both. If `--fastq` is used, you have to use `--shift` and `–threshold`.
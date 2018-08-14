************************************************************************
file with basedata            : cr460_.bas
initial value random generator: 1346908903
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  128
RESOURCES
  - renewable                 :  4   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       25        9       25
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           7   9  11
   3        3          3           6  12  15
   4        3          3           5   9  11
   5        3          1           8
   6        3          3           7  10  11
   7        3          1          17
   8        3          2          12  14
   9        3          3          12  15  16
  10        3          2          13  14
  11        3          2          13  14
  12        3          1          17
  13        3          2          16  17
  14        3          1          16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0
  2      1     1       7    1    7    0    5    7
         2     2       6    0    4    0    4    6
         3     8       6    0    0    0    3    4
  3      1     5       7    0    5    6    7    4
         2     8       0    0    3    5    6    4
         3    10       6    5    0    0    6    3
  4      1     4       0    0    0    8    9    8
         2     6       0    7    0    0    8    7
         3     7       4    7    0    0    7    6
  5      1     7       0    0    5    4    6    6
         2     9       0    8    0    2    3    1
         3     9       1    4    0    3    2    1
  6      1     3       0    0    0    5    6    8
         2     5       0    0    6    4    4    4
         3     8       0    8    1    0    3    4
  7      1     2       4    0    0    4    7    6
         2     4       0    2    0    0    4    2
         3     4       0    0    7    0    3    3
  8      1     1       8    3    7    0    6    9
         2     6       3    0    0    4    4    8
         3     8       0    0    6    0    3    6
  9      1     5       9    0    0    9    7    7
         2     5       0   10    0    0    7    6
         3     9      10    0    0    0    3    5
 10      1     2      10    8    0    6    2    9
         2     3       0    5    1    0    2    7
         3     5       2    4    0    0    1    6
 11      1     3       4    9    8    8    4    3
         2     9       4    0    0    4    3    3
         3     9       4    9    0    4    2    3
 12      1     1       9    7    0    0    6    3
         2     2       0    0    8    0    5    3
         3     7       0    0    0    5    4    3
 13      1     6       9    6    2    9    8    6
         2     7       9    0    0    9    7    6
         3     8       0    4    0    0    6    6
 14      1     2       8    3    6    0   10    6
         2     7       6    0    6    0    8    4
         3     7       6    0    6   10    7    5
 15      1     2       0    0    0    4    6    6
         2     7       0    0    9    4    6    4
         3    10       0    0    8    4    6    4
 16      1     8       0    0    0    8    6    7
         2     9       0    0    0    8    5    3
         3     9       6    3    4    8    5    2
 17      1     3       6    2    9    3    7    8
         2    10       5    2    0    0    2    8
         3    10       0    0    9    0    4    7
 18      1     0       0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  N 1  N 2
   29   26   33   31  102  103
************************************************************************

************************************************************************
file with basedata            : cm555_.bas
initial value random generator: 253644628
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  143
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       12       13       12
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        5          3           7   9  12
   3        5          3           5   8  15
   4        5          2          10  14
   5        5          3           6  11  14
   6        5          1          16
   7        5          3           8  11  16
   8        5          1          17
   9        5          3          13  14  15
  10        5          3          11  12  16
  11        5          1          13
  12        5          2          13  15
  13        5          1          17
  14        5          1          17
  15        5          1          18
  16        5          1          18
  17        5          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     2       8    5    9    8
         2     5       8    4    9    8
         3     6       5    4    6    7
         4     8       5    3    4    4
         5     9       2    3    4    3
  3      1     2       5    8    8    8
         2     5       3    8    7    7
         3     6       2    5    5    6
         4     6       3    4    6    5
         5     9       2    2    4    4
  4      1     2       4    7    9    7
         2     5       4    6    9    7
         3     5       4    7    9    6
         4     6       4    6    8    5
         5     7       3    5    5    3
  5      1     2       6    6    7    7
         2     5       6    4    6    7
         3     6       6    4    6    6
         4     8       5    3    5    6
         5    10       5    2    5    5
  6      1     2       7    3    8   10
         2     3       6    3    7    9
         3     8       5    2    6    9
         4    10       5    1    4    9
         5    10       5    1    5    8
  7      1     3       7    9    9   10
         2     3       8    9    8    9
         3     4       5    8    7    8
         4     9       5    8    5    7
         5    10       1    8    2    6
  8      1     2       3    5    8    5
         2     6       3    4    7    3
         3     6       2    5    6    4
         4     6       2    4    7    5
         5     8       2    4    5    3
  9      1     3       7    9    7    6
         2     7       6    9    6    6
         3     8       6    8    6    6
         4     8       6    9    4    6
         5     9       5    8    3    4
 10      1     4       5   10    7    7
         2     7       5    7    7    7
         3     8       4    6    6    6
         4     9       4    5    5    6
         5    10       4    3    4    6
 11      1     1       8    5   10    7
         2     5       8    4    9    7
         3     6       7    3    6    7
         4     7       7    2    6    6
         5    10       6    1    4    5
 12      1     3       7    9    9    9
         2     4       7    7    9    7
         3     7       7    5    8    7
         4     8       6    5    7    3
         5     8       5    3    8    5
 13      1     2       9    7   10    8
         2     4       7    6   10    8
         3     5       6    5   10    8
         4     5       6    6    9    8
         5     7       4    4    9    7
 14      1     1       7    7    8    9
         2     2       6    5    6    8
         3     2       5    6    5    8
         4     8       4    5    5    7
         5     8       5    3    3    5
 15      1     1       8    5    7    8
         2     1       9    6    6    8
         3     5       7    5    6    8
         4     6       6    4    5    8
         5    10       6    3    3    8
 16      1     1       6    9    7    6
         2     2       4    7    7    5
         3     3       4    6    6    5
         4     8       3    4    6    4
         5     9       3    1    6    4
 17      1     1       6    2    1    8
         2     2       6    2    1    7
         3     4       6    2    1    5
         4     8       6    2    1    4
         5     9       6    2    1    2
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   22   28  110  111
************************************************************************

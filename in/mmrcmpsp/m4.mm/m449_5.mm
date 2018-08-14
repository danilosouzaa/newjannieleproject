************************************************************************
file with basedata            : cm449_.bas
initial value random generator: 1061086277
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  132
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       15        3       15
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        4          3           6  12  17
   3        4          3           7   8  16
   4        4          3           5  12  13
   5        4          2          11  16
   6        4          1           9
   7        4          3          10  11  14
   8        4          3           9  10  11
   9        4          1          13
  10        4          2          15  17
  11        4          2          15  17
  12        4          2          14  16
  13        4          1          14
  14        4          1          15
  15        4          1          18
  16        4          1          18
  17        4          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       0    7    9    6
         2     1      10    0    6    5
         3     3       0    7    5    5
         4     7       0    6    1    4
  3      1     1       4    0    6    8
         2     1       4    0    5    9
         3     3       0    6    5    8
         4     5       4    0    4    8
  4      1     4       5    0    6    8
         2     5       1    0    6    7
         3     8       0    7    5    4
         4    10       0    7    4    2
  5      1     2       0    9    6    9
         2     3       0    6    6    9
         3     6       7    0    6    9
         4     9       0    4    5    9
  6      1     2       0    3    2    9
         2     4       3    0    2    8
         3     8       3    0    2    7
         4     9       0    3    1    6
  7      1     2       9    0    8    5
         2     5       0    9    6    4
         3    10       0    9    5    1
         4    10       0    9    4    2
  8      1     4       8    0    5    7
         2     5       7    0    4    7
         3     8       7    0    4    6
         4    10       0    4    4    5
  9      1     2       8    0    2    8
         2     2       7    0    3    8
         3     9       0    9    2    6
         4    10       0    9    1    6
 10      1     2       0   10    9    7
         2     3       0   10    9    5
         3     4       0   10    9    3
         4     7       0    9    9    2
 11      1     4       6    0    7    6
         2     4       0    8    7    6
         3     6       6    0    5    6
         4     7       4    0    1    5
 12      1     2       0    8    3    8
         2     6       6    0    2    7
         3     8       0    7    2    7
         4     9       0    7    2    5
 13      1     4       8    0    8    3
         2     6       0   10    5    3
         3     6       7    0    5    3
         4     8       6    0    4    2
 14      1     1       5    0    9    2
         2     5       0    8    7    2
         3     7       4    0    5    2
         4     7       0    5    6    2
 15      1     3       5    0    9    9
         2     4       0   10    7    9
         3     8       4    0    7    9
         4     9       0    8    6    8
 16      1     1       7    0    6    8
         2     2       4    0    5    8
         3     4       0   10    4    8
         4     9       3    0    3    6
 17      1     1       0   10    5    8
         2     5       8    0    4    8
         3     6       0    9    4    8
         4     6       6    0    4    8
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
    6   16   90  104
************************************************************************

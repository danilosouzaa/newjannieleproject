************************************************************************
file with basedata            : cm461_.bas
initial value random generator: 21653
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  140
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       33        1       33
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        4          3           5   6  12
   3        4          3           6   7  12
   4        4          2           8  14
   5        4          3           9  11  17
   6        4          3           8   9  11
   7        4          3          13  16  17
   8        4          2          15  17
   9        4          2          10  16
  10        4          1          13
  11        4          2          13  16
  12        4          1          14
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
  2      1     5       9    7   10    7
         2     9       9    7    6    5
         3    10       9    4    4    3
         4    10       9    2    6    3
  3      1     6       7    5    2    5
         2     6       7    6    2    4
         3     7       6    4    2    3
         4     9       6    3    1    1
  4      1     1      10    8    5    3
         2     4      10    8    4    3
         3     7      10    6    3    3
         4     9       9    5    3    2
  5      1     1       5    5    7    6
         2     2       4    5    7    5
         3     6       4    3    6    5
         4     8       2    1    4    4
  6      1     2      10    2    9    5
         2     4       9    2    9    4
         3     7       8    2    7    2
         4     9       6    2    5    2
  7      1     3       7    9    4   10
         2     4       7    8    4   10
         3     8       7    6    3    9
         4     9       6    5    3    9
  8      1     2      10    8    9    8
         2     6      10    7    9    8
         3     7      10    7    8    8
         4    10       9    7    7    8
  9      1     4       7    6    9    8
         2     5       7    6    8    8
         3     8       6    6    6    7
         4     9       6    6    5    7
 10      1     6      10    8    7    4
         2     7       9    7    7    3
         3    10       9    4    7    3
         4    10       9    5    7    2
 11      1     2       8    8    8    6
         2     4       5    5    7    5
         3     6       4    2    5    5
         4     6       5    4    5    4
 12      1     1       4    6    8    1
         2     4       3    5    8    1
         3     5       2    5    8    1
         4    10       1    4    8    1
 13      1     1       6    9    3    2
         2     1       7    8    3    2
         3     3       6    7    3    1
         4     8       6    5    2    1
 14      1     7       9    8    8    8
         2     9       8    8    7    7
         3     9       7    8    7    8
         4    10       7    8    6    7
 15      1     7      10    9    3    9
         2     7       6    9    4    9
         3     7       6   10    3    9
         4     9       6    9    1    7
 16      1     1       6    9    6    9
         2     3       6    6    5    8
         3     5       5    5    5    5
         4     8       5    4    4    5
 17      1     2       3    6    5    5
         2     2       3    6    6    4
         3     6       2    6    4    3
         4     6       3    5    4    4
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   15   15  105   96
************************************************************************

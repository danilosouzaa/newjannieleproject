************************************************************************
file with basedata            : md114_.bas
initial value random generator: 1541979461
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  14
horizon                       :  94
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     12      0       14        5       14
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           5  10
   3        3          2           7  10
   4        3          2           5  10
   5        3          3           6   7   8
   6        3          1           9
   7        3          3           9  11  12
   8        3          3           9  11  12
   9        3          1          13
  10        3          3          11  12  13
  11        3          1          14
  12        3          1          14
  13        3          1          14
  14        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       0    4    4    9
         2     8       0    3    4    9
         3     9       2    0    4    7
  3      1     4       7    0   10    8
         2     5       0    5    7    8
         3     6       2    0    6    7
  4      1     1       0    7    9    3
         2     6       8    0    7    2
         3     7       5    0    6    1
  5      1     4       7    0    9   10
         2     9       3    0    6    7
         3     9       4    0    4    7
  6      1     5       0    8    5    2
         2     7       3    0    4    2
         3     8       0    7    4    1
  7      1     1       9    0    3    7
         2     1       0    7    3    6
         3     9       8    0    3    5
  8      1     3       0    9    8   10
         2     5       2    0    8    7
         3     8       0    4    7    6
  9      1     1       0    7    3    3
         2     5       0    6    3    3
         3    10       4    0    2    2
 10      1     6       0    3    9    9
         2     7       0    2    7    6
         3     8       4    0    6    5
 11      1     4       0    6    4    9
         2     7       9    0    4    8
         3     8       8    0    1    7
 12      1     1       6    0    9    7
         2     3       0    3    9    5
         3     6       3    0    8    4
 13      1     1       0    5    9   10
         2     5       0    5    8    8
         3     6       0    5    8    6
 14      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   11   16   76   80
************************************************************************

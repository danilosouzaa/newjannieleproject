************************************************************************
file with basedata            : cm423_.bas
initial value random generator: 1146923310
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  112
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       23        2       23
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        4          3           7   9  14
   3        4          2          15  17
   4        4          3           5   6  10
   5        4          2           7   8
   6        4          3           9  12  16
   7        4          1          12
   8        4          3           9  11  16
   9        4          2          13  17
  10        4          2          11  16
  11        4          2          12  14
  12        4          1          13
  13        4          1          15
  14        4          2          15  17
  15        4          1          18
  16        4          1          18
  17        4          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       2    6    2    0
         2     2       1    4    0    4
         3     2       2    5    0    3
         4     3       1    2    0    1
  3      1     1       8    8    0    6
         2     3       4    5    9    0
         3     5       2    5    9    0
         4     5       1    4    0    4
  4      1     3       8    9    9    0
         2     3       7    9    0    4
         3     4       6    8    8    0
         4     8       5    8    7    0
  5      1     3       2    4    0    6
         2     8       1    3    7    0
         3     8       2    3    6    0
         4     9       1    2    5    0
  6      1     1      10    8   10    0
         2     4      10    6   10    0
         3     6      10    6    9    0
         4     8      10    2    9    0
  7      1     1      10   10   10    0
         2     2       8    8   10    0
         3     3       7    7    9    0
         4     5       6    4    0    9
  8      1     3       8    7    7    0
         2     4       7    7    0    9
         3     5       3    7    0    8
         4     6       2    5    0    7
  9      1     3       2    3    6    0
         2     5       2    3    0    8
         3     9       1    1    0    7
         4     9       1    1    5    0
 10      1     3      10    5    6    0
         2     8       8    3    4    0
         3     8       8    3    0    6
         4     9       7    2    0    4
 11      1     3       3    7    5    0
         2     3       3    7    0    3
         3     5       3    4    0    3
         4     6       2    2    0    2
 12      1     1       9    6    6    0
         2     1       7    6    0    7
         3     5       6    6    0    2
         4     5       5    5    6    0
 13      1     7       3    6    4    0
         2     9       3    6    2    0
         3     9       2    6    4    0
         4    10       1    5    0    3
 14      1     2       9    7    9    0
         2     5       9    7    0    7
         3     5       9    6    8    0
         4     7       8    6    6    0
 15      1     3       6    7    4    0
         2     3       7    7    0    2
         3     6       6    6    0    1
         4     9       5    6    5    0
 16      1     1       6   10    7    0
         2     2       5    8    0    2
         3     2       5    9    7    0
         4     4       5    7    2    0
 17      1     3       8    9    0    4
         2     5       7    8    9    0
         3     7       7    8    0    4
         4     9       5    7    0    4
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   19   19   86   60
************************************************************************

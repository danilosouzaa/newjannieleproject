************************************************************************
file with basedata            : cr129_.bas
initial value random generator: 114147013
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  126
RESOURCES
  - renewable                 :  1   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       22        8       22
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   9  14
   3        3          3           6  10  13
   4        3          3           5   9  14
   5        3          3           7   8  11
   6        3          3          11  12  16
   7        3          1          15
   8        3          1          10
   9        3          2          11  16
  10        3          2          12  17
  11        3          2          15  17
  12        3          1          15
  13        3          1          14
  14        3          2          16  17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0
  2      1     2      10    0    9
         2     7       8    5    0
         3     7       9    0    7
  3      1     2       4    2    0
         2     2       6    0    9
         3     4       4    0    5
  4      1     5       9    6    0
         2     9       5    3    0
         3    10       2    3    0
  5      1     8       2    0    5
         2     8       2    8    0
         3    10       2    5    0
  6      1     1       8    6    0
         2     1       6    8    0
         3     7       6    4    0
  7      1     3       8    9    0
         2     9       6    0    5
         3    10       3    0    4
  8      1     3       8    9    0
         2     5       7    9    0
         3     6       7    0    7
  9      1     5       4    8    0
         2     6       4    6    0
         3     9       2    2    0
 10      1     3       6    6    0
         2     9       4    0    8
         3    10       3    3    0
 11      1     3       8    0    3
         2     3       8    9    0
         3     5       3    0    2
 12      1     2       8    3    0
         2     2      10    0    5
         3     4       7    0    4
 13      1     4       9    0    7
         2     7       7    9    0
         3    10       5    8    0
 14      1     1       9    9    0
         2     6       9    6    0
         3     9       8    0    6
 15      1     1       7    7    0
         2     3       7    0    8
         3     5       5    0    7
 16      1     1      10    0    7
         2     3      10    2    0
         3    10       9    0    5
 17      1     2       6    0    6
         2     7       5    1    0
         3    10       4    1    0
 18      1     0       0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  N 1  N 2
   13  101   85
************************************************************************

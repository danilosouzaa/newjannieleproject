************************************************************************
file with basedata            : md137_.bas
initial value random generator: 1294092442
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  16
horizon                       :  111
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     14      0       18        8       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5  11  14
   3        3          3           7  10  12
   4        3          3           6   7  10
   5        3          1           8
   6        3          3           9  12  15
   7        3          3           9  14  15
   8        3          2          10  13
   9        3          1          11
  10        3          1          15
  11        3          1          13
  12        3          2          13  14
  13        3          1          16
  14        3          1          16
  15        3          1          16
  16        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       8    0    8    0
         2     3       1    0    0    4
         3     7       0    6    8    0
  3      1     4       0    7    6    0
         2    10       3    0    0   10
         3    10       0    2    4    0
  4      1     4       3    0    0   10
         2     6       3    0    4    0
         3     8       2    0    0    9
  5      1     2       9    0    0    2
         2     4       0    4    0    2
         3     5       0    3    0    2
  6      1     6       0    6    7    0
         2     6       0    4    0    6
         3     7       8    0    0    6
  7      1     6       4    0    5    0
         2     7       3    0    0    8
         3     8       0    5    0    8
  8      1     5       9    0    0    3
         2     5       0    5    9    0
         3     5       0    9    0    3
  9      1     3       6    0    0    4
         2     3       0    5    0    4
         3    10       0    4    0    4
 10      1     3       0    5    0   10
         2     6       9    0    0    7
         3     8       0    5    0    5
 11      1     2       5    0    5    0
         2     5       5    0    0    9
         3     8       3    0    0    8
 12      1     5       4    0    3    0
         2     7       3    0    0    8
         3     9       0    8    0    7
 13      1     3       7    0    0    5
         2     4       0    9    6    0
         3     6       0    8    5    0
 14      1     3       0    5    0    9
         2     6       5    0    8    0
         3    10       0    5    6    0
 15      1     5       9    0    0    8
         2     6       0    4    5    0
         3    10       0    4    0    8
 16      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
    8    6   33   54
************************************************************************

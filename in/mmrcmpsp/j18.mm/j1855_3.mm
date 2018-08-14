************************************************************************
file with basedata            : md311_.bas
initial value random generator: 1064483876
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  20
horizon                       :  146
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     18      0       24        1       24
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           7   9  15
   3        3          3           5  14  16
   4        3          3           6   8  18
   5        3          3           6  11  15
   6        3          1          10
   7        3          2          11  16
   8        3          3          10  14  15
   9        3          3          12  14  18
  10        3          2          12  19
  11        3          1          13
  12        3          1          13
  13        3          1          17
  14        3          2          17  19
  15        3          1          17
  16        3          1          18
  17        3          1          20
  18        3          1          20
  19        3          1          20
  20        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       8    8   10    4
         2     1       8    8    9    5
         3     8       7    7    8    3
  3      1     3       6    6    6    7
         2     5       4    5    5    7
         3    10       2    5    5    6
  4      1     1       8    7    4    8
         2     8       7    6    3    7
         3    10       5    6    1    5
  5      1     1      10   10    9    6
         2     3       9    4    9    5
         3     8       9    1    8    4
  6      1     3       2    9    6    9
         2     4       2    6    5    9
         3     6       2    3    3    8
  7      1     6       7    7    7    9
         2     8       2    7    7    7
         3     8       4    6    7    8
  8      1     6       5    6    2    4
         2     6       7    5    2    5
         3    10       3    4    2    1
  9      1     3       6    6    5    8
         2     9       6    5    5    5
         3    10       6    5    2    4
 10      1     1       2    6    6   10
         2     3       1    5    6    9
         3     5       1    3    3    8
 11      1     3       9    6    9    2
         2     5       8    5    8    2
         3     8       6    3    6    1
 12      1     4       5    6   10    8
         2     5       5    4    7    8
         3     7       4    1    6    5
 13      1     8       2    6    3   10
         2    10       1    5    1    7
         3    10       2    6    2    6
 14      1     1       5    7    7    9
         2     3       3    5    7    6
         3     7       2    3    7    6
 15      1     2       8   10    8    7
         2     2       8    9    9    8
         3     6       5    4    8    5
 16      1     3       8    4    7    6
         2     7       6    4    6    4
         3     8       3    3    6    3
 17      1     4       2   10    9    6
         2     5       2    9    8    5
         3     7       2    6    8    5
 18      1     6       6    5    6    8
         2     7       4    5    5    8
         3     8       2    4    4    7
 19      1     2       7    7    5    8
         2     5       3    6    5    7
         3    10       3    6    2    7
 20      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   30   27  112  122
************************************************************************

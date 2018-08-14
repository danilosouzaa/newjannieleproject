************************************************************************
file with basedata            : cn38_.bas
initial value random generator: 846396232
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  128
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  3   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       20        9       20
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           5  10
   3        3          3           7   8  11
   4        3          3           5   6   9
   5        3          1          12
   6        3          1           7
   7        3          3          14  15  16
   8        3          2          12  14
   9        3          3          11  15  17
  10        3          3          11  12  13
  11        3          1          16
  12        3          3          15  16  17
  13        3          1          14
  14        3          1          17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2  N 3
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     3       5   10    0    0    2
         2     4       4    9    0    0    1
         3     4       5    9    4    0    0
  3      1     1       5    8    6    0    9
         2     2       4    7    6    0    0
         3     9       4    4    5    0    9
  4      1     8       7    5    0    8    6
         2     8       7    6    5    0    0
         3     9       7    5    1    8    0
  5      1     7       6    6    5    0    6
         2    10       5    4    1    0    0
         3    10       4    2    0    0    3
  6      1     3       9    7    5    0    0
         2     8       8    6    3    6    0
         3     9       8    4    3    0    0
  7      1     1       6    8    6    7    0
         2     2       6    7    0    0    6
         3     8       5    4    0    0    3
  8      1     6       4   10    5   10    0
         2     9       3   10    0   10    4
         3    10       2    9    4    0    2
  9      1     1       5    8    0    3    0
         2     2       5    7    0    2    0
         3     3       4    7    0    2    6
 10      1     1       8    7    7    6    0
         2     4       8    7    7    0    0
         3    10       6    2    3    0    9
 11      1     1       7    8    0    0    8
         2     5       6    7    4    0    0
         3     8       6    6    0    8    0
 12      1     1       8    4    0    4    0
         2     3       7    4    0    0    8
         3     4       5    4    2    1    0
 13      1     8       9    9    6    8    0
         2     8       9    6    0    0    7
         3     9       7    3    0    0    5
 14      1     1       6    5   10    6    0
         2     6       6    4    0    0    4
         3     8       4    2    8    4    0
 15      1     1       3    6    0    6    4
         2     7       3    3    5    3    0
         3     7       1    5    0    4    4
 16      1     3       5   10    9    0    8
         2     9       5   10    0    4    7
         3    10       3   10    5    0    6
 17      1     4       7    8    0    8    0
         2     4       5    8    0    4    9
         3    10       3    8    3    0    6
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2  N 3
   31   30   29   25   29
************************************************************************

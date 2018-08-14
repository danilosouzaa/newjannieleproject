************************************************************************
file with basedata            : cr521_.bas
initial value random generator: 9261
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  132
RESOURCES
  - renewable                 :  5   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       25       10       25
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   8  14
   3        3          3           5   8  14
   4        3          2           8   9
   5        3          2           6  11
   6        3          2           7  13
   7        3          2           9  15
   8        3          1          10
   9        3          1          10
  10        3          1          12
  11        3          2          12  13
  12        3          2          16  17
  13        3          3          15  16  17
  14        3          3          15  16  17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  R 5  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0    0
  2      1     2       5   10    7    8    7    0    3
         2    10       5    2    4    4    6    0    3
         3    10       5    3    2    4    6    9    0
  3      1     2       5   10    2    5    4   10    0
         2     4       4    9    2    5    3   10    0
         3     6       2    9    1    4    3    0    2
  4      1     3       6    3   10    6    3    7    0
         2     6       5    3    9    5    1    0    7
         3     6       5    3    8    3    1    3    0
  5      1     1       8    4    7    2    7    8    0
         2     3       6    4    5    2    7    7    0
         3     4       4    4    3    2    7    7    0
  6      1     4       5    4    6    2    4    6    0
         2     8       3    4    5    1    3    5    0
         3     9       3    3    5    1    3    0    3
  7      1     6       7    9    7    7    5    8    0
         2     8       5    7    7    6    5    8    0
         3    10       1    6    5    3    2    7    0
  8      1     6       8    9    7    8    9    0    4
         2     6       9    7    9    7    8    0    4
         3     8       2    6    4    5    2    9    0
  9      1     2       9    7    6    6    5    1    0
         2     8       7    6    3    4    5    0   10
         3     9       5    3    3    1    4    0    9
 10      1     4       7    7    8    7    7    0    9
         2     5       5    6    8    7    7    0    9
         3     9       5    5    8    7    7    0    9
 11      1     2       5    6    8    9    8    5    0
         2     8       4    3    8    8    4    0   10
         3     9       4    3    7    6    4    5    0
 12      1     3       7    2    6    7   10    9    0
         2     6       7    2    5    4    9    7    0
         3     7       6    2    3    2    9    5    0
 13      1     1       8    8    7    4    5    6    0
         2     3       3    8    6    4    5    0    5
         3    10       3    7    3    3    4    5    0
 14      1     1       9    2    6    5    8    0    2
         2     2       9    1    5    3    7    4    0
         3     7       9    1    5    1    2    3    0
 15      1     4       7   10    5    3    6    9    0
         2     9       5   10    5    3    5    7    0
         3    10       2    9    4    1    5    4    0
 16      1     3       9    6    6    6    6    0    5
         2     6       6    6    6    5    6    8    0
         3     8       6    5    5    3    4    0    3
 17      1     3       2    5    6    5    7   10    0
         2     3       2    4    7    5    7    6    0
         3    10       1    4    6    3    7    0   10
 18      1     0       0    0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  R 5  N 1  N 2
   13   13   12   10   12   88   55
************************************************************************

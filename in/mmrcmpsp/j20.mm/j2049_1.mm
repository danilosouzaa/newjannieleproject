************************************************************************
file with basedata            : md369_.bas
initial value random generator: 5887
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  22
horizon                       :  140
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     20      0       16        4       16
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6   8
   3        3          3           8  10  17
   4        3          3           5   9  12
   5        3          2          11  19
   6        3          2           7  13
   7        3          1          10
   8        3          3          11  15  16
   9        3          2          11  15
  10        3          2          15  19
  11        3          1          18
  12        3          2          13  18
  13        3          1          14
  14        3          1          16
  15        3          2          20  21
  16        3          3          19  20  21
  17        3          1          18
  18        3          2          20  21
  19        3          1          22
  20        3          1          22
  21        3          1          22
  22        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  0      1     0       0    0    0    0
  1      1     1       0    3    6    9
         2     5       0    3    4    8
         3     6       9    0    4    7
  2      0     6       9    0    7    5
         1     9       9    0    7    4
         2     9       0    7    7    4
  4      1     1       0    9    6    6
         2     4       0    9    6    3
         3     4       5    0    5    6
  5      1     1       0    5   10    8
         2     6       0    4    7    6
         3     8       0    3    2    4
  6      1     1       0    8    5    8
         2     4       0    8    4    4
         3     5       4    0    2    3
  7      1     7       0    5    8    9
         2    10       0    5    5    7
         3    10       8    0    4    7
  8      1     1       0    6   10    6
         2     1       5    0    9    7
         3     5       4    0    8    3
  9      1     4       5    0    6    7
         2     4       0    3    6    8
         3    10       0    2    5    6
 10      1     3       6    0    8    7
         2     4       0    7    8    4
         3     8       0    4    8    2
 11      1     1       0    8    7   10
         2     5       1    0    7    9
         3    10       0    1    6    9
 12      1     2       0    7    4    9
         2     3       6    0    4    8
         3     5       5    0    3    7
 13      1     4       7    0    5   10
         2     5       0    7    3    9
         3     8       0    5    1    9
 14      1     2       6    0    6    7
         2     3       0    4    6    6
         3     6       6    0    6    4
 15      1     2       6    0    9    7
         2     6       0    7    7    7
         3     7       0    6    3    5
 16      1     1       7    0    7    1
         2     2       0    5    7    1
         3     4       5    0    5    1
 17      1     2       0    9    6    5
         2     2      10    0    6    4
         3     5       0    8    4    2
 18      1     2       0    4    4    7
         2     5       2    0    4    6
         3     7       0    3    3    4
 19      1     2       0   10    6    5
         2     6       4    0    6    3
         3     9       0    7    5    2
 20      1     1       7    0    6    8
         2     2       0    6    6    8
         3     4       7    0    5    8
 21      1     2       0    4    5    7
         2     4       9    0    5    5
         3    10       9    0    5    3
 22      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
    8   10  121  131
************************************************************************

************************************************************************
file with basedata            : cm457_.bas
initial value random generator: 1988388651
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  130
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       16        4       16
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        4          3           6  10  12
   3        4          3          13  14  17
   4        4          3           5  10  14
   5        4          3           7   8   9
   6        4          1           8
   7        4          3          12  13  16
   8        4          2          11  17
   9        4          2          12  15
  10        4          2          11  13
  11        4          2          15  16
  12        4          1          17
  13        4          1          15
  14        4          1          16
  15        4          1          18
  16        4          1          18
  17        4          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       0    8    2   10
         2     5       9    0    2    8
         3     6       0    7    1    6
         4     6       8    0    1    6
  3      1     3       0    7    2    9
         2     4       0    5    2    7
         3     6       0    2    2    7
         4     9       5    0    1    5
  4      1     1       9    0    9    9
         2     5       6    0    8    7
         3     6       0    1    8    5
         4     6       5    0    8    4
  5      1     3       3    0    4    4
         2     4       3    0    3    4
         3     6       3    0    3    3
         4    10       1    0    2    2
  6      1     2       2    0    9   10
         2     2       0    8    9   10
         3     4       0    7    9    8
         4    10       2    0    6    8
  7      1     1       0    7    3    7
         2     2       8    0    2    7
         3     4       8    0    2    6
         4     6       0    6    1    6
  8      1     2       0    9    7   10
         2     5       0    6    6    9
         3     7       5    0    6    8
         4     8       0    4    3    5
  9      1     5       6    0    4    6
         2     5       0    7    4    6
         3     8       6    0    3    6
         4    10       4    0    3    5
 10      1     2       0    3    6    6
         2     4       0    3    6    5
         3    10       5    0    5    2
         4    10       0    2    5    2
 11      1     4       6    0    9   10
         2     5       0    9    8    8
         3     6       6    0    8    7
         4     7       5    0    5    3
 12      1     6       7    0    5    8
         2     7       0    5    5    8
         3     8       6    0    2    5
         4     9       6    0    1    5
 13      1     2       9    0    9    9
         2     2       0    4    8    9
         3     3       0    2    8    8
         4     6       8    0    4    7
 14      1     2       0    8    7    5
         2     5       5    0    6    4
         3     5       0    6    7    4
         4     8       4    0    6    3
 15      1     1       6    0    6    3
         2     4       0    3    5    3
         3     6       5    0    5    2
         4     8       5    0    3    2
 16      1     6       0    6    9    7
         2     8       9    0    8    6
         3     9       6    0    8    5
         4    10       5    0    6    3
 17      1     1       0    8    5    8
         2     2       8    0    5    8
         3     5       7    0    5    7
         4     7       0    8    4    7
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
    8    6   96  121
************************************************************************

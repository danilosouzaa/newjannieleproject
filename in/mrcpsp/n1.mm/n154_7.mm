************************************************************************
file with basedata            : cn154_.bas
initial value random generator: 1939192279
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  132
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  1   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       17       12       17
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           5  13
   3        3          3           9  12  17
   4        3          2           5  13
   5        3          3           6   7  12
   6        3          2           8  10
   7        3          2          11  14
   8        3          2           9  11
   9        3          2          14  16
  10        3          2          11  16
  11        3          2          15  17
  12        3          2          15  16
  13        3          2          14  17
  14        3          1          15
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1
------------------------------------------------------------------------
  1      1     0       0    0    0
  2      1     1       7    9    6
         2     4       7    7    6
         3     5       3    3    6
  3      1     3       5    6    6
         2     4       4    5    4
         3     9       2    3    3
  4      1     2       5    9    8
         2     3       3    6    8
         3     7       1    6    8
  5      1     3       9    6    8
         2     3      10    5    7
         3    10       8    5    4
  6      1     1       4    9    3
         2     7       4    8    3
         3     9       2    5    2
  7      1     8       7    6    2
         2     9       5    4    2
         3     9       6    5    1
  8      1     6       6    4    3
         2     7       5    3    2
         3    10       4    1    2
  9      1     1       9    5   10
         2     2       8    4    6
         3     3       8    3    4
 10      1     5       5    3    4
         2     6       5    2    2
         3     8       4    2    2
 11      1     1      10    8   10
         2     6       7    4   10
         3     8       3    1   10
 12      1     4       8    6    7
         2    10       8    2    6
         3    10       7    3    6
 13      1     2       6    7    7
         2     2       5   10    7
         3    10       1    4    6
 14      1     1       9    6    2
         2     5       9    3    1
         3     8       6    2    1
 15      1     1       9    5    8
         2     4       6    4    8
         3     9       5    4    7
 16      1     3       7    7    5
         2     6       6    7    4
         3     8       6    6    4
 17      1     3       7   10    5
         2     3       4    9    6
         3     9       1    9    4
 18      1     0       0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1
   17   17   89
************************************************************************

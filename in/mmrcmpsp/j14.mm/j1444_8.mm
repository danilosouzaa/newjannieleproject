************************************************************************
file with basedata            : md172_.bas
initial value random generator: 62065475
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  16
horizon                       :  109
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     14      0       21        3       21
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2          12  13
   3        3          3           5   7   8
   4        3          2           6   7
   5        3          3           9  10  11
   6        3          3           9  11  14
   7        3          1          10
   8        3          3          11  13  14
   9        3          2          13  15
  10        3          2          12  14
  11        3          1          12
  12        3          1          15
  13        3          1          16
  14        3          1          16
  15        3          1          16
  16        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       8    0   10   10
         2     2       0    2   10    7
         3     6       6    0    9    2
  3      1     5       7    0    3    2
         2     9       0    5    3    2
         3    10       0    4    1    2
  4      1     2      10    0    6    3
         2     4       4    0    6    2
         3     6       3    0    5    2
  5      1     3       0    7    7   10
         2     4       0    5    4    8
         3     7       0    4    2    6
  6      1     1       2    0    8    6
         2     9       0    3    5    5
         3    10       2    0    2    5
  7      1     1       0    6    9    6
         2     4       0    6    5    5
         3     7       0    1    3    5
  8      1     4       6    0    2    4
         2     6       0    7    2    3
         3     8       0    4    1    1
  9      1     2       9    0    4    6
         2     2       0    1    6    7
         3     3       9    0    4    4
 10      1     2       6    0    8    5
         2     3       3    0    8    5
         3     8       0    6    5    3
 11      1     8      10    0    6    7
         2     9       0    6    6    4
         3     9       9    0    4    6
 12      1     2       8    0    6   10
         2     2       0    4    8   10
         3     7       8    0    3   10
 13      1     2       7    0    8    8
         2     7       0    6    5    8
         3    10       0    3    4    8
 14      1     5       0    9    7    6
         2     5       8    0    8    7
         3     9       8    0    5    3
 15      1     2       0    9    9    8
         2     8       9    0    8    8
         3     9       5    0    8    7
 16      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   25   23   77   78
************************************************************************

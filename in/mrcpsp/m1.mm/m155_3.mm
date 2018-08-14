************************************************************************
file with basedata            : cm155_.bas
initial value random generator: 1934807374
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  67
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       37        3       37
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        1          2           6  15
   3        1          3           5   7  12
   4        1          3          10  11  13
   5        1          3           6  11  15
   6        1          1           8
   7        1          3           8   9  15
   8        1          1          10
   9        1          3          10  11  13
  10        1          2          16  17
  11        1          2          14  16
  12        1          2          13  16
  13        1          1          14
  14        1          1          17
  15        1          1          18
  16        1          1          18
  17        1          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     2       7    4    5    5
  3      1     5       9    7    7    4
  4      1     9       8    3    9    9
  5      1     4      10    6    1    6
  6      1     1       5    6    8    1
  7      1     5       1    3    5    5
  8      1     1       7    4    5    4
  9      1     2       8    3    6    4
 10      1     1       9    2    2    6
 11      1     5       9    8    6    4
 12      1     3       7   10    3    5
 13      1     8       3    8    2    3
 14      1     7       6    3    2    2
 15      1     2       5    4    9    7
 16      1     2       5    9    7   10
 17      1    10       5    5    4    5
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   22   19   81   80
************************************************************************

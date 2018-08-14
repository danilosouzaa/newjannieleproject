************************************************************************
file with basedata            : cm123_.bas
initial value random generator: 515476848
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  79
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       42       13       42
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        1          3           5  12  13
   3        1          1           8
   4        1          3           6   7   9
   5        1          1          15
   6        1          2           8  13
   7        1          2          11  14
   8        1          2          10  16
   9        1          3          11  12  13
  10        1          2          12  17
  11        1          3          15  16  17
  12        1          1          14
  13        1          3          14  16  17
  14        1          1          15
  15        1          1          18
  16        1          1          18
  17        1          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     5       6    6    8    0
  3      1     9       5    5    7    0
  4      1     3       3    8    0    3
  5      1     6       2    7    0    5
  6      1     1       4    8    0    3
  7      1     8       2    2    6    0
  8      1    10      10    6    0    3
  9      1     1       4    1    0   10
 10      1     6       8    7    0    8
 11      1     4       6    8    8    0
 12      1     6       8    6    0    7
 13      1     5       7    6    9    0
 14      1     1       8    2    0    9
 15      1    10       4    2    2    0
 16      1     3       2    1    8    0
 17      1     1       9    2    0    5
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   18   19   48   53
************************************************************************

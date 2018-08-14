************************************************************************
file with basedata            : me39_.bas
initial value random generator: 11829
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  20
horizon                       :  135
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  0   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     18      0       23        0       23
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          1           5
   3        3          3           6   7  16
   4        3          3           7  13  17
   5        3          3           6   7   9
   6        3          2           8  19
   7        3          2          14  15
   8        3          3          11  13  15
   9        3          3          10  11  16
  10        3          3          17  18  19
  11        3          1          12
  12        3          1          14
  13        3          1          14
  14        3          1          18
  15        3          1          18
  16        3          2          17  19
  17        3          1          20
  18        3          1          20
  19        3          1          20
  20        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2
------------------------------------------------------------------------
  1      1     0       0    0
  2      1     3       7    5
         2     4       6    4
         3     5       3    3
  3      1     2       9    5
         2     5       5    4
         3    10       4    4
  4      1     1       9    8
         2     1       8   10
         3     9       6    7
  5      1     5       7    9
         2     7       6    8
         3     9       6    6
  6      1     1       7    3
         2     2       5    3
         3     4       4    3
  7      1     3       8    6
         2     4       7    6
         3     9       6    4
  8      1     1       5    7
         2     5       3    7
         3     6       2    7
  9      1     5       7    9
         2     6       5    7
         3     7       2    3
 10      1     3       6    9
         2     6       4    9
         3     6       5    8
 11      1     5       8    4
         2     9       6    4
         3     9       7    3
 12      1     3       6    6
         2     4       4    6
         3     9       1    5
 13      1     1       5    9
         2     7       4    8
         3     9       3    7
 14      1     1       6    6
         2     4       5    6
         3     7       4    5
 15      1     1      10    8
         2     2       9    6
         3     5       9    4
 16      1     3       2    5
         2     9       1    2
         3    10       1    1
 17      1     2       1    7
         2     3       1    5
         3     4       1    4
 18      1     1      10    9
         2     6       8    8
         3     7       7    5
 19      1     5       6    5
         2     9       6    4
         3    10       5    3
 20      1     0       0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2
   21   20
************************************************************************

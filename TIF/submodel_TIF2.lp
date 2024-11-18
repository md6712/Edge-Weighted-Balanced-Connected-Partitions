\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1: 85 x(0,0) + 65 x(1,0) + 96 x(2,0)
Subject To
 OrderTrees(0)#0:      85 x(0,0) + 65 x(1,0) + 96 x(2,0) - 85 x(0,1)
                       - 65 x(1,1) - 96 x(2,1) >= 0
 PicknkV#1:            x(0,0) + x(1,0) + x(2,0) + x(0,1) + x(1,1) + x(2,1)  = 3
 _EachEdgeToTree(0)#2: x(0,0) + x(1,0) + x(2,0) <= 1
 _EachEdgeToTree(1)#3: x(0,1) + x(1,1) + x(2,1) <= 1
Bounds
 0 <= x(0,0) <= 1
 0 <= x(1,0) <= 1
 0 <= x(2,0) <= 1
 0 <= x(0,1) <= 1
 0 <= x(1,1) <= 1
 0 <= x(2,1) <= 1
Binaries
 x(0,0)  x(1,0)  x(2,0)  x(0,1)  x(1,1)  x(2,1) 
End

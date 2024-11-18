\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1: x1
Subject To
 ConsEachVertexInATree(0): y(0,0) + y(0,1)  = 1
 ConsEachVertexInATree(1): y(1,0) + y(1,1)  = 1
 ConsEachVertexInATree(2): y(2,0) + y(2,1)  = 1
 ConsEachVertexInATree(3): y(3,0) + y(3,1)  = 1
 ConsEachVertexInATree(4): y(4,0) + y(4,1)  = 1
 ConsEachVertexInATree(5): y(5,0) + y(5,1)  = 1
 ConsEachVertexInATree(6): y(6,0) + y(6,1)  = 1
 ConsEachVertexInATree(7): y(7,0) + y(7,1)  = 1
 ConsEachVertexInATree(8): y(8,0) + y(8,1)  = 1
 ConsEachVertexInATree(9): y(9,0) + y(9,1)  = 1
Bounds
 0 <= y(0,0) <= 1
 0 <= y(0,1) <= 1
 0 <= y(1,0) <= 1
 0 <= y(1,1) <= 1
 0 <= y(2,0) <= 1
 0 <= y(2,1) <= 1
 0 <= y(3,0) <= 1
 0 <= y(3,1) <= 1
 0 <= y(4,0) <= 1
 0 <= y(4,1) <= 1
 0 <= y(5,0) <= 1
 0 <= y(5,1) <= 1
 0 <= y(6,0) <= 1
 0 <= y(6,1) <= 1
 0 <= y(7,0) <= 1
 0 <= y(7,1) <= 1
 0 <= y(8,0) <= 1
 0 <= y(8,1) <= 1
 0 <= y(9,0) <= 1
 0 <= y(9,1) <= 1
Binaries
 y(0,0)  y(0,1)  y(1,0)  y(1,1)  y(2,0)  y(2,1)  y(3,0)  y(3,1)  y(4,0) 
 y(4,1)  y(5,0)  y(5,1)  y(6,0)  y(6,1)  y(7,0)  y(7,1)  y(8,0)  y(8,1) 
 y(9,0)  y(9,1) 
End

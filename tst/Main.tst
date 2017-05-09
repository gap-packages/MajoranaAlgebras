gap> a := (1,2)(3,4)(5,6);; b:=(1,4)(2,3);; c:=(2,4);;
gap> G := Group(a,b,c);;
gap> T2 := [a,b,c,a*b,(1,4)(2,3)(5,6),(1,2)(4,3),(1,3),(1,3)(5,6),(2,4)(5,6)];;
gap> res := MajoranaRepresentation(G,T2);;
gap> res[1][4];
[ "1A", "2B", "4A", "4A", "2B", "2A", "2A", "1A", "4A", "4A", "2B", "2A", 
  "1A", "2B", "2A", "2B", "2A", "1A", "2A", "2B", "1A" ]
gap> res[1][1];
"Error"
gap> res[1][2];
"Fusion of 0,0 eigenvectors does not hold"
gap> res[1][3];
[ [ 9, 1, 3 ], [ 9, 2, 3 ], [ 9, 3, 1 ], [ 9, 3, 2 ], [ 9, 3, 4 ], 
  [ 9, 3, 5 ], [ 9, 4, 3 ], [ 9, 4, 5 ], [ 9, 5, 3 ], [ 9, 5, 4 ] ]
gap> T := [ (1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (1,2)(3,4), (1,3)(2,4), (1,4)(2,3) ];;
gap> G := Group(T);;
gap> res := MajoranaRepresentation(G,T);;
gap> res[1][4];
[ "1A", "3C", "2A", "2A", "4B", "1A", "2A" ]
gap> res[1][1];
"Success"
gap> Size(res[1][6][1]);
9
gap> res[2][4];
[ "1A", "3A", "2A", "2A", "4B", "1A", "2A" ]
gap> res[2][1];
"Success"
gap> dim := Size(res[2][6][1]);
13
gap> u := [1..dim]*0;; u[7] := 1;;
gap> v := [1..dim]*0;; v[10] := 1;;
gap> MAJORANA_AlgebraProduct(u,v,res[2][6],res[2][8]);
[ 0, 0, 0, 0, 0, 0, 1/9, 0, 0, 5/64, 3/64, -1/16, -1/16 ]
gap> u := [1..dim]*0;; u[11] := 1;;
gap> MAJORANA_AlgebraProduct(u,v,res[2][6],res[2][8]);
[ 0, 0, 0, 0, 0, 0, 128/2025, -64/675, -64/675, 1/5, 1/5, -1/18, -1/18 ]

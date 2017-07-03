gap> T := [ (1,2), (3,4), (1,2)(5,6), (3,4)(5,6), (1,3)(2,4)(5,6), (1,4)(2,3)(5,6), (5,6), (1,3)(2,4), (1,4)(2,3), (1,2)(3,4)(5,6) ];;
gap> G := Group(T);;
gap> input := ShapesOfMajoranaRepresentation(G,T);;
gap> res := MajoranaRepresentation(input,1);;
gap> res[4];
[ "1A", "2A", "2A", "2B", "2A", "2A", "1A", "2A", "2B", "2A", "2A", "4A", 
  "4A", "1A", "2B", "2A", "4A", "4A", "1A", "2A", "2A", "1A", "2A", "2B", 
  "2A", "1A", "2B" ]
gap> res[1];
"Success"
gap> res[8][6];
[ [ 12 ], 
  [ [ -2/3, -2/3, -2/3, -2/3, -2/3, -2/3, 2/3, -2/3, -2/3, 2/3, 1, 1 ] ] ]
gap> T := [ (1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (1,2)(3,4), (1,3)(2,4), (1,4)(2,3) ];; G := Group(T);;  
gap> input := ShapesOfMajoranaRepresentation(G,T);;
gap> res := MajoranaRepresentation(input,1);;
gap> res[4];
[ "1A", "3C", "2A", "2A", "4B", "1A", "2A" ]
gap> res[1];
"Success"
gap> Size(res[6][1]);
9
gap> res := MajoranaRepresentation(input,2);;
gap> res[4];
[ "1A", "3A", "2A", "2A", "4B", "1A", "2A" ]
gap> res[1];
"Success"
gap> dim := Size(res[6][1]);
13
gap> u := [1..dim]*0;; u[7] := 1;;
gap> v := [1..dim]*0;; v[10] := 1;;
gap> MAJORANA_AlgebraProduct(u,v,res[6],res[8]);
[ 0, 0, 0, 0, 0, 0, 1/9, 0, 0, 5/64, 3/64, -1/16, -1/16 ]
gap> u := [1..dim]*0;; u[11] := 1;;
gap> MAJORANA_AlgebraProduct(u,v,res[6],res[8]);
[ 0, 0, 0, 0, 0, 0, 128/2025, -64/675, -64/675, 1/5, 1/5, -1/18, -1/18 ]

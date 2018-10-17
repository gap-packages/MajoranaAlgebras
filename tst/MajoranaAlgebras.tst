##
## Test inner product solution functions
##

gap> mat := SparseMatrix( 1, 5, [ [ 1, 4 ] ], [ [ 1, -1 ] ], Rationals );;
gap> vec := SparseMatrix( 1, 1, [ [ 1 ] ], [ [ 7123/518400 ] ], Rationals );;
gap> unknowns := [1..5];;
gap> innerproducts := [false, 289/57600, 1321/518400, false, 23/5184 ];;
gap> MAJORANA_RemoveKnownInnProducts(mat, vec, unknowns, innerproducts);
rec( mat := <a 1 x 2 sparse matrix over Rationals>, unknowns := [ 1, 4 ], 
  vec := <a 1 x 1 sparse matrix over Rationals> )
gap> eq := [ SparseMatrix( 1, 3, [ [ 1 ] ], [ [ -1 ] ], Rationals ), SparseMatrix( 1, 1, [ [ 1 ] ], [ [ -1/8192 ] ], Rationals ) ];;
gap> mat := SparseMatrix( 0, 3, [  ], [  ], Rationals );;
gap> vec := SparseMatrix( 0, 1, [  ], [  ], Rationals );;
gap> unknowns := [ 1, 2, 3 ];;
gap> innerproducts := [ false, false, false ];;
gap> MAJORANA_SingleInnerSolution( eq, mat, vec, unknowns, innerproducts );;
gap> innerproducts;
[ 1/8192, false, false ]
gap> mat := SparseMatrix( 1, 1, [ [ 1 ] ], [ [ 1 ] ], Rationals );;
gap> vec := SparseMatrix( 1, 1, [ [ 1 ] ], [ [ 1/2 ] ], Rationals );;
gap> unknowns := [ 1 ];;
gap> innerproducts := [ false ];;
gap> MAJORANA_SolutionInnerProducts(mat, vec, unknowns, innerproducts);;
gap> innerproducts;
[ 1/2 ]

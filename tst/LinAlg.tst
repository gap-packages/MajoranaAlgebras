
##
## Test solve system
##
gap> system := rec();;
gap> system.mat := SparseMatrix( 5, 6, [ [ 1, 2 ], [ 3, 4, 5, 6 ], [ 3, 4, 5, 6 ], [ 5, 6 ], [ 3, 4 ] ], [ [ -1, 1 ], [ 1, 1, 1, 1 ], [ -1,-1, 1, 1 ], [ -1, 1 ], [ -1, 1 ] ], Rationals );;
gap> system.vec := SparseMatrix( 5, 13, [ [ 12, 13 ], [ 5 ], [ 10, 11, 12, 13 ], [ 12, 13 ], [ 10, 11 ] ], [ [ -1/32, 1/32 ], [ 4/9 ], [ -1/4, -1/4, 1/4, 1/4 ], [ -1/32, 1/32 ], [ -1/32, 1/32 ] ], Rationals );;
gap> system.unknowns := [ [ 1, 12 ], [ 1, 13 ], [ 5, 10 ], [ 5, 11 ], [ 5, 12 ], [ 5, 13 ] ];;
gap> MAJORANA_SolveSystem(system);;
gap> system.mat!.indices;
[ [ 1, 2 ], [ 3 ], [ 4 ], [ 5 ], [ 6 ] ]
gap> system.mat!.entries;
[ [ 1, -1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ] ]
gap> system.vec!.indices;
[ [ 12, 13 ], [ 5, 10, 11, 12, 13 ], [ 5, 10, 11, 12, 13 ], 
  [ 5, 10, 11, 12, 13 ], [ 5, 10, 11, 12, 13 ] ]
gap> system.vec!.entries;
[ [ 1/32, -1/32 ], [ 1/9, 5/64, 3/64, -1/16, -1/16 ], 
  [ 1/9, 3/64, 5/64, -1/16, -1/16 ], [ 1/9, -1/16, -1/16, 5/64, 3/64 ], 
  [ 1/9, -1/16, -1/16, 3/64, 5/64 ] ]
gap> system.solutions{[1,2]};
[ fail, fail ]
gap> system.solutions[3]!.indices;
[ [ 5, 10, 11, 12, 13 ] ]
gap> system.solutions[3]!.entries;
[ [ 1/9, 5/64, 3/64, -1/16, -1/16 ] ]

##
## Test LDLT decomposition and positive definite
##
gap> mat := SparseMatrix( 6, 6, [ [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 6 ] ], [ [ 1, 3/128, 3/128, 3/128, 3/128 ], [ 3/128, 1, 3/128, 3/128, 3/128 ], [ 3/128, 3/128, 1, 3/128, 3/128 ], [ 3/128, 3/128, 3/128, 1, 3/128 ], [ 3/128, 3/128, 3/128, 3/128, 1 ], [ 875/524288 ] ], Rationals );;
gap> mat := ConvertSparseMatrixToMatrix(mat);;
gap> MAJORANA_PositiveDefinite(mat);
1
gap> L := MAJORANA_LDLTDecomposition(mat);;
gap> L[1];
[ [ 1, 0, 0, 0, 0, 0 ], [ 3/128, 1, 0, 0, 0, 0 ], 
  [ 3/128, 3/131, 1, 0, 0, 0 ], [ 3/128, 3/131, 3/134, 1, 0, 0 ], 
  [ 3/128, 3/131, 3/134, 3/137, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ]
gap> L[2];
[ [ 1, 0, 0, 0, 0, 0 ], [ 0, 16375/16384, 0, 0, 0, 0 ], 
  [ 0, 0, 8375/8384, 0, 0, 0 ], [ 0, 0, 0, 17125/17152, 0, 0 ], 
  [ 0, 0, 0, 0, 4375/4384, 0 ], [ 0, 0, 0, 0, 0, 875/524288 ] ]
gap> mat := SparseIdentityMatrix(5, Rationals);;
gap> mat := ConvertSparseMatrixToMatrix(mat);;
gap> L := MAJORANA_LDLTDecomposition(mat);;
gap> L[1];
[ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
  [ 0, 0, 0, 0, 1 ] ]
gap> L[2];
[ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
  [ 0, 0, 0, 0, 1 ] ]
gap> MAJORANA_PositiveDefinite( NullMat(6, 6) );
0

##
## Test iterator for sparse matrices and _IsRowOfSparseMatrix
##
gap> mat := SparseMatrix( 6, 6, [ [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 1, 2, 3, 4, 5 ], [ 6 ] ], [ [ 1, 3/128, 3/128, 3/128, 3/128 ], [ 3/128, 1, 3/128, 3/128, 3/128 ], [ 3/128, 3/128, 1, 3/128, 3/128 ], [ 3/128, 3/128, 3/128, 1, 3/128 ], [ 3/128, 3/128, 3/128, 3/128, 1 ], [ 875/524288 ] ], Rationals );;
gap> for v in Iterator(mat) do if not _IsRowOfSparseMatrix(mat, v ) then Error(); fi; od;

##
## Test reversed echelon mat transformation
##
gap> mat := ReversedEchelonMatDestructive(mat);;
gap> mat.heads;
[ 6, 5, 4, 3, 2, 1 ]
gap> mat.vectors!.indices;
[ [ 6 ], [ 5 ], [ 4 ], [ 3 ], [ 2 ], [ 1 ] ]
gap> mat.vectors!.entries;
[ [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ] ]

##
## Test remove mat with heads
##
gap> mat := SparseMatrix( 1, 21, [ [ 1, 5, 9, 10, 15, 16 ] ], [ [ 3/128, 3/128, -1/128, -1/128, -1/128, 1 ] ], Rationals );;
gap> null := rec( heads := [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ], vectors := SparseMatrix( 1, 21, [ [ 16, 17, 18, 19, 20, 21 ] ], [ [ 1, 1, 1, 1, 1, 1 ] ], Rationals ) );
rec( heads := [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 
     ], vectors := <a 1 x 21 sparse matrix over Rationals> )
gap> mat := RemoveMatWithHeads(mat, null);;
gap> mat!.indices;
[ [ 1, 5, 9, 10, 15, 16 ] ]
gap> mat!.entries;
[ [ 3/128, 3/128, -1/128, -1/128, -1/128, 1 ] ]

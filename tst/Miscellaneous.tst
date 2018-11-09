gap> SetInfoLevel( InfoMajorana, 0);

##
## Set up infinite family example
##
gap> ex := min3gen9();;
gap> rep := MajoranaRepresentation(ex, 2);;
gap> MAJORANA_Dimension(rep);
7
gap> rep := MajoranaRepresentation(ex, 1);;
gap> rep.innerproducts[19] := 100;;
gap> MAJORANA_MainLoop(rep);;
gap> NClosedMajoranaRepresentation(rep);;
gap> MAJORANA_Dimension(rep);
12
gap> MajoranaAlgebraTest(rep);
true

##
## Test subalgebra
##
gap> vecs := SparseMatrix(2, 15, [[7], [1,2]], [[1], [1, 1]], Rationals);;
gap> subalg := MAJORANA_Subalgebra(vecs, rep);;
gap> subalg!.indices;
[ [ 1, 2 ], [ 3, 4 ], [ 7 ] ]
gap> subalg!.entries;
[ [ 1, 1 ], [ 1, 1 ], [ 1 ] ]

##
## Test IsJordanAlgebra
##
gap> MAJORANA_IsJordanAlgebra(subalg, rep);
true
gap> MAJORANA_IsJordanAlgebra(SparseIdentityMatrix(15, Rationals), rep);
false

##
## Test AdjointAction
##
gap> axis := SparseMatrix(1, 15, [[1]], [[1]], Rationals);;
gap> basis := SparseIdentityMatrix(15, Rationals);;
gap> adj := MAJORANA_AdjointAction(axis, basis, rep);;
gap> adj := ConvertSparseMatrixToMatrix(adj);;
gap> AsSet(Eigenvalues(Rationals, adj)) = AsSet([1, 0, 1/4, 1/32]);
true

##
## Test ConvertToBasis
##
gap> vec := SparseMatrix( 1, 15, [[1, 2, 7]], [[1, 1, -1]], Rationals);;
gap> vec := MAJORANA_ConvertToBasis(subalg, vec);;
gap> vec!.indices;
[ [ 1, 3 ] ]
gap> vec!.entries;
[ [ 1, -1 ] ]

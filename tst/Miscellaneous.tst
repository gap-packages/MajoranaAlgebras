gap> SetInfoLevel( InfoMajorana, 0);

##
## Set up infinite family example
##
gap> ex := MAJORANA_Example_min3gen9();;
gap> rep := MajoranaRepresentation(ex, 2);;
gap> MAJORANA_Dimension(rep);
7
gap> rep := MajoranaRepresentation(ex, 1);;
gap> rep.innerproducts[19] := 1/12;;
gap> MAJORANA_MainLoop(rep);;
gap> NClosedMajoranaRepresentation(rep);;
gap> MAJORANA_Dimension(rep);
12
gap> MajoranaAlgebraTest(rep);
true

##
## Test eigenvectors
##
gap> u := SparseMatrix( 1, 15, [[2]], [[1]], Rationals);;
gap> evecs := MAJORANA_Eigenvectors( 2, 0, rep);;
gap> for v in Iterator( evecs ) do prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup); if prod <> (0)*v then Error(); fi; od;
gap> evecs := MAJORANA_Eigenvectors( 2, 1/4, rep);;
gap> for v in Iterator( evecs ) do prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup); if prod <> (1/4)*v then Error(); fi; od;
gap> evecs := MAJORANA_Eigenvectors( 2, 1/32, rep);;
gap> for v in Iterator( evecs ) do prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup); if prod <> (1/32)*v then Error(); fi; od;

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
gap> basis := MAJORANA_Basis(rep);;
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

##
## Test eigenvectors (again)
##
gap> rep := MAJORANA_DihedralAlgebras("3C");;
gap> List( [1, 0, 1/4, 1/32], i -> Nrows(MAJORANA_Eigenvectors( 1, i, rep)) );
[ 1, 1, 0, 1 ]

##
## Test main funcs for A5
##
gap> a5 := MAJORANA_Example_A5();;
gap> Length(a5.shapes) = 4;
true
gap> setup := MAJORANA_SetUp(a5, 1, rec( axioms := "AllAxioms"));;
gap> Set([ "1A", "2B", "3C", "5A", "5A" ]) = Set(setup.shape);
true
gap> MajoranaAlgebraTest(setup);
true
gap> setup := MAJORANA_SetUp(a5, 2, rec( axioms := "NoAxioms"));;
gap> MajoranaAlgebraTest(setup);
true

##
##  Test main funcs for S5
##
gap> s5 := MAJORANA_Example_S5();;
gap> setup := MAJORANA_SetUp(s5, 1, rec( axioms := "AllAxioms"));;
gap> MAJORANA_TestSetup(setup);
true
gap> s5 := MAJORANA_Example_S5();;
gap> setup := MAJORANA_SetUp(s5, 1, rec( axioms := "NoAxioms"));;
gap> MAJORANA_TestSetup(setup);
true

##
## Test main funcs for A7
##
gap> a7 := MAJORANA_Example_A7();;
gap> rep := MAJORANA_SetUp(a7, 2, rec( axioms := "AllAxioms") );;
gap> Size(rep.setup.coords);
406
gap> MajoranaAlgebraTest(rep);
true
gap> gens := List(rep.setup.pairconjelts, x -> List(x, AbsInt));;
gap> gens := List(gens, PermList);;
gap> Size(Group(gens)) = Size(rep.group);
true

##
## Test FindEmbedding
##
gap> MAJORANA_FindEmbedding(rep, MAJORANA_DihedralAlgebras("5A"), [ 1, 5 ]);
[ 1, 5, 9, 15, 10, 141 ]
gap> MAJORANA_FindEmbedding(rep, MAJORANA_DihedralAlgebras("2A"), [ 1, 2 ]);
[ 1, 2, 3 ]

##
## Test ExtendPerm
##
gap> perm := ShallowCopy( rep.setup.pairconjelts[2]{[1..105]} );;
gap> MAJORANA_ExtendPerm(perm, rep);
gap> perm = rep.setup.pairconjelts[2];
true

##
## Test main funcs for 2^3
##
gap> T := [ (1,2), (3,4), (5,6) ];; G := Group(T);;
gap> shapes := ShapesOfMajoranaRepresentation(G, T);;
gap> Length(shapes.shapes) = 8;
true
gap> MAJORANA_RemoveDuplicateShapes(shapes);
gap> Length(shapes.shapes) = 4;
true

##
## Test main funcs for S4
##
gap> s4 := MAJORANA_Example_S4T1();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> Size(setup.shapes);
2
gap> s4 := MAJORANA_Example_S4T2();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> Size(setup.shapes);
2
gap> s5 := MAJORANA_Example_S5();;
gap> setup := ShapesOfMajoranaRepresentation(s5.group, s5.involutions);;
gap> Size(setup.shapes);
2
gap> ex := MAJORANA_Example_min3gen9();;
gap> shapes := ShapesOfMajoranaRepresentationAxiomM8(ex.group, ex.involutions);;
gap> Size(setup.shapes);
2

##
## Test shape for a group that isn't a permutation group
##
gap> G := SU(3,2);;
gap> cc := ConjugacyClasses(G);;
gap> cc := Filtered(cc, x -> Order(Representative(x)) = 2);;
gap> T := AsList(cc[1]);;
gap> MAJORANA_IsSixTranspositionGroup(G, T);
false
gap> G := Group(T);;
gap> MAJORANA_IsSixTranspositionGroup(G, T);
true
gap> ex := ShapesOfMajoranaRepresentation(G,T);;
gap> Size(ex.shapes);
16
gap> ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
gap> Size(ex.shapes);
16

##
## Test MAJORANA_IsSixTranspositionGroup
##
gap> G := AlternatingGroup(5);;
gap> T := AsList(ConjugacyClass(G, (1,2,3)));;
gap> MAJORANA_IsSixTranspositionGroup(G, T);
false
gap> T := [(1,2)(3,4), (1,2)(3,5), (1,3)(4,5)];;
gap> MAJORANA_IsSixTranspositionGroup(G, T);
false
gap> G := DihedralGroup(14);;
gap> T := AsList(ConjugacyClass(G, G.1));;
gap> MAJORANA_IsSixTranspositionGroup(G,T);
false

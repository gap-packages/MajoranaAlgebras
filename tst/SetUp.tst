##
## Test main funcs for A5
##
gap> a5 := A5();;
gap> Length(a5.shapes) = 4;
true
gap> setup := MAJORANA_SetUp(a5, 1, "AllAxioms");;
gap> Set([ "1A", "2B", "3C", "5A", "5A" ]) = Set(setup.shape);
true
gap> MajoranaAlgebraTest(setup);
true
gap> setup := MAJORANA_SetUp(a5, 2, "NoAxioms");;
gap> MajoranaAlgebraTest(setup);
true
##
##  Test main funcs for S5
##
gap> s5 := S5();;
gap> setup := MAJORANA_SetUp(s5, 1, "AllAxioms");;
gap> Set([ "1A", "3A", "2A", "2A", "4B", "6A", "1A", "2A", "3A", "5A" ]) = Set(setup.shape);
true
gap> s5 := S5();;
gap> setup := MAJORANA_SetUp(s5, 1, "NoAxioms");;
gap> Set([ "1A", "3A", "2A", "2A", "4B", "6A", "1A", "2A", "3A", "5A" ]) = Set(setup.shape);
true
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
gap> s4 := S4T1();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> s4 := S4T2();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> s5 := S5();;
gap> setup := ShapesOfMajoranaRepresentation(s5.group, s5.involutions);;
gap> ex := min3gen9();;
gap> shapes := ShapesOfMajoranaRepresentationAxiomM8(ex.group, ex.involutions);;
##
## Test main funcs for A7
##
gap> a7 := A7();;
gap> setup := MAJORANA_SetUp(a7, 2, "AllAxioms");;
gap> Size(rep.setup.coords);
406
gap> MajoranaAlgebraTest(setup);
true
##

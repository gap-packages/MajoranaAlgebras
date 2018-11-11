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
gap> Size(setup.shapes);
2
gap> s4 := S4T2();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> Size(setup.shapes);
2
gap> s5 := S5();;
gap> setup := ShapesOfMajoranaRepresentation(s5.group, s5.involutions);;
gap> Size(setup.shapes);
2
gap> ex := min3gen9();;
gap> shapes := ShapesOfMajoranaRepresentationAxiomM8(ex.group, ex.involutions);;
gap> Size(setup.shapes);
2

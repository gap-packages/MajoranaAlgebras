# 
gap> a5 := A5();;
gap> setup := MAJORANA_SetUp(a5, 1, "AllAxioms");;
gap> Set([ "1A", "2B", "3C", "5A", "5A" ]) = Set(setup.shape);
true
gap> MajoranaAlgebraTest(setup);
true
gap> setup := MAJORANA_SetUp(a5, 2, "NoAxioms");;
gap> MajoranaAlgebraTest(setup);
true
gap> s5 := S5();;
gap> setup := MAJORANA_SetUp(s5, 1, "AllAxioms");;
gap> Set([ "1A", "3A", "2A", "2A", "4B", "6A", "1A", "2A", "3A", "5A" ]) = Set(setup.shape);
true
gap> s5 := S5();;
gap> setup := MAJORANA_SetUp(s5, 1, "NoAxioms");;
gap> Set([ "1A", "3A", "2A", "2A", "4B", "6A", "1A", "2A", "3A", "5A" ]) = Set(setup.shape);
true
gap> T := [ (1,2), (3,4), (5,6) ];; G := Group(T);;
gap> shapes := ShapesOfMajoranaRepresentation(G, T);;
gap> origshapes := ShallowCopy(shapes);;
gap> MAJORANA_RemoveDuplicateShapes(shapes);
gap> IsSubset(Set(shapes.shapes), Set(origshapes.shapes));
true
gap> s4 := S4T1();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> s4 := S4T2();;
gap> setup := ShapesOfMajoranaRepresentationAxiomM8(s4.group, s4.involutions);;
gap> s5 := S5();;
gap> setup := ShapesOfMajoranaRepresentation(s5.group, s5.involutions);;
gap> ex := min3gen9();;                                                    
gap> shapes := ShapesOfMajoranaRepresentationAxiomM8(ex.group, ex.involutions);;

#
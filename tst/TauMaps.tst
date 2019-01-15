gap> SetInfoLevel(InfoMajorana, 0);
gap> taumaps := [(3,4)(5,6), (3,4)(5,6), (1,2)(5,6), (1,2)(5,6), (1,2)(3,4), (1,2)(3,4)];;
gap> ex := TauMapsShapesOfMajoranaRepresentation(taumaps);;
gap> rep := MajoranaRepresentation(ex, 1, rec(taumaps := true, axioms := "NoAxioms"));;
gap> rep.innerproducts[12] := 123456;;
gap> MAJORANA_MainLoop(rep);;
gap> NClosedMajoranaRepresentation(rep);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
12
gap> ex := MAJORANA_Example_A5();;
gap> taumaps := [];;
gap> hom := ActionHomomorphism(ex.group, ex.involutions);;
gap> for t in ex.involutions do Add(taumaps, Image(hom, t)); od;
gap> ex := TauMapsShapesOfMajoranaRepresentation(taumaps);;
gap> rep := TauMapsMajoranaRepresentation(ex, 1);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
21

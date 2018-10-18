gap> SetInfoLevel( InfoMajorana, 0);

##
## Test funcs for infinite family example
##
gap> tau := [(3,4)(5,6), (3,4)(5,6), (1,2)(5,6), (1,2)(5,6), (1,2)(3,4), (1,2)(3,4)];;
gap> ex := TauShapesOfMajoranaRepresentation(tau);;
gap> rep := TauMapMajoranaRepresentation(ex, 2);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
7
gap> rep := TauMapMajoranaRepresentation(ex, 1);;
gap> MAJORANA_IsComplete(rep);
false
gap> pos := rep.setup.coordmap[[3,5]];;
gap> k := rep.setup.pairorbit[1][pos];
19
gap> rep.innerproducts[k] := 100;;
gap> MAJORANA_MainLoop(rep);;
gap> NClosedMajoranaRepresentation(rep);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
0

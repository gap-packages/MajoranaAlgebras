gap> SetInfoLevel( InfoMajorana, 0 );

##
## Test main funcs for A5
##
gap> ex := A5();;
gap> rep := MajoranaRepresentation(ex, 1);;
gap> MAJORANA_IsComplete(rep);
true
gap> MajoranaAlgebraTest(rep);
true
gap> MAJORANA_Dimension(rep);
21
gap> rep := MajoranaRepresentation(ex, 1, "NoAxioms");;
gap> MajoranaAlgebraTest(rep);
true
gap> MAJORANA_Dimension(rep);
21
gap> MAJORANA_IsComplete(rep);
true
gap> rep := MajoranaRepresentation(ex, 4);;
gap> MAJORANA_IsComplete(rep);
true
gap> MajoranaAlgebraTest(rep);
true
gap> MAJORANA_Dimension(rep);
26
gap> MAJORANA_TestOrthogonality(rep);
true
gap> MAJORANA_TestEvecs(rep);
true
gap> rep := MajoranaRepresentation(ex, 3);;

##
## Test Axiom M2 and positive definiteness on S4
##
gap> ex := S4T1();;
gap> rep := MajoranaRepresentation(ex, 1);;
gap> MAJORANA_Dimension(rep);
12
gap> MAJORANA_TestAxiomM2(rep);
true
gap> MAJORANA_TestPositiveDefiniteForm(rep);
true

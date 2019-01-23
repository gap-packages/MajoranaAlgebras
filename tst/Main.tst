gap> SetInfoLevel( InfoMajorana, 0 );

##
## Test main funcs for A5
##
gap> ex := MAJORANA_Example_A5();;
gap> rep := MajoranaRepresentation(ex, 1);;
gap> MAJORANA_IsComplete(rep);
true
gap> MajoranaAlgebraTest(rep);
true
gap> MAJORANA_Dimension(rep);
21
gap> rep := MajoranaRepresentation(ex, 1, rec( axioms := "NoAxioms") );;
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
gap> MAJORANA_TestEvecs(rep);
true
gap> rep := MajoranaRepresentation(ex, 1, rec(form := false));;
gap> MAJORANA_Dimension(rep);
21
gap> MajoranaAlgebraTest(rep);
true

##
## Test an A6 example that gives no algebra
##
gap> ex := MAJORANA_Example_A6();;
gap> rep := MajoranaRepresentation(ex, 2);;
gap> MAJORANA_Dimension(rep);
0

##
## Test Axiom M2 and positive definiteness and fusion on S4
##
gap> ex := MAJORANA_Example_S4T1();;
gap> rep := MajoranaRepresentation(ex, 1);;
gap> MAJORANA_Dimension(rep);
12
gap> MAJORANA_TestAxiomM2(rep);
true
gap> MAJORANA_TestInnerProduct(rep);
true
gap> MAJORANA_TestFusion(rep);
true

##
## Test a 3-closed example with no form
##
gap> ex := MAJORANA_Example_S4T2();;
gap> rep := MajoranaRepresentation(ex, 3, rec( form := false ));;
gap> NClosedMajoranaRepresentation(rep);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
13

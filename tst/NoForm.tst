##
## Test a 3-closed example with no form
##
gap> ex := S4T2();;
gap> rep := MajoranaRepresentationNoForm(ex, 3);;
gap> NClosedMajoranaRepresentationNoForm(rep);;
gap> MAJORANA_IsComplete(rep);
true
gap> MAJORANA_Dimension(rep);
13

##
## Test IntersectEigenspaces
##
gap> ex := S4T1();;
gap> rep := MAJORANA_SetUp(ex, 2, "AllAxioms");;
gap> MAJORANA_IntersectEigenspaces(rep);;
gap> MAJORANA_Dimension(rep);
0

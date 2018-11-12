gap> SetInfoLevel(InfoMajorana, 0);
gap> SetInfoLevel(InfoPerformance, 0);

##
## Test MaximalSubgps func
##
gap> ex := A6();;
gap> rep := MAJORANA_SetUp(ex, 1, rec());;
gap> MAJORANA_MaximalSubgps(rep, rec());;
gap> MAJORANA_TestAxiomM1(rep);

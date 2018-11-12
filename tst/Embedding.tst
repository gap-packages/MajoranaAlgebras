gap> SetInfoLevel(InfoMajorana, 0);
gap> SetInfoLevel(InfoPerformance, 0);

##
## Test MaximalSubgps func
##
gap> ex := A6();;
gap> rep := MajoranaRepresentation(ex, 1, rec( embeddings := true));;
gap> MAJORANA_Dimension(rep);
70

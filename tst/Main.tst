gap> SetInfoLevel(InfoPackageLoading,0);
gap> LoadPackage("MajoranaAlgebras");;
gap> SetInfoLevel(InfoMajorana, 0);
gap> ex := A5();;
gap> rep := MajoranaRepresentation(ex,1);;
gap> MajoranaAlgebraTest(rep); 
true
gap> rep := MajoranaRepresentation(ex,2);;
gap> MajoranaAlgebraTest(rep);            
true

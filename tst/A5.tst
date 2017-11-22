# A5
gap> G := AlternatingGroup(5);;
gap> T := Filtered(Elements(G), x -> Order(x) = 2);;
gap> rep := MajoranaRepresentation(G, T);;
gap> for r in rep do Print(r[1],"\n",r[2], "\n"); od;
[ "1A", "2A", "2A", "3C", "5A", "5A" ]
Success
[ "1A", "2A", "2A", "3A", "5A", "5A" ]
Success

gap> for type in RecNames(MAJORANA_DihedralAlgebras) do MajoranaAlgebraTest(MAJORANA_DihedralAlgebras.(type)); od;
gap> List( RecNames(MAJORANA_DihedralAlgebras), type -> MAJORANA_Dimension(MAJORANA_DihedralAlgebras.(type)));
[ 3, 2, 4, 3, 5, 5, 6, 8 ]

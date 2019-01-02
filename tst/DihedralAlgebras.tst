gap> types := ["2A", "2B", "3A", "3C", "4A", "4B", "5A", "6A"];;
gap> List( types, type -> MAJORANA_Dimension(MAJORANA_DihedralAlgebras(type)));
[ 3, 2, 4, 3, 5, 5, 6, 8 ]

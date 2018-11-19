gap> examples := [ S4T1, S4T2, A5, S5, L32, A6, S6, 3A6, 3S6, S4S3A7, A7, S7, 3A7, 3S7, PSL211, L33, S5S3A8  ];;
gap> for func in examples do ex := func(); Display(Size(ex.shapes)); od;
8
4
4
1
2
4
2
4
2
8
2
2
2
2
1
2
2

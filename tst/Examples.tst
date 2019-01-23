gap> examples := [  MAJORANA_Example_S4T1,
>                   MAJORANA_Example_S4T2,
>                   MAJORANA_Example_A5,
>                   MAJORANA_Example_S5,
>                   MAJORANA_Example_L32,
>                   MAJORANA_Example_A6,
>                   MAJORANA_Example_S6,
>                   MAJORANA_Example_3A6,
>                   MAJORANA_Example_3S6,
>                   MAJORANA_Example_S4S3A7,
>                   MAJORANA_Example_A7,
>                   MAJORANA_Example_S7,
>                   MAJORANA_Example_3A7,
>                   MAJORANA_Example_3S7,
>                   MAJORANA_Example_PSL211,
>                   MAJORANA_Example_L33,
>                   MAJORANA_Example_S5S3A8  ];;
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

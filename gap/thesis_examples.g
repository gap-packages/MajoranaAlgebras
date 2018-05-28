InfoThesis := NewInfoClass("InfoThesis");
SetInfoLevel(InfoThesis, 100);
SetInfoLevel(InfoMajorana, 0);
SetInfoLevel(InfoPerformance, 0);

##### 2 x D8

a := (1,2); b := (1,2)(5,6); c := (1,3)(2,4)(5,6);
G := Group(a, b, c);;

T1 := [a, a^c, b, b^c, c, c^a, a*b, (a*c)^2];
T2 := [a, a^c, b, b^c, c, c^a, a*b, a*c*b*c];
T4 := [a, a^c, b, b^c, c, c^a, a*b, a*b*c, (a*b*c)^a, a*c*b*c];

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T1);
rep_2XD8_T1 := MajoranaRepresentation(ex, 1, "AxiomM8");;
MajoranaAlgebraTest(rep_2XD8_T1);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x D8, T1, dim ", MAJORANA_Dimension(rep_2XD8_T1) ) ); 

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T2);
rep_2XD8_T2 := MajoranaRepresentation(ex, 1, "AxiomM8");;
MajoranaAlgebraTest(rep_2XD8_T2);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x D8, T2, dim ", MAJORANA_Dimension(rep_2XD8_T2) ) ); 

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T4);
rep_2XD8_T4 := MajoranaRepresentation(ex, 1, "AxiomM8");;
MajoranaAlgebraTest(rep_2XD8_T4);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x D8, T4, dim ", MAJORANA_Dimension(rep_2XD8_T4) ) ); 

##### 2^2 x S3

a := (1,2); b := (3,4); c := (4,5)(6,7);
G := Group(a,b,c);

T := [a, b, b^c, b^(c*b), c, c^b, c^(b*c), a*b, (a*b)^c, (a*b)^(c*b), (b*c)^3, (a*b*c)^3];

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T);
rep_22xS3 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_22xS3);
Info(InfoThesis, 50, STRINGIFY("Group: 2^2 x S3, dim ", MAJORANA_Dimension(rep_22xS3) ) ); 

##### 242

a := (1,2)(3,4); b := (5,6)(7,8); c := (1,3)(5,7);
G := Group(a,b,c);

T1 := [a, a^c, b, b^c, c, c^a, c^b, c^(a*b), a*b, (a*b)^c];;
T2 := [a, a^c, b, b^c, c, c^a, c^b, c^(a*b), a*b, (a*b)^c, (a*c)^2, (b*c)^2, a*c*b*c, c*a*c*b];;

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T1);
rep_242_T1 := MajoranaRepresentation(ex, 1, "AxiomM8");;
NClosedMajoranaRepresentation(rep_242_T1);
MajoranaAlgebraTest(rep_242_T1);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:2, T1, dim ", MAJORANA_Dimension(rep_242_T1) ) ); 

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T2);
rep_242_T2 := MajoranaRepresentation(ex, 1, "AxiomM8");;
MajoranaAlgebraTest(rep_242_T2);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:2, T2, dim ", MAJORANA_Dimension(rep_242_T2) ) ); 

##### S3 x S3

G := Group((1,2),(1,3),(4,5),(4,6));
T := Filtered(G, x -> Order(x) = 2);

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
rep_S3S3 := MajoranaRepresentation(ex, 2, "AxiomM8");
MajoranaAlgebraTest(rep_S3S3);
Info(InfoThesis, 50, STRINGIFY("Group: S3 x S3, dim ", MAJORANA_Dimension(rep_S3S3) ) ); 

##### 2 x S4 

G := Group((1,2), (1,2,3,4), (5,6));;
T1 := [ (1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(1,2)(3,4),(1,3)(2,4),(1,4)(2,3),
        (1,2)(5,6),(1,3)(5,6),(1,4)(5,6),(2,3)(5,6),(2,4)(5,6),(3,4)(5,6), (5,6)];
T2 := [ (1,2),(1,3),(1,4),(2,3),(2,4),(3,4),
        (1,2)(5,6),(1,3)(5,6),(1,4)(5,6),(2,3)(5,6),(2,4)(5,6),(3,4)(5,6), 
        (1,2)(3,4)(5,6), (1,3)(2,4)(5,6), (1,4)(2,3)(5,6), (5,6) ];
        
ex := ShapesOfMajoranaRepresentationAxiomM8(G, T1);
rep_2S4_T1 := MajoranaRepresentation(ex, 1);
MajoranaAlgebraTest(rep_2S4_T1);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x S4, T1, dim ", MAJORANA_Dimension(rep_2S4_T1) ) ); 

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T2);
rep_2S4_T2 := MajoranaRepresentation(ex, 1);
MajoranaAlgebraTest(rep_2S4_T2);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x S4, T2, dim ", MAJORANA_Dimension(rep_2S4_T2) ) ); 
        
##### 2wr2

a := (1,2)(3,4); b := (1,3)(2,4)(5,6)(7,8); c := (1,5)(2,7);
G := Group(a,b,c);
T2 := [ a, a^c, a^(c*b), a^(c*b*c), b, b^c, b^(c*a), b^(c*a*c), c, c^a, c^b, c^(a*b), 
        a*b, (a*b)^c, (a*b)^(c*a), (a*b)^(c*a*c), (a*c)^2, ((a*c)^2)^b];
T6 := [a, a^c, a^(c*b), a^(c*b*c), b, b^c, b^(c*a), b^(c*a*c), c, c^a, c^b, c^(a*b), 
        a*b, (a*b)^c, (a*b)^(c*a), (a*b)^(c*a*c), (a*c)^2, ((a*c)^2)^b, 
        (b*c)^2, ((b*c)^2)^a, (a*b*c)^2, ((a*b*c)^2)^a ];
 
ex := ShapesOfMajoranaRepresentationAxiomM8(G, T2);
rep_2wr2_T2 := MajoranaRepresentation(ex, 1);
MajoranaAlgebraTest(rep_2wr2_T2);
Info(InfoThesis, 50, STRINGIFY("Group: 2wr2, T2, dim ", MAJORANA_Dimension(rep_2wr2_T2) ) ); 

ex := ShapesOfMajoranaRepresentationAxiomM8(G, T6);
rep_2wr2_T6 := MajoranaRepresentation(ex, 1);
MajoranaAlgebraTest(rep_2wr2_T6);
Info(InfoThesis, 50, STRINGIFY("Group: 2wr2, T6, dim ", MAJORANA_Dimension(rep_2wr2_T6) ) ); 

##### (S3 x S3):2

a := (1,2); b := (4,5); c := (1,6)(2,5)(3,4);
G := Group(a,b,c);;
T := Filtered(G, x -> Order(x) = 2);

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_S3S32 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_S3S32);
Info(InfoThesis, 50, STRINGIFY("Group: (S3 x S3):2, dim ", MAJORANA_Dimension(rep_S3S32) ) ); 

##### 2^2 x S4

a := (1,2)(3,4)(5,6); b := (1,2)(7,8); c := (2,3)(5,6);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,a*b)));
Append(T, AsList(ConjugacyClass(G,(b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*b*c*b*c)^3)));
Append(T, AsList(ConjugacyClass(G,a*b*c*b*c*b*c)));
Append(T, AsList(ConjugacyClass(G,b*c*b*a*c*a*b*c)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_22S4 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_22S4);
Info(InfoThesis, 50, STRINGIFY("Group: 2^2 x S4, dim ", MAJORANA_Dimension(rep_22S4) ) ); 

##### 3^{1+2}_+:2^2

a := (1,2)(3,4)(5,6); b := (1,3)(2,4)(7,8); c := (1,9)(3,8)(5,7);
G := Group(a,b,c);
T := Filtered(G, x -> Order(x) = 2);

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_exsp := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_exsp);
Info(InfoThesis, 50, STRINGIFY("Group: 3^{1+2}_+:2^2, dim ", MAJORANA_Dimension(rep_exsp) ) ); 

##### S5

G := SymmetricGroup(5);
T := Filtered(G, x -> Order(x) = 2);

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_S5 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_S5);
Info(InfoThesis, 50, STRINGIFY("Group: S5, dim ", MAJORANA_Dimension(rep_S5) ) ); 

##### (S3 x S3):2^2

a := (1,2)(3,4); b := (1,2)(5,6); c := (1,10)(2,9)(3,8)(4,5)(6,7);
G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,a*b)));
Append(T, AsList(ConjugacyClass(G,(a*b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*c*b*c)^3)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_S3S322 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_S3S322);
Info(InfoThesis, 50, STRINGIFY("Group: (S3 x S3):2^2, dim ", MAJORANA_Dimension(rep_S3S322) ) ); 

##### 2^4:D10

a := (1,2)(3,4); b := (1,3)(2,4)(5,6)(7,8)(9,10); c := (1,2)(3,5)(4,7)(6,9)(8,10);
G := Group(a,b,c);

T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,(a*c)^2)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_24D10 := MajoranaRepresentation(ex, 1, "AxiomM8" );
MajoranaAlgebraTest(rep_24D10);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:D10, dim ", MAJORANA_Dimension(rep_24D10) ) ); 

##### 2^4:D12

a := (1,2)(3,4); b := (1,2)(5,6); c := (1,5)(2,7)(3,6)(4,8);

G := Group(a,b,c);
T2 := [];
Append(T2, AsList(ConjugacyClass(G,a)));
Append(T2, AsList(ConjugacyClass(G,b)));
Append(T2, AsList(ConjugacyClass(G,c)));
Append(T2, AsList(ConjugacyClass(G,(b*c)^3)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T2);
rep_24D12 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_24D12);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:D12, dim ", MAJORANA_Dimension(rep_24D12) ) ); 

##### 2S5

a := (1,2); b := (1,2)(3,4)(6,7); c := (1,5)(2,3)(6,7);
G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,a*b)));
Append(T, AsList(ConjugacyClass(G,(a*b*c*a*c)^3)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_2S5 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_2S5);
Info(InfoThesis, 50, STRINGIFY("Group: 2 x S5, dim ", MAJORANA_Dimension(rep_2S5) ) ); 

##### 2^5:D12

a := (1,2)(3,4)(9,10); b := (1,2)(5,6); c := (1,5)(2,7)(3,8)(4,6)(9,10);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,a*b)));
Append(T, AsList(ConjugacyClass(G,(b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*b*c*b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(a*c)^2)));
Append(T, AsList(ConjugacyClass(G,a*b*c*a*c*b*c*a*c)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_25D12 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_25D12);
Info(InfoThesis, 50, STRINGIFY("Group: 2^5:D12, dim ", MAJORANA_Dimension(rep_25D12) ) ); 

##### 2^4:(S3 x S3)

a := (1,2)(3,4)(5,6)(7,8)(9,10)(11,12);
b := (1,3)(2,4)(5,6)(7,8)(9,10)(11,12);
c := (2,5)(6,7);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,(a*b*c)^3)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_24S3S3 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_24S3S3);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:(S3 x S3), dim ", MAJORANA_Dimension(rep_24S3S3) ) ); 

##### 2^4:A5

a := (1, 2)(3, 4)(5, 6)(7, 8)(9, 10)(11, 12);
b := (1, 11)(2, 12)(3, 9)(4, 10)(5, 6)(13, 14);
c := (1, 3)(2, 15)(4, 13)(6, 12)(7, 11)(14, 16);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,(a*c)^3)));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_24A5 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_24A5);
Info(InfoThesis, 50, STRINGIFY("Group: 2^4:A5, dim ", MAJORANA_Dimension(rep_24A5) ) );

##### 2 x S6

a := (1,2); b := (1,2)(3,4)(5,6)(7,8); c := (2,3)(4,5)(7,8);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,b)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,(a*c)^3)));
Append(T, AsList(ConjugacyClass(G,(b*c)^3)));
Append(T, AsList(ConjugacyClass(G,(7,8))));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_2S6 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_2S6);
Info(InfoThesis, 50, STRINGIFY("Group: 2:S6, dim ", MAJORANA_Dimension(rep_2S6) ) );

##### 2^5:S5

a := (1,2)(3,4)(5,6)(7,8)(9,10)(11,12);
b := (1,3)(2,4)(5,7)(6,8)(9,11)(10,12);
c := (1,8)(2,6)(3,9)(4,12)(5,10)(7,11);

G := Group(a,b,c);
T := [];
Append(T, AsList(ConjugacyClass(G,a)));
Append(T, AsList(ConjugacyClass(G,c)));
Append(T, AsList(ConjugacyClass(G,(a*c)^3)));
Append(T, AsList(ConjugacyClass(G,Center(G).1)));
Append(T, AsList(ConjugacyClass(G,(a*((a*c)^3)*(b*c))));

ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
rep_25S5 := MajoranaRepresentation(ex, 1, "AxiomM8");
MajoranaAlgebraTest(rep_25S5);
Info(InfoThesis, 50, STRINGIFY("Group: 2^5:S5, dim ", MAJORANA_Dimension(rep_25S5) ) );

##### 242 non-existence

a := (1,2)(3,4); b := (5,6)(7,8); c := (1,3)(5,7);
G := Group(a,b,c);

T5 := [a, a^c, b, b^c, c, c^a, c^b, c^(a*b), a*b, (a*b)^c, a*c*b*c, c*a*c*b, a*b*c*b*c, c*a*b*c*b, a*b*c*a*c, c*a*b*c*a];;










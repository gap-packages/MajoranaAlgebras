BindGlobal( "MAJORANA_AllShapes", function( arg )
    local ex, axioms, reps, i;

    ex := arg[1];

    if Length(arg) = 2 then
        axioms := arg[2];
    else
        axioms := "AllAxioms";
    fi;

    reps := [];

    for i in [1 .. Size(ex.shapes)] do
        Add(reps, MajoranaRepresentation(ex, i, axioms));
    od;

    return reps;

end );

BindGlobal("S3S3", function()
    local G, T, ex;
    G := Group((1,2),(1,3),(4,5),(4,6));;
    T := [(1,2),(1,3),(2,3),(4,5),(4,6),(5,6)];;
    ex := ShapesOfMajoranaRepresentation(G,T);;
    return ex;
    end
    );

BindGlobal("A8",function()
    local G, T,input;
    G:=AlternatingGroup(8);
    T:=ShallowCopy(AsList(ConjugacyClass(G, (1,2)(3,4))));
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );

BindGlobal("J2",function() # shape 2
    local G, T, t;
    G := AtlasGroup("J2");
    t := (1,7)(2,47)(3,79)(4,25)(5,97)(6,41)(8,86)(9,52)(10,28)(14,42)(15,70)(16,80)(17,72)(19,20)(21,65)(22,90)(24,96)(26,64)(27,36)(29,67)(31,92)(32,69)(33,49)(34,88)(35,83)(37,74)(40,99)(45,48)(50,93)(51,82)(53,84)(55,63)(58,60)(59,77)(61,87)(62,95)(68,100)(73,81)(76,85)(89,98);
    T := ShallowCopy(AsList(ConjugacyClass(G, t)));
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    end );

BindGlobal("L34",function()
    local G, T, ex;
    G := PSL(3,4);;
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return ex;
    end );

BindGlobal("L42",function()
    local G, T, C, ex;
    G := PSL(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 210);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("U33",function()
    local G, T, ex;
    G := PSU(3,3);
    T := Filtered(G, x -> Order(x) = 2);
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(ex);
    end);

BindGlobal("U42T1",function()
    local G, T, C, ex;
    G := PSU(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 45);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("U42T2",function()
    local G, T, C, ex;
    G := PSU(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 270);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("25S5", function()
    local a, b, c, C, G, T, ex;
    a := (1,2)(3,4)(5,6)(7,8)(9,10)(11,12) ;;
    b := (1,3)(2,4)(5,7)(6,8)(9,11)(10,12) ;;
    c := (1,8)(2,6)(3,9)(4,12)(5,10)(7,11) ;;
    G := Group(a,b,c);;
    C := ConjugacyClasses(G);
    C := Filtered(C,x -> Order(Representative(x)) = 2);
    T := [];
    Append(T,AsList(C[1]));
    Append(T,AsList(C[3]));
    Append(T,AsList(C[4]));
    Append(T,AsList(C[6]));
    Append(T,AsList(C[7]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end);

BindGlobal("24A5", function()
    local a, b, c, C, G, T, ex;
    a := (1, 2)(3, 4)(5, 6)(7, 8)(9, 10)(11, 12) ;;
    b := (1, 11)(2, 12)(3, 9)(4, 10)(5, 6)(13, 14) ;;
    c := (1, 3)(2, 15)(4, 13)(6, 12)(7, 11)(14, 16) ;;
    G := Group(a,b,c);
    T := [];
    Append(T, AsList(ConjugacyClass(G,a)));
    Append(T, AsList(ConjugacyClass(G,(a*c)^3)));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return  ex;
    end);

BindGlobal("2wr2", function()
    local a, b, c, G, T, ex;
    a := (1,2)(3,4);;
    b := (1,3)(2,4)(5,6)(7,8);;
    c := (1,5)(2,7);;
    G := Group(a,b,c);;
    T := [ (1,2)(3,4), (3,4)(5,7), (1,2)(6,8), (5,7)(6,8), (1,3)(2,4)(5,6)(7,8),
        (1,6)(2,8)(3,5)(4,7), (1,8)(2,6)(3,7)(4,5), (1,4)(2,3)(5,8)(6,7), (1,5)(2,7),
        (1,7)(2,5), (3,6)(4,8), (3,8)(4,6), (1,4)(2,3)(5,6)(7,8), (1,6)(2,8)(3,7)(4,5),
        (1,8)(2,6)(3,5)(4,7), (1,3)(2,4)(5,8)(6,7) ];;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end);

BindGlobal("thesis", function()
    local a, b, c, T, G, ex;
    a := (1,2)(3,4);;
    b := (5,6)(7,8);;
    c := (1,3)(5,7);;
    G := Group(a,b,c);;
    T := [a, a^c, b, b^c, c, c^a, c^b, c^(a*b), a*b, (a*b)^c];;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("M12", function()
    local G, T, ex;
    G := MathieuGroup(12);;
    T := ShallowCopy(AsList(ConjugacyClass(G, (1,11)(2,7)(3,5)(4,6)(8,9)(10,12))));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return ex;
    end );

BindGlobal( "min3gen9", function()
    local a, b, c, G, T;
    a := (1,3)(2,4);; b := (1,5)(2,6)(3,7)(4,8);; c := (1,8)(2,5)(3,6)(4,7);;
    G := Group(a, b, c);;
    T := [a, a^b, b, b^a, c, c^a];;
    return ShapesOfMajoranaRepresentation(G,T);;
    end );

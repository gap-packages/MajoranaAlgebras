
BindGlobal("MAJORANA_Example_S3S3", function()
    local G, T, ex;
    G := Group((1,2),(1,3),(4,5),(4,6));;
    T := [(1,2),(1,3),(2,3),(4,5),(4,6),(5,6)];;
    ex := ShapesOfMajoranaRepresentation(G,T);;
    return ex;
    end
    );

BindGlobal("MAJORANA_Example_A8",function()
    local G, T,input;
    G:=AlternatingGroup(8);
    T:=ShallowCopy(AsList(ConjugacyClass(G, (1,2)(3,4))));
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );

BindGlobal("MAJORANA_Example_J2",function() # shape 2
    local G, T, t, b11, b21;
    b11 := (1,84)(2,20)(3,48)(4,56)(5,82)(6,67)(7,55)(8,41)(9,35)(10,40)(11,78)(12,100)(13,49)(14,37)(15,94)(16,76)(17,19)(18,44)(21,34)(22,85)(23,92)(24,57)(25,75)(26,28)(27,64)(29,90)(30,97)(31,38)(32,68)(33,69)(36,53)(39,61)(42,73)(43,91)(45,86)(46,81)(47,89)(50,93)(51,96)(52,72)(54,74)(58,99)(59,95)(60,63)(62,83)(65,70)(66,88)(71,87)(77,98)(79,80);

    b21 := (1,80,22)(2,9,11)(3,53,87)(4,23,78)(5,51,18)(6,37,24)(8,27,60)(10,62,47)(12,65,31)(13,64,19)(14,61,52)(15,98,25)(16,73,32)(17,39,33)(20,97,58)(21,96,67)(26,93,99)(28,57,35)(29,71,55)(30,69,45)(34,86,82)(38,59,94)(40,43,91)(42,68,44)(46,85,89)(48,76,90)(49,92,77)(50,66,88)(54,95,56)(63,74,72)(70,81,75)(79,100,83);

    G := Group(b11, b21);

    t := (1,7)(2,47)(3,79)(4,25)(5,97)(6,41)(8,86)(9,52)(10,28)(14,42)(15,70)(16,80)(17,72)(19,20)(21,65)(22,90)(24,96)(26,64)(27,36)(29,67)(31,92)(32,69)(33,49)(34,88)(35,83)(37,74)(40,99)(45,48)(50,93)(51,82)(53,84)(55,63)(58,60)(59,77)(61,87)(62,95)(68,100)(73,81)(76,85)(89,98);
    T := ShallowCopy(AsList(ConjugacyClass(G, t)));
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    end );

BindGlobal("MAJORANA_Example_L34",function()
    local G, T, ex;
    G := PSL(3,4);;
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return ex;
    end );

BindGlobal("MAJORANA_Example_L42",function()
    local G, T, C, ex;
    G := PSL(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 210);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("MAJORANA_Example_U33",function()
    local G, T, ex;
    G := PSU(3,3);
    T := Filtered(G, x -> Order(x) = 2);
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(ex);
    end);

BindGlobal("MAJORANA_Example_U42T1",function()
    local G, T, C, ex;
    G := PSU(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 45);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("MAJORANA_Example_U42T2",function()
    local G, T, C, ex;
    G := PSU(4,2);
    C := ConjugacyClasses(G);
    C := Filtered(C, x -> Order(Representative(x)) = 2);
    C := Filtered(C, x -> Size(x) = 270);
    T := ShallowCopy(AsList(C[1]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("MAJORANA_Example_25S5", function()
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

BindGlobal("MAJORANA_Example_24A5", function()
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

BindGlobal("MAJORANA_Example_2wr2", function()
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

BindGlobal("MAJORANA_Example_thesis", function()
    local a, b, c, T, G, ex;
    a := (1,2)(3,4);;
    b := (5,6)(7,8);;
    c := (1,3)(5,7);;
    G := Group(a,b,c);;
    T := [a, a^c, b, b^c, c, c^a, c^b, c^(a*b), a*b, (a*b)^c];;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end );

BindGlobal("MAJORANA_Example_M12", function()
    local G, T, ex;
    G := MathieuGroup(12);;
    T := ShallowCopy(AsList(ConjugacyClass(G, (1,11)(2,7)(3,5)(4,6)(8,9)(10,12))));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return ex;
    end );

BindGlobal( "MAJORANA_Example_min3gen9", function()
    local a, b, c, G, T;
    a := (1,3)(2,4);; b := (1,5)(2,6)(3,7)(4,8);; c := (1,8)(2,5)(3,6)(4,7);;
    G := Group(a, b, c);;
    T := [a, a^b, b, b^a, c, c^a];;
    return ShapesOfMajoranaRepresentation(G,T);;
    end );

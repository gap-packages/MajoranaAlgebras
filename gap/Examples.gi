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

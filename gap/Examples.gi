BindGlobal("PSL211",function() 
    local G, T, input; 
    G:=PSL(2,11); 
    T:=Filtered(AsList(G),x->Order(x) = 2); 
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input); 
    end
    );
    
BindGlobal("2XD8",function()
    local G, T, input, res;
    T:=[ (1,2)(5,6), (3,4)(5,6), (1,3)(2,4)(5,6), (1,4)(2,3)(5,6), (1,3)(2,4), (1,4)(2,3) ];
    G := Group(T);
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    res := MajoranaRepresentation(input,1);
    return res;
    end
    );
    
BindGlobal("S5",function()
    local G, T, ex;
    G:=SymmetricGroup(5);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return  ex;
    end
    );
 
BindGlobal("A5",function()
    local G, T,input;
    G:=AlternatingGroup(5);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );
    
BindGlobal("A6",function()
    local G, T,input;
    G:=AlternatingGroup(6);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );
    
BindGlobal("A7",function()
    local G, T,input;
    G:=AlternatingGroup(7);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
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
    
BindGlobal("S4",function()
    local G, T, input, res;
    G:=SymmetricGroup(4);
    T := [ (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)];;
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return(input);
    end
    );
    
BindGlobal("S6", function()
    local G, T, input;
    G := SymmetricGroup(6);
    T := ShallowCopy(AsList(ConjugacyClass(G,(1,2))));
    Append(T,AsList(ConjugacyClass(G,(1,2)(3,4))));
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );
    
BindGlobal("S7", function()
    local G, T, input;
    G := SymmetricGroup(7);
    T := ShallowCopy(AsList(ConjugacyClass(G,(1,2))));
    Append(T,AsList(ConjugacyClass(G,(1,2)(3,4))));
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );

BindGlobal("L33", function()
    local G, T, ex;
    G := PSL(3,3);
    T := Filtered(G, x -> Order(x) = 2);
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end 
    );

BindGlobal("M11",function()
    local G, T,input;
    G:=MathieuGroup(11);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return(input);
    end
    );
    
BindGlobal("3A6",function()
    local z1,z2,G,T,ex;
    z1 := (2,6)(4,11)(7,9)(8,13)(10,14)(12,16);;
    z2 := (1,2,7,4)(3,8,6,10)(5,9,13,12)(11,15)(14,17)(16,18);;
    G := Group(z1,z2);;
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return  ex;
    end );
    
BindGlobal("3S6",function()
    local z1, z2, G, T, ex;
    z1 := (2,6)(3,5)(4,10)(8,14)(9,13)(11,17)(15,16);
    z2 := (1,2,7,11,4)(3,8,15,17,10)(5,9,16,18,12);
    G := Group(z1,z2);
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return  ex;
    end );

BindGlobal("3A7",function()
    local z1, z2, G, T, ex;
    z1 := PermList([2,4,8,1,13,15,17,10,21,3,5,25,11,28,16,6,18,7,33,27,22,9,14,35,26, 12,37,23,39,38,43,45,36,32,40,19,20,42,41,24,29,30,44,31,34]);
    z2 := PermList([3,6,9,11,1,13,2,10,14,23,24,4,19,5,29,31,33,35,7,8,38,18,
    28,27,36,40,12,20,26,15, 43,16,44,17,25,22,34,42,21,30,39,41,45,37,32]);
    G := Group(z1,z2);
    T := Filtered(G, x -> Order(x) = 2);;
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return  ex;
    end );
    
BindGlobal("3S7",function()
    local z1, z2, G, T, C, ex;
    z1 := PermList([ 1,6,5,9,3,2,11,8,4,14,7,17,20,10,23,25,12,30,29,13,34,36,15,39,16, 42,43,28,19,18,40,46,33,21,50,22,37,44,24,31,49,26,27,38,55,32,56,48,41,35, 58,52,53,54,45,47,62,51,63,60,61,57,59]);
    z2 := PermList([2,7,8,1,6,10,3,4,12,5,15,18,9,19,24,11,27,31,32,33,14,13,37,40,41, 16,20,17,43,42,22,48,45,49,21,51,25,23,36,26,52,54,34,29,28,50,30,35,55,58, 59,38,39,61,44,46,47,63,60,53,57,56,62]);
    G := Group(z1,z2);
    C := ConjugacyClasses(G);
    C := Filtered(C,x -> Order(Representative(x)) = 2);
    T := [];
    Append(T,AsList(C[1]));
    Append(T,AsList(C[3]));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);;
    return  ex;
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
    T := [];
    Append(T, AsList(ConjugacyClass(G,a)));
    Append(T, AsList(ConjugacyClass(G,b)));
    Append(T, AsList(ConjugacyClass(G,c)));
    Append(T, AsList(ConjugacyClass(G,a*b)));
    ex := ShapesOfMajoranaRepresentationAxiomM8(G,T);
    return ex;
    end);

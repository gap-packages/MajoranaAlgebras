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
    local G, T;
    G:=SymmetricGroup(5);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    return([G,T]);
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

#
# These examples give all the results from Table 4 of the paper 
# "Construcing Majorana Representations" - https://arxiv.org/abs/1803.10723.
#
# Where the index of a shape is given with an asterix, this implies that
# the associated representation is 3-closed
#

BindGlobal("S4T1", function() # Shapes 1, 2, 3, 4

    local G, T;
    
    G := SymmetricGroup(4);
    T := Filtered(AsList(G), x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentation(G,T);
    
    end );
    
BindGlobal("S4T2", function() # Shapes 1, 2, 3*

    local G, T;
    
    G := SymmetricGroup(4);
    T := [ (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)];;
    
    return ShapesOfMajoranaRepresentation(G,T);
    
    end );

BindGlobal("A5", function() # Shapes 1, 2, 3*, 4

    local G, T;
    
    G := AlternatingGroup(5);
    T := Filtered(AsList(G), x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentation(G,T);
    
    end );
    
BindGlobal("S5", function() # Shape 1

    local G, T;
    
    G := SymmetricGroup(5);
    T := Filtered(AsList(G), x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("L32", function() # Shapes 1 and 2

    local G, T;
    
    G := PSL(3,2);
    T := Filtered(AsList(G), x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("A6", function() # Shapes 1, 4

    local G, T;
    
    G := AlternatingGroup(6);
    T := Filtered(AsList(G),x->Order(x) = 2);

    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("S6", function() # Shape 2

    local G, T;
    
    G := SymmetricGroup(6);
    T := [];
    
    Append(T, AsList(ConjugacyClass(G,(1,2))));
    Append(T, AsList(ConjugacyClass(G,(1,2)(3,4))));

    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("3A6",function() # Shapes 1, 2, 4

    local z1, z2, G, T;
    
    z1 := (2,6)(4,11)(7,9)(8,13)(10,14)(12,16);;
    z2 := (1,2,7,4)(3,8,6,10)(5,9,13,12)(11,15)(14,17)(16,18);;
    
    G := Group(z1,z2);;
    T := Filtered(G, x -> Order(x) = 2);;

    return  ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("3S6",function() # Shape 1

    local z1, z2, G, T;
    
    z1 := (2,6)(3,5)(4,10)(8,14)(9,13)(11,17)(15,16);
    z2 := (1,2,7,11,4)(3,8,15,17,10)(5,9,16,18,12);
    
    G := Group(z1,z2);
    T := [];
    
    Append(T, AsList(ConjugacyClass(G, z1)));
    Append(T, AsList(ConjugacyClass(G, (4,14)(6,10)(7,17)(11,16)(12,13)(15,18))));
    
    return  ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );

BindGlobal("S4S3A7", function()  # Shape 8 

    local G, T;
    
    G := Group((1,2),(1,2,3,4),(5,6),(5,7));;
    G := Intersection(G, AlternatingGroup(7));;
    
    T := Filtered(G, x -> Order(x) = 2);;
    
    return  ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("A7",function() # Shape 2

    local G, T;
    
    G := AlternatingGroup(7);
    T := Filtered(AsList(G),x->Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("S7", function() # Shape 2

    local G, T;
    
    G := SymmetricGroup(7);
    T := [];
    
    Append(T, AsList(ConjugacyClass(G,(1,2))));
    Append(T, AsList(ConjugacyClass(G,(1,2)(3,4))));
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("3A7",function() # Shape 1, 2

    local z1, z2, G, T;
    
    z1 := PermList([2,4,8,1,13,15,17,10,21,3,5,25,11,28,16,6,18,7,33,27,22,9,14,35,26, 12,37,23,39,38,43,45,36,32,40,19,20,42,41,24,29,30,44,31,34]);
    z2 := PermList([3,6,9,11,1,13,2,10,14,23,24,4,19,5,29,31,33,35,7,8,38,18,
    28,27,36,40,12,20,26,15, 43,16,44,17,25,22,34,42,21,30,39,41,45,37,32]);
    
    G := Group(z1,z2);
    T := Filtered(G, x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("3S7",function() # Shape 1

    local z1, z2, G, T, C;
    
    z1 := PermList([ 1,6,5,9,3,2,11,8,4,14,7,17,20,10,23,25,12,30,29,13,34,36,15,39,16, 42,43,28,19,18,40,46,33,21,50,22,37,44,24,31,49,26,27,38,55,32,56,48,41,35, 58,52,53,54,45,47,62,51,63,60,61,57,59]);
    z2 := PermList([2,7,8,1,6,10,3,4,12,5,15,18,9,19,24,11,27,31,32,33,14,13,37,40,41, 16,20,17,43,42,22,48,45,49,21,51,25,23,36,26,52,54,34,29,28,50,30,35,55,58, 59,38,39,61,44,46,47,63,60,53,57,56,62]);
    
    G := Group(z1,z2);
    T := [];
    
    C := ConjugacyClasses(G);
    C := Filtered(C,x -> Order(Representative(x)) = 2);
    C := Filtered(C,x -> Size(x) in [63, 105]); 
    
    Append(T,AsList(C[1]));
    Append(T,AsList(C[2]));
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("PSL211", function() 
    
    local G, T, input; 
    
    G := PSL(2,11); 
    T := Filtered(AsList(G), x -> Order(x) = 2); 
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("L33", function()

    local G, T;
    
    G := PSL(3,3);
    T := Filtered(G, x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("M11",function() # Shape 1

    local G, T;
    
    G := MathieuGroup(11);
    T := Filtered(AsList(G), x -> Order(x) = 2);
    
    return ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );
    
BindGlobal("S5S3A8", function() # Shape 2*

    local G,T;
    
    G := Group((1,2),(1,2,3,4,5),(6,7),(6,8));;
    G := Intersection(G, AlternatingGroup(8));;
    
    T := Filtered(G, x -> Order(x) = 2);;
    
    return  ShapesOfMajoranaRepresentationAxiomM8(G,T);
    
    end );

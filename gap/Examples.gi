BindGlobal("PSL211",function() 
    local G, T, input; 
    G:=PSL(2,11); 
    T:=Filtered(AsList(G),x->Order(x) = 2); 
    input := ShapesOfMajoranaRepresentation(G,T);
    return(input); 
    end
    );
    
BindGlobal("2XD8",function()
    local G, T, input;
    T:=[ (1,2)(5,6), (3,4)(5,6), (1,3)(2,4)(5,6), (1,4)(2,3)(5,6), (1,3)(2,4), (1,4)(2,3) ];
    G := Group(T);
    input := ShapesOfMajoranaRepresentation(G,T);
    return(input);
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
    input := ShapesOfMajoranaRepresentation(G,T);
    return(input);
    end
    );
    
BindGlobal("A6",function()
    local G, T,input;
    G:=AlternatingGroup(6);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentation(G,T);
    return(input);
    end
    );
    
BindGlobal("A7",function()
    local G, T,input;
    G:=AlternatingGroup(7);
    T:=Filtered(AsList(G),x->Order(x) = 2);
    input := ShapesOfMajoranaRepresentation(G,T);
    return(input);
    end
    );
    
BindGlobal("S4",function()
    local G, T;
    G:=SymmetricGroup(4);
    T := [ (1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (1,2)(3,4), (1,3)(2,4), (1,4)(2,3) ];;
    return([G,T]);
    end
    );
    

    

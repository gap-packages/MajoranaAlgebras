BindGlobal("PSL211",function() 
    local G, T; 
    G:=PSL(2,11); 
    T:=Filtered(AsList(G),x->Order(x) = 2); 
    return([G,T]); 
    end
    );
    
BindGlobal("2XD8",function()
    local G, T;
    G:=Group((1,2)(3,4)(5,6),(1,4)(2,3),(2,4));
    T:=[(1,2)(3,4)(5,6),(1,4)(2,3),(2,4),(1,3)(2,4)(5,6),(1,4)(2,3)(5,6),(1,2)(3,4),(1,3)];
    return([G,T]);
    end
    );

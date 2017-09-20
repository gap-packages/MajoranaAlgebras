InstallGlobalFunction(ThreeClosedMajorana,

    function(GramMatrix, AlgebraProducts, EigenVectors, ProductList)
    
    local   unknowns,
            NewProductList,
            x,
            dim,
            new_dim,
            NewAlgebraProducts,
            NewGramMatrix,
            i,
            j;
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    dim := Size(ProductList[1]);
    new_dim := dim + Size(unknowns);
    
    # set up new product list   
    
    NewProductList := ShallowCopy(ProductList);
    NewProductList := List(NewProductList, x -> ShallowCopy(x));
    
    for x in unknowns do 
        Add(NewProductList[1], [ProductList[x[1]], ProductList[x[1]]]);
    od;
    
    # add new coords to existing inner and algebra products
    
    for i in [1..Size(AlgebraProducts)] do 
        if AlgebraProducts[i] <> false then 
            AlgebraProducts[i] := ShallowCopy(AlgebraProducts[i]);
            AlgebraProducts[i] := Concatenation(AlgebraProducts[i],[1..Size(unknowns)]*0);
        fi;
    od;
    
    for i in ProductList[10][1] do 
        for j in [1..3] do 
            EigenVectors[i][j] := List(EigenVectors[i][j], x -> ShallowCopy(x));
            EigenVectors[i][j] := List(EigenVectors[i][j], x -> Concatenation(x, [1..Size(unknowns)]*0));
        od;
    od;
    
    # set up new inner and algebra products
    
    NewGramMatrix := NullMat(new_dim, new_dim);
    NewAlgebraProducts := NullMat(new_dim, new_dim); 
    
    for i in [1 .. new_dim] do 
        for j in [1..new_dim] do 
            NewAlgebraProducts[i][j] := false; 
            NewGramMatrix[i][j] := false;
        od;
    od;
    
    return [NewAlgebraProducts, NewProductList];
    
    end );
    
# InstallGlobalFunction(MAJORANA_ThreeClosedAxiomM1,

    
InstallGlobalFunction(MAJORANA_ThreeClosedAlgebraProduct,

    function(u,v,AlgebraProducts,NewAlgebraProducts,ProductList,NewProductList,unknowns)
    
    local   dim,
            new_dim,
            x, 
            i, 
            j,
            prod;
    
    dim := Size(ProductList[1]);
    new_dim := Size(NewProductList[1]);
    
    x := MAJORANA_SeparateAlgebraProduct(u{[1..dim]}, v{[1..dim]}, unknowns, AlgebraProducts, ProductList);
    
    prod := Concatenation(-x[2]{[1..dim]},x[1]);
    
    for i in [1..new_dim] do 
        if u[i] <> 0 then 
            for j in [dim + 1..new_dim] do 
                if v[j] <> 0 then 
                    
                    x := NewAlgebraProducts[i][j];
                    
                    if x = false then 
                        return false;
                    fi;
                
                    prod := prod + NewAlgebraProducts[i][j];
                fi;
            od;
        fi;
    od;
    
    return prod;
    
    end );
            
    
InstallGlobalFunction(MAJORANA_ThreeClosedFusion,

    function(GramMatrix, AlgebraProducts, NewAlgebraProducts, EigenVectors, ProductList, NewProductList, unknowns)
    
    local   dim,
            new_dim,
            i,
            j,
            new,
            ev,
            new_ev,
            pos,
            a,
            b,
            x,
            y,
            z,
            u;
            
    dim := Size(ProductList[1]);
    new_dim := Size(NewProductList[1]);

    for i in ProductList[10][1] do 
        for j in [1..3] do 
            if EigenVectors[i][j] <> [] then
                MAJORANA_ReversedEchelonForm(EigenVectors[i][j]);                
            fi;
        od;
    od;

    for i in ProductList[10][1] do
        
        new := [ [], [], [] ];
        u := [1..new_dim]*0; u[i] := 1;
        
        for j in [[1,1],[1,2],[1,3],[2,3]] do
            
            new_ev := MAJORANA_FusionTable[j[1] + 1][j[2] + 1];
            pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
            
            for a in EigenVectors[i][j[1]] do 
                for b in EigenVectors[i][j[2]] do 
                    x := MAJORANA_ThreeClosedAlgebraProduct(a, b,  AlgebraProducts, NewAlgebraProducts, ProductList, NewProductList, unknowns);
                    
                    if x <> false then 
                        Add(new[pos], x);
                    fi;
                od;
            od;            
        od;
        
        if false then 
        for a in EigenVectors[i][2] do 
            for b in EigenVectors[i][2] do
                x := MAJORANA_ThreeClosedAlgebraProduct(a, b,  AlgebraProducts, NewAlgebraProducts, ProductList, NewProductList, unknowns);
                
                if x <> false then
                    y := MAJORANA_InnerProduct(a,b, GramMatrix, ProductList);
                    if y <> false then 
                        Add(new[1], x - (1/4)*u*y);
                    fi;
                fi;
            od;
        od;
        fi;
        
        for a in EigenVectors[i][3] do 
            for b in EigenVectors[i][3] do 
                x := MAJORANA_ThreeClosedAlgebraProduct(a, b,  AlgebraProducts, NewAlgebraProducts, ProductList, NewProductList, unknowns);
                
                if x <> false then
                    y := MAJORANA_InnerProduct(a, b, GramMatrix, ProductList);
                    
                    if y <> false then 
                        z := MAJORANA_ThreeClosedAlgebraProduct(u, x, AlgebraProducts,NewAlgebraProducts,ProductList,NewProductList,unknowns);
                        
                        if z <> false then  
                            Add(new[2], z - (1/32)*u*y); 
                        fi;
                    fi;
                fi;
            od;
        od;
        
        for j in [1..3] do 
            Append(EigenVectors[i][j], new[j]);
            
            if Size(EigenVectors[i][j]) > 0 then 
                EigenVectors[i][j] := ShallowCopy(BaseMat(EigenVectors[i][j]));
            fi;
        od;    
    od;
            
    end );
            
            
            
            

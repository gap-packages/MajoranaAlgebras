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
            j,
            u,
            v,
            pos;
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    dim := Size(ProductList[1]);
    new_dim := dim + Size(unknowns);
    
    # set up new product list   
    
    NewProductList := ShallowCopy(ProductList);
    NewProductList := List(NewProductList, x -> ShallowCopy(x));
    
    for x in unknowns do 
        Add(NewProductList[1], [ProductList[x[1]], ProductList[x[1]]]);
    od;
    
    # add new coords to existing eigenvectors
    
    for i in ProductList[10][1] do 
        for j in [1..3] do 
            EigenVectors[i][j] := List(EigenVectors[i][j], x -> ShallowCopy(x));
            EigenVectors[i][j] := List(EigenVectors[i][j], x -> Concatenation(x, [1..Size(unknowns)]*0));
        od;
    od;
    
    # set up new inner and algebra products
    
    NewGramMatrix := NullMat(new_dim, new_dim);
    NewAlgebraProducts := NullMat(new_dim, new_dim);
    
    for i in [1..new_dim] do 
        for j in [1..new_dim] do
            NewGramMatrix[i][j] := false;
            NewAlgebraProducts[i][j] := false;
        od;
    od;

    for i in [1..dim] do 
        for j in [i..dim] do 
            u := [1..dim]*0;; u[i] := 1;;
            v := [1..dim]*0;; v[j] := 1;;
            
            NewGramMatrix[i][j] := MAJORANA_InnerProduct(u,v,GramMatrix,ProductList);
            NewGramMatrix[j][i] := MAJORANA_InnerProduct(u,v,GramMatrix,ProductList);
            
            x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
            
            if x = false then  
                pos := Position(unknowns, [i,j]);
                
                NewAlgebraProducts[i][j] := [1..new_dim]*0;;
                NewAlgebraProducts[i][j][pos] := 1;
                
                NewAlgebraProducts[j][i] := [1..new_dim]*0;;
                NewAlgebraProducts[j][i][pos] := 1;
            else
                NewAlgebraProducts[i][j] := Concatenation(x, [1..Size(unknowns)]*0);
                NewAlgebraProducts[j][i] := Concatenation(x, [1..Size(unknowns)]*0);
            fi;
            
           
        od;    
    od;       
    
    return [NewGramMatrix, NewAlgebraProducts, NewProductList];
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedProduct,

    function(u,v,products)
    
    local   i,
            j,
            prod;
            
    prod := 0;
            
    for i in [1..Size(u)] do 
        if u[i] <> 0 then 
            for j in [1..Size(v)] do 
                if v[j] <> 0 then 
                    if products[i][j] = false then 
                        return false;
                    fi;
                    
                    prod := prod + u[i]*v[j]*products[i][j];
                fi;
            od;
        fi;
    od;
          
    return prod;
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedExtractUnknownProducts,

    function(products,dim)
    
    local   i,
            j,
            list,
            new_dim;
            
    list := [];
    new_dim := Size(products);
    
    for i in [1..new_dim] do 
        for j in [Maximum([i, dim + 1])..new_dim] do 
            if products[i][j] = false then 
               Add(list, [i,j]);
            fi;
        od;
    od;
    
    return list;
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedSeparateProduct,

    function(u, v, new_unknowns, products)
    
    local   new_dim,
            i,
            j,
            row,
            sum,
            pos;
            
    new_dim := Size(products);
    
    row := [1..Size(new_unknowns)]*0;
    sum := 0;
    
    for i in [1..new_dim] do 
        if u[i] <> 0 then 
            for j in [1..new_dim] do 
                if v[j] <> 0 then 
                    if products[i][j] <> false then 
                        sum := sum - u[i]*v[j]*products[i][j];
                    else
                        if i < j then 
                            pos := Position(new_unknowns, [i,j]);
                        else 
                            pos := Position(new_unknowns, [j,i]);
                        fi;
                        
                        row[pos] := row[pos] + u[i]*v[j];
                        
                    fi;
                fi;
            od;
        fi;
    od;
    
    return [row,sum];

    end);
    
                    
InstallGlobalFunction(MAJORANA_ThreeClosedAxiomM1,

    function(NewGramMatrix, NewAlgebraProducts, ProductList, NewProductList, unknowns)
    
    local   dim,
            new_dim,
            i,
            j,
            u,
            v,
            w,
            x;

    dim := Size(ProductList[1]);
    new_dim := Size(NewProductList[1]);
    
    for i in [1..dim] do 
        u := [1..dim]*0; u[i] := 1;
        
        for j in [1..Size(unknowns)] do 
            if NewGramMatrix[i][dim + j] = false then 
                v := [1..dim]*0; v[unknowns[j][1]] := 1;
                w := [1..dim]*0; w[unknowns[j][2]] := 1;
                
                x := MAJORANA_ThreeClosedProduct(u,v,NewAlgebraProducts);
                
                if x <> false then
                    NewGramMatrix[i][dim + j] := MAJORANA_ThreeClosedProduct(x,w,NewGramMatrix);
                    NewGramMatrix[dim + j][i] := MAJORANA_ThreeClosedProduct(x,w,NewGramMatrix);
                else 
                    x := MAJORANA_ThreeClosedProduct(u,w,NewAlgebraProducts);
                    
                    if x <> false then 
                        NewGramMatrix[i][dim + j] := MAJORANA_ThreeClosedProduct(x,v,NewGramMatrix);
                        NewGramMatrix[dim + j][i] := MAJORANA_ThreeClosedProduct(x,v,NewGramMatrix);
                    fi;
                fi;
            fi;
        od;
    od;
    
    end );
    
## We are getting extra long vectors again here :(
    
InstallGlobalFunction(MAJORANA_ThreeClosedFusion,

    function(NewGramMatrix, NewAlgebraProducts, EigenVectors, ProductList, NewProductList, unknowns)
    
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
                    x := MAJORANA_ThreeClosedProduct(a, b,  NewAlgebraProducts);
                    
                    if x <> false then 
                        Add(new[pos], x);
                    fi;
                od;
            od;            
        od;
        
        if false then 
        for a in EigenVectors[i][2] do 
            for b in EigenVectors[i][2] do
                x := MAJORANA_ThreeClosedProduct(a, b, NewAlgebraProducts);
                
                if x <> false then
                    y := MAJORANA_ThreeClosedProduct(a,b, NewGramMatrix);
                    if y <> false then 
                        Add(new[1], x - (1/4)*u*y);
                    fi;
                fi;
            od;
        od;
        fi;
        
        for a in EigenVectors[i][3] do 
            for b in EigenVectors[i][3] do 
                x := MAJORANA_ThreeClosedProduct(a, b,  NewAlgebraProducts);
                
                if x <> false then
                    y := MAJORANA_ThreeClosedProduct(a, b, NewGramMatrix);
                    
                    if y <> false then 
                        z := MAJORANA_ThreeClosedProduct(u, x, NewAlgebraProducts);
                        
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
            
    
            
            

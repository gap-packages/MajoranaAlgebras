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
            pos,
            falsecount,
            newfalsecount,
            maindimensions,
            newdimensions,
            switchmain;
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    dim := Size(ProductList[1]);
    new_dim := dim + Size(unknowns);
    
    # set up new product list   
    
    NewProductList := ShallowCopy(ProductList);
    NewProductList := List(NewProductList, x -> ShallowCopy(x));
    
    for x in unknowns do 
        Add(NewProductList[1], [ProductList[1][x[1]], ProductList[1][x[2]]]);
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
                NewAlgebraProducts[i][j][dim + pos] := 1;
                
                NewAlgebraProducts[j][i] := [1..new_dim]*0;;
                NewAlgebraProducts[j][i][dim + pos] := 1;
            else
                NewAlgebraProducts[i][j] := Concatenation(x, [1..Size(unknowns)]*0);
                NewAlgebraProducts[j][i] := Concatenation(x, [1..Size(unknowns)]*0);
            fi;
            
           
        od;    
    od;       
    
    falsecount := [,];
    
    falsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts,dim));
    falsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix,dim));
    
    maindimensions := [];
    
    for i in NewProductList[10][1] do 
        for j in [1..3] do 
            if Size(EigenVectors[i][j]) > 0 then 
                EigenVectors[i][j] := ShallowCopy(BaseMat(EigenVectors[i][j]));
            fi;
        od;
        
        Add(maindimensions, Sum(List(EigenVectors[i], x -> Size(x))) + 1);
    od;
    
    if falsecount = [0,0] and ForAll(maindimensions, x -> x = new_dim) then 
        switchmain := 1;
    else
        switchmain := 0;
    fi;
    
    while switchmain = 0 do 
    
        MAJORANA_ThreeClosedAxiomM1(NewGramMatrix, NewAlgebraProducts, NewProductList, unknowns);
        
        MAJORANA_ThreeClosedFusion(NewGramMatrix, NewAlgebraProducts, EigenVectors, NewProductList, unknowns);
        
        MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns(NewAlgebraProducts,EigenVectors,NewProductList, unknowns);
        
        MAJORANA_ThreeClosedResurrection(NewAlgebraProducts, EigenVectors, NewProductList, unknowns);
    
        newdimensions := [];
        
        for i in NewProductList[10][1] do
            Add(newdimensions, Sum(List(EigenVectors[i], x -> Size(x))) + 1);
        od;
    
        newfalsecount := [0,0];
        
        newfalsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts, dim));
        newfalsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix,dim));
    
        if newfalsecount = [0,0] then
            break;
        elif newdimensions = maindimensions and newfalsecount = falsecount then 
            Error( "fail" );
            
             return [NewGramMatrix, NewAlgebraProducts, NewProductList];
        else
            maindimensions := StructuralCopy(newdimensions);
            falsecount := StructuralCopy(newfalsecount);
        fi;
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
    
InstallGlobalFunction( MAJORANA_ThreeClosedConjugateVector,
    
    function(v, g, NewProductList, unknowns)
    
    local   dim,
            new_dim,
            vec,
            list,
            i,
            j,
            k,
            h,
            pos;
    
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
    vec := [1..new_dim]*0;;
    
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
    list := [[],[1],[1,2],[1,3],[1,2,3,4]];
    
    for i in [1..dim] do 
        if v[i] <> 0 then 
            
            h := NewProductList[1][i]^g;
        
            for j in list[Order(NewProductList[1][i])] do 
                pos := Position(NewProductList[1], h^j);
                
                if pos <> fail then
                    if Order(NewProductList[1][i]) = 5 and j in [2,3] then 
                        vec[pos] := -v[i];
                    else
                        vec[pos] := v[i];
                    fi;
                fi;
            od;
        fi;
    od;
    
    for i in [dim + 1 .. new_dim] do 
        
        h := List(NewProductList[1][i], x -> x^g);
        
        for j in list[Order(NewProductList[1][i][1])] do 
            for k in list[Order(NewProductList[1][i][2])] do 
                    
                pos := Position(NewProductList[1], [h[1]^j,h[2]^k]);
                    
                if pos <> fail then                     
                    vec[pos] := v[i];
                    
                    if Order(NewProductList[1][i][1]) = 5 and j in [2,3] then 
                        vec[pos] := vec[pos]*(-1);
                    fi;
                    
                    if Order(NewProductList[1][i][2]) = 5 and j in [2,3] then 
                        vec[pos] := vec[pos]*(-1);
                    fi;
                    
                fi; 
            od;
        od;
    od;
    
    return vec;
    
    end);
         
InstallGlobalFunction( MAJORANA_ThreeClosedAxiomM1,

    function(NewGramMatrix, NewAlgebraProducts, NewProductList, unknowns)
    
    local   dim,
            new_dim,
            i,
            j,
            k,
            u,
            v,
            w,
            x;

    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
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
                fi; 
                
                if NewGramMatrix[i][dim + j] = false then
                    x := MAJORANA_ThreeClosedProduct(u,w,NewAlgebraProducts);
                    
                    if x <> false then 
                        NewGramMatrix[i][dim + j] := MAJORANA_ThreeClosedProduct(x,v,NewGramMatrix);
                        NewGramMatrix[dim + j][i] := MAJORANA_ThreeClosedProduct(x,v,NewGramMatrix);
                    fi;
                fi;
            fi;
        od;
    od;
    
    for i in [1..Size(unknowns)] do 
        for j in [1..Size(unknowns)] do 
            for k in [[1,2],[2,1]] do 
                if NewGramMatrix[dim + i][dim + j] = false then 
                    u := [1..new_dim]*0; u[unknowns[i][k[1]]] := 1;
                    v := [1..new_dim]*0; v[unknowns[i][k[2]]] := 1;
                    w := [1..new_dim]*0; w[dim + j] := 1;
                    
                    x := MAJORANA_ThreeClosedProduct(v,w,NewAlgebraProducts);
                    
                    if x <> false then 
                       NewGramMatrix[dim + i][dim + j] := MAJORANA_ThreeClosedProduct(u,x,NewGramMatrix);
                       NewGramMatrix[dim + j][dim + i] := MAJORANA_ThreeClosedProduct(u,x,NewGramMatrix);
                    fi;
                fi;
            od;
        od;
    od;
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedFusion,

    function(NewGramMatrix, NewAlgebraProducts, EigenVectors, NewProductList, unknowns)
    
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
            
    # need to 
            
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);

    for i in NewProductList[10][1] do 
        for j in [1..3] do 
            if EigenVectors[i][j] <> [] then
                MAJORANA_ReversedEchelonForm(EigenVectors[i][j]);                
            fi;
        od;
    od;

    for i in NewProductList[10][1] do
        
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
    
    InstallGlobalFunction( MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns,

    function(NewAlgebraProducts,EigenVectors,NewProductList, unknowns)
    
    local   mat,
            vec,
            record,
            table,
            new_dim,
            dim,
            new_unknowns,
            i,
            j,
            u,
            v,
            x,
            sol;
    
    mat := [];
    vec := [];
    record := [];
    
    table := [0,1/4,1/32];
    
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts, dim);
    
    for i in NewProductList[10][1] do 
        u := [1..new_dim]*0;; u[i] := 1;;
    
        for j in [1..3] do 
            for v in EigenVectors[i][j] do 
                x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, NewAlgebraProducts);
            
                if ForAny(x[1] , y -> y <> 0) then 
                    Add(mat, x[1]);
                    Add(vec, x[2] + table[j]*v);
                    Add(record, [i,j, Position(EigenVectors[i][j],v)]);
                fi;
            od;
        od;
    od;
    
    MAJORANA_ThreeClosedSolutionAlgProducts(mat, vec, NewAlgebraProducts,NewProductList,new_unknowns,unknowns);
    
    end );
            
InstallGlobalFunction( MAJORANA_ThreeClosedSolutionAlgProducts,

    function(mat,vec,NewAlgebraProducts,NewProductList,new_unknowns,unknowns)
    
    local   sol,
            i,
            x;
    
    if mat <> [] then 
        
        sol := MAJORANA_SolutionMatVecs(mat,vec);
    
        if sol <> false then 
            for i in [1..Size(new_unknowns)] do 
                if not i in sol[2] then 
                    
                    x := new_unknowns[i];

                    MAJORANA_ThreeClosedAllConjugates(x,sol[1][i],NewAlgebraProducts,NewProductList,unknowns);
                    
                fi;
            od;
        fi;
    fi;   
    
    end);         
            
InstallGlobalFunction( MAJORANA_ThreeClosedAllConjugates,
            
    function(x,vec,NewAlgebraProducts,NewProductList,unknowns)
        
    local   dim,    
            new_dim,
            list,
            u,
            y,
            v,
            w,
            i,
            j,
            G,
            H,
            reps,
            g,
            pos;
    
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
    list := [[],[1],[1,2],[1,3],[1,2,3,4]];
    
    G := NewProductList[8];
    
    if x[1] <= dim then  
        u := NewProductList[1][x[1]];
        y := ();
    else
        u := NewProductList[1][x[1]][1];
        y := NewProductList[1][x[1]][2];
    fi; 
    
    v := NewProductList[1][x[2]][1];
    w := NewProductList[1][x[2]][2];
    
    H := Stabilizer(G, [u, v, w], OnTuples);
    
    reps := List(RightTransversal(G,H), i -> CanonicalRightCosetElement(H,i));
    
    for g in reps do 
        pos := [0,0];
    
        if x[1] <= dim then 
            for i in list[Order(u)] do 
                if Position(NewProductList[1], (u^i)^g) <> fail then 
                    pos[1] := Position(NewProductList[1], (u^i)^g);
                fi;
            od;
        else
            for i in list[Order(u)] do 
                for j in list[Order(y)] do
                    if Position(NewProductList[1], [(u^i)^g,(y^j)^g]) <> fail then 
                        pos[1] := Position(NewProductList[1], [(u^i)^g,(y^j)^g]);
                    fi;
                od;
            od;
            
            if pos[1] = fail then 
                for i in list[Order(u)] do 
                    for j in list[Order(y)] do
                        if Position(NewProductList[1], [(y^j)^g,(u^i)^g]) <> fail then 
                            pos[1] := Position(NewProductList[1], [(y^j)^g,(u^i)^g]);
                        fi;
                    od;
                od;
            fi;
        fi;
        
        for i in list[Order(v)] do
            for j in list[Order(w)] do
                if Position(NewProductList[1], [(v^i)^g,(w^j)^g]) <> fail then 
                    pos[2] := Position(NewProductList[1], [(v^i)^g,(w^j)^g]);
                fi;
            od;
        od;
        
        if pos[2] = fail then 
            for i in list[Order(v)] do
                for j in list[Order(w)] do
                    if Position(NewProductList[1], [(w^j)^g,(v^i)^g]) <> fail then 
                        pos[2] := Position(NewProductList[1], [(w^j)^g,(v^i)^g]);
                    fi;
                od;
            od;
        fi;
        
        if pos <> x then 
            NewAlgebraProducts[pos[1]][pos[2]] := MAJORANA_ThreeClosedConjugateVector(vec, g, NewProductList, unknowns);
            NewAlgebraProducts[pos[2]][pos[1]] := MAJORANA_ThreeClosedConjugateVector(vec, g, NewProductList, unknowns);
        fi;
    od;
            
    
    end );

InstallGlobalFunction( MAJORANA_ThreeClosedResurrection,

    function( NewAlgebraProducts, EigenVectors, NewProductList, unknowns)

    local   mat,
            vec,
            new_dim,
            dim,
            new_unknowns,
            evals,
            alpha,
            beta,
            gamma,
            i,
            u,
            x,
            y,
            z,
            row,
            sum;
    
    mat := [];
    vec := [];
    
    new_dim := Size(NewProductList[1]);
    dim := new_dim - Size(unknowns);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts,dim);
    
    # need to implement (2,1) resurrection
    
    for i in NewProductList[10][1] do 
        
        u := [1..new_dim]*0; u[i] := 1;
    
        for evals in [[1,2],[2,1],[1,3],[2,3]] do
            for beta in EigenVectors[i][evals[2]] do
                for gamma in EigenVectors[i][evals[1]] do 
                    if MAJORANA_ThreeClosedProduct(beta,gamma,NewAlgebraProducts) = false then
                
                    
                        x := (1/4)*MAJORANA_ThreeClosedSeparateProduct(beta, gamma, new_unknowns, NewAlgebraProducts);
                        
                        for alpha in EigenVectors[i][1] do 
                            y := MAJORANA_ThreeClosedProduct(alpha - beta, gamma, NewAlgebraProducts);
                            
                            if y <> false then 
                            
                                row := x[1]; sum := x[2];
                            
                                z := MAJORANA_ThreeClosedSeparateProduct(u, y, new_unknowns, NewAlgebraProducts);
                                
                                row := row + z[1]; sum := sum + z[2];
                            fi;
                        od;
                    fi;
                od;
            od;
        od;
    od;

    MAJORANA_ThreeClosedSolutionAlgProducts(mat, vec, NewAlgebraProducts,NewProductList,new_unknowns,unknowns);

    end );









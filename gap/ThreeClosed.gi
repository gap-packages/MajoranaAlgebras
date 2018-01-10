InstallGlobalFunction(ThreeClosedMajorana,

    function(innerproducts, algebraproducts, evecs, setup)
    
    local   unknowns,
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
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts, setup);
    
    dim := Size(setup.coords);
    new_dim := dim + Size(unknowns);
    
    # new setup 
    
    for x in unknowns do 
        Add(setup.coords, [setup[1][x[1]], setup[1][x[2]]]);
    od;
    
    # add new coords to existing eigenvectors
    
    for i in setup.orbitreps do 
        for j in [1..3] do 
            evecs[i][j] := SparseMatrix(    Nrows(evecs[i][j]), new_dim, 
                                            evecs[i][j]!.indices,
                                            evecs[i][j]!.entries,
                                            Rationals);
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
            u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
            v := SparseMatrix(1, new_dim, [[j]], [[1]], Rationals);
            
            NewGramMatrix[i][j] := MAJORANA_InnerProduct(u,v,innerproducts,setup);
            NewGramMatrix[j][i] := MAJORANA_InnerProduct(u,v,innerproducts,setup);
            
            x := MAJORANA_AlgebraProduct(u,v,algebraproducts,setup);
            
            if x = false then  
                pos := Position(unknowns, [i,j]);
                
                NewAlgebraProducts[i][j] := SparseMatrix(1, new_dim, [[dim + pos]], [[1]], Rationals);
                NewAlgebraProducts[j][i] := SparseMatrix(1, new_dim, [[dim + pos]], [[1]], Rationals);
            fi;
        od;    
    od;       
    
    falsecount := [,];
    
    falsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts,dim));
    falsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix,dim));
    
    while true do 
    
        MAJORANA_ThreeClosedAxiomM1(NewGramMatrix, NewAlgebraProducts, setup, unknowns);
        
        MAJORANA_ThreeClosedFusion(NewGramMatrix, NewAlgebraProducts, evecs, setup, unknowns);
        
        x := MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns(NewAlgebraProducts, evecs,  setup, unknowns);
        
        MAJORANA_ThreeClosedResurrection(x[1],x[2],NewGramMatrix, NewAlgebraProducts, evecs, setup, unknowns);
        
        MAJORANA_ThreeClosedOrthogonality(NewGramMatrix, evecs, setup, unknowns);
    
        newfalsecount := [0,0];
        
        newfalsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts, dim));
        newfalsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix,dim));
        
        if newfalsecount[2] = 0 then 
        #    MAJORANA_ThreeClosedCheckNullspace(NewGramMatrix,NewAlgebraProducts,evecs,setup);
        fi;
    
        Display([falsecount,newfalsecount]);
    
        if newfalsecount = [0,0] then
            Display("Success");
            return [NewGramMatrix, NewAlgebraProducts, setup];
        elif newfalsecount = falsecount then 
            Display( "fail" );
            return [NewGramMatrix, NewAlgebraProducts, setup];
        else
            falsecount := StructuralCopy(newfalsecount);
        fi;
    od;

    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedProduct,

    function(u,v,products)
    
    local   i,
            j,
            x,
            prod;
            
    prod := 0;
    
    if IsSparseMatrix(products[1][1]) then 
        prod := SparseZeroMatrix(1, Ncols(u), Rationals);
    fi;
            
    for i in [1..Size(u!.indices[1])] do 
        for j in [1..Size(v!.indices[1])] do 
            
            x := products[u!.indices[1][i]][v!.indices[1][j]];
        
            if x = false then 
                return false;
            fi;
            
            prod := prod + u!.entries[1][i]*v!.entries[1][j]*x;
        od;
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
            x,
            row,
            sum,
            pos;
            
    new_dim := Size(products);
    
    row := SparseZeroMatrix(1, Size(new_unknowns), Rationals);
    
    if IsSparseMatrix(products[1][1]) then 
        sum := SparseZeroMatrix(1, new_dim, Rationals);
    else
        sum := 0;
    fi;
    
    for i in [1..Size(u!.indices[1])] do 
        for j in [1..Size(v!.indices[1])] do 
        
            x := products[u!.indices[1][i]][v!.indices[1][j]];
        
            if x <> false then 
                sum := sum - u!.entries[1][i]*v!.entries[1][j]*x;
            else
                if i < j then 
                    pos := Position(new_unknowns, [i,j]);
                else 
                    pos := Position(new_unknowns, [j,i]);
                fi;
                
                AddToEntry(row, 1, pos, u!.entries[1][i]*v!.entries[1][j]);
            fi;
        od;
    od;
    
    if IsSparseMatrix(sum) then 
        return [row,sum];
    else 
        return [row, SparseMatrix(1, 1, [[1]], [[sum]], Rationals)];
    fi;
    
    end);
    
InstallGlobalFunction( MAJORANA_ThreeClosedAxiomM1,

    function(NewGramMatrix, NewAlgebraProducts, setup, unknowns)
    
    local   dim,
            new_dim,
            i,
            j,
            k,
            u,
            v,
            w,
            x,
            y;

    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);
    
    for i in [1..dim] do 
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        for j in [1..Size(unknowns)] do 
            if NewGramMatrix[i][dim + j] = false then 
                v := SparseMatrix(1, dim, [[unknowns[j][1]]], [[1]], Rationals);
                w := SparseMatrix(1, dim, [[unknowns[j][2]]], [[1]], Rationals);
                
                x := MAJORANA_ThreeClosedProduct(u, v, NewAlgebraProducts);
                
                if x <> false then
                    NewGramMatrix[i][dim + j] := MAJORANA_ThreeClosedProduct(x,w,NewGramMatrix);
                    NewGramMatrix[dim + j][i] := MAJORANA_ThreeClosedProduct(x,w,NewGramMatrix);
                fi; 
                
                if NewGramMatrix[i][dim + j] = false then
                    x := MAJORANA_ThreeClosedProduct(u,w,NewAlgebraProducts);
                    
                    if x <> false then 
                        
                        y := MAJORANA_ThreeClosedProduct(x,v,NewGramMatrix);
                        
                        if y <> false then 
                            
                          NewGramMatrix[i][dim + j] := y;
                          NewGramMatrix[dim + j][i] := y;
                          
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;
    
    for i in [1..Size(unknowns)] do 
        for j in [1..Size(unknowns)] do 
            for k in [[1,2],[2,1]] do 
                if NewGramMatrix[dim + i][dim + j] = false then 
                    
                    u := SparseMatrix(1, new_dim, [[unknowns[i][k[1]]]], [[1]], Rationals);
                    v := SparseMatrix(1, new_dim, [[unknowns[i][k[2]]]], [[1]], Rationals);
                    w := SparseMatrix(1, new_dim, [[dim + j]], [[1]], Rationals);
                    
                    x := MAJORANA_ThreeClosedProduct(v,w,NewAlgebraProducts);
                    
                    if x <> false then 
                        y := MAJORANA_ThreeClosedProduct(x,u,NewGramMatrix);
                       
                        if y <> false then 
                            NewGramMatrix[dim + i][dim + j] := y;
                            NewGramMatrix[dim + j][dim + i] := y;
                        fi;
                    fi;
                fi;
            od;
        od;
    od;
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedFuseEigenvectors,
    function(a, b, i, evals, new, NewGramMatrix, NewAlgebraProducts, setup)
    
    local   x,
            y,
            z,
            new_dim,
            u,
            new_ev,
            pos;
    
    new_dim := Size(setup.coords);
    
    u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
    new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
    
    x := MAJORANA_ThreeClosedProduct(a, b, NewAlgebraProducts);
    
    if x <> false then 
        if evals = [2,2] then 
            y := MAJORANA_ThreeClosedProduct(a,b,NewGramMatrix);
            
            if y <> false then 
                new[1] := UnionOfRows(new[1], x - (1/4)*u*y);
            fi;
        elif evals = [3,3] then 
            y := MAJORANA_ThreeClosedProduct(a,b,NewGramMatrix);
            z := MAJORANA_ThreeClosedProduct(u,x,NewAlgebraProducts);
            
            if y <> false and z <> false then 
                new[2] := UnionOfRows(new[2], z - (1/32)*u*y);
                new[1] := UnionOfRows(new[1], x + (3/32)*u*y - 4*z);            
            fi;
        else
            new[pos] := UnionOfRows(new[pos],x);
        fi;
    fi;    
    
    end);
    
InstallGlobalFunction(MAJORANA_ThreeClosedFusion,

    function(NewGramMatrix, NewAlgebraProducts, evecs, setup, unknowns)
    
    local   dim,
            new_dim,
            t,
            bad,
            null,
            evals,
            i,
            j,
            k,
            new,
            ev,
            new_ev,
            pos,
            a,
            b,
            x,
            y,
            z,
            u,
            w,
            mat,
            vecs;
            
    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);

    t := Size(evecs);

    for i in [1..t] do
        
        new := [ [], [], [] ];
        
        for j in [1..3] do 
            new[i] := SparseMatrix(0, new_dim, [], [], Rationals);
        od;
        
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        
        for evals in [[1,1],[1,2],[1,3],[2,1],[2,3],[3,1],[3,2]] do
    
            new_ev := MAJORANA_FusionTable[j[1] + 1][j[2] + 1];
            pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
            
            for j in [1..Nrows(evecs[i][evals[1]])] do 
            
                a := CertainRows(evecs[i][evals[1]], [j]);
                
                bad := MAJORANA_ThreeClosedFindBadIndices(a, NewAlgebraProducts, setup);
                
                if bad <> [] then 
                    null := KernelMat(CertainColumns(evecs[i][evals[2]], bad)).relations;
                else
                    null := SparseIdentityMatrix(Nrows(evecs[i][evals[2]]));
                fi;
                
                for k in [1..Nrows(null)] do 
                    b := CertainRows(null, [k])*evecs[i][evals[2]];
                    
                    MAJORANA_ThreeClosedFuseEigenvectors(a, b, i, evals, new, NewGramMatrix, NewAlgebraProducts, setup);
                    
                od;
            od;
        od;
        
        for j in [1..3] do 
            Append(evecs[i][j], new[j]);
            
            if Size(evecs[i][j]) > 0 then 
                evecs[i][j] := UnionOfRows(evecs[i][j], new[j]);
                evecs[i][j] := MAJORANA_BasisOfEvecs(evecs[i][j]);
            fi;
        od;    
    od;
            
    end );
    
    InstallGlobalFunction( MAJORANA_ThreeClosedFindBadIndices,
    
    function(v, NewAlgebraProducts, setup)
    
    local list, i, j;
    
    list := [];
    
    for i in v!.indices do 
        for j in [1..Size(setup.coords)] do 
            if NewAlgebraProducts[i][j] = false then 
                Add(list, j);
            fi;
        od;
    od;
    
    list := DuplicateFreeList(list);
    Sort(list);
    
    return list;
    
    end);
    
    
    
    InstallGlobalFunction( MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns,

    function(NewAlgebraProducts,evecs,setup, unknowns)
    
    local   mat,
            vec,
            table,
            new_dim,
            dim,
            new_unknowns,
            i,
            j,
            k,
            u,
            v,
            x,
            y,
            g,
            pos,
            sol;
    
    table := [0,1/4,1/32];
    
    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);
            
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts);
    
    mat := SparseMatrix(0, Size(new_unknowns), [], [], Rationals);
    vec := SparseMatrix(0, new_dim, [], [], Rationals);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts, dim);
    
    for i in [1..Size(evecs)] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
        for j in [1..3] do 
            for k in [1..Nrows(evecs[i][j])] do 
                
                v := CertainRows(evecs[i][j], [k]);
            
                x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, NewAlgebraProducts);
            
                if x[1]!.indices <> [] then 
                    mat := UnionOfRows(mat, x[1]);
                    vec := UnionOfRows(vec, x[2]);
                fi;
            od;
        od;
    od;
    
    return [mat,vec];
    
    end );
            
InstallGlobalFunction( MAJORANA_ThreeClosedSolutionProducts,

    function(mat,vec,products,setup,new_unknowns,unknowns)
    
    local   sol,
            i,
            x;
    
    if Nrows(mat) > 0 then 
    
        Display([Size(mat),Size(mat[1])]);
        
        sol := MAJORANA_SolutionMatVecs(mat,vec);
    
        for i in [1..Size(new_unknowns)] do 
            if sol.solutions[i] <> fail then 
                
                x := new_unknowns[i];
                
                products[x[1]][x[2]] := sol.solutions[i];
                products[x[2]][x[1]] := sol.solutions[i];
                
            fi;
        od;
    fi;   
    
    end);         

InstallGlobalFunction( MAJORANA_ThreeClosedResurrection,

    function( mat, vec, NewGramMatrix, NewAlgebraProducts, evecs, setup, unknowns)

    local   pos,
            new_dim,
            dim,
            new_unknowns,
            evals,
            alpha,
            beta,
            gamma,
            ev,
            i,
            j,
            k,
            l,
            u,
            x,
            y,
            z,
            w,
            g,
            row,
            sum;
    
    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts,dim);
    
    for i in [1..Size(evecs)] do 
        
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
        for evals in [[1,2],[2,1],[1,3],[2,3]] do
            for j in [1..Nrows(evecs[i][evals[2]])] do 
                beta := CertainRows(evecs[i][evals[2]], [j]);
                for k in [1..Nrows(evecs[i][evals[1]])] do 
                    gamma := CertainRows(evecs[i][evals[1]], [k]);
                     
                    if MAJORANA_ThreeClosedProduct(beta,gamma,NewAlgebraProducts) = false then
                
                        ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
                        x := ev*MAJORANA_ThreeClosedSeparateProduct(beta, gamma, new_unknowns, NewAlgebraProducts);
                        
                        for l in [1..Nrows(evecs[i][evals[1]])] do 
                        
                            alpha := CertainRows(evecs[i][evals[1]], [l]);
                        
                            row := SparseMatrix(0, Size(new_unknowns), [], [], Rationals); 
                            sum := SparseMatrix(0, new_dim, [], [], Rationals);
                        
                            y := MAJORANA_ThreeClosedProduct(alpha - beta, gamma, NewAlgebraProducts);
                            
                            if y <> false then 
                            
                                row := row + x[1]; sum := sum + x[2];
                            
                                z := MAJORANA_ThreeClosedSeparateProduct(u, y, new_unknowns, NewAlgebraProducts);
                                
                                row := row + z[1]; sum := sum + z[2];
                                
                                if evals[1] = 2 then 
                                    w := MAJORANA_ThreeClosedProduct(alpha,gamma,NewGramMatrix);
                                    
                                    if w <> false then 
                                        sum := sum + (1/4)*u*w;
                                    else 
                                        row := [];
                                    fi;
                                    
                                fi;
                                
                                if row <> [] and row!.indices <> [] then
                                    
                                    mat := UnionOfRows(mat, row);
                                    vec := UnionOfRows(vec, sum);
                                fi;
                            fi;
                        od;
                    fi;
                od;
            od;
        od;
    od;

    Display("Resurrection");

    MAJORANA_ThreeClosedSolutionProducts(mat, vec, NewAlgebraProducts,setup,new_unknowns,unknowns);

    end );

InstallGlobalFunction( MAJORANA_ThreeClosedOrthogonality,

    function(NewGramMatrix, evecs, setup, unknowns)
    
    local   new_dim,
            dim,
            new_unknowns,
            mat,
            vec,
            i,
            j,
            k,
            l,
            m,
            u,
            v,
            x;
            
    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix, dim);
    
    mat := SparseMatrix(0, Size(new_unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);
    
    for i in [1..Size(evecs)] do 
        for j in [1..3] do 
            for k in [j + 1.. 3] do 
                for l in [1..Nrows(evecs[i][j])] do
                    u := CertainRows(evecs[i][j],[l]);
                    for m in [1..Nrows(evecs[i][k])] do 
                        v := CertainRows(evecs[i][k], [m]); 
                        if MAJORANA_ThreeClosedProduct(u,v,NewGramMatrix) = false then
                            x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, NewGramMatrix);
                            
                            mat := UnionOfRows(mat, x[1]);
                            vec := UnionOfRows(vec, x[2]);
                            
                        fi;
                    od;
                od;
            od;
        od;
    od;

    MAJORANA_ThreeClosedSolutionProducts(mat,vec,NewGramMatrix,setup,new_unknowns,unknowns);
    
    end );        
    

InstallGlobalFunction(MAJORANA_ThreeClosedSetUp,

    function(rep)

    local   dim, new_dim,
            unknowns,
            x, i, j, u, v, pos,
            newalgebraproducts,
            newinnerproducts;
            
    MAJORANA_ThreeClosedAllEvecs(rep.evecs, rep.group, rep.setup);
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(rep.algebraproducts, rep.setup);
    
    dim := Size(rep.setup.coords);
    new_dim := dim + Size(unknowns);
    
    # new rep.setup 
    
    for x in unknowns do 
        Add(rep.setup.coords, [rep.setup.coords[x[1]], rep.setup.coords[x[2]]]);
    od;
    
    # add new coords to existing eigenvectors
    
    for i in [1..Size(rep.evecs)] do 
        for j in [1..3] do 
            rep.evecs[i][j] := SparseMatrix(    Nrows(rep.evecs[i][j]), new_dim, 
                                            rep.evecs[i][j]!.indices,
                                            rep.evecs[i][j]!.entries,
                                            Rationals);
        od;
    od;
    
    # add new coords to existing nullspace vectors
    
    rep.nullspace := SparseMatrix(  Nrows(rep.nullspace), new_dim, 
                                    rep.nullspace!.indices,
                                    rep.nullspace!.entries,
                                    Rationals);
    
    # set up new inner and algebra products
    
    newinnerproducts := NullMat(new_dim, new_dim);
    newalgebraproducts := NullMat(new_dim, new_dim);
    
    for i in [1..new_dim] do 
        for j in [1..new_dim] do
            newinnerproducts[i][j] := false;
            newalgebraproducts[i][j] := false;
        od;
    od;

    for i in [1..dim] do 
        for j in [i..dim] do 
            u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
            v := SparseMatrix(1, new_dim, [[j]], [[1]], Rationals);
            
            newinnerproducts[i][j] := MAJORANA_InnerProduct(u,v,rep.innerproducts,rep.setup);
            newinnerproducts[j][i] := MAJORANA_InnerProduct(u,v,rep.innerproducts,rep.setup);
            
            x := MAJORANA_AlgebraProduct(u,v,rep.algebraproducts,rep.setup);
            
            if x = false then  
                pos := Position(unknowns, [i,j]);
                
                newalgebraproducts[i][j] := SparseMatrix(1, new_dim, [[dim + pos]], [[1]], Rationals);
                newalgebraproducts[j][i] := SparseMatrix(1, new_dim, [[dim + pos]], [[1]], Rationals);
            else
                newalgebraproducts[i][j] := x;
                newalgebraproducts[j][i] := x;
            fi;
        od;    
    od;
    
    return rec( unknowns := unknowns,
                group := rep.group, 
                shape := rep.shape, 
                involutions := rep.involutions,
                algebraproducts := newalgebraproducts,
                innerproducts := newinnerproducts,
                evecs := rep.evecs, 
                nullspace := rep.nullspace );
    
    end);

InstallGlobalFunction(ThreeClosedMajorana,

    function(rep)
    
    local   new, x,
            falsecount,
            newfalsecount;
    
    new := MAJORANA_ThreeClosedSetUp(rep);

    falsecount := [,];
    
    falsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.algebraproducts));
    falsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.innerproducts));
    
    while true do 
    
        MAJORANA_ThreeClosedAxiomM1(new.innerproducts, new.algebraproducts, rep.setup, new.unknowns);
        
        MAJORANA_ThreeClosedFusion(new.innerproducts, new.algebraproducts, rep.evecs, rep.setup, new.unknowns);
        
        x := MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns(new.algebraproducts, rep.evecs,  rep.setup, new.unknowns, new.nullspace);
        
        MAJORANA_ThreeClosedResurrection(x[1],x[2],new.innerproducts, new.algebraproducts, rep.evecs, rep.setup, new.unknowns);
        
        MAJORANA_ThreeClosedOrthogonality(new.innerproducts, rep.evecs, rep.setup, new.unknowns);
    
        newfalsecount := [0,0];
        
        newfalsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.algebraproducts));
        newfalsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.innerproducts));
        
        if falsecount[1] <> 0 and newfalsecount[2] = 0 then 
            new.nullspace := NullspaceMat(new.innerproducts);
        fi;
    
        Display([falsecount,newfalsecount]);
    
        if newfalsecount = [0,0] then
            Display("Success");
            return [new.innerproducts, new.algebraproducts, rep.setup];
        elif newfalsecount = falsecount then 
            Display( "fail" );
            return [new.innerproducts, new.algebraproducts, rep.setup];
        else
            falsecount := StructuralCopy(newfalsecount);
        fi;
    od;

    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedAllEvecs,

    function(evecs, group, setup)
    
    local   i,
            j, 
            dim,
            reps, 
            elts,
            pos,
            g;
            
    dim := Size(setup.coords);
    
    for i in [1..Size(evecs)] do 
        if not i in setup.orbitreps then 
        
            evecs[i] := [[],[],[]];
            
            for j in [1..3] do 
                evecs[i][j] := SparseMatrix(0, dim, [], [], Rationals);
            od;
        
            reps := setup.coords{setup.orbitreps};
            elts := List(reps, x -> RepresentativeAction(group, x, setup.coords[i]));
            
            pos := PositionNot(elts, fail);
            
            for j in [1..3] do 
                g := MAJORANA_FindVectorPermutation(elts[pos], setup);
                evecs[i][j] := MAJORANA_ConjugateVector(evecs[pos][j], g, setup);
            od;
        fi;
    od;
    
    end);
            
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

    function(products)
    
    local   i,
            j,
            list,
            new_dim;
            
    list := [];
    new_dim := Size(products);
    
    for i in [1..new_dim] do 
        for j in [i..new_dim] do 
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
                if u!.indices[1][i] < v!.indices[1][j] then 
                    pos := Position(new_unknowns, [u!.indices[1][i],v!.indices[1][j]]);
                else 
                    pos := Position(new_unknowns, [v!.indices[1][j],u!.indices[1][i]]);
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
            new[j] := SparseMatrix(0, new_dim, [], [], Rationals);
        od;
        
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        
        for evals in [[1,1],[1,2],[1,3],[2,1],[2,3],[3,1],[3,2]] do
    
            new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
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
            evecs[i][j] := UnionOfRows(evecs[i][j], new[j]);
            evecs[i][j] := MAJORANA_BasisOfEvecs(evecs[i][j]);
        od;    
    od;
            
    end );
    
    InstallGlobalFunction( MAJORANA_ThreeClosedFindBadIndices,
    
    function(v, NewAlgebraProducts, setup)
    
    local list, i, j;
    
    list := [];
    
    for i in v!.indices[1] do 
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

    function(NewAlgebraProducts,evecs,setup, unknowns, nullspace)
    
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
    
    for i in [1..new_dim] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        for j in [1..Nrows(nullspace)] do 
            v := CertainRows(nullspace, [j]);
            
            x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, NewAlgebraProducts);
            
            if x[1]!.indices <> [] then 
                mat := UnionOfRows(mat, x[1]);
                vec := UnionOfRows(vec, x[2]);
            fi;
        od;
    od;
    
    return [mat, vec];
    
    end );
            
InstallGlobalFunction( MAJORANA_ThreeClosedSolutionProducts,

    function(mat,vec,products,setup,new_unknowns,unknowns)
    
    local   sol,
            i,
            x;
    
    if Nrows(mat) > 0 then 
        
        sol := MAJORANA_SolutionMatVecs(mat,vec);
    
        for i in [1..Size(new_unknowns)] do 
            if sol.solutions[i] <> fail then 
                
                x := new_unknowns[i];
                
                if IsSparseMatrix(products[1][1]) then 
                    products[x[1]][x[2]] := sol.solutions[i];
                    products[x[2]][x[1]] := sol.solutions[i];
                elif sol.solutions[i]!.entries[1] = [] then 
                    products[x[1]][x[2]] := 0;
                    products[x[2]][x[1]] := 0;
                else
                    products[x[1]][x[2]] := sol.solutions[i]!.entries[1][1];
                    products[x[2]][x[1]] := sol.solutions[i]!.entries[1][1];
                fi;
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
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewAlgebraProducts);
    
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
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(NewGramMatrix);
    
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

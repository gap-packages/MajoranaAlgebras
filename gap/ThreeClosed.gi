InstallGlobalFunction(MAJORANA_ThreeClosedSetUp,

    function(rep)
    
    local   orders, signs, unknowns, dim, new_dim, x, elts, o1, o2, k, i, j, perms, gp, pos, g;
    
    orders := [[], [1], [1,2], [1,3], [1,2,3,4]];
    signs := [[], [1], [1, 1], [1, 1], [1,-1,-1,1]];

    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(rep.algebraproducts, rep.setup);
    
    dim := Size(rep.setup.coords);
    new_dim := dim + Size(unknowns);

    for x in unknowns do 
        
        elts := rep.setup.coords{x};
        
        Add(rep.setup.coords, elts);
        
        o1 := Order(elts[1]); o2 := Order(elts[2]);
        
        k := Size(rep.setup.coords);
        
        for i in orders[o1] do 
            for j in orders[o2] do 
                Add(rep.setup.longcoords, [elts[1]^i, elts[2]^j]);
                Add(rep.setup.poslist, signs[o1][i]*signs[o2][j]*k);
                Add(rep.setup.longcoords, [elts[2]^i, elts[1]^j]);
                Add(rep.setup.poslist, signs[o1][i]*signs[o2][j]*k);
            od;
        od;
    od;
    
    for i in [1..dim] do 
        Append(rep.setup.pairorbit[i], [dim + 1 .. new_dim]*0);
        Append(rep.setup.pairconj[i], [dim + 1 .. new_dim]*0);
    od;
    
    Append(rep.setup.pairorbit, NullMat(new_dim - dim, new_dim));
    Append(rep.setup.pairconj, NullMat(new_dim - dim, new_dim));
    
    MAJORANA_Orbitals(rep.group, dim, rep.setup);
    
    gp      := AsSet(rep.group);
    perms   := List(gp, x -> MAJORANA_ThreeClosedFindVectorPermutation(x, unknowns, rep.setup));
    
    for i in [1..new_dim] do 
        for j in [Maximum(i, dim + 1) .. new_dim] do
            g := rep.setup.pairconj[i][j];
            
            pos := [Position(gp, g)];
            pos[2] := Position(gp, Inverse(g));
            
            rep.setup.pairconj[i][j] := perms{pos};
            rep.setup.pairconj[j][i] := perms{pos};
            
        od;
    od;
    
    for i in [1..Size(rep.algebraproducts)] do 
        if rep.algebraproducts[i] <> false then 
            rep.algebraproducts[i] := SparseMatrix( 1, new_dim, 
                                                    rep.algebraproducts[i]!.indices, 
                                                    rep.algebraproducts[i]!.entries, 
                                                    Rationals);
        else
            pos := Position(unknowns, rep.setup.pairreps[i]);
            rep.algebraproducts[i] := SparseMatrix( 1, new_dim, [[dim + pos]], [[1]], Rationals);
        fi;
    od;
    
    rep.algebraproducts := List(rep.algebraproducts, 
    x -> SparseMatrix(1, new_dim, x!.indices, x!.entries, Rationals));
    
    for i in [Size(rep.algebraproducts) + 1 .. Size(rep.setup.pairreps)] do 
        rep.algebraproducts[i] := false;
        rep.innerproducts[i] := false;
    od;
    
    for i in rep.setup.orbitreps do 
        rep.evecs[i] := List(rep.evecs[1], 
        x -> SparseMatrix(Nrows(x), new_dim, x!.indices, x!.entries, Rationals));
    od;
    
    rep.nullspace := SparseMatrix(  Nrows(rep.nullspace), new_dim, 
                                    rep.nullspace!.indices,
                                    rep.nullspace!.entries,
                                    Rationals);
    
    end );
    
InstallGlobalFunction( MAJORANA_ThreeClosedFindVectorPermutation, 
    
    function(g, unknowns, setup)
    
    local   new_dim,
            dim,        # size of coordinates
            j,          # loop over coordinates
            list,       # list to build permutation
            pos;
    
    new_dim := Size(setup.coords);
    dim := new_dim - Size(unknowns);
    
    list := [1..new_dim]*0;
    
    if g = () then 
        return ();
    else        
        for j in [1..dim] do 
            pos := Position(setup.longcoords,setup.coords[j]^g); 
            list[j] := setup.poslist[pos];
        od;
        
        for j in [dim + 1 .. new_dim] do 
            pos := Position(setup.longcoords, OnPairs(setup.coords[j], g));
            list[j] := setup.poslist[pos];
        od; 
        
        return list;
    fi;
    
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
    
        MAJORANA_ThreeClosedAxiomM1(new);
        
        MAJORANA_ThreeClosedFusion(new);
        
        x := MAJORANA_ThreeClosedEigenvectorsAlgebraUnknowns(new);
        
        MAJORANA_ThreeClosedResurrection(x[1],x[2],new);
        
        MAJORANA_ThreeClosedOrthogonality(new);
    
        newfalsecount := [0,0];
        
        newfalsecount[1] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.algebraproducts));
        newfalsecount[2] := Size(MAJORANA_ThreeClosedExtractUnknownProducts(new.innerproducts));
        
        if falsecount[2] <> 0 and newfalsecount[2] = 0 then 
            new.nullspace := SparseMatrix(NullspaceMat(new.innerproducts), Rationals);
        fi;
    
        Display([falsecount,newfalsecount]);
    
        if newfalsecount = [0,0] then
            Display("Success");
            return new;
        elif newfalsecount = falsecount then 
            Display( "fail" );
            return new;
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
                evecs[i][j] := MAJORANA_ConjugateVec(evecs[pos][j], g, setup);
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

    function(new)
    
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

    new_dim := Size(new.setup.coords);
    dim := new_dim - Size(new.unknowns);
    
    for i in [1..dim] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        
        for j in [1..Size(new.unknowns)] do 
            if new.innerproducts[i][dim + j] = false then 
                v := SparseMatrix(1, new_dim, [[new.unknowns[j][1]]], [[1]], Rationals);
                w := SparseMatrix(1, new_dim, [[new.unknowns[j][2]]], [[1]], Rationals);
                
                x := MAJORANA_ThreeClosedProduct(u, v, new.algebraproducts);
                
                if x <> false then
                    new.innerproducts[i][dim + j] := MAJORANA_ThreeClosedProduct(x,w,new.innerproducts);
                    new.innerproducts[dim + j][i] := MAJORANA_ThreeClosedProduct(x,w,new.innerproducts);
                fi; 
                
                if new.innerproducts[i][dim + j] = false then
                    x := MAJORANA_ThreeClosedProduct(u,w,new.algebraproducts);
                    
                    if x <> false then 
                        
                        y := MAJORANA_ThreeClosedProduct(x,v,new.innerproducts);
                        
                        if y <> false then 
                            
                          new.innerproducts[i][dim + j] := y;
                          new.innerproducts[dim + j][i] := y;
                          
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;
    
    for i in [1..Size(new.unknowns)] do 
        for j in [1..Size(new.unknowns)] do 
            for k in [[1,2],[2,1]] do 
                if new.innerproducts[dim + i][dim + j] = false then 
                    
                    u := SparseMatrix(1, new_dim, [[new.unknowns[i][k[1]]]], [[1]], Rationals);
                    v := SparseMatrix(1, new_dim, [[new.unknowns[i][k[2]]]], [[1]], Rationals);
                    w := SparseMatrix(1, new_dim, [[dim + j]], [[1]], Rationals);
                    
                    x := MAJORANA_ThreeClosedProduct(v,w,new.algebraproducts);
                    
                    if x <> false then 
                        y := MAJORANA_ThreeClosedProduct(x,u,new.innerproducts);
                       
                        if y <> false then 
                            new.innerproducts[dim + i][dim + j] := y;
                            new.innerproducts[dim + j][dim + i] := y;
                        fi;
                    fi;
                fi;
            od;
        od;
    od;
    
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedFuseEigenvectors,
    function(a, b, i, evals, fusion, innerproducts, algebraproducts, evecs)
    
    local   x,
            y,
            z,
            new_dim,
            u, v, k, j, l, m,
            new_ev,
            pos;
    
    new_dim := Size(algebraproducts);
    
    u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
    new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
    
    x := MAJORANA_ThreeClosedProduct(a, b, algebraproducts);
    
    if x <> false then 
        if evals = [2,2] then 
            y := MAJORANA_ThreeClosedProduct(a,b,innerproducts);
            
            if y <> false then 
                fusion[1] := UnionOfRows(fusion[1], x - (1/4)*u*y);
            fi;
        elif evals = [3,3] then 
            y := MAJORANA_ThreeClosedProduct(a,b,innerproducts);
            z := MAJORANA_ThreeClosedProduct(u,x,algebraproducts);
            
            if y <> false and z <> false then 
                fusion[2] := UnionOfRows(fusion[2], z - (1/32)*u*y);
                fusion[1] := UnionOfRows(fusion[1], x + (3/32)*u*y - 4*z);            
            fi;
        else
            fusion[pos] := UnionOfRows(fusion[pos],x);
        fi;
    fi;    
    
    if evecs <> false then 
    
    for j in [1..3] do 
        for k in [1..Nrows(fusion[j])] do 
            u := CertainRows(fusion[j], [k]);
            for l in Difference([1..3], [j]) do 
                for m in [1..Nrows(evecs[i][l])] do 
                    v := CertainRows(evecs[i][l], [m]);
                    
                    if not MAJORANA_ThreeClosedProduct(u, v, innerproducts) in [0, false] then 
                        Error("Orthogonality in fusion");
                    fi;
                od;
            od;
        od;
    od;
    
    fi;
    
    end);
    
InstallGlobalFunction(MAJORANA_ThreeClosedFusion,

    function(new)
    
    local   new_dim,
            bad,
            null,
            evals,
            i,
            j,
            k,
            fusion,
            a,
            b;
            
    new_dim := Size(new.setup.coords);

    for i in [1..Size(new.evecs)] do
        
        fusion := [ [], [], [] ];
        
        for j in [1..3] do 
            fusion[j] := SparseMatrix(0, new_dim, [], [], Rationals);
        od;
        
        for evals in [[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]] do
            
            for j in [1..Nrows(new.evecs[i][evals[1]])] do 
            
                a := CertainRows(new.evecs[i][evals[1]], [j]);
                
                bad := MAJORANA_ThreeClosedFindBadIndices(a, new.algebraproducts, new.setup);
                
                if bad <> [] then 
                    null := KernelMat(CertainColumns(new.evecs[i][evals[2]], bad)).relations;
                else
                    null := SparseIdentityMatrix(Nrows(new.evecs[i][evals[2]]));
                fi;
                
                for k in [1..Nrows(null)] do 
                    b := CertainRows(null, [k])*new.evecs[i][evals[2]];
                    
                    MAJORANA_ThreeClosedFuseEigenvectors(a, b, i, evals, fusion, new.innerproducts, new.algebraproducts, new.evecs);
                    
                od;
            od;
        od;
        
        for j in [1..3] do             
            new.evecs[i][j] := UnionOfRows(new.evecs[i][j], fusion[j]);
            new.evecs[i][j] := MAJORANA_BasisOfEvecs(new.evecs[i][j]);
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

    function(new)
    
    local   mat,
            vec,
            table,
            new_dim,
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
    
    new_dim := Size(new.setup.coords);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(new.algebraproducts);
    
    mat := SparseMatrix(0, Size(new_unknowns), [], [], Rationals);
    vec := SparseMatrix(0, new_dim, [], [], Rationals);
    
    for i in [1..Size(new.evecs)] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
        for j in [1..3] do 
            for k in [1..Nrows(new.evecs[i][j])] do 
                
                v := CertainRows(new.evecs[i][j], [k]);
            
                x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, new.algebraproducts);
            
                if x[1]!.indices <> [] then 
                    mat := UnionOfRows(mat, x[1]);
                    vec := UnionOfRows(vec, x[2]);
                fi;
            od;
        od;
    od;
    
    for i in [1..new_dim] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        for j in [1..Nrows(new.nullspace)] do 
            v := CertainRows(new.nullspace, [j]);
            
            x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, new.algebraproducts);
            
            if x[1]!.indices <> [] then 
                mat := UnionOfRows(mat, x[1]);
                vec := UnionOfRows(vec, x[2]);
            fi;
        od;
    od;
    
    return [mat, vec];
    
    end );
            
InstallGlobalFunction( MAJORANA_ThreeClosedSolutionProducts,

    function(mat,vec,products,new_unknowns)
    
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

    function( mat, vec, new)

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
    
    new_dim := Size(new.setup.coords);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(new.algebraproducts);
    
    for i in [1..Size(new.evecs)] do 
        
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
    
        for evals in [[1,2],[2,1],[1,3],[2,3]] do
            for j in [1..Nrows(new.evecs[i][evals[2]])] do 
                beta := CertainRows(new.evecs[i][evals[2]], [j]);
                for k in [1..Nrows(new.evecs[i][evals[1]])] do 
                    gamma := CertainRows(new.evecs[i][evals[1]], [k]);
                     
                    if MAJORANA_ThreeClosedProduct(beta,gamma,new.algebraproducts) = false then
                
                        ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
                        x := ev*MAJORANA_ThreeClosedSeparateProduct(beta, gamma, new_unknowns, new.algebraproducts);
                        
                        for l in [1..Nrows(new.evecs[i][evals[1]])] do 
                        
                            alpha := CertainRows(new.evecs[i][evals[1]], [l]);
                        
                            row := SparseMatrix(0, Size(new_unknowns), [], [], Rationals); 
                            sum := SparseMatrix(0, new_dim, [], [], Rationals);
                        
                            y := MAJORANA_ThreeClosedProduct(alpha - beta, gamma, new.algebraproducts);
                            
                            if y <> false then 
                            
                                row := row + x[1]; sum := sum + x[2];
                            
                                z := MAJORANA_ThreeClosedSeparateProduct(u, y, new_unknowns, new.algebraproducts);
                                
                                row := row + z[1]; sum := sum + z[2];
                                
                                if evals[1] = 2 then 
                                    w := MAJORANA_ThreeClosedProduct(alpha,gamma,new.innerproducts);
                                    
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

    MAJORANA_ThreeClosedSolutionProducts(mat, vec, new.algebraproducts, new_unknowns);

    end );

InstallGlobalFunction( MAJORANA_ThreeClosedOrthogonality,

    function(new)
    
    local   new_dim,
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
            
    new_dim := Size(new.setup.coords);
    
    new_unknowns := MAJORANA_ThreeClosedExtractUnknownProducts(new.innerproducts);
    
    mat := SparseMatrix(0, Size(new_unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);
    
    for i in [1..Size(new.evecs)] do 
        for j in [1..3] do 
            for k in [j + 1.. 3] do 
                for l in [1..Nrows(new.evecs[i][j])] do
                    u := CertainRows(new.evecs[i][j],[l]);
                    for m in [1..Nrows(new.evecs[i][k])] do 
                        v := CertainRows(new.evecs[i][k], [m]); 
                        if MAJORANA_ThreeClosedProduct(u,v,new.innerproducts) = false then
                            x := MAJORANA_ThreeClosedSeparateProduct(u, v, new_unknowns, new.innerproducts);
                            
                            mat := UnionOfRows(mat, x[1]);
                            vec := UnionOfRows(vec, x[2]);
                            
                        fi;
                    od;
                od;
            od;
        od;
    od;

    MAJORANA_ThreeClosedSolutionProducts(mat, vec, new.innerproducts, new_unknowns);
    
    end );        

InstallGlobalFunction(MAJORANA_ThreeClosedTestOrthogonality,

    function(new)
    
    local i, j, k, evals, u, v;
    
    for i in [1..Size(new.evecs)] do 
        for evals in [[1,2],[1,3],[2,3]]  do 
            for j in [1..Nrows(new.evecs[i][evals[1]])] do 
                u := CertainRows(new.evecs[i][evals[1]],[j]);
                for k in [1..Nrows(new.evecs[i][evals[2]])] do 
                    v := CertainRows(new.evecs[i][evals[2]],[k]);
                    
                    if not MAJORANA_ThreeClosedProduct(u, v, new.innerproducts) in [0, false] then 
                        Error("Orthogonality");
                    fi;
                od;
            od;
        od;
    od;
    
    end);

InstallGlobalFunction(MAJORANA_ThreeClosedTestEvecs,

    function(new) 
    
    local i, j, k, u, v, x, y, ev, new_dim;
    
    new_dim := Size(new.algebraproducts);
    
    for i in [1..Size(new.evecs)] do 
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        for j in [1..3] do 
            ev := MAJORANA_FusionTable[1][j + 1];
            for k in [1..Nrows(new.evecs[i][j])] do 
                v := CertainRows(new.evecs[i][j],[k]);
                x := MAJORANA_ThreeClosedProduct(u, v, new.algebraproducts);
                if x <> v*ev and x <> false then 
                    y := MAJORANA_ThreeClosedProduct(x - ev*v, x - ev*v, new.innerproducts);
                    if y <> 0 and y <> false then 
                        Error("Eigenvectors");
                    fi;
                fi;
            od;
        od;
    od;
    
    end );
                
InstallGlobalFunction(MAJORANA_ThreeClosedTestAxiomM1,

    function(new)
    
    local   new_dim, i, j, k, u, v, w, x, y, z, a;
    
    new_dim := Size(new.algebraproducts);
    
    for i in [1..new_dim] do 
        for j in [1..new_dim] do 
            for k in [1..new_dim] do 
                u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
                v := SparseMatrix(1, new_dim, [[j]], [[1]], Rationals);
                w := SparseMatrix(1, new_dim, [[k]], [[1]], Rationals);
                
                x := MAJORANA_ThreeClosedProduct(u, v, new.algebraproducts);
                
                if x <> false then 
                    y := MAJORANA_ThreeClosedProduct(x, w, new.innerproducts);
                else 
                    y := false;
                fi;
                
                z := MAJORANA_ThreeClosedProduct(v, w, new.algebraproducts);
                
                if z <> false then 
                    a := MAJORANA_ThreeClosedProduct(z, u, new.innerproducts);
                else 
                    a := false;
                fi;
                
                if a <> false and y <> false and a <> y then 
                    Error("Axiom M1");
                fi;
            od;
        od;
    od;
                
    end );
    
InstallGlobalFunction(MAJORANA_ThreeClosedTestFusion,

    function(new)
    
    local new_dim, evals, ev, u, i, j, k, fusion, a, b, x, y;
    
    new_dim := Size(new.setup.coords);

    for i in [1..Size(new.evecs)] do
    
        u := SparseMatrix(1, new_dim, [[i]], [[1]], Rationals);
        
        fusion := [ [], [], [] ];
        
        for j in [1..3] do 
            fusion[j] := SparseMatrix(0, new_dim, [], [], Rationals);
        od;
        
        for evals in [[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]] do
            for j in [1..Nrows(new.evecs[i][evals[1]])] do 
                a := CertainRows(new.evecs[i][evals[1]], [j]);
                for k in [1..Nrows(new.evecs[i][evals[2]])] do 
                    b := CertainRows(new.evecs[i][evals[2]], [k]);
                    
                    MAJORANA_ThreeClosedFuseEigenvectors(a, b, i, evals, fusion, new.innerproducts, new.algebraproducts, false);
                    
                od;
            od;
        od;  

        for j in [1..3] do 
            ev := MAJORANA_FusionTable[1][j + 1];
                
            fusion[j] := EchelonMatDestructive(fusion[j]).vectors;
            
            for k in [1..Nrows(fusion[j])] do 
                a := CertainRows(fusion[j], [k]);
                x := MAJORANA_ThreeClosedProduct(u, a, new.algebraproducts);
                if x <> ev*a and x <> false then 
                    y := MAJORANA_ThreeClosedProduct(x - ev*a, x - ev*a, new.innerproducts);
                    if y <> 0 and y <> false then 
                        Error("Algebra does not obey the fusion rules");
                    fi;
                fi;
            od;
        od;
    od;

    end );
        
InstallGlobalFunction(ThreeClosedAlgebraTest,

    function(new);
    
    MAJORANA_ThreeClosedTestEvecs(new);
    
    MAJORANA_ThreeClosedTestAxiomM1(new);
    
    MAJORANA_ThreeClosedTestOrthogonality(new);
    
    MAJORANA_ThreeClosedTestFusion(new);
    
    end );

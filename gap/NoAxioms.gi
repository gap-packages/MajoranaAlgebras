InstallGlobalFunction( MAJORANA_SetUpNoAxioms, 

    function(input, index)
    
    local   rep, t, T, i, j, k, pos, x, im, sign, gens, dim;
    
    rep         := rec( group       := input.group,
                        involutions := input.involutions,
                        shape       := input.shapes[index] );
                        
    T := rep.involutions;
    t := Size(T);
                        
    rep.setup   := rec( coords      := ShallowCopy(input.involutions),
                        longcoords  := ShallowCopy(input.involutions),
                        pairorbit   := input.pairorbit,
                        pairconj    := input.pairconj,
                        pairreps    := input.pairreps,
                        poslist     := [1..t]               );
    
    rep.setup.pairconjelts := List(input.pairconjelts, 
        x -> MAJORANA_FindVectorPermutation(x, rep.setup));
    
    # Add axes from 6A algebras    
    
    for i in [1..t] do 
        for j in [i + 1..t] do 
            k := rep.setup.pairorbit[i][j];
            if rep.shape[k] = ['6','A'] then 
                pos := [0, 0, 0, 0];
                pos[1] := Position(T, T[j]*T[i]*T[j]);
                pos[2] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                pos[3] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                pos[4] := Position(T, T[i]*T[j]*T[i]);
                pos[5] := Position(T, (T[i]*T[j])^3); 
                
                # Add the 3A axes
                
                if not SortedList([i, pos[1]]) in rep.setup.longcoords then 
                    Add(rep.setup.coords, SortedList([i, pos[1]]));
                    
                    x := [  [i, pos[1]], [i, pos[3]], [j, pos[2]], [j, pos[4]],
                            [pos[1], pos[3]], [pos[2], pos[4]]];
                    
                    Append(rep.setup.longcoords, List(x, SortedList));
                    Append(rep.setup.poslist, [1, 1, 1, 1, 1, 1]*Size(rep.setup.coords)); 
                fi;
                
                # Add the 2A axes
                
                if not SortedList([i, pos[2]]) in rep.setup.longcoords then
                    
                    x := [[i, pos[2]], [j, pos[3]], [pos[1], pos[4]]];
                    
                    Append(rep.setup.longcoords, List(x, SortedList));
                    
                    if pos[5] = fail then 
                        Add(rep.setup.coords, SortedList([i, pos[2]]));
                        Append(rep.setup.poslist, [1, 1, 1]*Size(rep.setup.coords));
                    else
                        Append(rep.setup.poslist, [1, 1, 1]*pos[5]);
                    fi;
                fi;
                
                if pos[5] <> fail and not SortedList([pos[5],i]) in rep.setup.longcoords then 
                    
                    x := Cartesian([pos[5]], [i, j, pos[1], pos[2], pos[3], pos[4]]);
                
                    Append(rep.setup.longcoords, List(x, SortedList)); 
                    Append(rep.setup.poslist, [pos[2], pos[3], pos[4], i, j, pos[1]] );
                fi;
                
            fi;
        od;
    od;
    
    # Add axes from 4B algebra
    
    for i in [1..t] do 
        for j in [i + 1..t] do 
            k := rep.setup.pairorbit[i][j];
            if rep.shape[k] = ['4','B'] then
                pos := [0, 0, 0];
                pos[1] := Position(T, T[j]*T[i]*T[j]);
                pos[2] := Position(T, T[i]*T[j]*T[i]);
                pos[3] := Position(T, (T[i]*T[j])^2);
                
                if not SortedList([i, pos[1]]) in rep.setup.longcoords then 
                    
                    Append(rep.setup.longcoords, List([[i, pos[1]], [j, pos[2]]], SortedList));

                    if pos[3] = fail then 
                        Add(rep.setup.coords, SortedList([i, pos[1]]));
                        Append(rep.setup.poslist, [1, 1]*Size(rep.setup.coords));
                    else
                        Append(rep.setup.poslist, [1, 1]*pos[3]);
                    fi;
                fi;
                
                if pos[3] <> fail and not SortedList([pos[3],i]) in rep.setup.longcoords then
                    
                    x := Cartesian([pos[3]], [i, j, pos[1], pos[2]]);
                 
                   Append(rep.setup.longcoords, List(x, SortedList));
                   Append(rep.setup.poslist, [pos[1], pos[2], i, j]);
                fi;
                
            fi;
        od;
    od;
 
    for i in [1..t] do 
        for j in [i + 1 .. t] do 
            if not [i,j] in rep.setup.longcoords then
            
                k := rep.setup.pairorbit[i][j];
                
                if rep.shape[k] = ['2','A'] then 
                
                    Add(rep.setup.coords, [i,j]);
                    Add(rep.setup.longcoords, [i,j]);
                    Add(rep.setup.poslist, Size(rep.setup.coords));
                    
                elif rep.shape[k] = ['3','A'] then
                
                    pos := Position(T, T[j]*T[i]*T[j]);
                    
                    Add(rep.setup.coords, [i,j]);
                    Append(rep.setup.longcoords, [[i,j], [i,pos], [j,pos]]);
                    Append(rep.setup.poslist, [1, 1, 1]*Size(rep.setup.coords));
                    
                elif rep.shape[k] = ['4', 'A'] then 
                
                    pos := [0,0];
                    
                    pos[1] := Position(T, T[j]*T[i]*T[j]);
                    pos[2] := Position(T, T[i]*T[j]*T[i]);
                    
                    Add(rep.setup.coords, [i,j]);
                    Add(rep.setup.longcoords, [[i, pos[1]], SortedList([i, pos[2]])]);
                    
                elif rep.shape[k] = ['5','A'] then 
                
                    pos := [0, 0, 0];
                    
                    pos[1] := Position(T, T[j]*T[i]*T[j]);
                    pos[2] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    
                    Add(rep.setup.coords, [i,j]);
                    Append(rep.setup.longcoords, [[i,j], [i, pos[1]], [i, pos[2]], [i, pos[3]]]);
                    Append(rep.setup.poslist, [1, -1, -1, 1]*Size(rep.setup.coords));
                    Append(rep.setup.longcoords, [[j, pos[1]], [j, pos[2]], [j, pos[3]]]);
                    Append(rep.setup.poslist, [1, -1, -1]*Size(rep.setup.coords));
                    Append(rep.setup.longcoords, List(Combinations(pos, 2), SortedList));
                    Append(rep.setup.poslist, [1, -1, 1]*Size(rep.setup.coords));
                fi;
            fi;
        od;
    od;
    
    dim := Size(rep.setup.coords);
    
    # Extend pairconjelts permutations
    
    for i in [1..Size(rep.setup.pairconjelts)] do 
        for j in [t + 1 .. dim] do 
            x := rep.setup.coords[j];
            im := rep.setup.pairconjelts[i]{x};
            
            sign := 1;
                
            if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
            if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;
            
            if im[1] > im[2] then im := im{[2,1]}; fi;
            
            pos := Position(rep.setup.longcoords, im);
            
            Add(rep.setup.pairconjelts[i], sign*rep.setup.poslist[pos]);
        od;
    od;
    
    for i in [1..t] do 
        Append(rep.setup.pairorbit[i], [t + 1 .. dim]*0);
        Append(rep.setup.pairconj[i], [t + 1 .. dim]*0);
    od;
    
    Append(rep.setup.pairorbit, NullMat(dim - t, dim));
    Append(rep.setup.pairconj, NullMat(dim - t, dim));
    
    gens := GeneratorsOfGroup(rep.group);
    gens := List(gens, x -> Position(AsList(rep.group), x));
    gens := rep.setup.pairconjelts{gens};
    
    MAJORANA_Orbitals(gens, t, rep.setup);
    
    return rep;
    
    end );
    
InstallGlobalFunction(MAJORANA_DihedralProductsNoAxioms,

    function(rep)
    
    local i, j, k;
    
    t := Size(rep.involutions);
    dim := Size(rep.setup.coords);

    for i in [1..Size(rep.setup.pairreps)] do 
    
        j := rep.setup.pairreps[i][1]; k := rep.setup.pairreps[i][2];
        
        if j <= t and k <= t then 
            
            if rep.shape[i] = ['1','A'] then 
                rep.algebraproducts[i] := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
                rep.innerproducts[i] := 1;
            if rep.shape[i] = ['2','A'] then 
                pos := [j, k, 0];
                pos[3] := Position(rep.setup.longcoords, [j,k]);
                pos[3] := rep.setup.poslist[pos[3]];
                
                vals := [1/8, 1/8, -1/8];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector( pos, vals, dim);

                rep.innerproducts[i] := 1/8;

            elif rep.shape[i] = ['2','B'] then

                rep.algebraproducts[i] := SparseZeroMatrix(1, dim, Rationals);

                rep.innerproducts[i] := 0;

            elif rep.shape[i] = ['3','A'] then
            
                pos := [j, k, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := rep.setup.poslist[Position(rep.setup.longcoords,[j,k])];

                vals := [1/16, 1/16, 1/32, -135/2048];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 13/256;

            elif rep.shape[i] = ['3','C'] then
            
                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                
                vals := [1/64, 1/64, -1/64];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=1/64;

            elif rep.shape[i] = ['4','A'] then

                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := rep.setup.poslist[Position(rep.setup.longcoords,T[j]*T[k])];

                vals := [3/64, 3/64, 1/64, 1/64, -3/64]; 

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 1/32;

            elif rep.shape[i] = ['4','B'] then
            
                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(rep.setup.coords,(T[j]*T[k])^2);

                vals := [1/64, 1/64, -1/64, -1/64, 1/64];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=1/64;

            elif  rep.shape[i] = ['5','A'] then
            
                pos := [j, k, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := rep.setup.poslist[Position(rep.setup.longcoords,T[j]*T[k])];

                if pos[6] < 0 then
                    pos[6] := -pos[6];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [3/128, 3/128, -1/128, -1/128, -1/128, sign];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=3/128;

            elif rep.shape[i] = ['6','A'] then
            
                pos := [j, k, 0, 0, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := Position(T,T[k]*T[j]*T[k]*T[j]*T[k]);
                pos[7] := Position(rep.setup.coords,(T[j]*T[k])^3);
                pos[8] := rep.setup.poslist[Position(rep.setup.longcoords,(T[j]*T[k])^2)];
                
                vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=5/256;
            fi;
            
            
        fi;
    
    od;
    
    end);

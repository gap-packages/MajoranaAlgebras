InstallGlobalFunction(MAJORANA_ThreeClosedSetUp,

    function(rep)
    
    local   orders, signs, unknowns, dim, new_dim, x, elts, o1, o2, k, i, j, gens, pos;
    
    orders := [[], [1], [1,2], [1,3], [1,2,3,4]];
    signs := [[], [1], [1, 1], [1,1,1], [1,-1,-1,1]];

    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(rep.algebraproducts, rep.setup);
    
    dim := Size(rep.setup.coords);
    new_dim := dim + Size(unknowns);

    MAJORANA_ThreeClosedExtendPerm(unknowns, rep.setup);

    for x in unknowns do 
        
        elts := rep.setup.coords{x};
        
        Add(rep.setup.coords, elts);
        
        o1 := Order(elts[1]); o2 := Order(elts[2]);
        
        k := Size(rep.setup.coords);
        
        # TODO fix this, no index 3 in signs[4] for example
        
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
    
    gens := GeneratorsOfGroup(rep.group);
    gens := List(gens, x -> MAJORANA_FindVectorPermutation(x, rep.setup));
    
    MAJORANA_Orbitals(gens, dim, rep.setup);
    
    for i in [1..Size(rep.algebraproducts)] do 
        if rep.algebraproducts[i] = false then 
            pos := Position(unknowns, rep.setup.pairreps[i]);
            rep.algebraproducts[i] := SparseMatrix(     1, new_dim, [[dim + pos]], 
                                                        [[1]], Rationals    );
        fi;
    od;
    
    rep.algebraproducts := List(rep.algebraproducts, x -> SparseMatrix( 1, new_dim, 
                                                                        x!.indices, 
                                                                        x!.entries, 
                                                                        Rationals)  );
    
    for i in [Size(rep.algebraproducts) + 1 .. Size(rep.setup.pairreps)] do 
        rep.algebraproducts[i] := false;
        rep.innerproducts[i] := false;
    od;
    
    for i in rep.setup.orbitreps do 
        rep.evecs[i] := List(rep.evecs[i], 
        x -> SparseMatrix(Nrows(x), new_dim, x!.indices, x!.entries, Rationals));
    od;
    
    rep.nullspace := SparseMatrix(  Nrows(rep.nullspace), new_dim, rep.nullspace!.indices,
                                    rep.nullspace!.entries, Rationals);
    end );
    
InstallGlobalFunction( MAJORANA_ThreeClosedExtendPerm,

    function(unknowns, setup)
    
    local x, im, sign, pos, i, j, dim;
    
    dim := Size(setup.coords);
    
    for i in [1..Size(setup.pairconjelts)] do 
        if Size(setup.pairconjelts[i]) <= dim + Size(unknowns) then 
            for x in unknowns do 
                im := setup.pairconjelts[i]{x};
                sign := 1;
                
                if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
                if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;
                
                if im[1] > im[2] then 
                    im := im{[2,1]};
                fi;
                
                pos := Position(unknowns, im);

                Add(setup.pairconjelts[i], sign*(dim + pos));
            od;
        fi;
    od;    
    
    for i in [1..Size(setup.conjelts)] do 
        for x in unknowns do 
            if setup.conjelts[i][1] <> () then 
                im := setup.conjelts[i]{x};
                sign := 1;
                
                if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
                if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;
                
                if im[1] > im[2] then 
                    im := im{[2,1]};
                fi;
                
                pos := Position(unknowns, im);

                Add(setup.conjelts[i], sign*(dim + pos));
            fi;
        od;
    od;
    
    end);
    
InstallGlobalFunction( ThreeClosedMajoranaRepresentation, 

    function(rep)
    
    local falsecount, newfalsecount;
    
    MAJORANA_ThreeClosedSetUp(rep);
    
    falsecount := [0,0];
    
    falsecount[1] := Size(Positions(rep.algebraproducts,false));
    falsecount[2] := Size(Positions(rep.innerproducts,false));
    
    while true do
                                
        MAJORANA_MainLoop(rep);
        
        newfalsecount := [0,0];

        newfalsecount[1] := Size(Positions(rep.algebraproducts,false));
        newfalsecount[2] := Size(Positions(rep.innerproducts,false));

        Info(InfoMajorana, 20,
            STRINGIFY( "There are ", newfalsecount[1], " unknown algebra products ") );
        Info(InfoMajorana, 20,
            STRINGIFY( "There are ", newfalsecount[2], " unknown inner products ") );

        if newfalsecount = [0,0] then
            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif newfalsecount = falsecount then
            Info( InfoMajorana, 10, "Fail" );
            return rep;
        else
            falsecount := StructuralCopy(newfalsecount);
        fi;
    od;
    
    end );

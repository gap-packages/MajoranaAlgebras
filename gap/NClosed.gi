InstallGlobalFunction(MAJORANA_NClosedSetUp,

    function(rep, index)
    
    local   unknowns, dim, new_dim, x, elts, k, i, j, gens, pos, sign, new;
    
    dim := Size(rep.setup.coords);
    
    for i in [1..dim] do 
        for j in [i..dim] do 
            if rep.setup.pairorbit[i, j] in [index, -index] then 
                Add(rep.setup.coords, [i,j]);
            fi;
        od;
    od;
    
    for x in rep.setup.pairconjelts do 
        MAJORANA_ExtendPerm(x, rep);
    od;
    
    for x in rep.setup.conjelts do 
        MAJORANA_ExtendPerm(x, rep);
    od;
    
    new_dim := Size(rep.setup.coords);
    
    for i in [1..dim] do 
        Append(rep.setup.pairorbit[i], [dim + 1 .. new_dim]*0);
        Append(rep.setup.pairconj[i], [dim + 1 .. new_dim]*0);
    od;
    
    Append(rep.setup.pairorbit, NullMat(new_dim - dim, new_dim));
    Append(rep.setup.pairconj, NullMat(new_dim - dim, new_dim));
    
    gens := GeneratorsOfGroup(rep.group);
    gens := List(gens, x -> Position(AsList(rep.group), x));
    gens := rep.setup.pairconjelts{gens};
    
    MAJORANA_Orbitals(gens, dim, rep.setup);
    
    pos := Position(rep.setup.coords, rep.setup.pairreps[index]);    
    rep.algebraproducts[index] := SparseMatrix(1, new_dim, [[pos]], [[1]], Rationals);
    
    for i in [1..Size(rep.algebraproducts)] do 
        if not rep.algebraproducts[i] in [false, fail] then 
            rep.algebraproducts[i]!.ncols := new_dim;
        fi;
    od;
    
    for i in [Size(rep.algebraproducts) + 1 .. Size(rep.setup.pairreps)] do 
        rep.algebraproducts[i] := false;
        if rep.innerproducts <> false then rep.innerproducts[i] := false; fi;
    od;
    
    for i in rep.setup.orbitreps do 
        for j in [1..3] do 
            rep.evecs[i, j]!.ncols := new_dim;
        od;
    od;
    
    # 1/32 evecs from conjugation
    
    for i in rep.setup.orbitreps do 
        pos := Position(AsList(rep.group), rep.involutions[i]);
        x := rep.setup.pairconjelts[pos];
        
        for j in [dim + 1 .. new_dim] do 
            if x[j] < 0 then sign := -1; else sign := 1; fi;
            
            if sign*x[j] <> j then 
                new := SparseMatrix(1, new_dim, [[j, sign*x[j]]], [[1, -sign]], Rationals);
                Sort(new!.indices[1]);
            
                rep.evecs[i, 3] := UnionOfRows(rep.evecs[i, 3], new); 
            fi;
        od;
    od;
    
    MAJORANA_NClosedNullspace(rep);
    
    end );

InstallGlobalFunction( MAJORANA_NClosedNullspace,

    function(rep)
    
    local i, j, v, x, pos;
    
    rep.vec!.ncols := Size(rep.setup.coords);
    rep.setup.nullspace.vectors!.ncols := Size(rep.setup.coords);
    
    for i in [1..Nrows(rep.mat)] do 
        if ForAll(rep.mat!.indices[i], x -> rep.unknowns[x] in rep.setup.coords) then
            v := CertainRows(rep.vec, [i]);
            for j in [1..Size(rep.mat!.indices[i])] do 
                x := rep.mat!.indices[i, j];
            
                pos := Position(rep.setup.coords, rep.unknowns[x]);
                
                SetEntry(v, 1, pos, -rep.mat!.entries[i, j]);
                
                rep.setup.nullspace.vectors := UnionOfRows(rep.setup.nullspace.vectors, v);
            od;
        fi;
    od;
    
    rep.setup.nullspace := ReversedEchelonMatDestructive(rep.setup.nullspace.vectors);
    
    end );
    
InstallGlobalFunction( NClosedMajoranaRepresentation, 

    function(rep)
    
    local products, unknowns;

    products := Positions(rep.algebraproducts, false);
    
    MAJORANA_NClosedSetUp(rep, products[1]);
    
    while true do
    
        unknowns := Positions(rep.algebraproducts, false);
                                
        MAJORANA_MainLoop(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );
        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.innerproducts, false)), " unknown inner products ") );

        if not false in rep.algebraproducts then 
            Info( InfoMajorana, 10, "Success" );
            return;
        fi;

        if ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            products := Filtered(products, x -> rep.algebraproducts[x] = false);
            
            if products = [] then 
                Info( InfoMajorana, 10, "Fail" );
                return;
            else
                MAJORANA_NClosedSetUp(rep, products[1]);
            fi;
        fi;
    od;
    
    end );

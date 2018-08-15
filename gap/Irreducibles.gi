InstallGlobalFunction( MAJORANA_FindMatrix, 
    
    function(p, null, dim)
    
    local   mat, i, sign;
    
    mat := NullMat(dim, dim);
    
    for i in [1..dim] do
        if p[i] > 0 then sign := 1; else sign := -1; fi;
       mat[i][sign*p[i]] := sign; 
    od;
    
    if null <> [] then 
        mat := ReduceMat(mat, null.vectors).reduced_matrix;;
        mat := mat{Positions(null.heads, 0)};
        mat := List(mat, x -> x{Positions(null.heads, 0)});
    fi;
    
    return mat;
    
    end );
    
InstallGlobalFunction( MAJORANA_Decomposition,

    function(rep)
    
    local gens, null, dim, M, D, p, o, i;
    
    dim  := Size(rep.setup.coords);;
    
    if Nrows(rep.nullspace) > 0 then 
        null := EchelonMat(ConvertSparseMatrixToMatrix(rep.nullspace));;
    else 
        null := [];
    fi;
    
    p := 1;
    o := Size(rep.group);
    
    for i in [1..1000] do 
        p := NextPrimeInt(p);
        if o mod p <> 0 then break; fi;
    od;

    gens := GeneratorsOfGroup(rep.group);;
    gens := List(gens, x -> Position(AsList(rep.group), x));; 
    gens := rep.setup.pairconjelts{gens};;
    gens := List(gens, x -> MAJORANA_FindMatrix(x, null, dim)*Z(p)^0);;
    
    M := GModuleByMats(gens, GF(p));

    D := MTX.HomogeneousComponents(M);;
    
    Display( STRINGIFY("Irreducible components have dimensions ", List(D, x -> x.component[2].dimension ), " with multiplicity ", List(D, x -> Size(x.indices) ) ) );
    
    return D;
    
    end );
    
InstallGlobalFunction( MAJORANA_MiyamotoGroup,

    function(rep)
    
    local null, gens, dim;
    
    dim  := Size(rep.setup.coords);;

    if Nrows(rep.nullspace) > 0 then 
        null := EchelonMat(ConvertSparseMatrixToMatrix(rep.nullspace));;
    else 
        null := [];
    fi;

    gens := GeneratorsOfGroup(rep.group);;
    gens := List(gens, x -> Position(AsList(rep.group), x));; 
    gens := rep.setup.pairconjelts{gens};;
    gens := List(gens, x -> MAJORANA_FindMatrix(x, null, dim));;
    
    return Group(gens);
    
    end);
    
InstallGlobalFunction( MAJORANA_Dimension, 

    function(rep)
    
    return Size(rep.setup.coords) - Nrows(rep.setup.nullspace.vectors);  
    
    end );
    

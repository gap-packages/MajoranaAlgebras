InstallGlobalFunction( MAJORANA_FindMatrix, 
    
    function(p, null, dim)
    
    local   mat, i;
    
    mat := NullMat(dim, dim);
    
    for i in [1..dim] do 
       mat[i][AbsInt(p[i])] := 1; 
    od;
    
    if null <> [] then 
        mat := ReduceMat(mat, null).reduced_matrix;;
        mat := mat{[Size(null) + 1 .. dim]};
        mat := List(mat, x -> x{[Size(null) + 1 .. dim]});
    fi;
    
    mat := mat*Lcm(List(mat, x -> _FoldList2(x, DenominatorRat, LcmInt)));
    
    return mat;
    
    end );
    
InstallGlobalFunction( MAJORANA_Decomposition,

    function(rep)
    
    local gens, null, dim, M, D, p, o, i;
    
    dim  := Size(rep.setup.coords);;
    null := ConvertSparseMatrixToMatrix(EchelonMat(rep.nullspace).vectors);;
    
    p := 1;
    o := Size(rep.group);
    
    for i in [1..1000] do 
        p := NextPrimeInt(p);
        if o mod p <> 0 then break; fi;
    od;

    gens := GeneratorsOfGroup(rep.group);;
    gens := List(gens, x -> Position(AsList(rep.group), x));; 
    gens := rep.setup.pairconjelts{gens};
    gens := List(gens, x -> MAJORANA_FindMatrix(x, null, dim)*Z(p)^0);;
    
    M := GModuleByMats(gens, GF(p));

    D := MTX.Indecomposition(M);;
    
    Info( InfoMajorana, 50, STRINGIFY("Irreducible components have dimensions ", 
        List(D, x -> x[2].dimension ) ) 
    );
    
    return D;
    
    end );
    
InstallGlobalFunction( MAJORANA_Dimension, 

    function(rep)
    
    if false in rep.innerproducts then return fail; fi;
    
    if false in rep.algebraproducts then 
        Info( InfoMajorana, "Warning: not all algebra products have been found!");
    fi;
    
    return Ncols(rep.nullspace) - Nrows(rep.nullspace);  
    
    end );
    

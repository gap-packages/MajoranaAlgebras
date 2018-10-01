InstallGlobalFunction( MAJORANA_Subalgebra,

    function( U, rep )
    
    local n, x, u, v, prod;
    
    while true do 
        n := Nrows(U);
        
        for x in Combinations( [1 .. n], 2 ) do 
            u := CertainRows( U, [x[1]]);    
            v := CertainRows( U, [x[2]]);
            
            prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
            U := UnionOfRows( U, prod );
        od;
        
        U := EchelonMatDestructive(U).vectors;
        
        if Nrows(U) = n then break; fi;
    od;
        
    return U;
    
    end );
    
InstallGlobalFunction( MAJORANA_IsJordanAlgebra,

    function( U, rep )
    
    local n, i, j, k, l, w, x, y, z, p1, p2, p3, p4, p5, p6, xz, yw, zw, xy, wx, yz, lhs, rhs;
    
    n := Nrows(U);
    
    for i in [1 .. n] do 
        w := CertainRows( U, [i] );
        for j in [1 .. n] do 
            x := CertainRows( U, [j] ); 
            for k in [1..n] do
                y := CertainRows( U, [k] );
                for k in [1..n] do
                    z := CertainRows( U, [k] );
                    
                    p1 := MAJORANA_AlgebraProduct( x, z, rep.algebraproducts, rep.setup);
                    p1 := MAJORANA_AlgebraProduct( p1, y, rep.algebraproducts, rep.setup);
                    p1 := MAJORANA_AlgebraProduct( p1, w, rep.algebraproducts, rep.setup);
                    
                    p2 := MAJORANA_AlgebraProduct( z, w, rep.algebraproducts, rep.setup);
                    p2 := MAJORANA_AlgebraProduct( p2, y, rep.algebraproducts, rep.setup);
                    p2 := MAJORANA_AlgebraProduct( p2, x, rep.algebraproducts, rep.setup);
                    
                    p3 := MAJORANA_AlgebraProduct( w, x, rep.algebraproducts, rep.setup);
                    p3 := MAJORANA_AlgebraProduct( p3, y, rep.algebraproducts, rep.setup);
                    p3 := MAJORANA_AlgebraProduct( p3, z, rep.algebraproducts, rep.setup);
                    
                    xz := MAJORANA_AlgebraProduct( x, z, rep.algebraproducts, rep.setup);
                    yw := MAJORANA_AlgebraProduct( y, w, rep.algebraproducts, rep.setup);
                    
                    p4 := MAJORANA_AlgebraProduct(xz, yw, rep.algebraproducts, rep.setup);
                    
                    zw := MAJORANA_AlgebraProduct( z, w, rep.algebraproducts, rep.setup);
                    xy := MAJORANA_AlgebraProduct( x, y, rep.algebraproducts, rep.setup);
                    
                    p5 := MAJORANA_AlgebraProduct(zw, xy, rep.algebraproducts, rep.setup);
                    
                    wx := MAJORANA_AlgebraProduct( w, x, rep.algebraproducts, rep.setup);
                    yz := MAJORANA_AlgebraProduct( y, z, rep.algebraproducts, rep.setup);
                    
                    p6 := MAJORANA_AlgebraProduct(wx, yz, rep.algebraproducts, rep.setup);
                    
                    lhs := p1 + p2 + p3;;
                    rhs := p4 + p5 + p6;
            
                    if lhs <> rhs then Error( "Algebra is not Jordan" ); fi;
                od;
            od;
        od;
    od;
    
    return true;
    
    end );

InstallGlobalFunction( MAJORANA_AdjointAction,

    function(rep, u) 
    
    local basis, adj, dim, field, v, prod, i;
    
    basis := Positions(rep.setup.nullspace.heads, 0);;
    dim := Size(rep.setup.coords);
    field := rep.field;
    
    adj := SparseMatrix(0, dim, [], [], field);
    
    for i in basis do 
        v := SparseMatrix(1, dim, [ [ i ] ], [ [ One(field) ] ], field);
        prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
        adj := UnionOfRows(adj, prod);
    od;
    
    return CertainColumns(adj, basis);
    
    end );
    
InstallGlobalFunction( MAJORANA_NaiveProduct,

    function(u, v, products)
    
    local vec, i, j;
    
    vec := SparseMatrix( 1, Ncols(u), [ [  ] ], [ [  ] ], u!.ring );
    
    for i in [1 .. Size(u!.indices[1])] do
        for j in [1 .. Size(v!.indices[1])] do
            vec := vec + u!.entries[1][i]*v!.entries[1][j]*products[u!.indices[1][i]][v!.indices[1][j]];
        od;
    od;
    
    return vec;
    
    end );

InstallGlobalFunction( MAJORANA_FindFusionTable,

    function(rep, axis, evals, basis)
    
    local field, products, dim, i, u, j, v, adj, decomp, e, table, a, b, pol, product, ev, vals, k, x;
    
    field := rep.field;;
    
    products := NullMat(Size(basis), Size(basis));;
    
    dim := Size(rep.setup.coords);;
    
    for i in [1..Nrows(basis)] do 
        
        u := CertainRows( basis, [i] );
        
        for j in [1..Size(basis)] do 
            
            v := CertainRows( basis, [j] );
            
            products[i][j] := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);;
            products[i][j] := CertainColumns(products[i][j], basis);;
        od;;
    od;;
            
    dim := Size(basis);;
    
    adj := SparseMatrix(0, dim, [], [], field);;
    
    for i in [1 .. dim] do 
        adj := UnionOfRows(adj, products[axis][i]);;
    od;;
    
    decomp := List(evals, ev -> KernelMat(adj - ev*SparseIdentityMatrix( dim, field )));;
    decomp := List(decomp, x -> x.relations);;
    
    e := Size(evals);
    
    table := NullMat( e, e);
    
    for i in [1 .. e] do 
        for j in [i .. e] do 
        
            Display( [i,j] );
        
            product := SparseMatrix( 0, dim, [], [], field);
            
            for a in [1 .. Nrows(decomp[i]) ] do 
                u := CertainRows( decomp[i], [ a ] );
                for b in [1 .. Nrows(decomp[j])] do
                    v := CertainRows( decomp[j], [ b ] );
                    
                    product := UnionOfRows(product, MAJORANA_NaiveProduct( u, v, products) );
                    
                    # if [i,j] = [1,4] then Error(); fi;
                    
                od;
            od;
            
            if ForAll(product!.entries, x -> x = []) then 
                table[i][j] := [];
                table[j][i] := [];
            else
                # product := CertainColumns( product, Positions(rep.setup.nullspace.heads, 0));
                for k in [1 .. e] do 
                    for x in Combinations([1 .. e], k) do 
                        
                        pol := One(field);
                        
                        for ev in evals{x} do 
                            pol := pol*(adj - ev*SparseIdentityMatrix(dim, field));
                        od;
                        
                        vals := product*pol;
                        
                        # if [i, j] = [1,4] and x = [4] then Error("pause"); fi;
                        
                        if ForAll(vals!.entries, x -> x = []) then 
                            table[i][j] := x;
                            table[j][i] := x;
                            break;                            
                        fi;
                    od;
                    
                    if table[i][j] <> 0 then break; fi;
                    
                od;
            fi;
        od;
    od;
    
    return table;
    
    end );
    
InstallGlobalFunction( MAJORANA_ConvertToBasis,
    
    function( basis, v )
    
    local mat, vec, x;
    
    mat := ConvertSparseMatrixToMatrix( basis );
    vec := ConvertSparseMatrixToMatrix( v );

    x := SolutionMat( mat, vec[1] );
    
    return SparseMatrix( [x], v!.ring );

    end );
    
            
            

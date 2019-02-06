##
## Functions for calculating with Majorana algebras
##

InstallGlobalFunction( MAJORANA_Dimension,

    function(rep)

    return Size(rep.setup.coords) - Nrows(rep.setup.nullspace.vectors);

    end );

InstallGlobalFunction( MAJORANA_IsComplete,

    function(rep)

    if false in rep.algebraproducts then
        return false;
    else
        return true;
    fi;

    end );

InstallGlobalFunction( MAJORANA_Eigenvectors,

    function( index, eval, rep)

    local g, i, evecs, v;

    # Find the element <g> that maps the orbit rep to index and use this to
    # find the orbit rep itself
    g := rep.setup.conjelts[index];
    i := Position(g, index);

    if not IsBound(rep.evecs[i].(String(eval))) then
        if eval = 1 then
            return SparseMatrix( 1, Size(rep.setup.coords), [ [ index ] ], [ [ 1 ] ], Rationals);
        else
            return SparseMatrix( 0, Size(rep.setup.coords), [], [], Rationals);
        fi;
    fi;

    # Conjugate the eigevectors of the orbit rep by g
    evecs := SparseMatrix( 0, Size(rep.setup.coords), [], [], Rationals);
    for v in Iterator( rep.evecs[i].(String(eval)) ) do
        evecs := UnionOfRows( evecs, MAJORANA_ConjugateVec( v, g) );
    od;

    return RemoveMatWithHeads( evecs, rep.setup.nullspace);

end );

InstallGlobalFunction( MAJORANA_Subalgebra,

    function( vecs, rep )

    local n, x, u, v, prod;

    vecs := EchelonMatDestructive(vecs).vectors;

    # Loop until subalgebra is closed under multiplication
    while true do
        n := Nrows(vecs);

        # For all pairs of vecs, find their prod and add it to the subalg
        for x in Combinations( [1 .. n], 2 ) do
            u := CertainRows( vecs, [x[1]]);
            v := CertainRows( vecs, [x[2]]);

            prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
            vecs := UnionOfRows( vecs, prod );
        od;

        vecs := EchelonMatDestructive(vecs).vectors;

        # If no more vectors have been found, the space is closed under multiplication
        if Nrows(vecs) = n then break; fi;
    od;

    return vecs;

    end );

##
## TODO comment this - what is the name of the test we use?
## Linearised Jordan something?
##

InstallGlobalFunction( MAJORANA_IsJordanAlgebra,

    function( subalg, rep )

    local n, i, j, k, l, w, x, y, z, p1, p2, p3, p4, p5, p6, xz, yw, zw, xy, wx, yz, lhs, rhs;

    n := Nrows(subalg);

    for i in [1 .. n] do
        w := CertainRows( subalg, [i] );
        for j in [1 .. n] do
            x := CertainRows( subalg, [j] );
            for k in [1..n] do
                y := CertainRows( subalg, [k] );
                for k in [1..n] do
                    z := CertainRows( subalg, [k] );

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

                    if lhs <> rhs then return false; fi;
                od;
            od;
        od;
    od;

    return true;

    end );

InstallGlobalFunction( MAJORANA_Basis,

    function(rep)

    local dim, basis, i;

    dim := Size(rep.setup.coords);

    basis := SparseMatrix( 0, dim, [], [], Rationals);

    for i in [1..dim] do
        if rep.setup.nullspace.heads[i] = 0 then
            basis := UnionOfRows(basis, SparseMatrix(1, dim, [[i]], [[1]], Rationals));
        fi;
    od;

    return basis;

end );

InstallGlobalFunction( MAJORANA_AdjointAction,

    function(axis, basis, rep)

    local adj, dim, field, v, prod, i;

    dim := Nrows(basis);
    field := Rationals;

    adj := SparseMatrix(0, dim, [], [], field);

    # For each basis elt, add its product with axis to the adjoint matrix
    # each time in terms of the basis itself
    for i in [1..dim] do
        v := CertainRows(basis, [i]);
        prod := MAJORANA_AlgebraProduct(axis, v, rep.algebraproducts, rep.setup);
        prod := MAJORANA_ConvertToBasis(basis, prod);
        adj := UnionOfRows(adj, prod);
    od;

    return adj;

    end );

InstallGlobalFunction( MAJORANA_ConvertToBasis,

    function( basis, v )

    local mat, vec, x;

    # Use the GAP function solution mat to find the vector v in terms of the basis
    mat := ConvertSparseMatrixToMatrix( basis );
    vec := ConvertSparseMatrixToMatrix( v );

    x := SolutionMat( mat, vec[1] );

    # If the vector is not in the linear span of the basis vectors then return fail
    if x = fail then
        return fail;
    else
        return SparseMatrix( [x], v!.ring );
    fi;

    end );

InstallGlobalFunction( MAJORANA_NaiveProduct,

    function(u, v, prod)

    local vec, i, j;

    vec := SparseMatrix(1, Ncols(u), [[]], [[]], u!.ring );

    for i in [ 1 .. Size(u!.indices[1]) ] do
        for j in [ 1 .. Size(v!.indices[1]) ] do
            if not prod[u!.indices[1][i], v!.indices[1][j]] in [fail, false] then
                vec := vec + u!.entries[1][i]*v!.entries[1][j]*prod[u!.indices[1][i], v!.indices[1][j]];
            else
                return prod[u!.indices[1][i], v!.indices[1][j]];
            fi;
        od;
    od;

    return vec;

end );

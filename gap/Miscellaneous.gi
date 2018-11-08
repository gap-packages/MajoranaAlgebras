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

InstallGlobalFunction( MAJORANA_Subalgebra,

    function( vecs, rep )

    local n, x, u, v, prod;

    vecs := EchelonMatDestructive(vecs).vectors;

    while true do
        n := Nrows(vecs);

        for x in Combinations( [1 .. n], 2 ) do
            u := CertainRows( vecs, [x[1]]);
            v := CertainRows( vecs, [x[2]]);

            prod := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
            vecs := UnionOfRows( vecs, prod );
        od;

        vecs := EchelonMatDestructive(vecs).vectors;

        if Nrows(vecs) = n then break; fi;
    od;

    return vecs;

    end );

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

InstallGlobalFunction( MAJORANA_AdjointAction,

    function(axis, basis, rep)

    local adj, dim, field, v, prod, i;

    dim := Nrows(basis);
    field := Rationals;

    adj := SparseMatrix(0, dim, [], [], field);

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

    mat := ConvertSparseMatrixToMatrix( basis );
    vec := ConvertSparseMatrixToMatrix( v );

    x := SolutionMat( mat, vec[1] );

    if x = fail then
        return fail;
    else
        return SparseMatrix( [x], v!.ring );
    fi;

    end );

BindGlobal( "MAJORANA_JordanFusionTable", function(eta)
            return [ [ 1, 0, eta], [0, 0, eta], [eta, eta, [1,0]] ]; end );

BindGlobal( "MAJORANA_DihedralAlgebraOfJordanType", rec() );

f := FreeGroup(2);

MAJORANA_DihedralAlgebraOfJordanType.2B := MAJORANA_DihedralAlgebras.2B;

MAJORANA_DihedralAlgebraOfJordanType.3X := function( eta )

    local rep;

    if not eta in [1/2, -1] then return false; fi;

    rep := rec(
        algebraproducts := [    SparseMatrix( 1, 2, [[1]], [[1]]),
                                SparseMatrix( 1, 2, [[2]], [[1]]),
                                SparseMatrix( 1, 3, [[1, 2]], [[eta, eta]]) ],
        innerproducts := false,
        evecs := [ [    SparseMatrix( 1, 2, [[1]], [[1]]),
                        SparseMatrix( 0, 2, [], []),
                        SparseMatrix( 1, 2, [[1, 2]], [[eta/(eta - 1), 1]] ) ],
                    [   SparseMatrix( 1, 2, [[2]], [[1]]),
                        SparseMatrix( 0, 2, [], []),
                        SparseMatrix( 1, 2, [[1, 2]], [[1, eta/(eta - 1)]] ) ] ],
        setup := rec(
            coords := [f.1, f.2],
            longcoords := [f.1, f.2],
            nullspace := rec(   vectors := SparseMatrix( 0, 3, [  ], [  ] ),
                                heads := [] ),
            orbitreps := [1, 2],
            pairreps := [[1, 1], [2, 2], [1, 2]],
            pairorbit := [[1, 3], [3, 2]],
            pairconj := [[1, 1], [1, 1]],
            pairconjelts := [[1,2]],
            poslist := [1,2] ) );

        return rep;

    end;

MAJORANA_DihedralAlgebraOfJordanType.3C := function( arg )

    local eta, phi, pi, rep;

    eta := arg[1];

    if eta = 1/2 then
        phi := arg[2];
    else
        phi := eta*(1/2);
    fi;

    pi := phi - eta - eta*phi;

    rep := rec(
        algebraproducts := [    SparseMatrix( 1, 3, [[1]], [[1]]),
                                SparseMatrix( 1, 3, [[1, 2, 3]], [[eta, eta, 1]]) ,
                                SparseMatrix( 1, 3, [[1]], [[pi]]),
                                SparseMatrix( 1, 3, [[2]], [[1]]),
                                SparseMatrix( 1, 3, [[2]], [[pi]]),
                                SparseMatrix( 1, 3, [[3]], [[pi]])],
        innerproducts := false,
        evecs := [ [    SparseMatrix( 1, 3, [[1]], [[1]]),
                        SparseMatrix( 1, 3, [[1, 3]], [[pi, -1 ]]),
                        SparseMatrix( 1, 3, [[1, 2, 3]], [[eta - phi, eta, 1]] ) ],
                    [   SparseMatrix( 1, 3, [[2]], [[1]]),
                        SparseMatrix( 1, 3, [[2, 3]], [[pi, -1 ]]),
                        SparseMatrix( 1, 3, [[1, 2, 3]], [[eta, eta - phi, 1]] ) ] ],
        setup := rec(
            coords := [f.1, f.2, [1,2]],
            longcoords := [f.1, f.2, [1,2]],
            nullspace := rec(   vectors := SparseMatrix( 0, 3, [  ], [  ] ),
                                heads := [] ),
            orbitreps := [1,2],
            pairreps := [[1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [3, 3]],
            pairorbit := [[1, 2, 3], [2, 4, 5], [3, 5, 6]],
            pairconj := [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            pairconjelts := [[1, 2, 3]],
            poslist := [1, 2, 3] ) );

        return rep;

    end;

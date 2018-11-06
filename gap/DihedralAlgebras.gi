BindGlobal( "MAJORANA_DihedralAlgebras", rec());

f := FreeGroup(2);
g := f/[f.1^2, f.2^2, (f.1*f.2)^2];

MAJORANA_DihedralAlgebras.2A :=

rec(    algebraproducts := [    SparseMatrix( 1, 3, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ 1/8, 1/8, -1/8 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ 1/8, -1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 2 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ -1/8, 1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 3 ] ], [ [ 1 ] ], Rationals ) ],
        evecs := [  rec(    ("0")     := SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ -1/4, 1, 1 ] ], Rationals ),
                            ("1/4")  := SparseMatrix( 1, 3, [ [ 2, 3 ] ], [ [ -1, 1 ] ], Rationals ) ),
                    rec(    ("0")     := SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ 1, -1/4, 1 ] ], Rationals ),
                            ("1/4")  := SparseMatrix( 1, 3, [ [ 1, 3 ] ], [ [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 1/8, 1/8, 1, 1/8, 1 ],
        involutions := [ g.1, g.2 ],
        nullspace := SparseMatrix( 0, 3, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 3 ] ],
                        coords := [ g.1, g.2, [1,2] ],
                        longcoords := [ g.1, g.2, [1,2] ],
                        orbitreps := [ 1, 2 ],
                        pairconj := [ [ 1, 1, 1 ], [ 1, 1, 1 ], [ 1, 1, 1 ] ],
                        pairconjelts := [ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ],
                        pairorbit := [ [ 1, 2, 3 ], [ 2, 4, 5 ], [ 3, 5, 6 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 2 ], [ 2, 3 ], [ 3, 3 ] ],
                        poslist := [ 1 .. 3 ] ),
        shape := [ "1A", "2A", "1A" ] );

MAJORANA_DihedralAlgebras.2B :=

rec(     algebraproducts := [   SparseMatrix( 1, 2, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 2, [ [  ] ], [ [  ] ], Rationals ),
                                SparseMatrix( 1, 2, [ [ 2 ] ], [ [ 1 ] ], Rationals ) ],
        evecs := [  rec( ("0") := SparseMatrix( 1, 2, [ [ 2 ] ], [ [ 1 ] ], Rationals ) ),
                    rec( ("0") := SparseMatrix( 1, 2, [ [ 1 ] ], [ [ 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 0, 1 ],
        involutions := [ g.1, g.2 ],
        nullspace := SparseMatrix( 0, 2, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 2 ] ],
                        coords := [ g.1, g.2 ],
                        longcoords := [ g.1, g.2 ],
                        orbitreps := [ 1, 2 ],
                        nullspace := rec(   vectors := SparseMatrix( 0, 2, [  ], [  ], Rationals ),
                                            heads := [] ),
                        pairconj := [ [ 1, 1 ], [ 1, 1 ] ],
                        pairconjelts := [ [ 1, 2 ], [ 1, 2 ], [ 1, 2 ], [ 1, 2 ] ],
                        pairorbit := [ [ 1, 2 ], [ 2, 3 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 2, 2 ] ],
                        poslist := [ 1 .. 2 ] ),
        shape := [ "1A", "2B", "1A" ] );

g := f/[f.1^2, f.2^2, (f.1*f.2)^3];

MAJORANA_DihedralAlgebras.3A :=

rec(    algebraproducts := [    SparseMatrix( 1, 4, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 4, [ [ 1, 2, 3, 4 ] ], [ [ 1/16, 1/16, 1/32, -135/2048 ] ], Rationals ),
                                SparseMatrix( 1, 4, [ [ 1, 2, 3, 4 ] ], [ [ 2/9, -1/9, -1/9, 5/32 ] ], Rationals ),
                                SparseMatrix( 1, 4, [ [ 4 ] ], [ [ 1 ] ], Rationals ) ],
        evecs := [  rec(    ("0")       := SparseMatrix( 1, 4, [ [ 1, 2, 3, 4 ] ], [ [ -10/27, 32/27, 32/27, 1 ] ], Rationals ),
                            ("1/4")    := SparseMatrix( 1, 4, [ [ 1, 2, 3, 4 ] ], [ [ -8/45, -32/45, -32/45, 1 ] ], Rationals ),
                            ("1/32")  := SparseMatrix( 1, 4, [ [ 2, 3 ] ], [ [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 13/256, 1/4, 8/5 ],
        involutions := [ g.1, g.2, g.1*g.2*g.1 ],
        nullspace := SparseMatrix( 0, 4, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 4 ], [ 2, 3, 1, 4 ], [ 3, 2, 1, 4 ] ],
                        coords := [ g.1, g.2, g.1*g.2*g.1, [1,2] ],
                        longcoords := [ g.1, g.2, g.1*g.2*g.1, [1,2], [1,3], [2,3] ],
                        orbitreps := [ 1 ],
                        pairconj := [ [ 1, 1, 3, 1 ], [ 1, 5, 6, 5 ], [ 3, 6, 6, 6 ], [ 1, 5, 6, 1 ] ],
                        pairconjelts := [ [ 1, 2, 3, 4 ], [ 2, 1, 3, 4 ], [ 1, 3, 2, 4 ], [ 3, 1, 2, 4 ], [ 2, 3, 1, 4 ], [ 3, 2, 1, 4 ] ],
                        pairorbit := [ [ 1, 2, 2, 3 ], [ 2, 1, 2, 3 ], [ 2, 2, 1, 3 ], [ 3, 3, 3, 4 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 4, 4 ] ],
                        poslist := [ 1, 2, 3, 4, 4, 4 ] ),
        shape := [ "1A", "3A" ] );

MAJORANA_DihedralAlgebras.3C :=

rec (   algebraproducts := [    SparseMatrix( 1, 3, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ 1/64, 1/64, -1/64 ] ], Rationals ) ],
        evecs := [  rec(    ("0")       := SparseMatrix( 1, 3, [ [ 1, 2, 3 ] ], [ [ -1/32, 1, 1 ] ], Rationals ),
                            ("1/4")    := SparseMatrix( 0, 3, [  ], [  ], Rationals ),
                            ("1/32")  := SparseMatrix( 1, 3, [ [ 2, 3 ] ], [ [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 1/64 ],
        involutions := [ g.1, g.2, g.1*g.2*g.1 ],
        nullspace := SparseMatrix( 0, 3, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 3 ], [ 2, 3, 1 ], [ 3, 2, 1 ] ],
                          coords := [ g.1, g.2, g.1*g.2*g.1 ],
                          longcoords := [ g.1, g.2, g.1*g.2*g.1 ],
                          nullspace := rec(   vectors := SparseMatrix( 0, 3, [  ], [  ], Rationals ),
                                            heads := [] ),
                          orbitreps := [ 1 ],
                          pairconj := [ [ 1, 1, 3 ], [ 1, 5, 6 ], [ 3, 6, 6 ] ],
                          pairconjelts := [ [ 1, 2, 3 ], [ 2, 1, 3 ], [ 1, 3, 2 ], [ 3, 1, 2 ], [ 2, 3, 1 ], [ 3, 2, 1 ] ],
                          pairorbit := [ [ 1, 2, 2 ], [ 2, 1, 2 ], [ 2, 2, 1 ] ],
                          pairreps := [ [ 1, 1 ], [ 1, 2 ] ],
                          poslist := [ 1 .. 3 ] ),
        shape := [ "1A", "3C" ] ) ;

g := f/[f.1^2, f.2^2, (f.1*f.2)^4];

MAJORANA_DihedralAlgebras.4A :=

rec (   algebraproducts := [    SparseMatrix( 1, 5, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ 3/64, 3/64, 1/64, 1/64, -3/64 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [  ] ], [ [  ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 2 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [  ] ], [ [  ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ 5/16, -1/8, -1/16, -1/8, 3/16 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ -1/8, 5/16, -1/8, -1/16, 3/16 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 5 ] ], [ [ 1 ] ], Rationals ) ],
        evecs := [  rec( ("0")    := SparseMatrix( 2, 5, [ [ 1, 2, 4, 5 ], [ 3 ] ], [ [ -1/2, 2, 2, 1 ], [ 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ -1/3, -2/3, -1/3, -2/3, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 1, 5, [ [ 2, 4 ] ], [ [ -1, 1 ] ], Rationals ) ),
                    rec( ("0")    := SparseMatrix( 2, 5, [ [ 1, 2, 3, 5 ], [ 4 ] ], [ [ 2, -1/2, 2, 1 ], [ 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ -2/3, -1/3, -2/3, -1/3, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 1, 5, [ [ 1, 3 ] ], [ [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 1/32, 0, 1, 0, 3/8, 3/8, 2 ],
        involutions := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1],
        nullspace := SparseMatrix( 0, 5, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 5 ], [ 3, 2, 1, 4, 5 ], [ 1, 4, 3, 2, 5 ] ],
                        coords := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1, [1,2] ],
                        longcoords := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1, [1,2], [1,4], [2,3], [3,4] ],
                        orbitreps := [ 1, 2 ],
                        pairconj := [ [ 1, 1, 1, 3, 1 ], [ 1, 1, 5, 1, 1 ], [ 1, 5, 5, 7, 5 ], [ 3, 1, 7, 3, 2 ], [ 1, 1, 5, 2, 1 ] ],
                        pairconjelts := [ [ 1, 2, 3, 4, 5 ], [ 1, 4, 3, 2, 5 ], [ 1, 4, 3, 2, 5 ], [ 1, 2, 3, 4, 5 ], [ 3, 2, 1, 4, 5 ], [ 3, 4, 1, 2, 5 ], [ 3, 4, 1, 2, 5 ], [ 3, 2, 1, 4, 5 ] ],
                        pairorbit := [ [ 1, 2, 3, 2, 6 ], [ 2, 4, 2, 5, 7 ], [ 3, 2, 1, 2, 6 ], [ 2, 5, 2, 4, 7 ], [ 6, 7, 6, 7, 8 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 2 ], [ 2, 4 ], [ 1, 5 ], [ 2, 5 ], [ 5, 5 ] ],
                        poslist := [ 1, 2, 3, 4, 5, 5, 5, 5 ] ),
        shape := [ "1A", "4A", "2B", "1A", "2B" ] ) ;

MAJORANA_DihedralAlgebras.4B :=

rec (   algebraproducts := [    SparseMatrix( 1, 5, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 2, 3, 4, 5 ] ], [ [ 1/64, 1/64, -1/64, -1/64, 1/64 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 3, 5 ] ], [ [ 1/8, 1/8, -1/8 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 1, 3, 5 ] ], [ [ 1/8, -1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 2 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 2, 4, 5 ] ], [ [ 1/8, 1/8, -1/8 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 2, 4, 5 ] ], [ [ 1/8, -1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 5, [ [ 5 ] ], [ [ 1 ] ], Rationals ) ],
        evecs :=  [ rec( ("0")    := SparseMatrix( 2, 5, [ [ 1, 3, 5 ], [ 1, 2, 3, 4 ] ], [ [ -1/4, 1, 1 ], [ -1/16, 1, 1/4, 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 1, 5, [ [ 3, 5 ] ], [ [ -1, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 1, 5, [ [ 2, 4 ] ], [ [ -1, 1 ] ], Rationals ) ),
                    rec( ("0")    := SparseMatrix( 2, 5, [ [ 1, 3, 5 ], [ 1, 2, 3, 4 ] ], [ [ -4, -4, 1 ], [ 4, -1/4, 4, 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 1, 5, [ [ 4, 5 ] ], [ [ -1, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 1, 5, [ [ 1, 3 ] ], [ [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 1/64, 1/8, 1/8, 1, 1/8, 1/8, 1 ],
        involutions := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1],
        nullspace := SparseMatrix( 0, 5, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 5 ], [ 3, 2, 1, 4, 5 ], [ 1, 4, 3, 2, 5 ] ],
                        coords := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1, [1,3]  ],
                      longcoords := [g.1, g.2, g.2*g.1*g.2, g.1*g.2*g.1, [1,3], [2,4] ],
                      orbitreps := [ 1, 2,],
                      pairconj := [ [ 1, 1, 1, 3, 1 ], [ 1, 1, 5, 1, 1 ], [ 1, 5, 5, 7, 5 ], [ 3, 1, 7, 3, 3 ], [ 1, 1, 5, 3, 1 ] ],
                      pairconjelts := [ [ 1, 2, 3, 4, 5 ], [ 1, 4, 3, 2, 5 ], [ 1, 4, 3, 2, 5 ], [ 1, 2, 3, 4, 5 ], [ 3, 2, 1, 4, 5 ], [ 3, 4, 1, 2, 5 ], [ 3, 4, 1, 2, 5 ], [ 3, 2, 1, 4, 5 ] ],
                      pairorbit := [ [ 1, 2, 3, 2, 4 ], [ 2, 5, 2, 6, 7 ], [ 3, 2, 1, 2, 4 ], [ 2, 6, 2, 5, 7 ], [ 4, 7, 4, 7, 8 ] ],
                      pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 1, 5 ], [ 2, 2 ], [ 2, 4 ], [ 2, 5 ], [ 5, 5 ] ],
                      poslist := [ 1, 2, 3, 4, 5, 5 ] ),
          shape := [ "1A", "4B", "2A", "1A", "2A" ] );

g := f/[f.1^2, f.2^2, (f.1*f.2)^5];

MAJORANA_DihedralAlgebras.5A :=

rec(    algebraproducts := [    SparseMatrix( 1, 6, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 6, [ [ 1, 2, 3, 4, 5, 6 ] ], [ [ 3/128, 3/128, -1/128, -1/128, -1/128, 1 ] ], Rationals ),
                                SparseMatrix( 1, 6, [ [ 1, 2, 3, 4, 5, 6 ] ], [ [ 3/128, -1/128, -1/128, 3/128, -1/128, -1 ] ], Rationals ),
                                SparseMatrix( 1, 6, [ [ 2, 3, 4, 5, 6 ] ], [ [ 7/4096, 7/4096, -7/4096, -7/4096, 7/32 ] ], Rationals ),
                                SparseMatrix( 1, 6, [ [ 1, 2, 3, 4, 5 ] ], [ [ 175/524288, 175/524288, 175/524288, 175/524288, 175/524288 ] ], Rationals ) ],
        evecs := [  rec( ("0")    := SparseMatrix( 2, 6, [ [ 1, 2, 3, 6 ], [ 1, 2, 3, 4, 5 ] ], [ [ 21/4096, -7/64, -7/64, 1 ], [ -3/32, 1, 1, 1, 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 1, 6, [ [ 2, 3, 4, 5, 6 ] ], [ [ 1/128, 1/128, -1/128, -1/128, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 2, 6, [ [ 4, 5 ], [ 2, 3 ] ], [ [ -1, 1 ], [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 3/128, 3/128, 0, 875/524288 ],
        involutions := [ g.1, g.2, g.1*g.2*g.1, g.1*g.2*g.1*g.2*g.1, g.2*g.1*g.2 ],
        nullspace := SparseMatrix( 0, 6, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 6 ], [ 2, 5, 1, 3, 4, 6 ], [ 3, 4, 1, 2, 5, 6 ], [ 4, 3, 5, 2, 1, 6 ], [ 5, 2, 4, 3, 1, 6 ] ],
                        coords := [ g.1, g.2, g.1*g.2*g.1, g.1*g.2*g.1*g.2*g.1, g.2*g.1*g.2, [1,2],  ],
                        longcoords := [ g.1, g.2, g.1*g.2*g.1, g.1*g.2*g.1*g.2*g.1, g.2*g.1*g.2, [1,2], [1,3], [1,4], [1,5], [2,3], [2,4], [2,5], [3,4], [3,5], [4,5]  ],
                        orbitreps := [ 1 ],
                        pairconj := [ [ 1, 1, 3, 1, 3, 1 ], [ 1, 8, 7, 9, 5, 8 ], [ 3, 7, 7, 9, 5, 7 ], [ 1, 9, 9, 9, 4, 9 ], [ 3, 5, 5, 4, 5, 5 ], [ 1, 8, 7, 9, 5, 1 ] ],
                        pairconjelts := [ [ 1, 2, 3, 4, 5, 6 ], [ 4, 5, 3, 1, 2, 6 ], [ 1, 3, 2, 5, 4, 6 ], [ 5, 4, 2, 1, 3, 6 ], [ 5, 2, 4, 3, 1, 6 ], [ 3, 1, 4, 5, 2, 6 ], [ 3, 4, 1, 2, 5, 6 ], [ 2, 5, 1, 3, 4, 6 ], [ 4, 3, 5, 2, 1, 6 ], [ 2, 1, 5, 4, 3, 6 ] ],
                        pairorbit := [ [ 1, 2, 2, 3, 3, 4 ], [ 2, 1, 3, 3, 2, 4 ], [ 2, 3, 1, 2, 3, 4 ], [ 3, 3, 2, 1, 2, 4 ], [ 3, 2, 3, 2, 1, 4 ], [ 4, 4, 4, 4, 4, 5 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 6 ], [ 6, 6 ] ],
                        poslist := [ 1, 2, 3, 4, 5, 6, 6, -6, -6, -6, -6, 6, 6, -6, 6] ),
        shape := [ "1A", "5A", "5A" ] ) ;

g := f/[f.1^2, f.2^2, (f.1*f.2)^6];

MAJORANA_DihedralAlgebras.6A :=

rec(    algebraproducts := [    SparseMatrix( 1, 8, [ [ 1 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 1, 2, 3, 4, 5, 6, 7, 8 ] ], [ [ 1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 1, 4, 7 ] ], [ [ 1/8, 1/8, -1/8 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 1, 5, 6, 8 ] ], [ [ 1/16, 1/16, 1/32, -135/2048 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 1, 4, 7 ] ], [ [ 1/8, -1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 2 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 2, 3, 4, 8 ] ], [ [ 1/16, 1/16, 1/32, -135/2048 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 2, 6, 7 ] ], [ [ 1/8, -1/8, 1/8 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 7 ] ], [ [ 1 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 1, 5, 6, 8 ] ], [ [ 2/9, -1/9, -1/9, 5/32 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [ 2, 3, 4, 8 ] ], [ [ 2/9, -1/9, -1/9, 5/32 ] ], Rationals ),
                                SparseMatrix( 1, 8, [ [  ] ], [ [  ] ], Rationals ), SparseMatrix( 1, 8, [ [ 8 ] ], [ [ 1 ] ], Rationals ) ],
        evecs := [  rec( ("0")    := SparseMatrix( 3, 8, [ [ 2, 3, 4, 8 ], [ 1, 4, 7 ], [ 1, 2, 3, 4, 5, 6 ] ], [ [ -32/9, -32/9, -8/9, 1 ], [ -1/4, 1, 1 ], [ -5/16, 3, 3, 3/4, 1, 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 2, 8, [ [ 1, 5, 6, 8 ], [ 4, 7 ] ], [ [ -8/45, -32/45, -32/45, 1 ], [ -1, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 2, 8, [ [ 5, 6 ], [ 2, 3 ] ], [ [ -1, 1 ], [ -1, 1 ] ], Rationals ) ),
                    rec( ("0")    := SparseMatrix( 3, 8, [ [ 2, 3, 4, 8 ], [ 1, 2, 3, 4, 5, 7 ], [ 1, 2, 3, 4, 5, 6 ] ], [ [ -10/27, 32/27, 32/27, 1 ], [ -4, 1/6, -4/3, -4/3, -4, 1 ], [ 4, -5/12, 4/3, 4/3, 4, 1 ] ], Rationals ),
                         ("1/4")  := SparseMatrix( 2, 8, [ [ 2, 3, 4, 8 ], [ 6, 7 ] ], [ [ -8/45, -32/45, -32/45, 1 ], [ -1, 1 ] ], Rationals ),
                         ("1/32") := SparseMatrix( 2, 8, [ [ 1, 5 ], [ 3, 4 ] ], [ [ -1, 1 ], [ -1, 1 ] ], Rationals ) ) ],
        group := g,
        innerproducts := [ 1, 5/256, 1/8, 13/256, 1/8, 1, 13/256, 1/8, 1, 1/4, 1/4, 0, 8/5 ],
        involutions := [ g.1, g.2, g.1*g.2*g.1, g.2*g.1*g.2*g.1*g.2, g.2*g.1*g.2, g.1*g.2*g.1*g.2*g.1],
        nullspace := SparseMatrix( 0, 8, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 8 ], [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ], [ 5, 2, 4, 3, 1, 6, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ] ],
                        coords := [ g.1, g.2, g.1*g.2*g.1, g.2*g.1*g.2*g.1*g.2, g.2*g.1*g.2, g.1*g.2*g.1*g.2*g.1, [1,4], [1,5]],
                        longcoords := [ g.1, g.2, g.1*g.2*g.1, g.2*g.1*g.2*g.1*g.2, g.2*g.1*g.2, g.1*g.2*g.1*g.2*g.1, [1,4], [1,5], [1,6], [2,3], [2,4], [2,6], [3,4], [3,5], [5,6] ],
                        orbitreps := [ 1, 2 ],
                        pairconj := [ [ 1, 1, 3, 1, 1, 3, 1, 1 ], [ 1, 1, 1, 5, 5, 11, 1, 1 ], [ 3, 1, 3, 11, 5, 11, 3, 3 ], [ 1, 5, 11, 4, 4, 2, 4, 4 ], [ 1, 5, 5, 4, 5, 4, 5, 5 ], [ 3, 11, 11, 2, 4, 11, 11, 6 ], [ 1, 1, 3, 4, 5, 11, 1, 1 ], [ 1, 1, 3, 4, 5, 6, 1, 1 ] ],
                        pairconjelts := [   [ 1, 2, 3, 4, 5, 6, 7, 8 ], [ 6, 4, 3, 2, 5, 1, 7, 8 ], [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ],
                                            [ 5, 2, 4, 3, 1, 6, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ], [ 6, 4, 3, 2, 5, 1, 7, 8 ], [ 1, 2, 3, 4, 5, 6, 7, 8 ],
                                            [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ], [ 5, 2, 4, 3, 1, 6, 7, 8 ] ],
                        pairorbit := [ [ 1, 2, 2, 3, 4, 4, 5, 10 ], [ 2, 6, 7, 7, 2, 3, 8, 11 ], [ 2, 7, 6, 7, 3, 2, 8, 11 ], [ 3, 7, 7, 6, 2, 2, 8, 11 ], [ 4, 2, 3, 2, 1, 4, 5, 10 ], [ 4, 3, 2, 2, 4, 1, 5, 10 ], [ 5, 8, 8, 8, 5, 5, 9, 12 ], [ 10, 11, 11, 11, 10, 10, 12, 13 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 5 ], [ 1, 7 ], [ 2, 2 ], [ 2, 3 ], [ 2, 7 ], [ 7, 7 ], [ 1, 8 ], [ 2, 8 ], [ 7, 8 ], [ 8, 8 ] ],
                        poslist := [ 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 7, 8, 7, 8] ),
        shape := [ "1A", "6A", "2A", "3A", "2A", "1A", "3A", "2A", "1A" ] );;

BindGlobal( "MAJORANA_DihedralAlgebrasTauMaps", rec());

for type in RecNames( MAJORANA_DihedralAlgebras ) do
    MAJORANA_DihedralAlgebrasTauMaps.(type) := StructuralCopy( MAJORANA_DihedralAlgebras.(type) );
od;

g := MAJORANA_DihedralAlgebrasTauMaps.2A.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.2A.setup.orbits := [ [ e, ], [ , e ] ];

g := MAJORANA_DihedralAlgebrasTauMaps.2B.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.2B.setup.orbits := [ [ e, ], [ , e ] ];

g := MAJORANA_DihedralAlgebrasTauMaps.3A.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.3A.setup.orbits := [ [ e, g.2*g.1, g.1 ] ];

g := MAJORANA_DihedralAlgebrasTauMaps.3C.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.3C.setup.orbits := [ [ e, g.2*g.1, g.1 ] ];

g := MAJORANA_DihedralAlgebrasTauMaps.4A.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.4A.setup.orbits := [ [ e, , g.2, ], [ , e, ,g.1] ];

g := MAJORANA_DihedralAlgebrasTauMaps.4B.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.4B.setup.orbits := [ [ e, , g.2, ], [ , e, ,g.1] ];

g := MAJORANA_DihedralAlgebrasTauMaps.5A.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.5A.setup.orbits := [ [ e, , , g.2*g.1, g.2 ], [ , e, g.1, , , ] ];

g := MAJORANA_DihedralAlgebrasTauMaps.6A.group; e := Identity(g);
MAJORANA_DihedralAlgebrasTauMaps.6A.setup.orbits := [ [ e, , , , g.2, g.2*g.1 ], [ , e, g.1, g.1*g.2, , ] ];

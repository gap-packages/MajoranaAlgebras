BindGlobal( "MAJORANA_DihedralAlgebrasAxiomM8", rec());

MAJORANA_DihedralAlgebrasAxiomM8.2A := MAJORANA_DihedralAlgebras.2A;

MAJORANA_DihedralAlgebrasAxiomM8.2B := MAJORANA_DihedralAlgebras.2B;

MAJORANA_DihedralAlgebrasAxiomM8.3A := MAJORANA_DihedralAlgebrasNoAxioms.3A;
                
MAJORANA_DihedralAlgebrasAxiomM8.3C := MAJORANA_DihedralAlgebras.3C;

MAJORANA_DihedralAlgebrasAxiomM8.4A := MAJORANA_DihedralAlgebrasNoAxioms.4A;

MAJORANA_DihedralAlgebrasAxiomM8.4B := MAJORANA_DihedralAlgebras.4B;

MAJORANA_DihedralAlgebrasAxiomM8.5A := MAJORANA_DihedralAlgebrasNoAxioms.5A;
        
MAJORANA_DihedralAlgebrasAxiomM8.6A := 

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
        evecs := [  [   SparseMatrix( 3, 8, [ [ 2, 3, 4, 8 ], [ 1, 4, 7 ], [ 1, 2, 3, 4, 5, 6 ] ], [ [ -32/9, -32/9, -8/9, 1 ], [ -1/4, 1, 1 ], [ -5/16, 3, 3, 3/4, 1, 1 ] ], Rationals ), 
                        SparseMatrix( 2, 8, [ [ 1, 5, 6, 8 ], [ 4, 7 ] ], [ [ -8/45, -32/45, -32/45, 1 ], [ -1, 1 ] ], Rationals ), 
                        SparseMatrix( 2, 8, [ [ 5, 6 ], [ 2, 3 ] ], [ [ -1, 1 ], [ -1, 1 ] ], Rationals ) ], 
                    [   SparseMatrix( 3, 8, [ [ 2, 3, 4, 8 ], [ 1, 2, 3, 4, 5, 7 ], [ 1, 2, 3, 4, 5, 6 ] ], [ [ -10/27, 32/27, 32/27, 1 ], [ -4, 1/6, -4/3, -4/3, -4, 1 ], [ 4, -5/12, 4/3, 4/3, 4, 1 ] ], Rationals ), 
                        SparseMatrix( 2, 8, [ [ 2, 3, 4, 8 ], [ 6, 7 ] ], [ [ -8/45, -32/45, -32/45, 1 ], [ -1, 1 ] ], Rationals ), 
                        SparseMatrix( 2, 8, [ [ 1, 5 ], [ 3, 4 ] ], [ [ -1, 1 ], [ -1, 1 ] ], Rationals ) ], [ false, false, false ], 
                    [ false, false, false ], [ false, false, false ], [ false, false, false ], 
                    [   SparseMatrix( 4, 8, [ [ 8 ], [ 1, 4, 7 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 5 ] ], [ [ 1 ], [ -4, -4, 1 ], [ -1, 1, -1, 1 ], [ -1, 1, -1, 1 ] ], Rationals ), 
                    SparseMatrix( 3, 8, [ [ 2, 6 ], [ 3, 5 ], [ 1, 4 ] ], [ [ -1, 1 ], [ -1, 1 ], [ -1, 1 ] ], Rationals ), 
                    SparseMatrix( 0, 8, [  ], [  ], Rationals ) ] ],
        group := Group( [ (1,2)(3,6)(4,5), (1,3)(4,6) ] ),
        innerproducts := [ 1, 5/256, 1/8, 13/256, 1/8, 1, 13/256, 1/8, 1, 1/4, 1/4, 0, 8/5 ],
        involutions := [ (1,2)(3,6)(4,5), (1,3)(4,6), (2,6)(3,5), (1,5)(2,4), (1,4)(2,3)(5,6), (1,6)(2,5)(3,4), (1,4)(2,5)(3,6) ],
        nullspace := SparseMatrix( 0, 8, [  ], [  ], Rationals ),
        setup := rec(   conjelts := [ [ 1 .. 8 ], [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ], [ 5, 2, 4, 3, 1, 6, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ] ],
                        coords := [ (1,2)(3,6)(4,5), (1,3)(4,6), (2,6)(3,5), (1,5)(2,4), (1,4)(2,3)(5,6), (1,6)(2,5)(3,4), (1,4)(2,5)(3,6), [1,5] ],
                        longcoords := [ (1,2)(3,6)(4,5), (1,3)(4,6), (2,6)(3,5), (1,5)(2,4), (1,4)(2,3)(5,6), (1,6)(2,5)(3,4), (1,4)(2,5)(3,6), [1,5], [1,6], [2,3], [2,4], [3,4], [5,6] ],
                        orbitreps := [ 1, 2, 7 ],
                        pairconj := [ [ 1, 1, 3, 1, 1, 3, 1, 1 ], [ 1, 1, 1, 5, 5, 11, 1, 1 ], [ 3, 1, 3, 11, 5, 11, 3, 3 ], [ 1, 5, 11, 4, 4, 2, 4, 4 ], [ 1, 5, 5, 4, 5, 4, 5, 5 ], [ 3, 11, 11, 2, 4, 11, 11, 6 ], [ 1, 1, 3, 4, 5, 11, 1, 1 ], [ 1, 1, 3, 4, 5, 6, 1, 1 ] ],
                        pairconjelts := [   [ 1, 2, 3, 4, 5, 6, 7, 8 ], [ 6, 4, 3, 2, 5, 1, 7, 8 ], [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ], 
                                            [ 5, 2, 4, 3, 1, 6, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ], [ 6, 4, 3, 2, 5, 1, 7, 8 ], [ 1, 2, 3, 4, 5, 6, 7, 8 ], 
                                            [ 1, 3, 2, 4, 6, 5, 7, 8 ], [ 5, 4, 2, 3, 6, 1, 7, 8 ], [ 6, 3, 4, 2, 1, 5, 7, 8 ], [ 5, 2, 4, 3, 1, 6, 7, 8 ] ],
                        pairorbit := [ [ 1, 2, 2, 3, 4, 4, 5, 10 ], [ 2, 6, 7, 7, 2, 3, 8, 11 ], [ 2, 7, 6, 7, 3, 2, 8, 11 ], [ 3, 7, 7, 6, 2, 2, 8, 11 ], [ 4, 2, 3, 2, 1, 4, 5, 10 ], [ 4, 3, 2, 2, 4, 1, 5, 10 ], [ 5, 8, 8, 8, 5, 5, 9, 12 ], [ 10, 11, 11, 11, 10, 10, 12, 13 ] ],
                        pairreps := [ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 5 ], [ 1, 7 ], [ 2, 2 ], [ 2, 3 ], [ 2, 7 ], [ 7, 7 ], [ 1, 8 ], [ 2, 8 ], [ 7, 8 ], [ 8, 8 ] ],
                        poslist := [ 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8 ] ),
        shape := [ "1A", "6A", "2A", "3A", "2A", "1A", "3A", "2A", "1A" ] );

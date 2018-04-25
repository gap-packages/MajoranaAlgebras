
InstallGlobalFunction( AXIAL_FindAll3DAlgebras,

    function(p)
    
    local algebras, eigenvalues, T, C, a12, a13, a23, a33, A, B, ev1, ev2;

    algebras := [];
    eigenvalues := [];

    T := EmptySCTable(3, Z(p)*0, "symmetric");;

    T[1][1] := [ [1], [Z(p)^0] ];;   
    T[2][2] := [ [2], [Z(p)^0] ];;

    C := Cartesian(GF(p), GF(p), GF(p));;

    for a12 in C do 
        T[1][2] := [    Filtered( [1 .. 3],   i -> a12[i] <> Z(p)*0 ),
                        Filtered( a12,        i -> i <> Z(p)*0 )        ];
        for a13 in C do 
            T[1][3] := [    Filtered( [1 .. 3],   i -> a13[i] <> Z(p)*0 ),
                            Filtered( a13,        i -> i <> Z(p)*0 )        ];    
            for a23 in C do 
                T[2][3] := [    Filtered( [1 .. 3],   i -> a23[i] <> Z(p)*0 ),
                                Filtered( a23,        i -> i <> Z(p)*0 )        ]; 
                for a33 in C do 
                    T[3][3] := [    Filtered( [1 .. 3],   i -> a33[i] <> Z(p)*0 ),
                                    Filtered( a33,        i -> i <> Z(p)*0 )        ];
                    
                    A := AlgebraByStructureConstants( GF(p), T);
                    
                    
                od;
            od;
        od;
    od;
    
    return rec( algebras := algebras, eigenvalues := eigenvalues );

    end );

InstallGlobalFunction( Is3DAxialAlgebra,
    
    function(A, p)
    
    local 
    
    B := AdjointBasis( Basis (A) );

    ev1 := Eigenvalues(GF(p), B[1]);
    ev2 := Eigenvalues(GF(p), B[2]);

    if AsSet(ev1) <> AsSet(ev2) then return false; fi;
    
    one1 := Position(ev1, Z(p)^0);
    one2 := Position(ev2, Z(p)^0);
    
    if Dimension(es1[one1]) <> 1 or  Dimension(es2[one2]) <> 1 then 
        return false;
    fi;
    
    es1 := AdjointEigenspaces(GF(p), A, B[1]);
    es2 := AdjointEigenspaces(GF(p), A, B[2]);
    
    if Size(ev1) = 2 then 
        minus1 := Difference([1, 2], one1);
        minus2 := Difference([1, 2], one2);
    else 
    
    for plus in Filtered(Combinations([1 .. 3], x -> one1 in x and x <> [1 .. 3]) do 
        
        minus := Difference( [1 .. 3], plus );
    
InstallGlobalFunction( AdjointEigenspaces,

    function( F, A, a);
    
    local es;
    
    es := Eigenspaces( F, a );
    
    for i in [1.. Size(es)] do 
        es[i] := List(Basis(es[i]), x -> List([1..3], i -> V.(i)*x[i]));
        es[i] := List(es[i], x -> Sum(x));
        es[i] := VectorSpace(F, es[i], "basis");
    od;
    
    end );
        
    
    
        
    

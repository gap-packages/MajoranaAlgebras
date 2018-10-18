
InstallGlobalFunction( AXIAL_FindAll3DAlgebras,

    function(p)
    
    local algebras, eigenvalues, T, C, a12, a13, a23, a33, A, B, ev;

    algebras := [];
    eigenvalues := [];

    T := EmptySCTable(3, Z(p)*0, "symmetric");;

    T[1, 1] := [ [1], [Z(p)^0] ];;   
    T[2, 2] := [ [2], [Z(p)^0] ];;

    C := Cartesian(GF(p), GF(p), GF(p));;

    for a12 in C do 
        T[1, 2] := [    Filtered( [1 .. 3],   i -> a12[i] <> Z(p)*0 ),
                        Filtered( a12,        i -> i <> Z(p)*0 )        ];
        for a13 in C do 
            T[1, 3] := [    Filtered( [1 .. 3],   i -> a13[i] <> Z(p)*0 ),
                            Filtered( a13,        i -> i <> Z(p)*0 )        ];    
            for a23 in C do 
                T[2, 3] := [    Filtered( [1 .. 3],   i -> a23[i] <> Z(p)*0 ),
                                Filtered( a23,        i -> i <> Z(p)*0 )        ]; 
                for a33 in C do 
                    T[3, 3] := [    Filtered( [1 .. 3],   i -> a33[i] <> Z(p)*0 ),
                                    Filtered( a33,        i -> i <> Z(p)*0 )        ];
                    
                    A := AlgebraByStructureConstants( GF(p), T);
                    
                    ev := Is3DAxialAlgebra(A, p);
                    
                    if ev <> false then 
                        Add(algebras, A);
                        Add(eigenvalues, ev);
                    fi;
                od;
            od;
        od;
    od;
    
    return rec( algebras := algebras, eigenvalues := eigenvalues );

    end );

InstallGlobalFunction( Is3DAxialAlgebra,
    
    function(A, p)
    
    local B, ev, one, es, minus, plus, x, plus_ev;
    
    B := AdjointBasis( Basis (A) );

    ev := List(B, x -> Eigenvalues(GF(p), x));
    plus_ev :=  [];
    
    if ev[1] <> ev[2] then 
        if SortedList(ev[1]) <> SortedList(ev[2]) then 
            return false;
        else
            Error("Eigenvalues are in different order");
        fi;
    fi;
    
    one := Positions(ev[1], Z(p)^0);
    
    if Size(one) <> 1 then return false; fi;
    
    es := List(B{[1,2]}, x -> AdjointEigenspaces(GF(p), A, x));
    
    if Size(ev[1]) = 2 then 
        plus  := List(es, x -> x[one[1]]);
        minus := List(es, x -> x[Difference([1, 2], one)[1]]);
        
        return  AXIAL_TestFusion(GF(p), A, plus[1], minus[1]) and 
                AXIAL_TestFusion(GF(p), A, plus[2], minus[2]);
    fi;
        
    for x in Filtered(Combinations([1 .. 3]), x -> one[1] in x and not x in [[1,2,3],[]]) do 
        plus  := List(es, y -> Union(y{x}));
        minus := List(es, y -> Union(y{Difference([1,2,3], x)}));
        
        if AXIAL_TestFusion(GF(p), A, plus[1], minus[1]) and AXIAL_TestFusion(GF(p), A, plus[2], minus[2]) then 
            Add(plus_ev, x);
        fi;
    od;
    
    if plus_ev = [] then return false; fi;
    
    return [ev, plus_ev];
    
    end );
    
InstallGlobalFunction( AXIAL_TestFusion,

    function(F, V, plus, minus)
    
    local u, v;
    
    for u in plus do 
        for v in plus do 
            if not u*v in VectorSpace(F, plus, "basis") then 
                return false;
            fi;
        od;
    od;
    
    for u in plus do 
        for v in minus do 
            if not u*v in VectorSpace(F, minus, "basis") then 
                return false;
            fi;
        od;
    od;
    
    for u in minus do 
        for v in minus do 
            if not u*v in VectorSpace(F, plus, "basis") then 
                return false;
            fi;
        od;
    od;
    
    return true;
    
    end );
    
InstallGlobalFunction( AdjointEigenspaces,

    function( F, V, a)
    
    local es, i;
    
    es := Eigenspaces( F, a );
    
    for i in [1.. Size(es)] do 
        es[i] := List(Basis(es[i]), x -> List([1..3], i -> V.(i)*x[i]));
        es[i] := List(es[i], x -> Sum(x));
        # es[i] := VectorSpace(F, es[i], "basis");
    od;
    
    return es;
    
    end );
        

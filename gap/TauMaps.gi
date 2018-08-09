InstallGlobalFunction( MAJORANA_TauShapes,

    function(tau)
    
    local G, d, orbs, shape, i, x, D, o1, o2, n, gph, comps, binaries, unknowns3X, j, k, shapeslist;
    
    G := Group(tau);
    d := Size(tau);
    
    orbs := MAJORANA_OrbitalsT(G, [1..d]);
    
    shape := [1..Size(orbs.pairreps)]*0;
    
    for i in [1..Size(orbs.pairreps)] do 
    
        x := orbs.pairreps[i];
    
        D := Group(tau{x});
        
        o1 := Orbit(D, x[1]);
        o2 := Orbit(D, x[2]);
        
        if x[1] in o2 then 
            n := Size(o1);
            
            if n = 1 then 
                shape[i] := "1A";
            elif n = 3 then 
                shape[i] := "3X";
            elif n = 5 then 
                shape[i] := "5A";
            else
                Error("This is not a valid tau map");
            fi;
        else
            n := Size(o1) + Size(o2);
            
            if n = 2 then   
                shape[i] := "2X";
            elif n = 4 then 
                shape[i] := "4X";
            elif n = 6 then 
                shape[i] := "6A";                                
            else
                Error("This is not a valid tau map");
            fi;
        fi;
    od;
    
    # Inclusions of 2A and 3A in 6A algebras
    
    for i in [1..Size(orbs.pairreps)] do
        if shape[i][1] = '6' then 
        
            for x in [orbs.pairreps[i], Reversed(orbs.pairreps[i])] do 
        
                k := orbs.pairorbit[x[1]][x[1]^tau[x[2]]];
                
                shape[k] := "3A";
            
                k := orbs.pairorbit[x[1]][x[2]^(tau[x[1]]*tau[x[2]])];
                
                shape[k] := "2A";
            od;
        fi;
    od;
    
    # Inclusions of 2X into 4X
    
    gph := NullMat(Size(orbs.pairreps), 0);
    
    for i in [1..Size(orbs.pairreps)] do 
        if shape[i][1] = '4' then 
        
            for x in [orbs.pairreps[i], Reversed(orbs.pairreps[i])] do
                Add(gph[i], orbs.pairorbit[x[1]][x[1]^tau[x[2]]]);
            od;
        fi;
    od;
    
    comps := AutoConnectedComponents(gph);
    
    comps := Filtered(comps, x -> shape[x[1]][2] = 'X' and shape[x[1]][1] in ['2','4']);
    
    # Put in any known (4B, 2A) pairs
    
    for x in comps do 
        if ForAny(shape{x}, y -> y[2] = 'A') then 
            for i in x do 
                if shape[i][1] = '2' then 
                    shape[i] := "2A";
                else
                    shape[i] := "4B";
                fi;
            od;
        fi;
    od;
    
    unknowns3X := Filtered([1..Size(shape)], i -> shape[i] = "3X");
    
    binaries := AsList(FullRowSpace(GF(2), Size(comps) + Size(unknowns3X) ));
    
    shapeslist := [];
    
    for i in [1..Size(binaries)] do
    
        for j in [1 .. Size(comps)] do 
            
            if binaries[i][j] = 1*Z(2) then
                for k in comps[j] do 
                    if shape[k][1] = '2' then 
                        shape[k] := "2A";
                    else
                        shape[k] := "4B";
                    fi;
                od;
            else
                for k in comps[j] do 
                    if shape[k][1] = '2' then 
                        shape[k] := "2B";
                    else
                        shape[k] := "4A";
                    fi;
                od;
            fi;
        od;
        
        for j in [1 .. Size(unknowns3X)] do
            
            if binaries[i][Size(comps) + j] = 1*Z(2) then
                shape[unknowns3X[j]] := "3A";
            else
                shape[unknowns3X[j]] := "3C";
            fi;
        od;
        
        Add(shapeslist,ShallowCopy(shape));
        
    od;
    
    return rec( group       := Group(tau),
                tau         := tau,
                shapes      := shapeslist,
                pairreps    := orbs.pairreps,
                pairorbit   := orbs.pairorbit,
                pairconj    := orbs.pairconj,
                pairconjelts := orbs.pairconjelts     );
    
    end );

        

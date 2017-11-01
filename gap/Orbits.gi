 
InstallGlobalFunction(MAJORANA_Orbits,

    function(G, t, SetUp)
    
    local   gens,
            dim,
            i,
            pnt,
            d,
            orb,
            gen,
            elts,
            count,
            p,
            q,
            h,
            g,
            pos;
    
    gens := GeneratorsOfGroup(G);
    dim := Size(SetUp.coords);
    SetUp.conjelts := [1..dim]*0;
    SetUp.orbitreps := [];
    
    for i in [1..dim] do 
        if SetUp.conjelts[i] = 0 then 
            
            Add(SetUp.orbitreps, i);
            SetUp.conjelts[i] := ();
        
            pnt := Immutable(SetUp.coords[i]);
            
            d := NewDictionary(pnt, false);
            
            orb := [pnt];
            elts := [()];
            
            count := 0;
            
            AddDictionary(d, pnt);
            
            for p in orb do 
            
                count := count + 1;
                h := elts[count];
                
                for gen in gens do 
                
                    q := OnPoints(p, gen);
                    g := h*gen;
                    
                    MakeImmutable(q);
                    
                    if not KnowsDictionary(d,q) then 
                        Add(orb, q);
                        AddDictionary(d,q);
                        Add(elts, g);
                        
                        pos := Position(SetUp.longcoords, q);
                        pos := SetUp.poslist[pos];
                        
                        if pos < 0 then pos := -pos; fi;
                        
                        SetUp.conjelts[pos] := g;
                    fi;
                od;
            od;
        fi;
    od;
    
    SetUp.conjelts := DuplicateFreeList(SetUp.conjelts);
    
    SetUp.orbitreps := [  Filtered(SetUp.orbitreps, x -> x <= t),
                          Filtered(SetUp.orbitreps, x -> x > t)   ];
                        
    end ); 
   
InstallGlobalFunction(MAJORANA_Orbitals,

    function(G,t,SetUp)
    
    local   dim,
            gens,
            i,j,k,l,
            pnt,
            d,
            orb,
            elts,
            count,
            o,
            p,
            h,
            g,
            q,
            x,
            y,
            table,
            gen,
            pos_1,
            pos_2,
            sign;
    
    dim := Size(SetUp.coords);
    
    table := [[], [1], [1,2], [1,3], [1,2,3,4]];

    gens := GeneratorsOfGroup(G);

    for i in [1..dim] do 
        for j in [Maximum(i,t + 1)..dim] do 

            if SetUp.pairorbit[i][j] = 0 then 
                
                Add(SetUp.pairreps, [i,j]);
                
                pnt := Immutable(SetUp.coords{[i,j]});
                
                d := NewDictionary(pnt, false);
                
                orb := [pnt];
                elts := [()];
                
                AddDictionary(d,pnt);
                
                count := 0;
                
                x := List(pnt, Order);
                y := Size(SetUp.pairreps);
                
                MAJORANA_AddPowers(table{x}, pnt, y, (), d, SetUp);
                
                for p in orb do 
                    
                    count := count + 1;
                    h := elts[count];
                    
                    for gen in gens do 
                    
                        q := OnPairs(p,gen);
                        g := h*gen;
                        
                        MakeImmutable(q);
                        
                        if not KnowsDictionary(d,q) then 
                        
                            Add( orb, q );
                            Add( elts, g);
                            
                            MAJORANA_AddPowers(table{x}, q, y, g, d, SetUp);
                        fi;
                    od;
                od;                
            fi;
        od;
    od;
    
    for i in [1 .. dim] do 
    
        o := Order(SetUp.coords[i]);
        
        for j in [Maximum(t, i) + 1 .. dim] do 
            if Order(SetUp.coords[j]) = 5 then 
            
                sign := 1;
                
                g := SetUp.pairconj[i][j];
                x := SetUp.pairorbit[i][j];
                
                p := SetUp.coords{SetUp.pairreps[x]};
                q := OnPairs(p, g);
                
                if o = 5 and SetUp.coords[i] in [q[1]^2,q[1]^3] then 
                    sign := -sign;
                fi;
                
                if SetUp.coords[j] in [q[2]^2,q[2]^3] then 
                    sign := -sign;
                fi;  
            
                SetUp.pairorbit[i][j] := sign*SetUp.pairorbit[i][j];
                SetUp.pairorbit[j][i] := sign*SetUp.pairorbit[j][i];
                              
            fi;
            
        od;
    od;
    
    end );
    
InstallGlobalFunction( MAJORANA_AddPowers,

    function(table, pnt, k, g, d, SetUp)
    
    local   i,
            j,
            pos_1,
            pos_2;
    
    for i in table[1] do 
        for j in table[2] do
        
            AddDictionary(d,[pnt[1]^i,pnt[2]^j]);
        
            pos_1 := Position(SetUp.longcoords, pnt[1]^i);
            pos_1 := SetUp.poslist[pos_1];
            
            if pos_1 < 0 then pos_1 := -pos_1; fi;
            
            pos_2 := Position(SetUp.longcoords, pnt[2]^j);
            pos_2 := SetUp.poslist[pos_2];
            
            if pos_2 < 0 then pos_2 := -pos_2; fi;
            
            SetUp.pairorbit[pos_1][pos_2] := k;
            SetUp.pairorbit[pos_2][pos_1] := k;
    
            SetUp.pairconj[pos_1][pos_2] := g;
            SetUp.pairconj[pos_2][pos_1] := g;
        
        od;
    od;
    
    end );
            
InstallGlobalFunction( MAJORANA_OrbitalsT,

    function(G, T)
    
    local   gens,
            t, 
            i,
            j,
            k,
            pairorbit,
            pairconj,
            pairreps,
            pnt,
            d,
            gen,
            orb,
            orbs,
            elts,
            count,
            p,
            q,
            h,
            g,
            res,
            pos_1,
            pos_2;
            
    gens := GeneratorsOfGroup(G);
    t := Size(T);
    
    pairorbit := NullMat(t,t);
    pairconj  := NullMat(t,t);
    pairreps  := [];
    orbs      := [];
    
    
    for i in [1..t] do 
        for j in [i..t] do 
            if pairorbit[i][j] = 0 then 
                
                Add(pairreps, [i,j]);
                
                k := Size(pairreps);
                
                pairorbit[i][j] := k;
                pairorbit[j][i] := k;
                
                pairconj[i][j] := ();
                pairconj[j][i] := ();
                
                pnt := Immutable(T{[i,j]});
                
                d := NewDictionary(pnt, false);
                
                orb := [pnt];
                elts := [()];
                
                AddDictionary(d,pnt);
                
                count := 0;
                
                for p in orb do 
                    
                    count := count + 1;
                    h := elts[count];
                    
                    for gen in gens do 
                    
                        q := OnPairs(p,gen);
                        g := h*gen;
                        
                        MakeImmutable(q);
                        
                        if not KnowsDictionary(d,q) then 
                        
                            Add( orb, q );
                            AddDictionary(d,q);
                            Add( elts, g);
                
                            pos_1 := Position(T,q[1]);
                            pos_2 := Position(T,q[2]);
                                
                            pairorbit[pos_1][pos_2] := k;
                            pairorbit[pos_2][pos_1] := k;
                            
                            pairconj[pos_1][pos_2] := g;
                            pairconj[pos_2][pos_1] := g;
                            
                        fi;
                    od;
                od;
                
                Add(orbs, Immutable(orb));
                
            fi;
        od;
    od;
    
    res := rec( pairorbit := pairorbit,
                pairconj  := pairconj,
                pairreps  := pairreps,
                orbitals  := orbs   );
                
    return res;
    
    end );

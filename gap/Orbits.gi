 
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
                
                SetUp.pairorbit[i][j] := y;
                SetUp.pairorbit[j][i] := y;
        
                SetUp.pairconj[i][j] := ();
                SetUp.pairconj[j][i] := ();
                
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
                            
                            for k in table[x[1]] do 
                                for l in table[x[2]] do
                                    sign := 1;
                                
                                    pos_1 := Position(SetUp.longcoords, q[1]^k);
                                    pos_1 := SetUp.poslist[pos_1];
                                    
                                    if pos_1 < 0 then 
                                        pos_1 := -pos_1;
                                        sign := - sign;
                                    fi;
                                    
                                    pos_2 := Position(SetUp.longcoords, q[2]^l);
                                    pos_2 := SetUp.poslist[pos_2];
                                    
                                    if pos_2 < 0 then 
                                        pos_2 := -pos_2;
                                        sign := - sign;
                                    fi;
                         
                                    if x = [5,5] then
                                        if k in [2,3] then  
                                            sign := -sign;
                                        fi;
                                        if l in [2,3] then 
                                            sign := -sign;
                                        fi;
                                    elif x[2] = 5 and k in [2,3] then 
                                        sign := -sign;
                                    fi;
                                    
                                    SetUp.pairorbit[pos_1][pos_2] := sign*y;
                                    SetUp.pairorbit[pos_2][pos_1] := sign*y;
                            
                                    SetUp.pairconj[pos_1][pos_2] := g;
                                    SetUp.pairconj[pos_2][pos_1] := g;
                                
                                od;
                            od;
                        fi;
                    od;
                od;                
            fi;
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

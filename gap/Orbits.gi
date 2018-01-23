 
InstallGlobalFunction(MAJORANA_Orbits,

    function(G, t, setup)
    
    local   gens,
            conjelts,
            orbitreps,
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
    conjelts := [1..t]*0;
    orbitreps := [];
    
    for i in [1..t] do 
        if conjelts[i] = 0 then 
            
            Add(orbitreps, i);
            conjelts[i] := ();
        
            pnt := Immutable(setup.coords[i]);
            
            d := NewDictionary(pnt, false);
            
            orb := [pnt];
            elts := [()];
            
            count := 0;
            
            AddDictionary(d, pnt);
            
            for p in orb do 
            
                count := count + 1;
                h := elts[count];
                
                for gen in gens do 
                
                    if IsRowVector(p) then 
                        q := OnTuples(p, gen);
                    else
                        q := OnPoints(p, gen);
                    fi;
                    g := h*gen;
                    
                    MakeImmutable(q);
                    
                    if not KnowsDictionary(d,q) then 
                        Add(orb, q);
                        AddDictionary(d,q);
                        Add(elts, g);
                        
                        pos := Position(setup.longcoords, q);
                        pos := setup.poslist[pos];
                        
                        if pos < 0 then pos := -pos; fi;
                        
                        conjelts[pos] := g;
                    fi;
                od;
            od;
        fi;
    od;
    
    conjelts := DuplicateFreeList(conjelts);
    
    return rec( conjelts := conjelts,
                orbitreps := orbitreps  );
                        
    end ); 
   
InstallGlobalFunction(MAJORANA_Orbitals,

    function(G,t,setup)
    
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
            y,
            table,
            gen,
            pos_1,
            pos_2,
            pos,
            sign;
    
    dim := Size(setup.coords);

    gens := GeneratorsOfGroup(G);

    for i in [1..dim] do 
        for j in [Maximum(i,t + 1)..dim] do 

            if setup.pairorbit[i][j] = 0 then 
                
                Add(setup.pairreps, [i,j]);
                
                pnt := Immutable(Flat(setup.coords{[i,j]}));
                
                orb := [pnt];
                elts := [()];
                
                count := 0;
                
                y := Size(setup.pairreps);
                
                setup.pairorbit[i][j] := y;
                setup.pairorbit[j][i] := y;
                
                setup.pairconj[i][j] := 1;
                setup.pairconj[j][i] := 1;
                
                for p in orb do 
                    
                    count := count + 1;
                    h := elts[count];
                    
                    for gen in gens do 
                    
                        q := OnTuples(p,gen);
                        g := h*gen;
                        
                        if Size(q) = 2 then 
                            pos_1 := Position(setup.longcoords,q[1]);
                            pos_2 := Position(setup.longcoords,q[2]);
                        elif Size(q) = 3 then 
                            pos_1 := Position(setup.longcoords, q[1]);
                            pos_2 := Position(setup.longcoords, q{[2,3]});
                        else
                            pos_1 := Position(setup.longcoords, q{[1,2]});
                            pos_2 := Position(setup.longcoords, q{[3,4]});
                        fi;
                        
                        pos_1 := setup.poslist[pos_1];
                        pos_2 := setup.poslist[pos_2];
                        
                        sign := 1;
                        
                        if pos_1 < 0 then pos_1 := -pos_1; sign := -sign; fi;
                        if pos_2 < 0 then pos_2 := -pos_2; sign := -sign; fi;
                        
                        if setup.pairorbit[pos_1][pos_2] = 0 then 
                        
                            Add( orb, q );
                            Add( elts, g);
                            
                            setup.pairorbit[pos_1][pos_2] := sign*y;
                            setup.pairorbit[pos_2][pos_1] := sign*y;
                            
                            pos := Position(setup.pairconjelts, g);
                            
                            if pos = fail then 
                                Add(setup.pairconjelts, g);
                                pos := Size(setup.pairconjelts);
                            fi;
                            
                            setup.pairconj[pos_1][pos_2] := pos;
                            setup.pairconj[pos_2][pos_1] := pos;
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
            setup,
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
            pos,
            pos_1,
            pos_2;
            
    gens := GeneratorsOfGroup(G);
    t := Size(T);
    
    setup := rec();
    
    setup.pairorbit := NullMat(t,t);
    setup.pairconj  := NullMat(t,t);
    setup.pairreps  := [];
    setup.orbitals  := [];
    setup.pairconjelts := [()];
    
    for i in [1..t] do 
        for j in [i..t] do 
            if setup.pairorbit[i][j] = 0 then 
                
                Add(setup.pairreps, [i,j]);
                
                k := Size(setup.pairreps);
                
                setup.pairorbit[i][j] := k;
                setup.pairorbit[j][i] := k;
                
                setup.pairconj[i][j] := 1;
                setup.pairconj[j][i] := 1;
                
                pnt := Immutable(T{[i,j]});
                
                orb := [pnt];
                elts := [()];
                
                count := 0;
                
                for p in orb do 
                    
                    count := count + 1;
                    h := elts[count];
                    
                    for gen in gens do 
                    
                        q := OnPairs(p,gen);
                        g := h*gen;
                        
                        MakeImmutable(q);

                        pos_1 := Position(T,q[1]);
                        pos_2 := Position(T,q[2]);
                        
                        if setup.pairorbit[pos_1][pos_2] = 0 then 
                        
                            Add( orb, q );
                            Add( elts, g);
                                
                            setup.pairorbit[pos_1][pos_2] := k;
                            setup.pairorbit[pos_2][pos_1] := k;
                            
                            pos := Position(setup.pairconjelts, g);
                            
                            if pos = fail then 
                                Add(setup.pairconjelts, g);
                                pos := Size(setup.pairconjelts);
                            fi;
                            
                            setup.pairconj[pos_1][pos_2] := pos;
                            setup.pairconj[pos_2][pos_1] := pos;
                            
                        fi;
                    od;
                od;
                
                Add(setup.orbitals, Immutable(orb));
                
            fi;
        od;
    od; 
                
    return setup;
    
    end );

 
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
                
                    q := OnPoints(p, gen);
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
                
                setup.pairconj[i][j] := ();
                setup.pairconj[j][i] := ();
                
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
                            pos_1 := Position(setup.longcoordinates, q[1]);
                            pos_2 := Position(setup.longcoordinates, q{[2,3]});
                        else
                            pos_1 := Position(setup.longcoordinates, q{[1,2]});
                            pos_2 := Position(setup.longcoordinates, q{[3,4]});
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
                            
                            setup.pairconj[pos_1][pos_2] := g;
                            setup.pairconj[pos_2][pos_1] := g;
                            
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
                        
                        if pairorbit[pos_1][pos_2] = 0 then 
                        
                            Add( orb, q );
                            Add( elts, g);
                                
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
    
    for i in [1..t] do 
        for j in [1..t] do
            k := pairorbit[i][j];
            g := pairconj[i][j];
            p := T{pairreps[k]}; 
            q := OnPairs(p,g);
            if not T{[i,j]} in [q,Reversed(q)] 
                then Error("pause");
            fi; 
        od; 
    od; 
    
    res := rec( pairorbit := pairorbit,
                pairconj  := pairconj,
                pairreps  := pairreps,
                orbitals  := orbs   );
                
    return res;
    
    end );

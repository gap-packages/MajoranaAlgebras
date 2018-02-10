InstallGlobalFunction(MAJORANA_Orbits,

    function(gens, t, setup)
    
    local   conjelts,
            orbitreps,
            i,
            orb,
            gen,
            elts,
            count,
            dim,
            p,
            q,
            h,
            g;

    conjelts := [1..t]*0;
    orbitreps := [];
    
    dim := Size(setup.coords);
    
    for i in [1..t] do 
        if conjelts[i] = 0 then 
            
            Add(orbitreps, i);
            conjelts[i] := [1..dim];
            
            orb := [i];
            elts := [[1..dim]];
            
            count := 0;
            
            for p in orb do 
            
                count := count + 1;
                h := elts[count];
                
                for gen in gens do 
                    q := gen[p];
                    g := SP_Product(h,gen);
                    
                    if q < 0 then q := -q; fi;
                    
                    if conjelts[q] = 0 then 
                        Add(orb, q);
                        Add(elts, g);
                        conjelts[q] := g;
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

    function(gens,t,setup)
    
    local   dim,
            i,j,
            orb,
            elts,
            count,
            p,
            h,
            g,
            q,
            y,
            gen,
            pos,
            sign;
    
    dim := Size(setup.coords);

    for i in [1..dim] do 
        for j in [Maximum(i,t + 1)..dim] do 

            if setup.pairorbit[i][j] = 0 then 
                
                Add(setup.pairreps, [i,j]);
                
                orb := [[i,j]];
                elts := [[1..dim]];
                
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
                    
                        q := gen{p};
                        
                        if q[1] < 0 then q[1] := -q[1]; fi;
                        if q[2] < 0 then q[2] := -q[2]; fi;
                        
                        if setup.pairorbit[q[1]][q[2]] = 0 then 
                        
                            g := SP_Product(h,gen);
                        
                            Add( orb, q );
                            Add( elts, g);
                            
                            if Product(g{orb[1]}) < 0 then 
                                sign := -1;
                            else
                                sign := 1;
                            fi;
                            
                            setup.pairorbit[q[1]][q[2]] := sign*y;
                            setup.pairorbit[q[2]][q[1]] := sign*y;
                            
                            pos := Position(setup.pairconjelts, g);
                            
                            if pos = fail then 
                                Add(setup.pairconjelts, g);
                                pos := Size(setup.pairconjelts);
                            fi;
                            
                            setup.pairconj[q[1]][q[2]] := pos;
                            setup.pairconj[q[2]][q[1]] := pos;
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
    setup.pairconjelts := [Identity(G)];
    
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
                elts := [Identity(G)];
                
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

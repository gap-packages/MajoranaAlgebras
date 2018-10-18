InstallGlobalFunction(MAJORANA_Orbits,

    function(gens, t, setup)
    
    local   conjelts,
            orbitreps,
            i, j,
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
                    
                    g := [];
    
                    for j in h do 
                        if j > 0 then 
                            Add(g, gen[j]);
                        else
                            Add(g, -gen[-j]);
                        fi;
                    od;
                    
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
    
    local   dim, i, j, orb;
    
    dim := Size(setup.coords);

    for i in [1..dim] do 
        for j in [Maximum(i,t + 1)..dim] do 

            if setup.pairorbit[i, j] = 0 then 
                
                orb := MAJORANA_NewOrbital([i,j], gens, setup);
                
                if IsBound(setup.orbitals) then 
                    Add(setup.orbitals, Immutable(orb));
                fi;
            fi;
        od;
    od;
    
    end );
    
InstallGlobalFunction( MAJORANA_NewOrbital,

    function( pnt, gens, setup)
    
    local orb, elts, count, y, p, h, q, z, g, i, j, k, pos, im, new, old, dim, sign, gen;

    Add(setup.pairreps, pnt);
    
    dim := Size(setup.coords);
                
    orb := [ pnt ];
    elts := [ [1..dim] ];
    
    count := 0;
    
    y := Size(setup.pairreps);
    
    setup.pairorbit[pnt[1], pnt[2]] := y;
    setup.pairorbit[pnt[2], pnt[1]] := y;
    
    setup.pairconj[pnt[1], pnt[2]] := 1;
    setup.pairconj[pnt[2], pnt[1]] := 1;
    
    for p in orb do 
        
        count := count + 1;
        h := elts[count];
        
        for gen in gens do 
        
            q := gen{p};
            
            if q[1] < 0 then q[1] := -q[1]; fi;
            if q[2] < 0 then q[2] := -q[2]; fi;
            
            z := setup.pairorbit[q[1], q[2]];
            
            if z = 0 then 

                g := [];

                for i in h do 
                    if i > 0 then 
                        Add(g, gen[i]);
                    else
                        Add(g, -gen[-i]);
                    fi;
                od;
            
                Add( orb, q );
                Add( elts, g);
                
                if Product(g{orb[1]}) < 0 then 
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                setup.pairorbit[q[1], q[2]] := sign*y;
                setup.pairorbit[q[2], q[1]] := sign*y;
                
                pos := Position(setup.pairconjelts, g);
                
                if pos = fail then 
                    Add(setup.pairconjelts, g);
                    pos := Size(setup.pairconjelts);
                fi;
                
                setup.pairconj[q[1], q[2]] := pos;
                setup.pairconj[q[2], q[1]] := pos;
            fi;
        od;
    od;

    return orb;

    end );

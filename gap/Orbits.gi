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
    
    local   dim,
            i,j,k,
            orb,
            elts,
            count,
            p,
            h,
            g,
            q,
            y,z,
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
                        
                        z := setup.pairorbit[q[1]][q[2]];
                        
                        if z = 0 then 

                            g := [];
    
                            for k in h do 
                                if k > 0 then 
                                    Add(g, gen[k]);
                                else
                                    Add(g, -gen[-k]);
                                fi;
                            od;
                        
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
                            
                        elif z <> y then # This orbit has already been found

                            Remove(setup.pairreps);
                            
                            g := setup.pairconj[q[1]][q[2]];
                            g := setup.pairconjelts[g];
                            
                            g := SP_Product( g, SP_Inverse(gen));
                            g := SP_Product( g, SP_Inverse(h));
                            
                            # Fix existing orbit info 
                            
                            elts := List( elts, x -> SP_Product( g, x) );
                            
                            for k in [1 .. Size(orb)] do 
                                
                                setup.pairorbit[orb[k][1]][orb[k][2]] := z;
                                setup.pairorbit[orb[k][2]][orb[k][1]] := z;
                                
                                pos := Position(setup.pairconjelts, elts[k]);
                            
                                if pos = fail then 
                                    Add(setup.pairconjelts, elts[k]);
                                    pos := Size(setup.pairconjelts);
                                fi;
                                
                                setup.pairconj[orb[k][1]][orb[k][2]] := pos;
                                setup.pairconj[orb[k][2]][orb[k][1]] := pos;
                            od;
                            
                            y := z;
                            h := elts[count];
                            
                        fi;
                    od;
                od; 
                
                if IsBound(setup.orbitals) then 
                    Add(setup.orbitals, Immutable(orb));
                fi;
            fi;
        od;
    od;
    
    end );

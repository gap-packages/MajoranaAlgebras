InstallGlobalFunction( "MAJORANA_AllEmbeddings",

    function(rep)
    
    local   unknowns,
            x,
            y,
            gens,
            subgp,
            T,
            subrep,
            ex,
            i,
            embs,
            emb,
            g;
            
    unknowns  := Positions(rep.algebraproducts, false);
    
    for x in unknowns do
        if rep.algebraproducts[x] = false then 
    
            y := rep.setup.pairreps[x];
            gens := rep.setup.coords{y};
            
            subgp := Subgroup(rep.group,gens);
            
            if Size(subgp) < Size(rep.group) then 
                T := Intersection(subgp, rep.involutions);

                if T <> [] and Size(Group(T)) = Size(subgp) then 
                    embs := IsomorphicSubgroups(rep.group,subgp);
                    ex := ShapesOfMajoranaRepresentationAxiomM8(subgp,T);
                    
                    for i in [1..Size(ex.shapes)] do 
                        if IsSubsetSet(AsSet(rep.shape), AsSet(ex.shapes[i])) then 
                            
                            Display(StructureDescription(subgp));
                            
                            subrep := MajoranaRepresentation(ex,i);
                            
                            for emb in embs do 
                                g := MAJORANA_CheckEmbedding(rep, subrep, emb);
                                
                                if g <> false then 
                                
                                    g := CompositionMapping2(emb, g);
                                    
                                    MAJORANA_Embed(rep,subrep,g);
                                fi;
                            od;    
                        fi;    
                    od;                
                fi;
            fi;  
        fi;
    od;
        
    end );

InstallGlobalFunction( "MAJORANA_CheckEmbedding",

    function(rep, subrep,emb)
    
    local   check,
            conj,
            g,
            aut_emb,
            list,
            aut;
    
    aut := AutomorphismGroup(subrep.group);
    list := [1..Size(subrep.shape)];
    
    for g in aut do 
    
        aut_emb := CompositionMapping2(emb, g);
    
        check := IsSubsetSet(AsSet(rep.involutions), AsSet(Image(aut_emb, subrep.involutions)));
        
        if check and ForAll(list, i -> MAJORANA_CheckShape(rep, subrep, aut_emb, i)) then 
            return g;
        fi;
    od;

    return false;
        
    end);
    
InstallGlobalFunction( "MAJORANA_CheckShape",

    function(rep, subrep, aut_emb, i)
    
    local   x,
            im1,
            im2,
            pos1,
            pos2,
            k;
    
    if subrep.shape[i][1] in [2,3,4] then 
        x := subrep.setup.pairreps[i];
        
        im1 := Image(aut_emb, subrep.setup.coords[x[1]]);
        im2 := Image(aut_emb, subrep.setup.coords[x[2]]);
        
        pos1 := Position(rep.setup.coords, im1);
        pos2 := Position(rep.setup.coords, im2);
        
        k := rep.setup.pairorbit[pos1][pos2];
        
        if subrep.shape[i] <> rep.shape[k] then 
            return false;
        else
            return true;
        fi; 
    fi;
    
    return true;
    
    end );


InstallGlobalFunction( "MAJORANA_Embed",
    
    function(rep, subrep, emb)
    
    local   i,
            j,
            x,
            im1,
            im2,
            pos1,
            pos2,
            k,
            g,
            perm,
            v,
            sign;
    
    for i in [1..Size(subrep.algebraproducts)] do 
    
        sign := 1;
    
        x := subrep.setup.pairreps[i]; 
        
        im1 := Image(emb, subrep.setup.coords[x[1]]);
        im2 := Image(emb, subrep.setup.coords[x[2]]);
        
        pos1 := Position(rep.setup.longcoords, im1);
        pos2 := Position(rep.setup.longcoords, im2);
        
        pos1 := rep.setup.poslist[pos1];
        pos2 := rep.setup.poslist[pos2];
        
        if pos1 < 0 then 
            sign := -sign;
            pos1 := -pos1;
        fi;        
        
        if pos2 < 0 then 
            sign := -sign;
            pos2 := -pos2;
        fi;
        
        k := rep.setup.pairorbit[pos1][pos2];
        
        if k < 0 then 
            sign := -sign;
            k := -k;
        fi;
        
        g := rep.setup.pairconj[pos1][pos2][2];
        
        if rep.algebraproducts[k] = false then 
            v := MAJORANA_ImageVector(subrep.algebraproducts[i], emb, rep, subrep);
            rep.algebraproducts[k] := sign*MAJORANA_ConjugateVector(v,g,rep.setup);
        fi;
        
        if rep.innerproducts[k] = false then 
            rep.innerproducts[k] := sign*subrep.innerproducts[i];
        fi;
        
    od;
    
    for i in subrep.setup.orbitreps do 
        
        im1 := Image(emb, subrep.setup.coords[i]);
        
        for g in List(rep.setup.conjelts, x -> Inverse(x[1])) do 
        
            pos1 := Position(rep.setup.coords, im1^g);
            
            if pos1 in rep.setup.orbitreps then 
                break; 
            fi;
        od;
        
        perm := MAJORANA_FindVectorPermutation(g, rep.setup); 
        
        for j in [1..3] do 
            for v in subrep.evecs[i][j] do 
                
                im2 := MAJORANA_ImageVector(v, emb, rep, subrep);
                im2 := MAJORANA_ConjugateVector(im2, perm, rep.setup);
                
                Add(rep.evecs[pos1][j], im2);
            od;
            
            rep.evecs[pos1][j] := ShallowCopy(BaseMat(rep.evecs[pos1][j]));
            
        od;
    od;    
    
    end );
    
InstallGlobalFunction( "MAJORANA_ImageVector",

    function(v, emb, rep, subrep)
    
    local   dim,
            i,
            im,
            res,
            sign,
            pos;
    
    dim := Size(rep.setup.coords);
    
    res := [1..dim]*0;
    
    for i in [1..Size(v)] do 
    
        if v[i] <> 0 then 
    
            sign := 1;;
            
            im := Image(emb, subrep.setup.coords[i]);
            
            pos := Position(rep.setup.longcoords, im);
            
            pos := rep.setup.poslist[pos];
            
            if pos < 0 then 
                sign := -sign;
                pos := -pos;
            fi;
            
            res[pos] := sign*v[i];
        fi;
    od;
    
    return res;
    
    end );

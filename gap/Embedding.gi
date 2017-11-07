InstallGlobalFunction( "MAJORANA_FindEmbeddings",

    function(rep, subrep)
    
    local   embs;
    
    embs := IsomorphicSubgroups(rep.group, subrep.group);

    end);
    

InstallGlobalFunction( "MAJORANA_Embed",
    
    function(rep, subrep, emb)
    
    local   i,
            x,
            im1,
            im2,
            pos1,
            pos2,
            k,
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
        
        #if rep.algebraproducts[k] = false then 
        #    rep.algebraproducts[k] := 
        #        sign*MAJORANA_ImageVector(subrep.algebraproducts[i], emb, rep, subrep);
        #fi;
        
        if rep.innerproducts[k] = false then 
            rep.innerproducts[k] := sign*subrep.innerproducts[i];
        fi;
        
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

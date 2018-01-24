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
            j,
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
                            
                            Info(   InfoMajorana, 10, 
                                    STRINGIFY("Constructing subrep of ", StructureDescription(subgp) ) );
                            
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
    
    for i in rep.setup.orbitreps do 
        for j in [1..3] do 
            if Nrows(rep.evecs[i][j]) > 0 then 
                rep.evecs[i][j] := MAJORANA_BasisOfEvecs(rep.evecs[i][j]);
            fi;
        od;
    od;
    
    if Nrows(rep.nullspace) > 0 then 
        rep.nullspace := ShallowCopy(MAJORANA_BasisOfEvecs(rep.nullspace));
    fi;
    
    end );

InstallGlobalFunction( "MAJORANA_CheckEmbedding",

    function(rep, subrep,emb)
    
    local   check,
            conj,
            g,
            aut_emb,
            im_gp,
            im_inv,
            list,
            aut;
    
    aut := AutomorphismGroup(subrep.group);
    list := [1..Size(subrep.shape)];
    
    for g in aut do 
    
        aut_emb := CompositionMapping2(emb, g);
        
        im_inv := AsSet(Image(aut_emb, subrep.involutions));
        im_gp  := Image(aut_emb, subrep.group);
    
        check := IsSubsetSet(AsSet(rep.involutions), im_inv); 
        
        #and Size(Intersection(rep.involutions, im_gp)) = Size(im_inv);
        
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
    
    if subrep.shape[i][1] in ['2','3','4'] then 
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
    
InstallGlobalFunction( "MAJORANA_EmbedKnownRep",

    function( rep, subrep)
    
    local   embs,
            i,
            g,
            emb;
            
    
    embs := IsomorphicSubgroups(rep.group, subrep.group);
    
    for i in [1..Size(embs)] do 
    
        g := MAJORANA_CheckEmbedding(rep, subrep, embs[i]);
        
        if g <> false then 
            emb := CompositionMapping2(embs[i], g);

            MAJORANA_Embed(rep, subrep,  emb);
        fi;
    od;
    
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
        
        g := SP_Inverse(rep.setup.pairconjelts[rep.setup.pairconj[pos1][pos2]]);
        
        if rep.algebraproducts[k] = false then 
            if subrep.algebraproducts[i] <> false then 
                v := MAJORANA_ImageVector(subrep.algebraproducts[i], emb, rep, subrep);
                rep.algebraproducts[k] := sign*MAJORANA_ConjugateVec(v,g,rep.setup);
            fi;
        fi;
        
        if rep.innerproducts[k] = false then 
            if subrep.innerproducts[i] <> false then 
                rep.innerproducts[k] := sign*subrep.innerproducts[i];
            fi;
        fi;
        
    od;
    
    im1 := MAJORANA_ImageVector( subrep.nullspace, emb, rep, subrep);

    UnionOfRows(rep.nullspace, im1);
    
    for i in subrep.setup.orbitreps do 
        
        im1 := Image(emb, subrep.setup.coords[i]);
        pos1 := Position(rep.setup.coords, im1);
        
        for g in List(rep.setup.conjelts, x -> SP_Inverse(x)) do
            if g[pos1] in rep.setup.orbitreps then 
                break; 
            fi;
        od;
        
        for j in [1..3] do 

            im2 := MAJORANA_ImageVector(subrep.evecs[i][j], emb, rep, subrep);
            im2 := MAJORANA_ConjugateVec(im2, g, rep.setup);
            
            UnionOfRows(rep.evecs[g[pos1]][j], im2);
        od;
    od;    
    
    end );
    
InstallGlobalFunction( "MAJORANA_ImageVector",

    function(mat, emb, rep, subrep)
    
    local   i,
            j,
            im,
            res,
            sign,
            pos,
            nrows,
            ncols,
            indices,
            entries;
    
    nrows := Nrows(mat);
    ncols := Ncols(mat);
    
    res := SparseZeroMatrix(nrows, ncols, Rationals);
    
    indices := IndicesOfSparseMatrix(mat);
    entries := EntriesOfSparseMatrix(mat);
    
    for i in [1..nrows] do 
        for j in [1..Size(indices[i])] do 
    
            sign := 1;;
            
            im := Image(emb, subrep.setup.coords[indices[i][j]]);
            
            pos := Position(rep.setup.longcoords, im);
            
            pos := rep.setup.poslist[pos];
            
            if pos < 0 then 
                sign := -sign;
                pos := -pos;
            fi;
            
            SetEntry(res, i, pos, sign*entries[i][j]);
        od;
    od;
    
    return res;
    
    end );

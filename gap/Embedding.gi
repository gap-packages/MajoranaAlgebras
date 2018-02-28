InstallGlobalFunction( "MAJORANA_AllEmbeddings",

    function(rep, axioms)
    
    local   x, gens, subgp, T, subrep, ex, i, j;
    
    for x in Positions(rep.algebraproducts, false) do
        if rep.algebraproducts[x] = false then 
    
            gens := ShallowCopy(rep.setup.coords{rep.setup.pairreps[x]});
            
            if IsRowVector(gens[1]) then gens[1] := rep.setup.coords{gens[1]}; fi;
            if IsRowVector(gens[2]) then gens[2] := rep.setup.coords{gens[2]}; fi;
            
            gens := Flat(gens);
            
            subgp := Subgroup(rep.group,gens);
            
            if Size(subgp) < Size(rep.group) then 
                T := Intersection(subgp, rep.involutions);

                if T <> [] and Size(Group(T)) = Size(subgp) then 
                
                    ex := ShapesOfMajoranaRepresentationAxiomM8(subgp,T);
                    
                    for i in [1..Size(ex.shapes)] do 
                        if IsSubsetSet(AsSet(rep.shape), AsSet(ex.shapes[i])) then 
                            
                            Info(   InfoMajorana, 10, 
                                    STRINGIFY("Constructing subrep of ", StructureDescription(subgp) ) );
                            
                            subrep := MajoranaRepresentation(ex,i,axioms);
                            
                            MAJORANA_EmbedKnownRep(rep, subrep);    
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
    
InstallGlobalFunction( "MAJORANA_MaximalSubgps",

    function(rep,axioms)
    
    local   max, inv, i, j, ex, subrep;
    
    max := ConjugacyClassesMaximalSubgroups(rep.group);
    max := List(max, Representative);
    max := Filtered(max, x -> Size(x) > 12);
    
    inv := List(max, x -> Intersection(AsList(x), rep.involutions));
    inv := Filtered(inv, x -> x <> []);
    
    max := List(inv, Group);
    
    for i in [1..Size(max)] do 
        ex := ShapesOfMajoranaRepresentationAxiomM8(max[i], inv[i]);
        for j in [1..Size(ex.shapes)] do 
            if IsSubsetSet(AsSet(rep.shape), AsSet(ex.shapes[j])) then 
                            
                Info(   InfoMajorana, 10, 
                        STRINGIFY("Constructing subrep of ", StructureDescription(ex.group) ) );
                
                subrep := MajoranaRepresentation(ex,j,axioms);
                
                MAJORANA_EmbedKnownRep(rep, subrep);    
            fi;    
        od;
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

    function(rep, subrep, emb)
    
    local   aut, g, aut_emb, im, i, x, pos1, pos2, k;
    
    aut := AutomorphismGroup(subrep.group);
    
    for g in aut do 
    
        aut_emb := CompositionMapping2(emb, g);
        
        im := AsSet(Image(aut_emb, subrep.involutions));
        
        if not IsSubsetSet(AsSet(rep.involutions), im) then 
            return false; 
        fi;
        
        for i in [1..Size(subrep.shape)] do 
            if subrep.shape[i][1] in ['2','3','4'] then 
                
                x := subrep.setup.pairreps[i];
        
                im := OnPairs(subrep.setup.coords{x}, aut_emb);
                
                pos1 := Position(rep.setup.coords, im[1]);
                pos2 := Position(rep.setup.coords, im[2]);
                
                k := rep.setup.pairorbit[pos1][pos2];
                
                if subrep.shape[i] <> rep.shape[k] then 
                    return false;
                fi; 
            fi;
        od;
        
        return g;
    od;

    return false;
        
    end);
    
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
    
    local   i, im, j, k, g, v, sign;
    
    emb := MAJORANA_FindPerm(emb, rep, subrep);
    
    for i in [1..Size(subrep.algebraproducts)] do 
    
        sign := 1;
    
        im := emb{subrep.setup.pairreps[i]};
        
        if im[1] < 0 then 
            sign := -sign; im[1] := -im[1];
        fi;        
        
        if im[2] < 0 then 
            sign := -sign; im[2] := -im[2];
        fi;
        
        k := rep.setup.pairorbit[im[1]][im[2]];
        
        if k < 0 then 
            sign := -sign; k := -k;
        fi;
        
        g := SP_Inverse(rep.setup.pairconjelts[rep.setup.pairconj[im[1]][im[2]]]);
        
        if rep.algebraproducts[k] = false then 
            if subrep.algebraproducts[i] <> false then 
                v := MAJORANA_ImageVector(subrep.algebraproducts[i], emb, rep, subrep);
                rep.algebraproducts[k] := sign*MAJORANA_ConjugateVec(v,g);
            fi;
        fi;
        
        if rep.innerproducts[k] = false then 
            if subrep.innerproducts[i] <> false then 
                rep.innerproducts[k] := sign*subrep.innerproducts[i];
            fi;
        fi;
        
    od;
    
    im := MAJORANA_ImageVector( subrep.nullspace, emb, rep, subrep);

    rep.nullspace := UnionOfRows(rep.nullspace, im);
    
    for i in subrep.setup.orbitreps do 
        
        k := emb[i];
        
        for g in List(rep.setup.conjelts, x -> SP_Inverse(x)) do
            if g[k] in rep.setup.orbitreps then 
                break; 
            fi;
        od;

        for j in [1..3] do 
            if Nrows(subrep.evecs[i][j]) > 0 then 
                im := MAJORANA_ImageVector(subrep.evecs[i][j], emb, rep, subrep);
                im := MAJORANA_ConjugateVec(im, g);
                rep.evecs[g[k]][j] := UnionOfRows(rep.evecs[g[k]][j], im);
            fi;
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
            
            im := emb[indices[i][j]];
            
            if im < 0 then 
                sign := -sign;
                im := -im;
            fi;
            
            SetEntry(res, i, im, sign*entries[i][j]);
        od;
    od;
    
    return res;
    
    end );
    
InstallGlobalFunction( MAJORANA_Image,

    function(rep, subrep, emb, x)
    
    local y, pos;
    
    if IsRowVector(x) then
        y := subrep.setup.coords{x};
        
        pos := List(y, i -> Position(rep.setup.longcoords, Image(emb, i))); # TODO fix signs 
        
        return SortedList(rep.setup.poslist{pos});
    else 
        return Image(emb, x);
    fi;
    
    end );

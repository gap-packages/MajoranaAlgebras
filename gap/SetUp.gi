InstallGlobalFunction(ShapesOfMajoranaRepresentationAxiomM8,
    
    function(G,T)
    
    local   t,              # size of T
            i,              # indices
            j,
            k,
            x,              # result of orbitals
            orbs,           # orbitals on T
            shape,          # one shape
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            unknowns,       # indices of 3X axes
            pos,            # positions
            Binaries,       # used to loop through options for shapes
            shapeslist,     # final list of shapes
            input;          #
    
    t := Size(T);

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;
    
    # Construct orbitals of  on T x T
    
    orbs := MAJORANA_OrbitalsT(G,T);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(orbs.pairreps))[1];

    RepsSquares6A := [];
    unknowns := [];;

    for i in [1..Size(orbs.pairreps)] do
    
        x := T{orbs.pairreps[i]};
        
        if Order(x[1]*x[2]) = 1 then 
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            shape[i]:="2A";
        elif Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            shape[i]:="2B";
        elif Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
            Add(unknowns,i);
        elif Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            shape[i]:="4A";
        elif Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            shape[i]:="4B";
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(RepsSquares6A,(x[1]*x[2])^2);
        else 
            Error("This is not a 6-transposition group");
        fi;
    od;

    # Check for inclusions of 2A and 3A in 6A

    for i in unknowns do
        if ForAny(orbs.orbitals[i], x -> x[1]*x[2] in RepsSquares6A) then
            shape[i]:="3A";;
            unknowns := Difference(unknowns, [i]);            
        fi;
    od;
    
    Binaries := AsList(FullRowSpace(GF(2),Size(unknowns)));
    
    shapeslist := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do
        
        for j in [1..Size(unknowns)] do
            k := unknowns[j];
            if Binaries[i][j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;            
        od;
        
        Add(shapeslist,ShallowCopy(shape));
        
    od;
    
    input  := rec(  group       := G,
                    involutions := T,
                    orbitals    := orbs.orbitals,
                    shapes      := shapeslist,
                    pairreps    := orbs.pairreps,
                    pairorbit   := orbs.pairorbit,
                    pairconj    := orbs.pairconj,
                    pairconjelts := orbs.pairconjelts     );
    
    return input;

    end );
    
InstallGlobalFunction(ShapesOfMajoranaRepresentation,
    
    function(G,T)
    
    local   t,              # size of T
            i,              # indices
            j,
            k,
            x,              # result of orbitals
            ind,            # list of indices  
            orbs,           # orbitals on T
            shape,          # one shape
            RepsSquares4X,  # (ts)^2 where o(ts) = 4
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            RepsCubes6A,    # (ts)^3 where o(ts) = 6
            gph,            # digraph of 2X, 4X inclusions
            cc,             # connected components of gph
            pos,            # positions
            Binaries,       # used to loop through options for shapes
            shapeslist,     # final list of shapes
            input;          #
    
    t := Size(T);
    
    # Construct orbitals of  on T x T
    
    orbs := MAJORANA_OrbitalsT(G,T);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(orbs.pairreps))[1];

    RepsSquares4X := [];
    RepsSquares6A := [];
    RepsCubes6A := [];
    
    ind := NullMat(6,0);;

    for i in [1..Size(orbs.pairreps)] do
    
        x := T{orbs.pairreps[i]};
        
        if Order(x[1]*x[2]) = 1 then 
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 then
            shape[i]:="2X";
            Add(ind[2],i);
        elif Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
            Add(ind[3],i);
        elif Order(x[1]*x[2]) = 4 then
            shape[i]:="4X";
            Add(ind[4],i);
            Add(RepsSquares4X, (x[1]*x[2])^2);
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(ind[6],i);
            Add(RepsSquares6A,(x[1]*x[2])^2);
            Add(RepsCubes6A,(x[1]*x[2])^3);
        fi;
    od;
    
    # Check for inclusions of 2X in 4X
    
    gph := NullMat(Size(orbs.orbitals), 0);
    
    for i in ind[2] do 
        for x in orbs.orbitals[i] do
            pos := Positions(RepsSquares4X, x[1]*x[2]);
            
            if pos <> [] then
                Append(gph[i], ind[4]{pos} + Size(ind[2]));
            fi;
        od;
    od;
    
    gph := List(gph, DuplicateFreeList);
    
    cc := AutoConnectedComponents(gph);

    # Check for inclusions of 2A and 3A in 6A

    for i in ind[3] do
        if ForAny(orbs.orbitals[i], x -> x[1]*x[2] in RepsSquares6A) then
            shape[i]:="3A";;
            ind[3] := Difference(ind[3], [i]);            
        fi;
    od;
    
    for i in ind[2] do 
        if ForAny(orbs.orbitals[i], x -> x[1]*x[2] in RepsCubes6A) then 
        
            shape[i]:="2A";;
            
            for x in cc do 
                if i in x then 
                    for j in Intersection(ind[2],x) do 
                        shape[j] := "2A";
                    od;
                    for j in Intersection(ind[4],x - Size(ind[2])) do 
                        shape[j] := "4B";
                    od;
                fi;
            od; 
             
            cc := Difference(cc, [x]);
        fi;
    od;
    
    cc := Filtered(cc, x -> Size(Intersection(x,ind[2])) > 0);

    Binaries := AsList(FullRowSpace(GF(2),Size(ind[3]) + Size(cc)));
    
    shapeslist := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do
        
        for j in [1..Size(ind[3])] do
            k:=ind[3][j];
            if Binaries[i][j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;            
        od;
        
        for j in [1 .. Size(cc)] do 
            
            if Binaries[i][j + Size(ind[3])] = 1*Z(2) then 
                for k in Intersection(ind[2],cc[j]) do 
                    shape[k] := "2A";
                od;
                for k in Intersection(ind[4],cc[j] - Size(ind[2])) do 
                    shape[k] := "4B";
                od;
            else
                for k in Intersection(ind[2],cc[j]) do 
                    shape[k] := "2B";
                od;
                for k in Intersection(ind[4],cc[j] - Size(ind[2])) do 
                    shape[k] := "4A";
                od;
            fi;
        od;            
        
        Add(shapeslist,ShallowCopy(shape));
        
    od;
    
    input  := rec(  group       := G,
                    involutions := T,
                    orbitals    := orbs.orbitals,
                    shapes      := shapeslist,
                    pairreps    := orbs.pairreps,
                    pairorbit   := orbs.pairorbit,
                    pairconj    := orbs.pairconj,
                    pairconjelts := orbs.pairconjelts     );
    
    return input;

    end );

InstallGlobalFunction( MAJORANA_SetUp,
    
    function(input, index, algebras)
    
    local   rep,
            t,              # size of T
            x,              # temporary variables
            i,              # index variables
            j,
            k,
            g,
            y,
            res,
            subrep,
            emb,
            gens,
            dim,            # size of coordinates
            s;              # number of orbits of G on T x T
            
    t := Size(input.involutions);
    
    rep         := rec( group       := input.group,
                        involutions := input.involutions,
                        shape       := input.shapes[index] );
                        
    rep.setup   := rec( coords      := ShallowCopy(input.involutions),
                        longcoords  := ShallowCopy(input.involutions),
                        poslist     := [1..t]               ); 

    gens := GeneratorsOfGroup(input.group);
    gens := List(gens, x -> MAJORANA_FindPerm(x, rep, rep));

    x := MAJORANA_Orbits(gens, t, rep.setup);

    rep.setup.conjelts := x.conjelts;
    rep.setup.orbitreps := x.orbitreps;    

    for i in [1..t] do
        for j in [i + 1 .. t] do 
            k := input.pairorbit[i][j];
            if rep.shape[k] in ["4B", "6A"] then 
                MAJORANA_RecordCoords(rep.involutions{[i,j]}, rep.shape[k], rep, algebras);
            fi;
         od;
    od;
    
    for i in [1..t] do
        for j in [i + 1 .. t] do 
            k := input.pairorbit[i][j];
            if not rep.shape[k] in ["4B", "6A"] then 
                MAJORANA_RecordCoords(rep.involutions{[i,j]}, rep.shape[k], rep, algebras);
            fi;
         od;
    od;
    
    dim := Size(rep.setup.coords);
    
    rep.setup.nullspace := rec( heads := [1..dim]*0, vectors := SparseMatrix(0, 0, [], [], Rationals));
    
    rep.setup.pairorbit := NullMat(dim,dim);
    rep.setup.pairconj  := NullMat(dim,dim);
    
    for i in [1..t] do 
        for j in [1..t] do 
            rep.setup.pairorbit[i][j] := input.pairorbit[i][j];
            rep.setup.pairconj[i][j]  := input.pairconj[i][j];
        od;
    od;

    rep.setup.pairreps  := ShallowCopy(input.pairreps);
    rep.setup.pairconjelts := List(input.pairconjelts, x -> MAJORANA_FindPerm(x, rep, rep));

    gens := GeneratorsOfGroup(input.group);
    gens := List(gens, x -> MAJORANA_FindPerm(x, rep, rep));

    x := MAJORANA_Orbits(gens, t, rep.setup);

    rep.setup.conjelts := x.conjelts;
    rep.setup.orbitreps := x.orbitreps;

    MAJORANA_Orbitals(gens, t, rep.setup);
    
    s := Size(rep.setup.pairreps);
    
    rep.algebraproducts := List([1..s], x -> false);
    rep.innerproducts   := List([1..s], x -> false);
    rep.evecs           := NullMat(t,3);
    rep.nullspace := SparseMatrix(0, dim, [], [], Rationals);

    # Set up eigenvector matrix

    for j in [1..t] do
        if j in rep.setup.orbitreps then
            for k in [1..3] do
                rep.evecs[j][k] := SparseMatrix(0, dim, [], [], Rationals);
            od;
        else
            for k in [1..3] do
                rep.evecs[j][k] := false;
            od;
        fi;
    od;

    # Embed dihedral algebras
    
    for i in [1..t] do 
        for j in [i + 1..t] do 
            
            k := rep.setup.pairorbit[i][j];
            
            subrep := algebras.(rep.shape[k]);
            gens := GeneratorsOfGroup(subrep.group);
            
            emb := GroupHomomorphismByImages(subrep.group, rep.group, gens, rep.setup.coords{[i,j]});
            
            MAJORANA_Embed(rep, subrep, emb);
        od;
    od;

    for i in rep.setup.orbitreps do
        for j in [1..3] do 
            rep.evecs[i][j] := MAJORANA_BasisOfEvecs(rep.evecs[i][j]);
        od; 
    od;
    
    for x in rep.algebraproducts do 
        if x <> false then x!.ncols := dim; fi;
    od;
    
    return rep;
    
    end );
    
InstallGlobalFunction( MAJORANA_FindPerm, 
    
    function(g, rep, subrep)
    
    local   dim, j, list, pos, im, sign;
    
    dim := Size(subrep.setup.coords);
    list := [1..dim]*0;
        
    for j in [1..dim] do 
        if IsRowVector(subrep.setup.coords[j]) then 
        
            im := list{subrep.setup.coords[j]};
            
            sign := 1;
            
            if im[1] < 0 then sign := -sign; im[1] := -im[1]; fi;
            if im[2] < 0 then sign := -sign; im[2] := -im[2]; fi;
            
            if im[1] > im[2] then im := im{[2,1]}; fi;
            
            pos := Position(rep.setup.longcoords,im); 

            list[j] := rep.setup.poslist[pos];
        else 
            pos := Position(rep.setup.longcoords,OnPoints(subrep.setup.coords[j],g)); 
            list[j] := rep.setup.poslist[pos];
        fi;
    od;

    return list;
    
    end);
    
InstallGlobalFunction( MAJORANA_RecordCoords,

    function(involutions, shape, rep, algebras)
    
    local subrep, gens, emb, t, i, k, x, im, list, pos;

    if not shape in ["6A","4B"] and Product(involutions) in rep.setup.longcoords then 
        return;
    fi;

    subrep := algebras.(shape);
    
    gens := GeneratorsOfGroup(subrep.group);
    
    # Add extra basis vectors
    
    for i in [1.. Size(subrep.setup.coords)] do 
        
        list := Positions(subrep.setup.poslist, i);
        im := List(subrep.setup.longcoords{list}, 
                y -> MAJORANA_MappedWord(rep, subrep, y, gens, involutions));
        
        x := First(im, y -> y in rep.setup.longcoords);
            
        if x = fail then 
            Add(rep.setup.coords, im[1]);
            k := Size(rep.setup.coords);
            
            Append(rep.setup.longcoords, im);
            Append(rep.setup.poslist, List(im, y -> k));
            
            if shape = "5A" then 
                list := Positions(subrep.setup.poslist, -i);
                x := List(subrep.setup.longcoords{list}, 
                        y -> MAJORANA_MappedWord(rep, subrep, y, gens, involutions));
                
                Append(rep.setup.longcoords, x);
                Append(rep.setup.poslist, List(list, y -> -k));
            fi;
        else
            im := Filtered(im, y -> not y in rep.setup.longcoords);
            pos := rep.setup.poslist[Position(rep.setup.longcoords, x)];
                        
            Append(rep.setup.longcoords, im); 
            Append(rep.setup.poslist, List(im, y -> pos));
        fi;
    od;
    
    end );
    
InstallGlobalFunction( MAJORANA_RemoveDuplicateShapes, 

    function(input)
    
    local autgp, inner_autgp, outer_auts, perm, g, i, pos, im;
    
    autgp := AutomorphismGroup(input.group);
    inner_autgp := InnerAutomorphismsAutomorphismGroup(autgp);
    outer_auts := [];
    
    for g in RightTransversal(autgp, inner_autgp) do 
        if AsSet(OnTuples(input.involutions, g)) = AsSet(input.involutions) then 
            perm := [];
            
            for i in [1..Size(input.orbitals)] do 
            
                im := OnPairs(Representative(input.orbitals[i]), g);
                
                pos := PositionProperty(input.orbitals, x -> im in x);
                
                if pos = fail then pos := PositionProperty(input.orbitals, x -> Reversed(im) in x); fi;
                
                if pos = fail then Error(); fi;
                
                Add(perm, pos);
            od;
        
            Add(outer_auts, perm);
        fi;
    od;
    
    for i in [1..Size(input.shapes)] do 
        if IsBound(input.shapes[i]) then 
            for g in outer_auts do 
                pos := Position(input.shapes, input.shapes[i]{g});
                if pos <> fail and pos <> i then 
                    Unbind(input.shapes[pos]);
                fi;
            od;
        fi;
    od;
    
    return Compacted(input.shapes);
    
    end );
    
InstallGlobalFunction(MAJORANA_MappedWord,

    function(rep, subrep, w, gens, imgs)
    
    local im;
    
    if IsRowVector(w) then 
        im := List(w, i -> MappedWord(subrep.setup.coords[i], gens, imgs));
    
        return List(im, x -> Position(rep.setup.coords, x ));
    else
        return MappedWord(w, gens, imgs);
    fi;
    
    end );
    
InstallGlobalFunction(SP_Inverse,

    function(perm)
    
    local l, inv, i;
    
    if perm = [] then return []; fi;
    
    l := Length(perm);
    
    inv := [1..l];
    
    for i in [1..l] do 
        if perm[i] > 0 then 
            inv[perm[i]] := i;
        else
            inv[-perm[i]] := -i;
        fi;
    od;
    
    return inv;
    
    end);

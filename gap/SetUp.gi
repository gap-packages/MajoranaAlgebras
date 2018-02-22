
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
    
InstallGlobalFunction( MAJORANA_SetUpNoAxioms, 

    function(input, index)
    
    local   rep, t, T, i, j, k, pos, x, im, sign, gens, dim;
    
    rep         := rec( group       := input.group,
                        involutions := input.involutions,
                        shape       := input.shapes[index] );
                        
    T := rep.involutions;
    t := Size(T);
                        
    rep.setup   := rec( coords      := ShallowCopy(input.involutions),
                        longcoords  := ShallowCopy(input.involutions),
                        pairorbit   := input.pairorbit,
                        pairconj    := input.pairconj,
                        pairreps    := input.pairreps,
                        poslist     := [1..t]               );
    
    rep.setup.pairconjelts := List(input.pairconjelts, 
        x -> MAJORANA_FindVectorPermutation(x, rep.setup));
    
    # Add axes from 6A algebras    
    
    for i in [1..t] do 
        for j in [i + 1..t] do 
            k := rep.setup.pairorbit[i][j];
            if rep.shape[k] = ['6','A'] then 
                pos := [0, 0, 0, 0];
                pos[1] := Position(T, T[j]*T[i]*T[j]);
                pos[2] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                pos[3] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                pos[4] := Position(T, T[i]*T[j]*T[i]);
                
                # Add the 3A axes
                
                if not SortedList([i, pos[1]]) in rep.setup.longcoords then 
                    Add(rep.setup.coords, SortedList([i, pos[1]]));
                    
                    x := [  [i, pos[1]], [i, pos[3]], [j, pos[2]], [j, pos[4]],
                            [pos[1], pos[3]], [pos[2], pos[4]]];
                    
                    Append(rep.setup.longcoords, List(x, SortedList));
                    Append(rep.setup.poslist, [1, 1, 1, 1, 1, 1]*Size(rep.setup.coords)); 
                fi;
                
                # Add the 2A axes
                
                if not SortedList([i, pos[2]]) in rep.setup.longcoords then
                    
                    x := [[i, pos[2]], [j, pos[3]], [pos[1], pos[4]]];
                    
                    Append(rep.setup.longcoords, List(x, SortedList));
                    
                    x := Position(T, (T[i]*T[j])^3); 
                    
                    if x = fail then 
                        Add(rep.setup.coords, SortedList([i, pos[2]]));
                        Append(rep.setup.poslist, [1, 1, 1]*Size(rep.setup.coords));
                    else
                        Append(rep.setup.poslist, [1, 1, 1]*x);
                    fi;
                fi;
                
                x := Position(T, (T[i]*T[j])^3); 
                
                if x <> fail and not SortedList([x,i]) in rep.setup.longcoords then 
                    Append(rep.setup.longcoords, List([[x, i], [x,j]], SortedList)); 
                    Append(rep.setup.poslist, pos{[2,3]});
                    Append(rep.setup.longcoords, List(Cartesian( [x], pos), SortedList)); 
                    Append(rep.setup.poslist, [pos[4], i, j, pos[1]]);
                fi;
                
            fi;
        od;
    od;
    
    # Add axes from 4B algebra
    
    for i in [1..t] do 
        for j in [i + 1..t] do 
            k := rep.setup.pairorbit[i][j];
            if rep.shape[k] = ['4','B'] then
                pos := [0, 0];
                pos[1] := Position(T, T[j]*T[i]*T[j]);
                pos[2] := Position(T, T[i]*T[j]*T[i]);
                
                x := Position(T, (T[i]*T[j])^2);
                
                if not SortedList([i, pos[1]]) in rep.setup.longcoords then 
                    
                    Append(rep.setup.longcoords, List([[i, pos[1]], [j, pos[2]]], SortedList));

                    if x = fail then 
                        Add(rep.setup.coords, SortedList([i, pos[1]]));
                        Append(rep.setup.poslist, [1, 1]*Size(rep.setup.coords));
                    else
                        Append(rep.setup.poslist, [1, 1]*x);
                    fi;
                fi;
                
                # TODO finish this 
                
                if x <> fail and not SortedList([x,i]) in rep.setup.longcoords then 
                    
                fi;
                
            fi;
        od;
    od;
 
    for i in [1..t] do 
        for j in [i + 1 .. t] do 
            if not [i,j] in rep.setup.longcoords then
            
                k := rep.setup.pairorbit[i][j];
                
                if rep.shape[k] = ['2','A'] then 
                    Add(rep.setup.coords, [i,j]);
                    Add(rep.setup.longcoords, [i,j]);
                    Add(rep.setup.poslist, Size(rep.setup.coords));
                elif rep.shape[k] = ['3','A'] then
                    pos := Position(T, T[j]*T[i]*T[j]);
                    Add(rep.setup.coords, [i,j]);
                    Append(rep.setup.longcoords, [[i,j], [i,pos], [j,pos]]);
                    Append(rep.setup.poslist, [1, 1, 1]*Size(rep.setup.coords));
                elif rep.shape[k] = ['4', 'A'] then 
                    pos := [0,0];
                    pos[1] := Position(T, T[j]*T[i]*T[j]);
                    pos[2] := Position(T, T[i]*T[j]*T[i]);
                    Add(rep.setup.coords, [i,j]);
                    Add(rep.setup.longcoords, [[i, pos[1]], SortedList([i, pos[2]])]);
                elif rep.shape[k] = ['5','A'] then 
                    pos := [0, 0, 0];
                    pos[1] := Position(T, T[j]*T[i]*T[j]);
                    pos[2] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    Add(rep.setup.coords, [i,j]);
                    Append(rep.setup.longcoords, [[i,j], [i, pos[1]], [i, pos[2]], [i, pos[3]]]);
                    Append(rep.setup.poslist, [1, -1, -1, 1]*Size(rep.setup.coords));
                    Append(rep.setup.longcoords, [[j, pos[1]], [j, pos[2]], [j, pos[3]]]);
                    Append(rep.setup.poslist, [1, -1, -1]*Size(rep.setup.coords));
                    Append(rep.setup.longcoords, List([[pos[1], pos[2]], [pos[1], pos[3]], [pos[2], pos[3]]], SortedList));
                    Append(rep.setup.poslist, [1, -1, 1]*Size(rep.setup.coords));
                fi;
            fi;
        od;
    od;
    
    dim := Size(rep.setup.coords);
    
    # Extend pairconjelts permutations
    
    for i in [1..Size(rep.setup.pairconjelts)] do 
        for j in [t + 1 .. dim] do 
            x := rep.setup.coords[j];
            im := rep.setup.pairconjelts[i]{x};
            
            sign := 1;
                
            if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
            if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;
            
            if im[1] > im[2] then im := im{[2,1]}; fi;
            
            pos := Position(rep.setup.coords, im);
            
            Add(rep.setup.pairconjelts[i], pos);
        od;
    od;
    
    for i in [1..t] do 
        Append(rep.setup.pairorbit[i], [t + 1 .. dim]*0);
        Append(rep.setup.pairconj[i], [t + 1 .. dim]*0);
    od;
    
    Append(rep.setup.pairorbit, NullMat(dim - t, dim));
    Append(rep.setup.pairconj, NullMat(dim - t, dim));
    
    gens := GeneratorsOfGroup(rep.group);
    gens := List(gens, x -> Position(AsList(rep.group), x));
    gens := rep.setup.pairconjelts{gens};
    
    MAJORANA_Orbitals(gens, dim, rep.setup);
    
    return rep;
    
    end );

InstallGlobalFunction( MAJORANA_SetUp,
    
    function(input, index)
    
    local   rep,
            t,              # size of T
            x,              # temporary variables
            i,              # index variables
            j,
            k,
            g,
            y,
            res,
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
    
    # Create coordinates lists

    for i in [1..Size(input.orbitals)] do
        
        if rep.shape[i] = ['2','A'] then 
            for j in [1..Size(input.orbitals[i])] do
            
                x := input.orbitals[i][j][1]*input.orbitals[i][j][2];
                
                if x in input.involutions then break; fi;
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                
                    k := Size(rep.setup.coords);
                
                    Append(rep.setup.poslist, [k]);
                    Append(rep.setup.longcoords, [x]);
                fi;
            od;
        elif rep.shape[i]=['3','A'] then
            for j in [1..Size(input.orbitals[i])] do
            
                x := input.orbitals[i][j][1]*input.orbitals[i][j][2];
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                
                    k := Size(rep.setup.coords);
                
                    Append(rep.setup.poslist, [k,k]);
                    Append(rep.setup.longcoords, [x, x^2]);
                fi;
            od;
        elif rep.shape[i]=['4','A'] then
            for j in [1..Size(input.orbitals[i])] do
            
                x := input.orbitals[i][j][1]*input.orbitals[i][j][2];
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                    
                    k := Size(rep.setup.coords);
                    
                    Append(rep.setup.poslist, [k,k]);
                    Append(rep.setup.longcoords, [x, x^3]);
                fi;
                
            od;
        elif rep.shape[i] = ['4','B'] then 
            for j in [1..Size(input.orbitals[i])] do
            
                x := (input.orbitals[i][j][1]*input.orbitals[i][j][2])^2;
                
                if x in input.involutions then break; fi;
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                
                    k := Size(rep.setup.coords);
                
                    Append(rep.setup.poslist, [k]);
                    Append(rep.setup.longcoords, [x]);
                fi;
            od;
        elif rep.shape[i]=['5','A'] then
            for j in [1..Size(input.orbitals[i])] do
            
                x := input.orbitals[i][j][1]*input.orbitals[i][j][2];
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                    
                    k := Size(rep.setup.coords);
                    
                    Append(rep.setup.poslist, [k,-k,-k,k]);
                    Append(rep.setup.longcoords, [x, x^2, x^3, x^4]);
                fi;    
            od;
        elif rep.shape[i] = ['6','A'] then 
            for j in [1..Size(input.orbitals[i])] do
            
                x := (input.orbitals[i][j][1]*input.orbitals[i][j][2])^3;
                
                if x in input.involutions then break; fi;
                
                if not x in rep.setup.longcoords then 
                    Add(rep.setup.coords,x);
                
                    k := Size(rep.setup.coords);
                
                    Append(rep.setup.poslist, [k]);
                    Append(rep.setup.longcoords, [x]);
                fi;
            od;
        fi;
    od;
    
    dim := Size(rep.setup.coords);

    rep.setup.pairorbit := NullMat(dim,dim);
    rep.setup.pairconj  := NullMat(dim,dim);
    
    rep.setup.pairreps  := ShallowCopy(input.pairreps);
    rep.setup.pairconjelts := List(input.pairconjelts, 
        x -> MAJORANA_FindVectorPermutation(x, rep.setup));
    
    for i in [1..t] do 
        for j in [1..t] do 
            rep.setup.pairorbit[i][j] := input.pairorbit[i][j];
            rep.setup.pairconj[i][j]  := input.pairconj[i][j];
        od;
    od;    

    gens := GeneratorsOfGroup(input.group);
    gens := List(gens, x -> MAJORANA_FindVectorPermutation(x, rep.setup));

    x := MAJORANA_Orbits(gens, t, rep.setup);

    rep.setup.conjelts := x.conjelts;
    rep.setup.orbitreps := x.orbitreps;

    MAJORANA_Orbitals(gens, t, rep.setup);
    
                                ## STEP 3: PRODUCTS AND EVECS I ##
                                
    s := Size(rep.setup.pairreps);

    # Set up rep record

    rep.algebraproducts := List([1..s], x -> false);
    rep.innerproducts   := List([1..s], x -> false);
    rep.evecs           := NullMat(t,3);
    rep.nullspace       := SparseMatrix(0, dim, [], [], Rationals);               

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
    
    # Start filling in values and products!

    MAJORANA_DihedralProducts(input.involutions, rep);

    for i in rep.setup.orbitreps do
        for j in [1..3] do 
            if rep.evecs[i][j] <> [] then 
                rep.evecs[i][j] := MAJORANA_BasisOfEvecs(rep.evecs[i][j]);
            fi;
        od; 
    od;
    
    return rep;
    
    end );
    
InstallGlobalFunction(MAJORANA_MakeVector,

    function(pos,vals,dim)
    
    local   i,
            vec;
    
    vec := [1..dim]*0;
    
    for i in [1..Size(pos)] do 
        vec[pos[i]] := vals[i];
    od; 
    
    return SparseMatrix([vec], Rationals);
    
    end);
    
InstallGlobalFunction( MAJORANA_FindVectorPermutation, 
    
    function(g,setup)
    
    local   dim, j, list, pos;
    
    dim := Size(setup.coords);
    list := [1..dim]*0;
        
    for j in [1..dim] do 
        if IsRowVector(setup.coords[j]) then 
            pos := Position(setup.longcoords,OnTuples(setup.coords[j],g)); 
            list[j] := setup.poslist[pos];
        else 
            pos := Position(setup.longcoords,OnPoints(setup.coords[j],g)); 
            list[j] := setup.poslist[pos];
        fi;
    od;

    return list;
    
    end);
    
InstallGlobalFunction(MAJORANA_DihedralProducts,

    function(T, rep)

    local   i,
            j,      # loop over reps of T
            k,      # loop over T
            l,
            x,      # orbital of i,j
            y,      #
            pos,    # position 
            vals,   # values    
            dim,    # size of coordinates
            s,      # elt of T
            h,      # axis
            t,      # size of T
            sign;   # corrects sign of 5A axes
            
    dim := Size(rep.setup.coords);
    t   := Size(T);
    
    # (2,2) products and rep.evecs from IPSS10

    # Add rep.evecs from IPSS10

    for i in rep.setup.orbitreps do
        for j in [1..t] do
        
            if i <> j then 

                x := rep.setup.pairorbit[i][j];

                if rep.shape[x] = ['2','A'] then
                
                    pos := [i, j, 0];
                    pos[3] := Position(rep.setup.coords,T[i]*T[j]);
                    
                    vals := [-1/4, 1, 1];
                    
                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));
                    
                    vals := [0, 1, -1];

                    rep.evecs[i][2] := UnionOfRows(rep.evecs[i][2], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['2','B'] then
                
                    pos := [j];
                    vals := [1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['3','A'] then
                
                    pos := [i, j, 0, 0];

                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := rep.setup.poslist[Position(rep.setup.longcoords,T[i]*T[j])];

                    vals := [-10/27, 32/27, 32/27, 1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));
                    
                    vals := [-8/45, -32/45, -32/45, 1];

                    rep.evecs[i][2] := UnionOfRows(rep.evecs[i][2], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['3','C'] then
                
                    pos := [i, j, 0];

                    pos[3] := Position(T, T[i]*T[j]*T[i]);

                    vals := [-1/32, 1, 1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['4','A'] then
                    
                    pos := [i, j, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := rep.setup.poslist[Position(rep.setup.longcoords,T[i]*T[j])];

                    vals := [-1/2, 2, 2, 1, 1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [-1/3, -2/3, -2/3, -1/3, 1]; 

                    rep.evecs[i][2] := UnionOfRows(rep.evecs[i][2], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['4','B'] then
                    
                    pos := [i, j, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(rep.setup.coords, (T[i]*T[j])^2);
                    
                    vals := [-1/32, 1, 1, 1/8, -1/8];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif rep.shape[x] = ['5','A'] then
                
                    x := Position(rep.setup.longcoords,T[i]*T[j]);
                
                    pos := [i, j, 0, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                    pos[6] := rep.setup.poslist[x]; 

                    if pos[6] < 0 then
                        pos[6] := -pos[6];
                        sign := -1;
                    else
                        sign := 1;
                    fi;

                    vals := [3/512, -15/128, -15/128, -1/128, -1/128, sign*1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [-3/512, 1/128, 1/128, 15/128, 15/128, sign*1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos, vals, dim));
                    
                    vals := [0, 1/128, 1/128, -1/128, -1/128, sign*1];

                    rep.evecs[i][2] := UnionOfRows(rep.evecs[i][2], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 1, -1, 0, 0, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 0, 0, 1, -1, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos, vals, dim));

                elif rep.shape[x] = ['6','A'] then

                    pos := [i, j, 0, 0, 0, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                    pos[6] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                    pos[7] := Position(rep.setup.coords, (T[i]*T[j])^3);
                    pos[8] := rep.setup.poslist[Position(rep.setup.longcoords,(T[i]*T[j])^2)];

                    vals := [2/45, -256/45, -256/45, -32/45, -32/45, -32/45, 32/45, 1];

                    rep.evecs[i][1] := UnionOfRows(rep.evecs[i][1], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [-8/45, 0, 0, -32/45, -32/45, -32/45, 32/45, 1];

                    rep.evecs[i][2] := UnionOfRows(rep.evecs[i][2], MAJORANA_MakeVector(pos, vals, dim));
                    
                    vals := [0, 1, -1, 0, 0, 0, 0, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 0, 0, 1, -1, 0, 0, 0];

                    rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos, vals, dim));
                    
                    # put in products of 2A and 3A axes
                    
                    x := rep.setup.pairorbit[pos[7]][pos[8]];
                    
                    rep.algebraproducts[x] := SparseZeroMatrix(1, dim, Rationals);
                    
                    rep.innerproducts[x] := 0;
                fi;
            fi;
        od;
        
        # 1/32 rep.evecs from conjugation
        
        for j in [t+1 .. dim] do 
        
            h := rep.setup.coords[j];
            
            if not h^T[i] in [h,h^2,h^3,h^4] then 
            
                pos := [j, 0];
                pos[2] := rep.setup.poslist[Position(rep.setup.longcoords,h^T[i])];
                
                if pos[2] < 0 then
                    pos[2] := -pos[2];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [1,-sign*1];
                
                rep.evecs[i][3] := UnionOfRows(rep.evecs[i][3], MAJORANA_MakeVector(pos, vals, dim));
            fi;
        od;
    od;
        
    # Products from IPSS10

    for i in [1..Size(rep.setup.pairreps)] do

        j := rep.setup.pairreps[i][1];
        k := rep.setup.pairreps[i][2];
        
        if rep.setup.coords[j] in T and rep.setup.coords[k] in T then 
            if rep.shape[i] = ['1','A'] then
            
                pos := [j];
                vals := [1];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector( pos, vals, dim);
                
                rep.innerproducts[i] := 1;
            
            elif rep.shape[i] = ['2','A'] then

                pos := [j, k, 0];
                pos[3] := Position(rep.setup.coords,T[j]*T[k]);

                vals := [1/8, 1/8, -1/8];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector( pos, vals, dim);

                rep.innerproducts[i] := 1/8;

            elif rep.shape[i] = ['2','B'] then

                rep.algebraproducts[i] := SparseZeroMatrix(1, dim, Rationals);

                rep.innerproducts[i] := 0;

            elif rep.shape[i] = ['3','A'] then
            
                pos := [j, k, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := rep.setup.poslist[Position(rep.setup.longcoords,T[j]*T[k])];

                vals := [1/16, 1/16, 1/32, -135/2048];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 13/256;

            elif rep.shape[i] = ['3','C'] then
            
                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                
                vals := [1/64, 1/64, -1/64];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=1/64;

            elif rep.shape[i] = ['4','A'] then

                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := rep.setup.poslist[Position(rep.setup.longcoords,T[j]*T[k])];

                vals := [3/64, 3/64, 1/64, 1/64, -3/64]; 

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 1/32;

            elif rep.shape[i] = ['4','B'] then
            
                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(rep.setup.coords,(T[j]*T[k])^2);

                vals := [1/64, 1/64, -1/64, -1/64, 1/64];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=1/64;

            elif  rep.shape[i] = ['5','A'] then
            
                pos := [j, k, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := rep.setup.poslist[Position(rep.setup.longcoords,T[j]*T[k])];

                if pos[6] < 0 then
                    pos[6] := -pos[6];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [3/128, 3/128, -1/128, -1/128, -1/128, sign];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=3/128;

            elif rep.shape[i] = ['6','A'] then
            
                pos := [j, k, 0, 0, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := Position(T,T[k]*T[j]*T[k]*T[j]*T[k]);
                pos[7] := Position(rep.setup.coords,(T[j]*T[k])^3);
                pos[8] := rep.setup.poslist[Position(rep.setup.longcoords,(T[j]*T[k])^2)];
                
                vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=5/256;
            fi;

        # 2,2 products

        elif rep.setup.coords[j] in T and Order(rep.setup.coords[k]) = 2 then
            if rep.setup.coords[j]*rep.setup.coords[k] in T then 

                s := rep.setup.coords[j]; h := rep.setup.coords[k];

                # Inside a 2A algebra
                
                pos := [j, 0, k];
                pos[2] := Position(T,s*h);
                
                vals := [1/8, -1/8, 1/8];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 1/8;
            fi;

        # 2,3 products

        elif rep.setup.coords[j] in T and Order(rep.setup.coords[k]) = 3 then
            if rep.setup.coords[j]*rep.setup.coords[k] in T then 

                s := rep.setup.coords[j]; h := rep.setup.coords[k];

                # Inside a 3A algebra
                
                pos := [j, 0, 0, k];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h);
                
                vals := [2/9, -1/9, -1/9, 5/32];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=1/4;
            fi;

        # 2,4 products

        elif rep.setup.coords[j] in T and Order(rep.setup.coords[k]) = 4 then

            s := rep.setup.coords[j];
            h := rep.setup.coords[k];

            if s*h in T then

                # Inside a 4A algebra
                
                pos := [j, 0, 0, 0, k];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h*h);
                pos[4] := Position(T,s*h*h);
                
                vals := [5/16, -1/8, -1/8, -1/16, 3/16];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 3/8;
            fi;

        # (2,5) values

        elif rep.setup.coords[j] in T and Order(rep.setup.coords[k]) = 5 then

            s := rep.setup.coords[j];
            h := rep.setup.coords[k];

            if s*h in T then

                # Inside a 5A algebra
                
                pos := [0, 0, 0, 0, k];
                pos[1] := Position(T,s*h);
                pos[2] := Position(T,s*h^4);
                pos[3] := Position(T,s*h^2);
                pos[4] := Position(T,s*h^3);
                
                vals := [7/4096, 7/4096, -7/4096, -7/4096, 7/32];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i] := 0;
            fi;
            
        elif j = k then  
        
            h := rep.setup.coords[j];
            
            if Order(h) = 3 then    # (3,3) values
                
                rep.algebraproducts[i] := SparseMatrix( 1, dim, [[j]], [[1]], Rationals);
                rep.innerproducts[i] := 8/5;
                
            elif Order(h) = 4 then  # (4,4) values
            
                rep.algebraproducts[i] := SparseMatrix( 1, dim, [[j]], [[1]], Rationals);
                rep.innerproducts[i] := 2;
                
            elif Order(h) = 5 then  # (5,5) values
            
                l := 1;

                while l < t+1 do

                    if T[l]*h in T then

                        s:=T[l]; 
                        
                        pos := [l, 0, 0, 0, 0];
                        pos[2] := Position(T,s*h); 
                        pos[3] := Position(T,s*h^2); 
                        pos[4] := Position(T,s*h^3); 
                        pos[5] := Position(T,s*h^4);
                        
                        vals := [1,1,1,1,1]*(175/524288);

                        rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                        l:=t+1;
                    else
                        l:=l+1;
                    fi;
                od;

                rep.innerproducts[i] := 875/2^(19);
            fi;
        fi;
    od;
    
    end );
    
InstallGlobalFunction(SP_Product,    
    
    function( perm1, perm2)
    
    local l1, l2, prod, i;
    
    l1 := Length(perm1);
    l2 := Length(perm2);
    
    # Make perms the same length
    
    if l1 < l2 then
        Append(perm1, [l1 + 1 .. l2]);
    else
        Append(perm2, [l2 + 1 .. l1]);
    fi;
    
    prod := [];
    
    for i in perm1 do 
        if i > 0 then 
            Add(prod, perm2[i]);
        else
            Add(prod, -perm2[-i]);
        fi;
    od;
    
    return prod;
    
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
    

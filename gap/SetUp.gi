# TODO check with Sasha what happens when we have g,t \in T st 
# \ll a_t, a_g \rr is of type 2A but tg \notin T and implement

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
    
        if rep.shape[i]=['3','A'] then
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
        fi;
    od;
    
    dim := Size(rep.setup.coords);

    rep.setup.pairorbit := NullMat(dim,dim);
    rep.setup.pairconj  := NullMat(dim,dim);
    
    rep.setup.pairconjelts := ShallowCopy(input.pairconjelts);
    rep.setup.pairreps  := ShallowCopy(input.pairreps);
    
    for i in [1..t] do 
        for j in [1..t] do 
            rep.setup.pairorbit[i][j] := input.pairorbit[i][j];
            rep.setup.pairconj[i][j]  := input.pairconj[i][j];
        od;
    od;    

    x := MAJORANA_Orbits(input.group, t, rep.setup);

    rep.setup.conjelts := x.conjelts;
    rep.setup.orbitreps := x.orbitreps;

    MAJORANA_Orbitals(input.group, t, rep.setup);

    MAJORANA_FindAllPermutations(input.group, rep.setup);
    
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
    
    if g = () then 
        return ();
    else        
        for j in [1..dim] do 
            pos := Position(setup.longcoords,setup.coords[j]^g); 
            list[j] := setup.poslist[pos];
        od;
    
        return list;
    fi;
    
    end);
    
InstallGlobalFunction( MAJORANA_FindAllPermutations,

    function(G, setup)
    
    local   i, g;
            
    setup.pairconjelts := List(setup.pairconjelts, x -> [   MAJORANA_FindVectorPermutation(x, setup),
                                                            MAJORANA_FindVectorPermutation(Inverse(x), setup) ]); 
    
    for i in [1..Size(setup.conjelts)] do 
        g := setup.conjelts[i];
        setup.conjelts[i] := [g, MAJORANA_FindVectorPermutation(g, setup)];
    od;
    
    end );
    
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
                    pos[3] := Position(T,T[i]*T[j]);
                    
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
                    pos[5] := Position(T, (T[i]*T[j])^2);
                    
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
                    pos[7] := Position(T, (T[i]*T[j])^3);
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
        
        if Order(rep.setup.coords[j]) = 2 and Order(rep.setup.coords[k]) = 2 then 
            if rep.shape[i] = ['1','A'] then
            
                pos := [j];
                vals := [1];
                
                rep.algebraproducts[i] := MAJORANA_MakeVector( pos, vals, dim);
                
                rep.innerproducts[i] := 1;
            
            elif rep.shape[i] = ['2','A'] then

                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]);

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
                pos[5] := Position(T,(T[j]*T[k])^2);

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
                pos[7] := Position(T,(T[j]*T[k])^3);
                pos[8] := rep.setup.poslist[Position(rep.setup.longcoords,(T[j]*T[k])^2)];
                
                vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                rep.algebraproducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                rep.innerproducts[i]:=5/256;
            fi;

        # 2,3 products

        elif Order(rep.setup.coords[j]) = 2 and Order(rep.setup.coords[k]) = 3 then
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

        elif Order(rep.setup.coords[j]) = 2 and Order(rep.setup.coords[k]) = 4 then

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

        elif Order(rep.setup.coords[j]) = 2 and Order(rep.setup.coords[k]) = 5 then

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
    
#InstallGlobalFunction( MAJORANA_CheckSetup,

#    function(rep)
    
#    local   dim,
#            i,
#            j,
#            k,
#            g;
            
#    table := [[], [1], [1,2], [1,3], [1,2,3,4]];
    
#    dim := Size(rep.setup.coords);
    
#    for i in [1..dim] do 
#        for j in [1..dim] do 
            
#            k := rep.setup.pairorbit[i][j];
#            g := rep.setup.pairconj[i][j][1];
            
#            x := rep.setup.pairreps[k];
            
#            im1 := rep.setup.coords[k[1]^g];
#            im2 := k[2]^g;
            
#            o := List(rep.setup.coords{x}, y -> Order(y));
            
            
            
#        od;
#    od;
    
#    end );
            

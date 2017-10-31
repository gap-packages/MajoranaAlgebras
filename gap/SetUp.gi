# TODO change instances of gp elements to permutations
# TODO check with Sasha what happens when we have g,t \in T st 
# \ll a_t, a_g \rr is of type 2A but tg \notin T and implement

InstallGlobalFunction(ShapesOfMajoranaRepresentationAxiomM8,
    
    function(G,T)
    
    local   t,              # size of T
            i,              # indexes
            j,              #
            k,              #
            x,              # result of orbitals
            OrbitalsT,      # orbitals on T
            shape,          # one shape
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            Unknowns3X,     # orbitals which may be 3A or 3C
            Binaries,       # used to loop through options for shapes
            shapeslist,     # final list of shapes
            result;         #
    
    t := Size(T);

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;
    
    # Construct orbitals of G on T x T
    
    if IsAbelian(G) then 
        x := Cartesian(T,T);
    
        x := List(x,y -> [y]);

    else
        x := OrbitsDomain(G,Cartesian(T,T),OnPairs);
    fi;

    
    OrbitalsT := [];
    
    for i in [1..Size(x)] do
        Add(OrbitalsT, ShallowCopy(x[i]));
    od;
    
    i := 1;
    
    while i < Size(OrbitalsT) do 
    
        if not [OrbitalsT[i][1][2],OrbitalsT[i][1][1]] in OrbitalsT[i] then
        
            j := i + 1;
            
            while j < Size(OrbitalsT) + 1 do
            
                if  [OrbitalsT[i][1][2],OrbitalsT[i][1][1]]  in OrbitalsT[j] then
                
                    Append(OrbitalsT[i],OrbitalsT[j]);
                    Remove(OrbitalsT,j);
                        
                    j := Size(OrbitalsT) + 1;
                    
                else
                    
                    j := j + 1;
                fi;
            od;
        fi;
        
        i := i + 1;
        
    od;

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(OrbitalsT))[1];

    RepsSquares6A := [];

    for i in [1..Size(OrbitalsT)] do
    
        x := Representative(OrbitalsT[i]);
        
        if Order(x[1]*x[2]) = 1 then
            shape[i]:="1A";
        fi;
        if Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            shape[i]:="2A";
        fi;
        if Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            shape[i]:="2B";
        fi;
        if Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
        fi;
        if Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            shape[i]:="4A";
        fi;
        if Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            shape[i]:="4B";
        fi;
        if Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        fi;
        if Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(RepsSquares6A,(x[1]*x[2])^2);
        fi;
    od;

    # Check for inclusions of 3A in 6A

    for i in [1..Size(OrbitalsT)] do
        if shape[i][1] = '3' then
            j := 0;
            while j < Size(OrbitalsT[i]) do
                j := j + 1;
                if OrbitalsT[i][j][1]*OrbitalsT[i][j][2] in RepsSquares6A then
                    shape[i]:="3A";;
                    j:=Size(OrbitalsT[i])+1;;
                fi;
            od;
        fi;
    od;

    Unknowns3X:=[];

    for i in [1..Size(OrbitalsT)] do
        if shape[i] = ['3','X'] then
            Add(Unknowns3X,i);
        fi;
    od;
    
    Binaries:=AsList(FullRowSpace(GF(2),Size(Unknowns3X)));
    
    shapeslist := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do
        
        for j in [1..Size(Unknowns3X)] do
            k:=Unknowns3X[j];
            if Binaries[i][j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;            
        od;
        
        Add(shapeslist,ShallowCopy(shape));
        
    od;
    
    result := rec(  group := G,
                    involutions := T,
                    orbitals := OrbitalsT,
                    shapes := shapeslist);
    
    return result;

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
            result;         #
    
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
    
    result := rec(  group       := G,
                    involutions := T,
                    orbitals    := orbs.orbitals,
                    shapes      := shapeslist,
                    pairreps    := orbs.pairreps,
                    pairorbit   := orbs.pairorbit,
                    pairconj    := orbs.pairconj     );
    
    return result;

    end );

InstallGlobalFunction( MAJORANA_SetUp,
    
    function(input, index)
    
    local   G,              # group
            T,              # involutions
            Shape,          # chosen shape
            OrbitalsT,      # orbits of G on T
            t,              # size of T
            x,              # temporary variables
            i,              # index variables
            j,
            k,
            g,
            y,
            res,
            dim,            # size of coordinates
            s,              # number of orbits of G on T x T
            GramMatrix,
            AlgebraProducts,
            EigenVectors,
            SetUp;

    G           := input.group;
    T           := input.involutions;
    Shape       := input.shapes[index];
    OrbitalsT   := input.orbitals;  
    
    t := Size(T);  
    
    SetUp := rec(   coords      := ShallowCopy(T),
                    longcoords  := ShallowCopy(T),
                    poslist     := [1..t]     );

    # Create coordinates lists

    for i in [1..Size(OrbitalsT)] do
    
        if Shape[i]=['3','A'] then
            for j in [1..Size(OrbitalsT[i])] do
            
                x := OrbitalsT[i][j][1]*OrbitalsT[i][j][2];
                
                if not x in SetUp.longcoords then 
                    Add(SetUp.coords,x);
                
                    k := Size(SetUp.coords);
                
                    Append(SetUp.poslist, [k,k]);
                    Append(SetUp.longcoords, [x, x^2]);
                fi;
            od;
        elif Shape[i]=['4','A'] then
            for j in [1..Size(OrbitalsT[i])] do
            
                x := OrbitalsT[i][j][1]*OrbitalsT[i][j][2];
                
                if not x in SetUp.longcoords then 
                    Add(SetUp.coords,x);
                    
                    k := Size(SetUp.coords);
                    
                    Append(SetUp.poslist, [k,k]);
                    Append(SetUp.longcoords, [x, x^3]);
                fi;
                
            od;
        elif Shape[i]=['5','A'] then
            for j in [1..Size(OrbitalsT[i])] do
            
                x := OrbitalsT[i][j][1]*OrbitalsT[i][j][2];
                
                if not x in SetUp.longcoords then 
                    Add(SetUp.coords,x);
                    
                    k := Size(SetUp.coords);
                    
                    Append(SetUp.poslist, [k,-k,-k,k]);
                    Append(SetUp.longcoords, [x, x^2, x^3, x^4]);
                fi;    
            od;
        fi;
    od;
    
    dim := Size(SetUp.coords);

    SetUp.pairorbit := NullMat(dim,dim);
    SetUp.pairconj  := NullMat(dim,dim);
    SetUp.nullspace := false;
    
    SetUp.pairreps  := ShallowCopy(input.pairreps);
    
    for i in [1..t] do 
        for j in [1..t] do 
            SetUp.pairorbit[i][j] := input.pairorbit[i][j];
            SetUp.pairconj[i][j]  := input.pairconj[i][j];
        od;
    od;    

    MAJORANA_Orbits(G, t, SetUp);

    MAJORANA_Orbitals(G, t, SetUp);

    MAJORANA_FindAllPermutations(G, SetUp);
    
                                ## STEP 3: PRODUCTS AND EVECS I ##
                                
    s := Size(SetUp.pairreps);

    # Set up algebra product and gram matrices

    AlgebraProducts := List([1..s], x -> false);
    GramMatrix := List([1..s], x -> false);

    # Set up eigenvector matrix

    EigenVectors := NullMat(t,3);

    for j in [1..t] do
        if j in SetUp.orbitreps[1] then
            for k in [1..3] do
                EigenVectors[j][k] := [];
            od;
        else
            for k in [1..3] do
                EigenVectors[j][k] := false;
            od;
        fi;
    od;
    
    # Start filling in values and products!

    MAJORANA_DihedralProducts(T,OrbitalsT, Shape, GramMatrix, AlgebraProducts, EigenVectors, SetUp);
    
    # Check algebra obeys axiom M1 at this stage

    #x := MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,ProductList);;
    
    if false then 
        return MAJORANA_OutputError( "Algebra does not obey axiom M1"
                            , x
                            , [Shape, GramMatrix, AlgebraProducts, EigenVectors, SetUp]);
    fi;
    
    return [Shape, GramMatrix, AlgebraProducts, EigenVectors, SetUp];
    
    end );
    
InstallGlobalFunction(MAJORANA_MakeVector,

    function(pos,vals,dim)
    
    local   vec,
            i;
    
    vec := [1..dim]*0;
    
    for i in [1..Size(pos)] do 
        vec[pos[i]] := vals[i];
    od;
    
    return vec;
    
    end);
    
InstallGlobalFunction( MAJORANA_FindVectorPermutation, 
    
    function(g,SetUp)
    
    local   dim,        # size of coordinates
            i,          # loop over group elements
            j,          # loop over coordinates
            list,       # list to build permutation
            perm,       # the permutation
            signlist,   # corrects signs of 5A axes
            pos_1,      # position of conjugated element in longcoordinates
            pos_2;      # corresponding position in coordinates
    
    dim := Size(SetUp.coords);
    
    signlist := ListWithIdenticalEntries(dim,1);
    
    if g = () then 
        return [(),signlist];
    else
        list := [1..dim]*0;
        for j in [1..dim] do 
        
            pos_1 := Position(SetUp.longcoords,SetUp.coords[j]^g);
            pos_2 := SetUp.poslist[pos_1];
            
            if pos_2 > 0 then 
                list[j] := pos_2;
            else
                list[j] := -pos_2;
                signlist[-pos_2] := -1;
            fi;
        od;
        
        perm := PermList(list);
    
        return [perm,signlist];
    fi;
    
    end);
    
InstallGlobalFunction( MAJORANA_FindAllPermutations,

    function(G, SetUp)
    
    local   gp,
            perms,
            dim,
            i,
            j,
            g,
            pos;
            
    gp := AsSet(G);
    perms := List(gp, x -> MAJORANA_FindVectorPermutation(x, SetUp)); 
    
    dim := Size(SetUp.coords);
    
    for i in [1..dim] do
        for j in [i..dim] do 
    
            g := SetUp.pairconj[i][j];
            
            pos := [Position(gp, g)];
            pos[2] := Position(gp, Inverse(g));
            
            SetUp.pairconj[i][j] := perms{pos};
            SetUp.pairconj[j][i] := perms{pos};
        od;
    od;
    
    for i in [1..Size(SetUp.conjelts)] do 
        
        g := SetUp.conjelts[i];
        
        pos := Position(gp, g);
        
        SetUp.conjelts[i] := [g, perms[pos]];
    od;
    
    end );
    
InstallGlobalFunction(MAJORANA_DihedralProducts,

    function(T,OrbitalsT, Shape, GramMatrix, AlgebraProducts, EigenVectors, SetUp)

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
            
    dim := Size(SetUp.coords);
    t   := Size(T);
    
    # (2,2) products and eigenvectors from IPSS10

    # Add eigenvectors from IPSS10

    for i in SetUp.orbitreps[1] do
        for j in [1..t] do
        
            if i <> j then 

                x := SetUp.pairorbit[i][j];

                if Shape[x] = ['2','A'] then
                
                    pos := [i, j, 0];
                    pos[3] := Position(T,T[i]*T[j]);
                    
                    vals := [-1/4, 1, 1];
                    
                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));
                    
                    vals := [0, 1, -1];

                    Add(EigenVectors[i][2], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['2','B'] then
                
                    pos := [j];
                    vals := [1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['3','A'] then
                
                    pos := [i, j, 0, 0];

                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := SetUp.poslist[Position(SetUp.longcoords,T[i]*T[j])];

                    vals := [-10/27, 32/27, 32/27, 1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));
                    
                    vals := [-8/45, -32/45, -32/45, 1];

                    Add(EigenVectors[i][2], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['3','C'] then
                
                    pos := [i, j, 0];

                    pos[3] := Position(T, T[i]*T[j]*T[i]);

                    vals := [-1/32, 1, 1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['4','A'] then
                    
                    pos := [i, j, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := SetUp.poslist[Position(SetUp.longcoords,T[i]*T[j])];

                    vals := [-1/2, 2, 2, 1, 1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [-1/3, -2/3, -2/3, -1/3, 1]; 

                    Add(EigenVectors[i][2], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['4','B'] then
                    
                    pos := [i, j, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(T, (T[i]*T[j])^2);
                    
                    vals := [-1/32, 1, 1, 1/8, -1/8];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos,vals,dim));

                    vals := [0, 1, -1, 0, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos,vals,dim));

                elif Shape[x] = ['5','A'] then
                
                    x := Position(SetUp.longcoords,T[i]*T[j]);
                
                    pos := [i, j, 0, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                    pos[6] := SetUp.poslist[x]; 

                    if pos[6] < 0 then
                        pos[6] := -pos[6];
                        sign := -1;
                    else
                        sign := 1;
                    fi;

                    vals := [3/512, -15/128, -15/128, -1/128, -1/128, sign*1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [-3/512, 1/128, 1/128, 15/128, 15/128, sign*1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos, vals, dim));
                    
                    vals := [0, 1/128, 1/128, -1/128, -1/128, sign*1];

                    Add(EigenVectors[i][2], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 1, -1, 0, 0, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 0, 0, 1, -1, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos, vals, dim));

                elif Shape[x] = ['6','A'] then

                    pos := [i, j, 0, 0, 0, 0, 0, 0];
                    pos[3] := Position(T, T[i]*T[j]*T[i]);
                    pos[4] := Position(T, T[j]*T[i]*T[j]);
                    pos[5] := Position(T, T[i]*T[j]*T[i]*T[j]*T[i]);
                    pos[6] := Position(T, T[j]*T[i]*T[j]*T[i]*T[j]);
                    pos[7] := Position(T, (T[i]*T[j])^3);
                    pos[8] := SetUp.poslist[Position(SetUp.longcoords,(T[i]*T[j])^2)];

                    vals := [2/45, -256/45, -256/45, -32/45, -32/45, -32/45, 32/45, 1];

                    Add(EigenVectors[i][1], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [-8/45, 0, 0, -32/45, -32/45, -32/45, 32/45, 1];

                    Add(EigenVectors[i][2], MAJORANA_MakeVector(pos, vals, dim));
                    
                    vals := [0, 1, -1, 0, 0, 0, 0, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos, vals, dim));

                    vals := [0, 0, 0, 1, -1, 0, 0, 0];

                    Add(EigenVectors[i][3], MAJORANA_MakeVector(pos, vals, dim));
                    
                    # put in products of 2A and 3A axes
                    
                    x := SetUp.pairorbit[pos[7]][pos[8]];
                    
                    AlgebraProducts[x] := [1..dim]*0;
                    
                    GramMatrix[x] := 0;
                fi;
            fi;
        od;
        
        # 1/32 eigenvectors from conjugation
        
        for j in [t+1 .. dim] do 
        
            h := SetUp.coords[j];
            
            if not h^T[i] in [h,h^2,h^3,h^4] then 
            
                pos := [j, 0];
                pos[2] := SetUp.poslist[Position(SetUp.longcoords,h^T[i])];
                
                if pos[2] < 0 then
                    pos[2] := -pos[2];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [1,-sign*1];
                
                Add(EigenVectors[i][3], MAJORANA_MakeVector(pos, vals, dim));
            fi;
        od;
    od;
        
    # Products from IPSS10

    for i in [1..Size(SetUp.pairreps)] do

        j := SetUp.pairreps[i][1];
        k := SetUp.pairreps[i][2];
        
        if Order(SetUp.coords[j]) = 2 and Order(SetUp.coords[k]) = 2 then 
            if Shape[i] = ['1','A'] then
            
                pos := [j];
                vals := [1];
                
                AlgebraProducts[i] := MAJORANA_MakeVector( pos, vals, dim);
                
                GramMatrix[i] := 1;
            
            elif Shape[i] = ['2','A'] then

                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]);

                vals := [1/8, 1/8, -1/8];
                
                AlgebraProducts[i] := MAJORANA_MakeVector( pos, vals, dim);

                GramMatrix[i] := 1/8;

            elif Shape[i] = ['2','B'] then

                AlgebraProducts[i] := NullMat(1,dim)[1];

                GramMatrix[i] := 0;

            elif Shape[i] = ['3','A'] then
            
                pos := [j, k, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := SetUp.poslist[Position(SetUp.longcoords,T[j]*T[k])];

                vals := [1/16, 1/16, 1/32, -135/2048];
                
                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i] := 13/256;

            elif Shape[i] = ['3','C'] then
            
                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                
                vals := [1/64, 1/64, -1/64];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i]:=1/64;

            elif Shape[i] = ['4','A'] then

                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := SetUp.poslist[Position(SetUp.longcoords,T[j]*T[k])];

                vals := [3/64, 3/64, 1/64, 1/64, -3/64]; 

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i] := 1/32;

            elif Shape[i] = ['4','B'] then
            
                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,(T[j]*T[k])^2);

                vals := [1/64, 1/64, -1/64, -1/64, 1/64];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i]:=1/64;

            elif  Shape[i] = ['5','A'] then
            
                pos := [j, k, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := SetUp.poslist[Position(SetUp.longcoords,T[j]*T[k])];

                if pos[6] < 0 then
                    pos[6] := -pos[6];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [3/128, 3/128, -1/128, -1/128, -1/128, sign];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i]:=3/128;

            elif Shape[i] = ['6','A'] then
            
                pos := [j, k, 0, 0, 0, 0, 0, 0];
                pos[3] := Position(T,T[j]*T[k]*T[j]);
                pos[4] := Position(T,T[k]*T[j]*T[k]);
                pos[5] := Position(T,T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := Position(T,T[k]*T[j]*T[k]*T[j]*T[k]);
                pos[7] := Position(T,(T[j]*T[k])^3);
                pos[8] := SetUp.poslist[Position(SetUp.longcoords,(T[j]*T[k])^2)];
                
                vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i]:=5/256;
            fi;

        # 2,3 products

        elif Order(SetUp.coords[j]) = 2 and Order(SetUp.coords[k]) = 3 then
            if SetUp.coords[j]*SetUp.coords[k] in T then 

                s := SetUp.coords[j]; h := SetUp.coords[k];

                # Inside a 3A algebra
                
                pos := [j, 0, 0, k];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h);
                
                vals := [2/9, -1/9, -1/9, 5/32];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i]:=1/4;
            fi;

        # 2,4 products

        elif Order(SetUp.coords[j]) = 2 and Order(SetUp.coords[k]) = 4 then

            s := SetUp.coords[j];
            h := SetUp.coords[k];

            if s*h in T then

                # Inside a 4A algebra
                
                pos := [j, 0, 0, 0, k];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h*h);
                pos[4] := Position(T,s*h*h);
                
                vals := [5/16, -1/8, -1/8, -1/16, 3/16];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i] := 3/8;
            fi;

        # (2,5) values

        elif Order(SetUp.coords[j]) = 2 and Order(SetUp.coords[k]) = 5 then

            s := SetUp.coords[j];
            h := SetUp.coords[k];

            if s*h in T then

                # Inside a 5A algebra
                
                pos := [0, 0, 0, 0, k];
                pos[1] := Position(T,s*h);
                pos[2] := Position(T,s*h^4);
                pos[3] := Position(T,s*h^2);
                pos[4] := Position(T,s*h^3);
                
                vals := [7/4096, 7/4096, -7/4096, -7/4096, 7/32];

                AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[i] := 0;
            fi;
            
        elif j = k then  
        
            h := SetUp.coords[j];
            
            if Order(h) = 3 then    # (3,3) values
                
                AlgebraProducts[i] := NullMat(1,dim)[1];
                AlgebraProducts[i][j] := 1;

                GramMatrix[i] := 8/5;
                
            elif Order(h) = 4 then  # (4,4) values
            
                AlgebraProducts[i] := NullMat(1,dim)[1];
                AlgebraProducts[i][j] := 1;

                GramMatrix[i] := 2;
                
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

                        AlgebraProducts[i] := MAJORANA_MakeVector(pos, vals, dim);

                        l:=t+1;
                    else
                        l:=l+1;
                    fi;
                od;

                GramMatrix[i] := 875/2^(19);
            fi;
        fi;
    od;
    
    end );

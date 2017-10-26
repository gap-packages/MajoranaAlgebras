InstallGlobalFunction(ShapesOfMajoranaRepresentation,
    
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

InstallGlobalFunction( MAJORANA_SetUp,
    
    function(input, index)
    
    local   G,              # group
            T,              # involutions
            Shape,          # chosen shape
            OrbitalsT,      # orbits of G on T
            t,              # size of T
            ProductList,    # lists needed to calculate products
            x,              # temporary variables
            i,              # index variables
            j,
            k,
            dim,            # size of coordinates
            SizeOrbitals,   # number of orbits of G on T x T
            GramMatrix,
            AlgebraProducts,
            EigenVectors;

    G           := input.group;
    T           := input.involutions;
    Shape       := input.shapes[index];
    OrbitalsT   := input.orbitals;  
    
    t := Size(T);  
    
    # Set up ProductList
    
    ProductList := [1..13]*0;
    
    ProductList[1] := StructuralCopy(T);

    # Create coordinates list

    for j in [1..Size(OrbitalsT)] do
        if Shape[j]=['3','A'] then
            for k in [1..Size(OrbitalsT[j])] do
                x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                Add(ProductList[1],Set([x,x^2]));
            od;
        fi;
    od;
    
    for j in [1..Size(OrbitalsT)] do
        if Shape[j]=['4','A'] then
            for k in [1..Size(OrbitalsT[j])] do
                x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                Add(ProductList[1],Set([x,x^3]));
            od;
        fi;
    od;
    
    for j in [1..Size(OrbitalsT)] do 
        if Shape[j]=['5','A'] then
            for k in [1..Size(OrbitalsT[j])] do
                x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                Add(ProductList[1],Set([x,x^2,x^3,x^4]));
            od;
        fi;
    od;

    ProductList[1] := DuplicateFreeList(ProductList[1]);
    
    dim := Size(ProductList[1]);

    for j in [t+1..dim] do
        ProductList[1][j] := ProductList[1][j][1];
    od;
    
    ProductList[2]  := StructuralCopy(T);
    ProductList[3]  := NullMat(dim,dim);
    ProductList[4]  := NullMat(dim,dim);
    ProductList[5]  := [1..t];
    ProductList[6]  := false;
    ProductList[7]  := [];    
    ProductList[8]  := G;
    ProductList[10] := [[],[]];
    ProductList[12] := [];
    ProductList[13] := [1..t]*0;
    ProductList[14] := [];

    MAJORANA_LongCoordinates(t, ProductList);
    
    MAJORANA_SetupOrbits(T, ProductList);
    
    MAJORANA_SetupOrbitals(ProductList, OrbitalsT);

    MAJORANA_PairRepresentatives(ProductList);

    MAJORANA_PairOrbits(ProductList);

    MAJORANA_PairConjElements(ProductList);
    
                                ## STEP 3: PRODUCTS AND EVECS I ##
                                
    SizeOrbitals := Size(ProductList[9]);

    # Set up algebra product and gram matrices

    AlgebraProducts := NullMat(1,SizeOrbitals)[1];
    GramMatrix := NullMat(1,SizeOrbitals)[1];

    for j in [1..SizeOrbitals] do
        AlgebraProducts[j]:=false;
        GramMatrix[j]:=false;
    od;

    # Set up eigenvector matrix

    EigenVectors := NullMat(t,3);

    for j in [1..t] do
        if j in ProductList[10][1] then
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

    MAJORANA_DihedralProducts(T,Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList);
    
    # Check algebra obeys axiom M1 at this stage

    x := MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,ProductList);;
    
    if Size(x) > 0 then 
        return MAJORANA_OutputError( "Algebra does not obey axiom M1"
                            , x
                            , [Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList]);
    fi;
    
    return [Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList];
    
    end );

InstallGlobalFunction( MAJORANA_PairConjElements,

    function(ProductList)
    
    local   x,      # input elements
            y,      # representative elements
            z,      # index of elements which will also have g
            list,   # list of indices z
            i,      # first basis element
            j,      # second basis element
            k,      # orbital of elements
            l,
            table,
            pos_1,
            pos_2,
            perm,
            gp,
            perms,
            dim,    # size of coordinates
            g;      # conjugating element
            
    table := [[],[1],[1,2],[1,3],[1,2,3,4]];
    
    dim := Size(ProductList[1]);
    
    gp := AsList(ProductList[8]);
    perms := List(gp, x -> MAJORANA_FindVectorPermutation(x,ProductList));
    
    for i in [1..dim] do
        for j in [1..dim] do
            if ProductList[4][i][j] = 0 then 
            
                x := [ProductList[1][i],ProductList[1][j]];
                
                # create list of other indices which are going to have this elt
                
                list := [];
                
                for k in table[Order(x[1])] do 
                    for l in table[Order(x[2])] do
                        Add(list,[x[1]^k,x[2]^l]);
                    od;
                od;
                
                # find orbit
                
                k := ProductList[3][i][j];
            
                if k < 0 then 
                    k := -k;
                fi; 
            
                # find rep of orbit
            
                y := ProductList[9][k][1];
                
                l := 1;
                
                while l < Size(list) + 1 do 
            
                    g := RepresentativeAction(ProductList[8],y,list[l],OnPairs);
                
                    if g <> fail then
                    
                        perm := [perms[Position(gp,g)],perms[Position(gp,Inverse(g))]];
                    
                        for z in list do
                            
                            pos_1 := ProductList[5][Position(ProductList[2],z[1])];
                            pos_2 := ProductList[5][Position(ProductList[2],z[2])];
                            
                            if pos_1 < 0 then 
                                pos_1 := -pos_1;
                            fi;
                            
                            if pos_2 < 0 then 
                                pos_2 := -pos_2;
                            fi;
                            
                            ProductList[4][pos_1][pos_2] := perm;
                            ProductList[4][pos_2][pos_1] := perm;
                        od;
                        
                        l := Size(list) + 1;
                        
                    else 
                             
                        g := RepresentativeAction(ProductList[8],y,Reversed(list[l]),OnPairs);
                    
                        if g <> fail then 
                        
                            perm := [perms[Position(gp,g)],perms[Position(gp,Inverse(g))]];
                        
                            for z in list do
                            
                                pos_1 := ProductList[5][Position(ProductList[2],z[1])];
                                pos_2 := ProductList[5][Position(ProductList[2],z[2])];
                                
                                if pos_1 < 0 then 
                                    pos_1 := -pos_1;
                                fi;
                                
                                if pos_2 < 0 then 
                                    pos_2 := -pos_2;
                                fi;
                                
                                ProductList[4][pos_1][pos_2] := perm;
                                ProductList[4][pos_2][pos_1] := perm;
                            od;
                            
                            l := Size(list) + 1;
                            
                        else
                            l := l + 1;
                        fi;
                    fi;
                od;
            fi;
        od;
    od;
    
    end );


    
InstallGlobalFunction( MAJORANA_PairOrbits,

    function(ProductList)
    
    local   i,    # first basis element
            j,    # second basis element
            dim,  # size of coordinates  
            k;    # loop through orbitals
            
    dim := Size(ProductList[1]);
    
    for i in [1..dim] do
        for j in [1..dim] do
        
            if ProductList[3][i][j] = 0 then 
                
                k := 1;
                
                while k < Size(ProductList[9]) + 1 do
                
                    if [ProductList[1][i],ProductList[1][j]] in ProductList[9][k] then
                    
                        if [ProductList[1][i],ProductList[1][j]] in ProductList[14] then
                        
                            ProductList[3][i][j] := -k;
                            ProductList[3][j][i] := -k;
                        
                        else
                    
                            ProductList[3][i][j] := k;
                            ProductList[3][j][i] := k;
                        fi;
                        
                        k := Size(ProductList[9]) + 1;
                        
                    else
                        k := k + 1;
                    fi;
                od;
            fi;
        od;
    od;
    
        
    end );
    
InstallGlobalFunction(MAJORANA_PairRepresentatives,

    function(ProductList)
    
    local   i,          # loop through orbitals
            x,          # representative elements
            y,          # positions of representatives
            list,       # list of equivalent pairs of elts
            j,          # orders of first element
            k,          # orders of second element
            pos_1,      # position of first element
            pos_2,      # position of second element
            table,      # table of orders of basis elements
            z;          # an equivalent pair of elements       
            
    table := [[],[1],[1,2],[1,3],[1,2,3,4]];

    for i in [1..Size(ProductList[9])] do
            
        x := ProductList[9][i][1];
        y := [Position(ProductList[1],x[1]), Position(ProductList[1],x[2])];
        
        Add(ProductList[7], y);
        
        # create list of other indices which are going to have this elt
                
        list := [];
        
        for j in table[Order(x[1])] do 
            for k in table[Order(x[2])] do
                Add(list,[x[1]^j,x[2]^k]);
            od;
        od;
        
        # put conj elt and pair orbit for representatives 
        
        for z in list do
                            
            pos_1 := ProductList[5][Position(ProductList[2],z[1])];
            pos_2 := ProductList[5][Position(ProductList[2],z[2])];
            
            if pos_1 < 0 then 
                pos_1 := -pos_1;
            fi;
            
            if pos_2 < 0 then 
                pos_2 := -pos_2;
            fi;
                            
            
            ProductList[4][pos_1][pos_2] := [[(),[]],[(),[]]];
            ProductList[4][pos_2][pos_1] := [[(),[]],[(),[]]];
            
            ProductList[3][pos_1][pos_2] := i;
            ProductList[3][pos_2][pos_1] := i;
        od;
    od;
    
    end);
    
InstallGlobalFunction( MAJORANA_LongCoordinates,
    
    function(t,ProductList)
    
    local   i,      # loop over coordinates
            dim,    # size of coordinates
            x;      # coordinate;
            
    dim := Size(ProductList[1]);

    for i in [t+1..dim] do
        
        x := ProductList[1][i];
    
        if Order(x) = 3 then 
    
            Append(ProductList[5],[i,i]);
            Append(ProductList[2],[x,x^2]);
            
        elif Order(x) = 4 then 

            Append(ProductList[5],[i,i]);
            Append(ProductList[2],[x,x^3]);
            
        elif Order(x) = 5 then 
            Append(ProductList[5],[i,-i,-i,i]);
            Append(ProductList[2],[x,x^2,x^3,x^4]); 
        fi;
    od;
    
    ProductList[2] := Flat(ProductList[2]);
    
    end );
    
InstallGlobalFunction( MAJORANA_MakeVector,

    function( pos, vals, dim)
    
    local   vec,    # output vector
            i;      # loop over input
            
    vec := [1..dim]*0;;
    
    for i in [1..Size(pos)] do
        vec[pos[i]] := vec[pos[i]] + vals[i];
    od;
    
    return vec;
    
    end );
    
InstallGlobalFunction( MAJORANA_SetupOrbitals,

    function(ProductList, OrbitalsT)
    
    local   G,      # the group
            i,      # loop over orbitals
            x,      # representative of orbital
            y,      # orders of representatives
            table,  # table of orders of axes
            j,      # loop over orders of first axis
            k,      # loop over orders of second axis
            l;      # loop over orbitals
            
   G := ProductList[8];
            
    x := Cartesian(ProductList[1],ProductList[1]);

    ProductList[9] := List(Orbits(G,x,OnPairs), y -> ShallowCopy(y));

    # This is a bit of a patch, ask Markus tomorrow

    j:=1;

    while j < Size(ProductList[9]) + 1 do
        if Order(ProductList[9][j][1][1]) = 2 and Order(ProductList[9][j][1][2]) = 2 then
            Remove(ProductList[9],j);
        else
            j := j+1;
        fi;
    od;

    ProductList[9] := Concatenation(OrbitalsT,ProductList[9]);
    
    j := Size(OrbitalsT) + 1;
    
    while j < Size(ProductList[9]) + 1 do 
        if not [ProductList[9][j][1][2],ProductList[9][j][1][1]] in ProductList[9][j] then
            k := j + 1;
            
            while k < Size(ProductList[9]) +1 do
            
                if  [ProductList[9][j][1][2],ProductList[9][j][1][1]]  in ProductList[9][k] then
                
                    if Order(ProductList[9][j][1][1]) < Order(ProductList[9][j][1][2]) then 
                
                        Append(ProductList[9][j],ProductList[9][k]);
                        Remove(ProductList[9],k);
                    
                    else 
                    
                        Append(ProductList[9][k],ProductList[9][j]);
                        Remove(ProductList[9],j);
                        
                        j := j - 1;
                        
                    fi;
                    
                    k := Size(ProductList[9]) + 1;
                    
                else
                    
                    k := k + 1;
                fi;
            od;
                                
        fi;
        
        j := j + 1;
        
    od;
    
    i := Size(OrbitalsT) + 1;
            
    while i < Size(ProductList[9]) + 1 do
        
        x := ProductList[9][i][1];

        y := [Order(x[1]),Order(x[2])];
        
        table  := [[],[1],[1,2],[1,3],[1,2,3,4]];
        
        for j in table[y[1]] do
            for k in table[y[2]] do  
                      
                if not [x[1]^j, x[2]^k] in ProductList[9][i] then
                    
                    l := i + 1;
                    
                    while l < Size(ProductList[9]) + 1 do
                    
                        if [x[1]^j, x[2]^k] in ProductList[9][l] then 
                        
                            Append(ProductList[9][i],ProductList[9][l]);
                           
                            if y = [5,5] then
                                if [j,k] in [[1,2],[2,1],[1,3],[3,1],[2,4],[4,2],[3,4],[4,3]] then 
                                    Append(ProductList[14], ProductList[9][l]);
                                fi;
                            elif y[2] = 5 and k in [2,3] then 
                                Append(ProductList[14], ProductList[9][l]);
                            fi;
                           
                            Remove(ProductList[9],l);
                           
                            l := Size(ProductList[9]) + 1;
                           
                        else
                            l := l + 1;
                        fi;
                        
                    od;
                fi;
            od;
        od;
        
        i := i + 1;       
        
    od;
        
end );     
    
InstallGlobalFunction(MAJORANA_SetupOrbits,

    function(T,ProductList)
    
    local   G,      # the group
            g,      # conjugating element
            x,      # orbit    
            i,      # loop over T
            j,      # loop over orbits 
            k,      # representative
            t,      # size of T
            orbits, # orbits of G on C
            pos;    # position in longcoordinates
    
    # Construct orbits of G on T
    
    G := ProductList[8];
    t := Size(T);
    
    ProductList[11] := OrbitsDomain(G,T);
    
    for x in ProductList[11] do
        Add(ProductList[10][1], Position( T, Representative(x)));
    od;
    
    for i in [1..t] do
        
        j := 1;
        
        while j < Size(ProductList[11]) + 1 do 
            if T[i] in ProductList[11][j] then 
            
                ProductList[13][i] := j;
                
                k := ProductList[10][1][j];
                
                g := [RepresentativeAction(G,T[k],T[i]),];
                
                if ForAll(ProductList[12], x -> not g[1] in x) then
                    g[2] := MAJORANA_FindVectorPermutation(g[1],ProductList);
                    
                    Add(ProductList[12],g);
                fi;
                
                j := Size(ProductList[11]) + 1;;
                
            else
                j := j + 1;
            fi;
        od;
    od;
    
    orbits := Orbits(G, ProductList[1]{[t+1 .. Size(ProductList[1])]});
    
    for i in [1..Size(orbits)] do 
        pos := Position(ProductList[2], Representative(orbits[i]));
        
        x := ProductList[5][pos];
        
        if x < 0 then 
            x := -x;
        fi;
        
        Add(ProductList[10][2],x);
    od;
    
    end );
    
    InstallGlobalFunction(MAJORANA_ConnectedComponents,

    function(T)
    
    local   graph,      # graph with vertices elements of T
            t,          # size of T
            i,          # loop over T
            j;          # loop over T
            
    LoadPackage("automata");
    
    t := Size(T);
    graph := NullMat(t,0);
    
    for i in [1..t] do 
        for j in [1..t] do 
            if Order(T[i]*T[j]) <> 2 or not T[i]*T[j] in T then 
                Add(graph[i],j);
            fi;
        od;
    od;
    
    return AutoConnectedComponents(graph);
    
    end);
    
    InstallGlobalFunction(MAJORANA_DihedralProducts,

function(T,Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList)

    local   j,      # loop over reps of T
            k,      # loop over T
            x,      # orbital of i,j
            y,      #
            pos,    # position 
            vals,   # values    
            dim,    # size of coordinates
            s,      # elt of T
            h,      # axis
            t,      # size of T
            sign;   # corrects sign of 5A axes
            
    dim := Size(ProductList[1]);
    t   := Size(T);
    
    # (2,2) products and eigenvectors from IPSS10

    # Add eigenvectors from IPSS10

    for j in ProductList[10][1] do
        for k in [1..t] do

            x := ProductList[3][j][k];

            if Shape[x] = ['2','A'] then
            
                pos := [j, k, 0];
                pos[3] := Position(T,T[j]*T[k]);
                
                vals := [-1/4, 1, 1];
                
                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));
                
                vals := [0, 1, -1];

                Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['2','B'] then
            
                pos := [k];
                vals := [1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['3','A'] then
            
                pos := [j, k, 0, 0];

                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := ProductList[5][Position(ProductList[2],T[j]*T[k])];

                vals := [-10/27, 32/27, 32/27, 1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));
                
                vals := [-8/45, -32/45, -32/45, 1];

                Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

                vals := [0, 1, -1, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['3','C'] then
            
                pos := [j, k, 0];

                pos[3] := Position(T, T[j]*T[k]*T[j]);

                vals := [-1/32, 1, 1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                vals := [0, 1, -1];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['4','A'] then
                
                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := Position(T, T[k]*T[j]*T[k]);
                pos[5] := ProductList[5][Position(ProductList[2],T[j]*T[k])];

                vals := [-1/2, 2, 2, 1, 1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                vals := [-1/3, -2/3, -2/3, -1/3, 1]; 

                Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

                vals := [0, 1, -1, 0, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['4','B'] then
                
                pos := [j, k, 0, 0, 0];
                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := Position(T, T[k]*T[j]*T[k]);
                pos[5] := Position(T, (T[j]*T[k])^2);
                
                vals := [-1/32, 1, 1, 1/8, -1/8];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                vals := [0, 1, -1, 0, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

            elif Shape[x] = ['5','A'] then
            
                x := Position(ProductList[2],T[j]*T[k]);
            
                pos := [j, k, 0, 0, 0, 0];
                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := Position(T, T[k]*T[j]*T[k]);
                pos[5] := Position(T, T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := ProductList[5][x]; 

                if pos[6] < 0 then
                    pos[6] := -pos[6];
                    sign := -1;
                else
                    sign := 1;
                fi;

                vals := [3/512, -15/128, -15/128, -1/128, -1/128, sign*1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));

                vals := [-3/512, 1/128, 1/128, 15/128, 15/128, sign*1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));
                
                vals := [0, 1/128, 1/128, -1/128, -1/128, sign*1];

                Add(EigenVectors[j][2], MAJORANA_MakeVector(pos, vals, dim));

                vals := [0, 1, -1, 0, 0, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

                vals := [0, 0, 0, 1, -1, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

            elif Shape[x] = ['6','A'] then

                pos := [j, k, 0, 0, 0, 0, 0, 0];
                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := Position(T, T[k]*T[j]*T[k]);
                pos[5] := Position(T, T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := Position(T, T[k]*T[j]*T[k]*T[j]*T[k]);
                pos[7] := Position(T, (T[j]*T[k])^3);
                pos[8] := ProductList[5][Position(ProductList[2],(T[j]*T[k])^2)];

                vals := [2/45, -256/45, -256/45, -32/45, -32/45, -32/45, 32/45, 1];

                Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));

                vals := [-8/45, 0, 0, -32/45, -32/45, -32/45, 32/45, 1];

                Add(EigenVectors[j][2], MAJORANA_MakeVector(pos, vals, dim));
                
                vals := [0, 1, -1, 0, 0, 0, 0, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

                vals := [0, 0, 0, 1, -1, 0, 0, 0];

                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));
                
                # put in products of 2A and 3A axes
                
                x := ProductList[3][pos[7]][pos[8]];
                
                AlgebraProducts[x] := [1..dim]*0;
                
                GramMatrix[x] := 0;
            fi;
        od;
        
        # 1/32 eigenvectors from conjugation
        
        for k in [t+1..dim] do 
        
            h := ProductList[1][k];
            
            if not h^T[j] in [h,h^2] then 
            
                pos := [k, 0];
                pos[2] := ProductList[5][Position(ProductList[2],h^T[j])];
                
                if pos[2] < 0 then
                    pos[2] := -pos[2];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [1,-sign*1];
                
                Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));
            fi;
        od;
    od;
        
    # Products from IPSS10

    for j in [1..Size(ProductList[9])] do

        x := ProductList[7][j][1];
        y := ProductList[7][j][2];

        if Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 2 then

            if Shape[j] = ['1','A'] then

                AlgebraProducts[j] := NullMat(1,dim)[1];

                AlgebraProducts[j][x] := 1;

                GramMatrix[j] := 1;

            elif Shape[j] = ['2','A'] then

                pos := [x, y, 0];
                pos[3] := Position(T,T[x]*T[y]);

                vals := [1/8, 1/8, -1/8];
                
                AlgebraProducts[j] := MAJORANA_MakeVector( pos, vals, dim);

                GramMatrix[j] := 1/8;

            elif Shape[j] = ['2','B'] then

                AlgebraProducts[j] := NullMat(1,dim)[1];

                GramMatrix[j] := 0;

            elif Shape[j] = ['3','A'] then
            
                pos := [x, y, 0, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                pos[4] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                vals := [1/16, 1/16, 1/32, -135/2048];
                
                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j] := 13/256;

            elif Shape[j] = ['3','C'] then
            
                pos := [x, y, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                
                vals := [1/64, 1/64, -1/64];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j]:=1/64;

            elif Shape[j] = ['4','A'] then

                pos := [x, y, 0, 0, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                pos[4] := Position(T,T[y]*T[x]*T[y]);
                pos[5] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                vals := [3/64, 3/64, 1/64, 1/64, -3/64]; 

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j] := 1/32;

            elif Shape[j] = ['4','B'] then
            
                pos := [x, y, 0, 0, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                pos[4] := Position(T,T[y]*T[x]*T[y]);
                pos[5] := Position(T,(T[x]*T[y])^2);

                vals := [1/64, 1/64, -1/64, -1/64, 1/64];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j]:=1/64;

            elif  Shape[j] = ['5','A'] then
            
                pos := [x, y, 0, 0, 0, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                pos[4] := Position(T,T[y]*T[x]*T[y]);
                pos[5] := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                pos[6] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                if pos[6] < 0 then
                    pos[6] := -pos[6];
                    sign := -1;
                else
                    sign := 1;
                fi;
                
                vals := [3/128, 3/128, -1/128, -1/128, -1/128, sign];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j]:=3/128;

            elif Shape[j] = ['6','A'] then
            
                pos := [x, y, 0, 0, 0, 0, 0, 0];
                pos[3] := Position(T,T[x]*T[y]*T[x]);
                pos[4] := Position(T,T[y]*T[x]*T[y]);
                pos[5] := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                pos[6] := Position(T,T[y]*T[x]*T[y]*T[x]*T[y]);
                pos[7] := Position(T,(T[x]*T[y])^3);
                pos[8] := ProductList[5][Position(ProductList[2],(T[x]*T[y])^2)];
                
                vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j]:=5/256;
            fi;

        # 2,3 products

        elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 3 then
            if ProductList[1][x]*ProductList[1][y] in T then 

                s := ProductList[1][x]; h := ProductList[1][y];

                # Inside a 3A algebra
                
                pos := [x, 0, 0, y];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h);
                
                vals := [2/9, -1/9, -1/9, 5/32];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j]:=1/4;
            fi;

        # 2,4 products

        elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 4 then

            s := ProductList[1][x];
            h := ProductList[1][y];

            if s*h in T then

                # Inside a 4A algebra
                
                pos := [x, 0, 0, 0, y];
                pos[2] := Position(T,s*h);
                pos[3] := Position(T,s*h*h*h);
                pos[4] := Position(T,s*h*h);
                
                vals := [5/16, -1/8, -1/8, -1/16, 3/16];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j] := 3/8;
            fi;

        # (2,5) values

        elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 5 then

            s := ProductList[1][x];
            h := ProductList[1][y];

            if s*h in T then

                # Inside a 5A algebra
                
                pos := [0, 0, 0, 0, y];
                pos[1] := Position(T,s*h);
                pos[2] := Position(T,s*h^4);
                pos[3] := Position(T,s*h^2);
                pos[4] := Position(T,s*h^3);
                
                vals := [7/4096, 7/4096, -7/4096, -7/4096, 7/32];

                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                GramMatrix[j] := 0;
            fi;
            
        elif x = y then 
        
            h := ProductList[1][x];
            
            if Order(h) = 3 then    # (3,3) values
                
                AlgebraProducts[j] := NullMat(1,dim)[1];
                AlgebraProducts[j][x] := 1;

                GramMatrix[j] := 8/5;
                
            elif Order(h) = 4 then  # (4,4) values
            
                AlgebraProducts[j] := NullMat(1,dim)[1];
                AlgebraProducts[j][x] := 1;

                GramMatrix[j] := 2;
                
            elif Order(h) = 5 then  # (5,5) values
            
                k := 1;

                while k < t+1 do

                    if T[k]*h in T then

                        s:=T[k]; 
                        
                        pos := [k, 0, 0, 0, 0];
                        pos[2] := Position(T,s*h); 
                        pos[3] := Position(T,s*h^2); 
                        pos[4] := Position(T,s*h^3); 
                        pos[5] := Position(T,s*h^4);
                        
                        vals := [1,1,1,1,1]*(175/524288);

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        k:=t+1;
                    else
                        k:=k+1;
                    fi;
                od;

                GramMatrix[j] := 875/2^(19);
            fi;
        fi;
    od;
    
    end );

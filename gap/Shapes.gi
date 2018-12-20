
##
## Finds all possible shapes of a Majorana representation of the form (G,T,V)
## that obey axiom M8.
##

InstallGlobalFunction(ShapesOfMajoranaRepresentationAxiomM8,

    function(G,T)

    local   g, perm, pos_1, pos_2, hom,
            t,              # size of T
            i,              # indices
            j,
            k,
            x,              # result of orbitals
            shape,          # one shape
            unknowns,       # indices of 3X axes
            pos,            # positions
            binaries,       # used to loop through options for shapes
            input;          #

    t := Size(T);

    # If G in not a permutation group then convert it to permutations on T

    if not IsPermGroup(G) then
        hom := ActionHomomorphism(G,T);
        T := Image(hom, T);
        G := Group(T);
    fi;

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;

    input := rec();

    input.setup := rec( pairrepsmap := HashMap( t*t ),pairreps := [], coords := [1..t] );
    input.involutions := T;
    input.group       := G;
    input.generators := [];

    for g in GeneratorsOfGroup(G) do
        perm := [];
        for i in [1..t] do
            Add(perm, Position(T, T[i]^g));
        od;
        Add(input.generators, perm);
    od;

    # Construct orbitals of  on T x T

    MAJORANA_FindOrbitals(input, input.generators, [1..t]);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := [1 .. Size(input.setup.pairreps)]*0;
    unknowns := [];;

    for i in [1..Size(input.setup.pairreps)] do

        x := T{input.setup.pairreps[i]};

        if Order(x[1]*x[2]) = 1 then
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            shape[i]:="2A";
        elif Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            shape[i]:="2B";
        elif Order(x[1]*x[2]) = 3 and shape[i] = 0 then
            shape[i]:="3X";
        elif Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            shape[i]:="4A";
        elif Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            shape[i]:="4B";
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2]) = 6 then
            shape[i]:="6A";

            MAJORANA_RecordSubalgebras( i, shape, input);

        elif Order(x[1]*x[2]) > 6 then
            Error("This is not a 6-transposition group");
        fi;
    od;

    unknowns := Positions(shape, "3X");

    binaries := AsList(FullRowSpace(GF(2),Size(unknowns)));

    input.shapes := [];

    # Add new values in the shape

    for i in [1..Size(binaries)] do

        for j in [1..Size(unknowns)] do
            k := unknowns[j];
            if binaries[i, j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;
        od;

        Add(input.shapes,ShallowCopy(shape));
    od;

    return input;

    end );

##
## Finds all possible shapes of a Majorana representation of the form (G,T,V)
##

InstallGlobalFunction(ShapesOfMajoranaRepresentation,

    function(G,T)

    local   hom,
            t,              # size of T
            i,              # indices
            j,
            k,
            g,
            perm,
            x,              # result of orbitals
            shape,          # one shape
            gph,            # digraph of 2X, 4X inclusions
            cc,             # connected components of gph
            binaries,       # used to loop through options for shapes
            input;          #

    t := Size(T);

    # If G in not a permutation group then convert it to permutations on T

    if not IsPermGroup(G) then
        hom := ActionHomomorphism(G,T);
        T := Image(hom, T);
        G := Group(T);
    fi;

    # Construct orbitals of  on T x T

    input := rec();

    input.setup := rec( pairrepsmap := HashMap( t*t ), pairreps := [], coords := [1..t] );
    input.involutions := T;
    input.group       := G;
    input.generators  := [];

    for g in GeneratorsOfGroup(G) do
        perm := [];
        for i in [1..t] do
            Add(perm, Position(T, T[i]^g));
        od;
        Add(input.generators, perm);
    od;

    MAJORANA_FindOrbitals(input, input.generators, [1..t]);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(input.setup.pairreps))[1];

    for i in [1..Size(input.setup.pairreps)] do

        x := T{input.setup.pairreps[i]};

        if Order(x[1]*x[2]) = 1 then
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 and shape[i] = 0 then
            shape[i]:="2X";
        elif Order(x[1]*x[2]) = 3 and shape[i] = 0 then
            shape[i]:="3X";
        elif Order(x[1]*x[2]) = 4 then
            shape[i]:="4X";
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2]) = 6 then
            shape[i]:="6A";
            MAJORANA_RecordSubalgebras( i, shape, input);
        elif Order(x[1]*x[2]) > 6 then
            Error("This is not a 6-transposition group");
        fi;
    od;

    # Check for inclusions of 2X in 4X

    gph := List( [1 .. Size(input.setup.pairreps)], x -> [] );

    for i in Positions(shape, "4X") do
        gph[i] := MAJORANA_RecordSubalgebras(i, shape, input);
    od;

    cc := AutoConnectedComponents(gph);
    cc := Filtered( cc, comp -> ForAny(comp, i -> shape[i][2] = 'X' )  );

    for i in [1 .. Size(cc)] do
        if ForAny( cc[i], j -> shape[j] = "2A") then
            for j in cc[i] do
                if shape[j] = "4X" then
                    shape[j] := "4B";
                elif shape[j] = "2X" then
                    shape[j] := "2A";
                fi;
            od;
        fi;
    od;

    cc := Filtered( cc, comp -> ForAny(comp, i -> shape[i][2] = 'X' )  );

    binaries := AsList( FullRowSpace(GF(2), Size(cc)) );

    input.shapes := [];

    # Add new values in the shape

    for i in [1 .. Size(binaries)] do
        for j in [1 .. Size(cc)] do
            if binaries[i, j] = 1*Z(2) then
                for k in cc[j] do
                    if shape[k][1] in [ '2', '3' ] then
                        shape[k][2] := 'A';
                    elif shape[k][1] = '4' then
                        shape[k][2] := 'B';
                    fi;
                od;
            else
                for k in cc[j] do
                    if shape[k][1] = '2' then
                        shape[k][2] := 'B';
                    elif shape[k][1] = '3' then
                        shape[k][2] := 'C';
                    elif shape[k][1] = '4' then
                        shape[k][2] := 'A';
                    fi;
                od;
            fi;
        od;

        Add(input.shapes,StructuralCopy(shape));
    od;

    return input;

    end );

##
## The index <i> should point to a position in <shape> that is equal to "4X"
#! or "6A". If it is equal to "6A" then it sets the positions in shape
## of the 2A and 3A algebras. If it is equal to "4X" then it returns a list
## giving the positions in shape of the two 2A or 2B algebras that it contains.
##

InstallGlobalFunction( MAJORANA_RecordSubalgebras,

    function( i, shape, input )

        local output, x, inv, pos, k;

        output := [];

        # Do this for both orderings of the pair representative in order to find all
        # subalgebras
        for x in [input.setup.pairreps[i], Reversed(input.setup.pairreps[i])] do

            inv := input.involutions{x};

            if shape[i] = "6A" then

                # Record the position of the 3A subalgebra
                pos := Position( input.involutions, inv[1]^inv[2] );
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                shape[k] := "3A";

                # Record the position of the 2A subalgebra
                pos := Position( input.involutions, inv[2]^Product(inv) );
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                shape[k] := "2A";

            elif shape[i][1] = '4' then

                # Add the position of the 2X subalgebra to the list <output>
                pos := Position( input.involutions, inv[1]^inv[2] );
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                Add( output, k );

            fi;
        od;

        return DuplicateFreeList(output);

    end );

##
## Optional function for use by the user after calling <ShapesOfMajoranaRepresentation>
##

InstallGlobalFunction( MAJORANA_RemoveDuplicateShapes,

    function(input)

    local autgp, inner_autgp, outer_auts, perm, g, i, im, pos, rep, k;

    autgp := AutomorphismGroup(input.group);
    inner_autgp := InnerAutomorphismsAutomorphismGroup(autgp);
    outer_auts := [];

    for g in RightTransversal(autgp, inner_autgp) do
        if AsSet(OnTuples(input.involutions, g)) = AsSet(input.involutions) then
            perm := [];

            for rep in input.setup.pairreps do

                im := OnPairs( input.involutions{rep}, g );
                im := List(im, x -> Position(input.involutions, x));

                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, im);
                Add(perm, input.setup.pairrepsmap[k]);
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

    input.shapes := Compacted(input.shapes);

    end );

InstallGlobalFunction( MAJORANA_IsSixTranspositionGroup,

    function(G,T)

    local t, s, g;

    # Test if T is a set of involutions
    if ForAny(T, t -> Order(t) <> 2) then
        return false;
    fi;

    # Test if T generates G
    if not Group(T) = G then
        return false;
    fi;

    # Test if T is closed under conjugation by G
    for g in GeneratorsOfGroup(G) do
        for t in T do
            if not t^g in T then
                return false;
            fi;
        od;
    od;

    # Test if the product of any two elements of T has order at most 6

    for t in T do
        for s in T do
            if Order(t*s) > 6 then
                return false;
            fi;
        od;
    od;

    return true;

end );

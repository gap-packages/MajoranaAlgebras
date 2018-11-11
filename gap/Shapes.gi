
##
## Finds all possible shapes of a Majorana representation of the form (G,T,V)
## that obey axiom M8.
##

InstallGlobalFunction(ShapesOfMajoranaRepresentationAxiomM8,

    function(G,T)

    local   g, perm, pos_1, pos_2,
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

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;

    input := rec();

    input.pairorbit := NullMat(t,t);
    input.pairconj  := NullMat(t,t);
    input.pairreps  := [];
    input.pairconjelts := [ [1..t] ];
    input.coords := [1..t];
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

    MAJORANA_Orbitals(input.generators, 0, input);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := [1 .. Size(input.pairreps)]*0;
    unknowns := [];;

    for i in [1..Size(input.pairreps)] do

        x := T{input.pairreps[i]};

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

    local   t,              # size of T
            i,              # indices
            j,
            k,
            perm,
            x,              # result of orbitals
            shape,          # one shape
            gph,            # digraph of 2X, 4X inclusions
            cc,             # connected components of gph
            binaries,       # used to loop through options for shapes
            input;          #

    t := Size(T);

    # Construct orbitals of  on T x T

    input := rec();

    input.pairorbit := NullMat(t,t);
    input.pairconj  := NullMat(t,t);
    input.pairreps  := [];
    input.pairconjelts := [ [1..t] ];
    input.coords := [1..t];
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

    MAJORANA_Orbitals(input.generators, 0, input);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(input.pairreps))[1];

    for i in [1..Size(input.pairreps)] do

        x := T{input.pairreps[i]};

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
        fi;
    od;

    # Check for inclusions of 2X in 4X

    gph := List( [1 .. Size(input.pairreps)], x -> [] );

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
## Optional function for use by the user after calling <ShapesOfMajoranaRepresentation>
##

InstallGlobalFunction( MAJORANA_RemoveDuplicateShapes,

    function(input)

    local autgp, inner_autgp, outer_auts, perm, g, i, im, pos, rep;

    autgp := AutomorphismGroup(input.group);
    inner_autgp := InnerAutomorphismsAutomorphismGroup(autgp);
    outer_auts := [];

    for g in RightTransversal(autgp, inner_autgp) do
        if AsSet(OnTuples(input.involutions, g)) = AsSet(input.involutions) then
            perm := [];

            for rep in input.pairreps do

                im := OnPairs( input.involutions{rep}, g );
                im := List(im, x -> Position(input.involutions, x));

                Add(perm, input.pairorbit[im[1], im[2]]);
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

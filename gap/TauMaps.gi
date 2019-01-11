
InstallGlobalFunction(TauMapsMajoranaRepresentation,

    function(arg)

    local input, index, options;

    input := arg[1];
    index := arg[2];

    if IsBound(arg[3]) then
        options := arg[3];
    else
        options := rec();
    fi;

    options.axioms := "NoAxioms";
    options.taumaps := true;

    return MajoranaRepresentation( input, index, options );

end );

InstallGlobalFunction(TauMapsShapesOfMajoranaRepresentation,

    function(taumaps)

    local t, input, shape, i, j, k, x, D, orbit, n, gph, cc, binaries;

    t := Size(taumaps);

    input := rec();

    input.setup := rec( pairrepsmap := HashMap( t*t ), pairreps := [], coords := [1..t] );
    input.involutions := taumaps;
    input.group       := Group(taumaps);
    input.generators  := GeneratorsOfGroup(input.group);
    input.generators  := List(input.generators, x -> ListPerm(x, t));

    MAJORANA_FindOrbitals(input, input.generators, [1..t]);

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(input.setup.pairreps))[1];

    for i in [1..Size(input.setup.pairreps)] do

        x := input.setup.pairreps[i];

        if x[1] = x[2] then
            shape[i] := "1A";
        else
            D := Group(taumaps{x});

            orbit := Orbit(D, x[1]);
            n := Length(orbit);

            if x[2] in orbit then
                if n = 3 and shape[i] = 0 then
                    shape[i] := "3X";
                elif n = 5 then
                    shape[i] := "5A";
                else
                    Error("This tau map is not admissible");
                fi;
            else
                if n = 1 and shape[i] = 0 then
                    shape[i] := "2X";
                elif n = 2 then
                    shape[i] := "4X";
                elif n = 3 then
                    shape[i] := "6A";
                    MAJORANA_TauMapsRecordSubalgebras(i, shape, input);
                else
                    Error("This tau map is not admissible");
                fi;
            fi;
        fi;
    od;

    # Check for inclusions of 2X in 4X

    gph := List( [1 .. Size(input.setup.pairreps)], x -> [] );

    for i in Positions(shape, "4X") do
        gph[i] := MAJORANA_TauMapsRecordSubalgebras(i, shape, input);
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

InstallGlobalFunction( MAJORANA_TauMapsRecordSubalgebras,

    function( i, shape, input )

        local output, x, inv, pos, k;

        output := [];

        # Do this for both orderings of the pair representative in order to find all
        # subalgebras
        for x in [input.setup.pairreps[i], Reversed(input.setup.pairreps[i])] do

            inv := input.involutions{x};

            if shape[i] = "6A" then

                # Record the position of the 3A subalgebra
                pos := x[1]^inv[2];
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                shape[k] := "3A";

                # Record the position of the 2A subalgebra
                pos := x[2]^Product(inv);
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                shape[k] := "2A";

            elif shape[i][1] = '4' then

                # Add the position of the 2X subalgebra to the list <output>
                pos := x[1]^inv[2];
                k := UnorderedOrbitalRepresentative(input.setup.orbitalstruct, [x[1], pos]);
                k := input.setup.pairrepsmap[k];
                Add( output, k );

            fi;
        od;

        return DuplicateFreeList(output);

    end );

##
##
##

InstallGlobalFunction( MAJORANA_TauMapsMappedWord,

    function(rep, subrep, w, gens, inv)

    local im, i, imgs;

    imgs := rep.involutions{inv};

    if IsRowVector(w) then
        im := List(w, i -> MAJORANA_TauMapsMappedWord( rep, subrep, subrep.setup.coords[i], gens, inv) );

        if im[1] = im[2] then Error(); fi;

        return SortedList( im );
    else
        i := Position(subrep.setup.coords, w);

        if IsBound( subrep.setup.orbits[1][i] ) then
            im := MappedWord( subrep.setup.orbits[1, i], gens, imgs );
            return inv[1]^im;
        else
            im := MappedWord( subrep.setup.orbits[2, i], gens, imgs );
            return inv[2]^im;
        fi;
    fi;

    end );

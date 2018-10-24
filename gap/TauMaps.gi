InstallGlobalFunction( TauShapesOfMajoranaRepresentation,

    function(tau)

    local input, d, shape, i, x, D, o1, o2, n, gph, comps, binaries, unknowns3X, j, k;

    d := Size(tau);

    input := rec();

    input.involutions := [ 1 .. d ];
    input.tau       := tau;
    input.pairorbit := NullMat(d,d);
    input.pairconj  := NullMat(d,d);
    input.pairreps  := [];
    input.orbitals  := [];
    input.pairconjelts := [ [1..d] ];
    input.coords := [1..d];
    input.group := Group(tau);
    input.generators := List(GeneratorsOfGroup(input.group), x -> ListPerm(x,d));

    MAJORANA_Orbitals(input.generators, 0, input);

    shape := [];

    for i in [1 .. Size(input.pairreps)] do

        x := input.pairreps[i];

        D := Group(tau{x});

        o1 := Orbit(D, x[1]);
        o2 := Orbit(D, x[2]);

        if x[1] in o2 then
            n := Size(o1);

            if n = 1 then
                Add( shape,  "1A" );
            elif n = 3 then
                Add( shape,  "3X" );
            elif n = 5 then
                Add( shape,  "5A" );
            else
                Error("This is not a valid tau map");
            fi;
        else
            n := Size(o1) + Size(o2);

            if n = 2 then
                Add( shape,  "2X" );
            elif n = 4 then
                Add( shape,  "4X" );
            elif n = 6 then
                Add( shape,  "6A" );
            else
                Error("This is not a valid tau map");
            fi;
        fi;
    od;

    # Inclusions of 2A and 3A in 6A algebras

    for i in [1..Size(input.pairreps)] do
        if shape[i, 1] = '6' then

            for x in [input.pairreps[i], Reversed(input.pairreps[i])] do

                k := input.pairorbit[x[1], x[1]^tau[x[2]]];

                shape[k] := "3A";

                k := input.pairorbit[x[1], x[2]^(tau[x[1]]*tau[x[2]])];

                shape[k] := "2A";
            od;
        fi;
    od;

    # Inclusions of 2X into 4X

    gph := NullMat(Size(input.pairreps), 0);

    for i in [1..Size(input.pairreps)] do
        if shape[i, 1] = '4' then

            for x in [input.pairreps[i], Reversed(input.pairreps[i])] do
                Add(gph[i], input.pairorbit[x[1], x[1]^tau[x[2]]]);
            od;
        fi;
    od;

    gph := List(gph, DuplicateFreeList);

    comps := AutoConnectedComponents(gph);

    comps := Filtered(comps, x -> ForAny(shape{x}, y -> y in ["2X", "4X"] ) );

    # Put in any known (4B, 2A) pairs

    for x in comps do
        if ForAny(shape{x}, y -> y[2] = 'A') then
            for i in x do
                if shape[i, 1] = '2' then
                    shape[i] := "2A";
                else
                    shape[i] := "4B";
                fi;
            od;
        fi;
    od;

    unknowns3X := Filtered([1..Size(shape)], i -> shape[i] = "3X");

    binaries := AsList(FullRowSpace(GF(2), Size(comps) + Size(unknowns3X) ));

    input.shapes := [];

    for i in [1..Size(binaries)] do

        for j in [1 .. Size(comps)] do

            if binaries[i, j] = 1*Z(2) then
                for k in comps[j] do
                    if shape[k, 1] = '2' then
                        shape[k] := "2A";
                    else
                        shape[k] := "4B";
                    fi;
                od;
            else
                for k in comps[j] do
                    if shape[k, 1] = '2' then
                        shape[k] := "2B";
                    else
                        shape[k] := "4A";
                    fi;
                od;
            fi;
        od;

        for j in [1 .. Size(unknowns3X)] do

            if binaries[i, Size(comps) + j] = 1*Z(2) then
                shape[unknowns3X[j]] := "3A";
            else
                shape[unknowns3X[j]] := "3C";
            fi;
        od;

        Add(input.shapes, StructuralCopy(shape));

    od;

    return input;

    end );

InstallGlobalFunction( TauMapMajoranaRepresentation,

    function(input, index)

    local rep, unknowns, main;

    rep := MAJORANA_SetUp(input, index, "NoAxioms", true);

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        main := MAJORANA_MainLoop(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );
        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.innerproducts, false)), " unknown inner products ") );

        if not false in rep.algebraproducts then
            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            Info( InfoMajorana, 10, "Fail" );
            rep.mat := main.mat; rep.vec := main.vec; rep.unknowns := main.unknowns;
            return rep;
        fi;
    od;

    end );

InstallGlobalFunction( MAJORANA_TauMappedWord,

    function(rep, subrep, w, gens, inv)

    local im, i, imgs;

    imgs := rep.tau{inv};

    if IsRowVector(w) then
        im := List(w, i -> MAJORANA_TauMappedWord( rep, subrep, subrep.setup.coords[i], gens, inv) );

        # if im[1] = im[2] then Error(); fi;

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

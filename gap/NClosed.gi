
##
## Takes an Majorana representation <rep> and an integer <index> that corresponds to
## an orbit of algebra products (given by rep.setup.pairreps[ index ])
##

InstallGlobalFunction(MAJORANA_NClosedSetUp,

    function(rep, index)

    local   unknowns, dim, new_dim, x, elts, k, i, j, ev, gens, pos, sign, new;

    dim := Size(rep.setup.coords);

    # Add all products that are in the chosen orbit
    for i in [1..dim] do
        for j in [i..dim] do
            if rep.setup.pairorbit[i, j] in [index, -index] then
                Add(rep.setup.coords, [i,j]);
                rep.setup.coordmap[[i,j]] := Size(rep.setup.coords);
            fi;
        od;
    od;

    # Extend all permutations that are stored in the setup record
    for x in rep.setup.pairconjelts do MAJORANA_ExtendPerm(x, rep); od;
    for x in rep.setup.conjelts do MAJORANA_ExtendPerm(x, rep); od;
    for x in rep.generators do MAJORANA_ExtendPerm(x, rep); od;

    new_dim := Size(rep.setup.coords);

    # Extend the matrices of pairorbit and pairconj
    for i in [1..dim] do
        Append(rep.setup.pairorbit[i], [dim + 1 .. new_dim]*0);
        Append(rep.setup.pairconj[i], [dim + 1 .. new_dim]*0);
    od;

    Append(rep.setup.pairorbit, NullMat(new_dim - dim, new_dim));
    Append(rep.setup.pairconj, NullMat(new_dim - dim, new_dim));

    # Calculate the new orbital
    MAJORANA_Orbitals(rep.generators, dim, rep.setup);

    # Add the new algebra product for the chosen orbit
    pos := rep.setup.coordmap[rep.setup.pairreps[index]];
    rep.algebraproducts[index] := SparseMatrix(1, new_dim, [[pos]], [[One(rep.field)]], rep.field);

    # Adjust the existing algebra products and eigenvectors
    for i in [1..Size(rep.algebraproducts)] do
        if not rep.algebraproducts[i] in [false, fail] then
            rep.algebraproducts[i]!.ncols := new_dim;
        fi;
    od;

    for i in [Size(rep.algebraproducts) + 1 .. Size(rep.setup.pairreps)] do
        rep.algebraproducts[i] := false;
        if IsBound(rep.innerproducts) then rep.innerproducts[i] := false; fi;
    od;

    for i in rep.setup.orbitreps do
        for ev in RecNames(rep.evecs[i]) do
            rep.evecs[i].(ev)!.ncols := new_dim;
        od;
    od;

    MAJORANA_NClosedNullspace(rep);

    end );

##
## Takes the last linear system outputted by the main algorithm and uses it
## to create new nullspace vectors
##

InstallGlobalFunction( MAJORANA_NClosedNullspace,

    function(rep)

    local i, j, v, x, pos;

    # Adjust the matrices of the system to the new spanning set
    rep.system.vec!.ncols := Size(rep.setup.coords);
    rep.setup.nullspace.vectors!.ncols := Size(rep.setup.coords);

    for i in [1..Nrows(rep.system.mat)] do
        # If the row of the matrix involves only products that are in the new spanning set
        if ForAll(rep.system.mat!.indices[i], x -> rep.system.unknowns[x] in rep.setup.coords) then
            # The vector <v> will be the new nullspace vector
            v := CertainRows(rep.system.vec, [i]);
            for j in [1..Size(rep.system.mat!.indices[i])] do
                # Add any coefficients coming from the matrix of the system
                x := rep.system.mat!.indices[i, j];
                pos := Position(rep.setup.coords, rep.system.unknowns[x]);
                SetEntry(v, 1, pos, -rep.system.mat!.entries[i, j]);

                rep.setup.nullspace.vectors := UnionOfRows(rep.setup.nullspace.vectors, v);
            od;
        fi;
    od;

    rep.setup.nullspace := ReversedEchelonMatDestructive(rep.setup.nullspace.vectors);

    end );

InstallGlobalFunction( NClosedMajoranaRepresentation,

    function(rep)

    local products, unknowns;

    # Find the positions of the unknown algebra products
    products := Positions(rep.algebraproducts, false);
    if products = [] then return; fi;

    # Added the first of these products to the spanning set of the algebra
    MAJORANA_NClosedSetUp(rep, products[1]);

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        MAJORANA_MainLoop(rep);

        MAJORANA_Fusion(rep, false);

        # MajoranaAlgebraTest(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );
        if IsBound(rep.innerproducts) then
            Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.innerproducts, false)), " unknown inner products ") );
        fi;

        if not false in rep.algebraproducts then
            Info( InfoMajorana, 10, "Success" );
            return;
        fi;

        # If no more algebra products have been found then add the next unknown product to the
        # spanning set of the algebra and run the algorithm again. If all products have been
        # used then quit with a fail.
        if ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            products := Filtered(products, x -> rep.algebraproducts[x] = false);

            if products = [] then
                Info( InfoMajorana, 10, "Fail" );
                return;
            else
                MAJORANA_NClosedSetUp(rep, products[1]);
            fi;
        fi;
    od;

    end );

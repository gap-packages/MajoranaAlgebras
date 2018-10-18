InstallGlobalFunction( MAJORANA_IntersectEigenspaces,

    function(rep)

    local dim, null, i, j, k, evecs_a, evecs_b, z, x, u, v, g, ev, gens;

    dim := Size(rep.setup.coords);

    null := [];

    # Add intersection of eigenspaces to nullspace

    for i in rep.setup.orbitreps do
        for j in [1 .. 3] do
            for k in [j + 1 .. 3] do
                evecs_a := ConvertSparseMatrixToMatrix(rep.evecs[i, j]);
                evecs_b := ConvertSparseMatrixToMatrix(rep.evecs[i, k]);

                z := SumIntersectionMat(evecs_a, evecs_b);

                Append(null, z[2]);
            od;
        od;
    od;

    null := SparseMatrix(null,rep.field);

    null!.ncols := dim;

    # Use eigenvectors to find potentially more nullspace vectors

    for i in Filtered(rep.setup.orbitreps, x -> rep.setup.nullspace.heads[x] = 0) do

        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        for j in [1 .. 3] do

            ev := MAJORANA_FusionTable[1, j + 1];

            for k in [1..Nrows(rep.evecs[i, j])] do
                v := CertainRows(rep.evecs[i, j], [k]);
                x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

                if not x in [fail, false] then
                    null := UnionOfRows(null, x - ev*v);
                fi;
            od;
        od;
    od;

    if Nrows(null) = 0 then return; fi;

    # Find basis and remove null vecs from products and evecs

    null!.ncols := dim;

    null := UnionOfRows(null, rep.setup.nullspace.vectors);

    rep.setup.nullspace := ReversedEchelonMatDestructive(null);

    if rep.setup.nullspace.heads = [] then Error(); fi;

    MAJORANA_RemoveNullspaceNoForm(rep);

    end );

InstallGlobalFunction( MAJORANA_RemoveNullspaceNoForm,

    function(rep)

    local null, g, i, j, x, dim;

    dim := Size(rep.setup.coords);

    if Nrows(rep.setup.nullspace.vectors) = 0 then return; fi;

    null := SparseMatrix(0, dim, [], [], rep.field);

    for g in rep.setup.conjelts do
        for i in [1..Nrows(rep.setup.nullspace.vectors)] do
            x := CertainRows(rep.setup.nullspace.vectors, [i]);
            null := UnionOfRows(null, MAJORANA_ConjugateVec(x, g));
        od;
    od;

    rep.setup.nullspace := ReversedEchelonMatDestructive(null);

    for i in [1..Size(rep.setup.pairreps)] do
        x := Filtered([1..dim], j -> i in rep.setup.pairorbit[j]);
        if ForAll(x, j -> rep.setup.nullspace.heads[j] <> 0) then
            rep.setup.pairreps[i] := fail;
            rep.algebraproducts[i] := fail;
        fi;
    od;

    for i in [1..Size(rep.algebraproducts)] do
        if not rep.algebraproducts[i] in [false, fail] then
            rep.algebraproducts[i] := RemoveMatWithHeads(rep.algebraproducts[i], rep.setup.nullspace);
        fi;
    od;

    for i in rep.setup.orbitreps do
        for j in [1..3] do
            rep.evecs[i, j] := RemoveMatWithHeads(rep.evecs[i, j], rep.setup.nullspace);
            rep.evecs[i, j] := MAJORANA_BasisOfEvecs(rep.evecs[i, j]);
        od;
    od;

    end );

InstallGlobalFunction( MAJORANA_FuseEigenvectorsNoForm,

    function(a, b, u, evals, new, innerproducts, algebraproducts, setup)

    local   dim,
            test,
            new_ev,
            pos,
            x,
            y,
            z;

    dim := Size(setup.coords);

    new_ev := MAJORANA_FusionTable[evals[1] + 1, evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;

    x := MAJORANA_AlgebraProduct(a,b,algebraproducts,setup);

    if x in [false, fail] then return; fi;

    if evals = [2,2] then
        y := MAJORANA_AlgebraProduct(u, x, algebraproducts, setup);

        if y in [fail, false] then return; fi;

        new[1] := MAJORANA_AddEvec(new[1], x - y);
    elif evals = [3,3] then
        y := MAJORANA_AlgebraProduct(u, x, algebraproducts, setup);

        if y in [fail, false] then return; fi;

        z := MAJORANA_AlgebraProduct(u,y,algebraproducts, setup);

        if z in [fail, false] then return; fi;

        new[2] := MAJORANA_AddEvec(new[2], y - z);
        new[1] := MAJORANA_AddEvec(new[1], x - 5*y + 4*z);
    else
        new[pos] := MAJORANA_AddEvec(new[pos],x);
    fi;

    end );

InstallGlobalFunction( MAJORANA_MainLoopNoForm,

    function(rep)

    MAJORANA_Fusion(rep, false);

    MAJORANA_EigenvectorsAlgebraUnknowns(rep);

    MAJORANA_Fusion(rep, false);

    return MAJORANA_UnknownAlgebraProducts(rep, false);

    end);

InstallGlobalFunction( MajoranaRepresentationNoForm,

    function(arg)

    local   rep, unknowns, algebras, main;

    if Size(arg) = 2 then
        arg[3] := "AllAxioms";
    fi;

    rep :=  MAJORANA_SetUp(arg[1], arg[2], arg[3]);

    rep.innerproducts := false;

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        main := MAJORANA_MainLoopNoForm(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );

        if not false in rep.algebraproducts then

            MAJORANA_Fusion(rep, false);
            MAJORANA_IntersectEigenspaces(rep);

            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            Info( InfoMajorana, 10, "Fail" );
            rep.mat := main.mat; rep.vec := main.vec; rep.unknowns := main.unknowns;
            return rep;
        fi;
    od;

    end );

InstallGlobalFunction( NClosedMajoranaRepresentationNoForm,

    function(rep)

    local products, unknowns;

    products := Positions(rep.algebraproducts, false);

    MAJORANA_NClosedSetUp(rep, products[1]);

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        MAJORANA_MainLoopNoForm(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );

        if not false in rep.algebraproducts then
            Info( InfoMajorana, 10, "Success" );
            return;
        fi;

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

InstallGlobalFunction( MajoranaAlgebraTestNoForm,

    function(rep)

    MAJORANA_TestPrimitivity(rep);

    MAJORANA_TestFusion(rep);

    end );

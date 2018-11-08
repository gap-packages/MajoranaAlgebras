InstallGlobalFunction( MAJORANA_IntersectEigenspaces,

    function(rep)

    local dim, null, i, ev, evecs_a, evecs_b, Z, x, u, v, g, conj;

    dim := Size(rep.setup.coords);

    null := [];

    # Add intersection of eigenspaces to nullspace

    for i in rep.setup.orbitreps do
        for ev in Combinations(RecNames(rep.evecs[i]), 2) do
            evecs_a := ConvertSparseMatrixToMatrix(rep.evecs[i].(ev[1]));
            evecs_b := ConvertSparseMatrixToMatrix(rep.evecs[i].(ev[2]));

            Z := SumIntersectionMat(evecs_a, evecs_b);

            Append(null, Z[2]);
        od;
    od;

    null := SparseMatrix(null, Rationals);

    null!.ncols := dim;

    # Use eigenvectors to find potentially more nullspace vectors

    for i in Filtered(rep.setup.orbitreps, x -> rep.setup.nullspace.heads[x] = 0) do

        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        for ev in rep.eigenvalues do
            for v in Iterator( rep.evecs[i].(String(ev)) ) do
                x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

                if not x in [fail, false] then
                    null := UnionOfRows(null, x - ev*v);
                fi;
            od;
        od;
    od;

    # Find basis and remove null vecs from products and evecs

    null := ReversedEchelonMatDestructive(null).vectors;

    if Nrows(null) = 0 then return; fi;

    conj := SparseMatrix(0, dim, [], [], Rationals);
    conj := UnionOfRows(conj, null);

    for g in rep.setup.conjelts do
        for i in [1..Nrows(null)] do
            x := CertainRows(null, [i]);
            conj := UnionOfRows(conj, MAJORANA_ConjugateVec(x, g));
        od;
    od;

    null := UnionOfRows(rep.setup.nullspace.vectors, conj);

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
        for ev in RecNames(rep.evecs[i]) do
            rep.evecs[i].(ev) := RemoveMatWithHeads(rep.evecs[i].(ev), rep.setup.nullspace);
            rep.evecs[i].(ev) := MAJORANA_BasisOfEvecs(rep.evecs[i].(ev));
        od;
    od;

    MAJORANA_IntersectEigenspaces(rep);

    end );

InstallGlobalFunction( MAJORANA_FuseEigenvectorsNoForm,

    function(a, b, u, evals, new, rep)

    local   ev, x, y, z;

    ev := MAJORANA_FusionTable[ evals ][1];

    x := MAJORANA_AlgebraProduct(a, b, rep.algebraproducts, rep.setup);

    if x in [false, fail] then return; fi;

    if evals = ["1/4", "1/4"] then
        y := MAJORANA_AlgebraProduct(u, x, rep.algebraproducts, rep.setup);

        if y in [fail, false] then return; fi;

        new.("0") := MAJORANA_AddEvec(new.("0"), x - y);
    elif evals = ["1/32", "1/32"] then
        y := MAJORANA_AlgebraProduct(u, x, rep.algebraproducts, rep.setup);

        if y in [fail, false] then return; fi;

        z := MAJORANA_AlgebraProduct(u, y, rep.algebraproducts, rep.setup);

        if z in [fail, false] then return; fi;

        new.("1/4") := MAJORANA_AddEvec(new.("1/4"), y - z);
        new.("0") := MAJORANA_AddEvec(new.("0"), x - 5*y + 4*z);
    else
        new.(String(ev)) := MAJORANA_AddEvec(new.(String(ev)), x);
    fi;

    end );

InstallGlobalFunction( MAJORANA_MainLoopNoForm,

    function(rep)

    MAJORANA_Fusion(rep, false);

    MAJORANA_EigenvectorsAlgebraUnknowns(false, rep);

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

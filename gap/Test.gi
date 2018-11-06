##
## The Majorana fusion table, implemented as a hashmap.
##

BindGlobal("MAJORANA_FusionTable", HashMap( 16 ) );

MAJORANA_FusionTable[ [ "0", "0" ] ]        := [ 0 ];
MAJORANA_FusionTable[ [ "0", "1/4" ] ]      := [ 1/4 ];
MAJORANA_FusionTable[ [ "0", "1/32" ] ]     := [ 1/32 ];

MAJORANA_FusionTable[ [ "1/4", "0" ] ]      := [ 1/4 ];
MAJORANA_FusionTable[ [ "1/4", "1/4" ] ]    := [ 1, 0 ];
MAJORANA_FusionTable[ [ "1/4", "1/32" ] ]   := [ 1/32 ];

MAJORANA_FusionTable[ [ "1/32", "0" ] ]     := [ 1/32 ];
MAJORANA_FusionTable[ [ "1/32", "1/4" ] ]   := [ 1/32 ];
MAJORANA_FusionTable[ [ "1/32", "1/32" ] ]  := [ 1, 0, 1/4 ];


InstallGlobalFunction(MAJORANA_TestEvecs,

    function(rep)

    local   u, v, i, ev, x, y;

    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

        for ev in rep.eigenvalues do
            for v in Iterator( rep.evecs[i].(String(ev)) ) do

                x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

                if not x in [ev*v, false, fail] then
                    Error("evecs");
                fi;
            od;
        od;
    od;

    return true;

    end );

InstallGlobalFunction( MAJORANA_TestFusionAxis,

    function(u, evecs, rep)

    local   dim, evals, new, ev, a, b, x;

    dim := Size(rep.setup.coords);

    for evals in UnorderedTuples( RecNames(evecs), 2 ) do

        new := rec();

        for ev in RecNames(evecs) do
            new.(ev) := SparseMatrix(0, dim, [], [], Rationals);
        od;

        for a in Iterator( evecs.(evals[1]) ) do
            for b in Iterator( evecs.(evals[2]) ) do
                MAJORANA_FuseEigenvectorsNoForm(  a, b, u, evals, new,
                                            rep.innerproducts,
                                            rep.algebraproducts,
                                            rep.setup );
            od;
        od;

        for ev in rep.eigenvalues do
            new.(String(ev)) := EchelonMatDestructive(new.(String(ev))).vectors;

            for a in Iterator( new.(String(ev)) ) do
                x := MAJORANA_AlgebraProduct(u, a, rep.algebraproducts, rep.setup);
                if not x in [ev*a, false, fail] then
                    Error("The algebra does not obey the fusion rules");
                fi;
            od;
        od;
    od;

    end );

# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(rep)

    local   dim, u, i;

    dim := Size(rep.setup.coords);

    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        MAJORANA_TestFusionAxis(u, rep.evecs[i], rep);
    od;

    return true;

    end );

InstallGlobalFunction(MajoranaAlgebraTest,

    function(rep)

    MAJORANA_TestAxiomM1(rep);

    MAJORANA_TestFusion(rep);

    MAJORANA_TestPrimitivity(rep);

    return true;

    end );

# Checks if bilinear and algebra products obey Fusion, outputs a list which is empty if they do obey the axiom

InstallGlobalFunction(MAJORANA_TestAxiomM1,

    function(rep)

    local   ErrorM1,    # setup of indices which do not obey Fusion
            j,          # loop over algebra products
            k,          # loop over setup.coords
            l,
            p,          # second product
            dim,        # size of setup.coords
            x,          # first inner product
            y,          # second inner product
            u,          # vectors
            w,          #
            v;          #

    dim := Size(rep.setup.coords);

    ErrorM1:=[];

    for j in Filtered([1..Size(rep.algebraproducts)], i -> rep.algebraproducts[i] <> fail) do
        if not rep.algebraproducts[j] in [false, fail] then
            for k in Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0) do
                for l in [rep.setup.pairreps[j], Reversed(rep.setup.pairreps[j])] do

                    u := SparseMatrix(1, dim, [[l[1]]], [[1]], Rationals);
                    v := SparseMatrix(1, dim, [[l[2]]], [[1]], Rationals);
                    w := SparseMatrix(1, dim, [[k]], [[1]], Rationals);

                    p := MAJORANA_AlgebraProduct(v,w,rep.algebraproducts,rep.setup);

                    if not p in [fail, false] then
                        x := MAJORANA_InnerProduct(u,p,rep.innerproducts, rep.setup);
                        y := MAJORANA_InnerProduct(rep.algebraproducts[j],w,rep.innerproducts, rep.setup);

                        if x <> false and y <> false and x <> y then
                            # return false;
                            Error("Axiom M1");
                            Add(ErrorM1,[l[1], l[2] ,k]);
                        fi;

                    fi;
                od;
            od;
        fi;
    od;

    if Size(ErrorM1) > 0 then Error("Axiom M1"); fi;

    return true;

    end );

InstallGlobalFunction( MAJORANA_TestPrimitivity,

    function(rep)

    local i, dim, u, v, j, x, mat, espace, basis;

    if false in rep.algebraproducts then return fail; fi;

    dim := Size(rep.setup.coords);

    basis := Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0);

    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        mat := SparseMatrix(0, Size(basis), [], [], Rationals);

        for j in basis do
            v := SparseMatrix(1, dim, [[j]], [[1]], Rationals);

            x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

            mat := UnionOfRows(mat, x);
        od;

        espace := KernelMat(mat - SparseIdentityMatrix(dim, Rationals));

        if Nrows(espace.relations) <> 1 then
            Error("Primitivity");
        fi;
    od;

    return true;

    end );

InstallGlobalFunction( MAJORANA_IsComplete,

    function(rep)

    if false in rep.algebraproducts then
        return false;
    else
        return true;
    fi;

    end );

InstallGlobalFunction( MAJORANA_TestPositiveDefiniteForm,

    function(rep)

    local dim, gram;

    dim := Size(rep.setup.coords);

    gram := MAJORANA_FillGramMatrix(Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0), rep.innerproducts, rep.setup);

    if MAJORANA_PositiveDefinite(ConvertSparseMatrixToMatrix(gram)) < 0 then
        return false;
    else
        return true;
    fi;

    end );

InstallGlobalFunction(MAJORANA_TestAxiomM2,

    function(rep) # Tests that the algebra obeys axiom M2

    local   B,      # matrix of inner products
            dim,    # size of setup.coords
            j,      # loop through setup.coords
            k,      #
            l,      #
            m,      #
            a,      # vectors
            b,      #
            c,      #
            d,      #
            x0,     # products
            x1,     #
            x2,     #
            x3,     #
            basis;

    dim:=Size(rep.setup.coords);

    basis := Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0);

    B:=NullMat(dim^2,dim^2);

    for j in basis do
        for k in basis do
            for l in basis do
                for m in basis do

                    a := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
                    b := SparseMatrix(1, dim, [[k]], [[1]], Rationals);
                    c := SparseMatrix(1, dim, [[l]], [[1]], Rationals);
                    d := SparseMatrix(1, dim, [[m]], [[1]], Rationals);

                    x0 := MAJORANA_AlgebraProduct(a,c,rep.algebraproducts,rep.setup);
                    x1 := MAJORANA_AlgebraProduct(b,d,rep.algebraproducts,rep.setup);
                    x2 := MAJORANA_AlgebraProduct(b,c,rep.algebraproducts,rep.setup);
                    x3 := MAJORANA_AlgebraProduct(a,d,rep.algebraproducts,rep.setup);

                    B[dim*(j-1) + k, dim*(l-1) +m]:=
                          MAJORANA_InnerProduct(x0,x1,rep.innerproducts, rep.setup)
                        - MAJORANA_InnerProduct(x2,x3,rep.innerproducts, rep.setup);
                od;
            od;
        od;
    od;

    if MAJORANA_PositiveDefinite(B) < 0 then
        return false;
    else
        return true;
    fi;

    end );

InstallGlobalFunction( MAJORANA_TestSetup,

    function(rep)

    local dim, i, j, k, g, sign, sign_k, im;

    dim := Size(rep.setup.coords);

    for i in [1 .. dim] do
        for j in [i .. dim] do
            k := rep.setup.pairorbit[i, j];
            g := rep.setup.pairconj[i, j];
            g := rep.setup.pairconjelts[g];

            sign_k := 1;

            if k < 0 then k := -k; sign_k := -1; fi;

            im := g{rep.setup.pairreps[k]};

            sign := 1;

            if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
            if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;

            if SortedList(im) <> [i,j] then Error("Does not conjugate to correct pair"); fi;

            if sign <> sign_k then Error("Sign error"); fi;
        od;
    od;

    end );

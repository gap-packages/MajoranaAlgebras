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

##
## The main test func, returns true if the algebra is a Majorana algebra and
## enters a break loop if not.
##

InstallGlobalFunction(MajoranaAlgebraTest,

    function(rep)

    if IsBound(rep.innerproducts) then
        if MAJORANA_TestFrobeniusForm(rep) = false then
            return false;
        elif MAJORANA_TestInnerProduct(rep) = false then
            return false;
        fi;
    fi;

    if MAJORANA_TestPrimitivity(rep) = false then
        return false;
    fi;

    return true;

    end );

##
## Tests that the vectors stored in rep.evecs are indeed eigenvectors
##

InstallGlobalFunction(MAJORANA_TestEvecs,

    function(rep)

    local   u, v, i, ev, x, y;

    # Loop over the representatives of the orbits of G on T
    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], rep.field);

        # For each of the three eigenvalues 0, 1/4, 1/32, check that the eqn u*v = ev*v holds
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

# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(rep)

    local   dim, u, i, evals, new, a, b, ev, x;

    dim := Size(rep.setup.coords);

    # Loop over the representatives of the orbits of G on T
    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, dim, [[i]], [[1]], rep.field);

        # Loop over all pairs of eigenvalues and perform fusion, storing new vecs in <new>
        for evals in UnorderedTuples( RecNames(rep.evecs[i]), 2 ) do

            new := rec();

            for ev in RecNames(rep.evecs[i]) do
                new.(ev) := SparseMatrix(0, dim, [], [], rep.field);
            od;

            for a in Iterator( rep.evecs[i].(evals[1]) ) do
                for b in Iterator( rep.evecs[i].(evals[2]) ) do
                    MAJORANA_FuseEigenvectorsNoForm(  a, b, u, evals, new, rep );
                od;
            od;

            # Check whether the new eigenvectors found by fusion are indeed eigevectors
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
    od;

    return true;

    end );

#
## Checks if algebra obeys axiom M1
##

InstallGlobalFunction(MAJORANA_TestFrobeniusForm,

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

                    u := SparseMatrix(1, dim, [[l[1]]], [[1]], rep.field);
                    v := SparseMatrix(1, dim, [[l[2]]], [[1]], rep.field);
                    w := SparseMatrix(1, dim, [[k]], [[1]], rep.field);

                    p := MAJORANA_AlgebraProduct(v,w,rep.algebraproducts,rep.setup);

                    if not p in [fail, false] then
                        x := MAJORANA_InnerProduct(u,p,rep.innerproducts, rep.setup);
                        y := MAJORANA_InnerProduct(rep.algebraproducts[j],w,rep.innerproducts, rep.setup);

                        if x <> false and y <> false and x <> y then
                            return false;
                            # Error("Axiom M1");
                            Add(ErrorM1,[l[1], l[2] ,k]);
                        fi;

                    fi;
                od;
            od;
        fi;
    od;

    # if Size(ErrorM1) > 0 then Error("Axiom M1"); fi;

    return true;

    end );

InstallGlobalFunction( MAJORANA_TestPrimitivity,

    function(rep)

    local i, dim, u, v, j, x, mat, espace, basis;

    if false in rep.algebraproducts then return fail; fi;

    if MAJORANA_Dimension(rep) = 0 then return true; fi;

    dim := Size(rep.setup.coords);

    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], rep.field);

        mat := SparseMatrix(0, dim, [], [], rep.field);

        for j in [1 .. dim] do
            v := SparseMatrix(1, dim, [[j]], [[1]], rep.field);

            x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

            mat := UnionOfRows(mat, x);
        od;

        espace := KernelMat(mat - SparseIdentityMatrix(dim, rep.field));

        if Nrows(espace.relations) <> 1 then
            Error("Primitivity");
        fi;
    od;

    return true;

    end );

InstallGlobalFunction( MAJORANA_TestInnerProduct,

    function(rep)

    local dim, gram;

    if IsBound(rep.innerproducts) and not false in rep.innerproducts then
        dim := Size(rep.setup.coords);

    gram := MAJORANA_FillGramMatrix(Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0), rep);

    if MAJORANA_PositiveDefinite(ConvertSparseMatrixToMatrix(gram), rep.field) < 0 then
        return false;
    else
        return fail;
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

    dim:=Size(rep.setup.coords);;

    B:=NullMat(dim^2,dim^2)*One(rep.field);;

    for j in [1..dim] do
        for k in [1..dim] do
            for l in [1..dim] do
                for m in [1..dim] do

                    a := SparseMatrix(1, dim, [[j]], [[One(rep.field)]], rep.field);
                    b := SparseMatrix(1, dim, [[k]], [[One(rep.field)]], rep.field);
                    c := SparseMatrix(1, dim, [[l]], [[One(rep.field)]], rep.field);
                    d := SparseMatrix(1, dim, [[m]], [[One(rep.field)]], rep.field);

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

    if MAJORANA_PositiveDefinite(B, rep.field) < 0 then
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

    return true;

    end );

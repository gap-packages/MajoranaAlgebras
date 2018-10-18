# Table of fusion rules

BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4,     0, 1/32]
            ,[ 1/32, 1/32,  1/32, 1/4 ] ]);

InstallGlobalFunction( MAJORANA_TestFusionTable,

    function(adjoint, evals, table, rep)

    local n, field, eigenvectors, ev, eigenspaces, i, j, k, l, u, v, x, space;

    n := Nrows(adjoint);

    field := adjoint!.ring;

    eigenvectors := [];

    for ev in evals do
        Add(eigenvectors, KernelMat(adjoint - ev*SparseIdentityMatrix(15, field)).relations);
    od;

    eigenspaces := List(eigenvectors, x -> ConvertSparseMatrixToMatrix(x));;
    eigenspaces := List(eigenspaces, x -> Subspace(field^n, x));;

    for i in [1..Size(evals)] do
        for j in [1..Size(evals)] do
            space := Subspace(field^n, Union(List(eigenspaces, Basis){table[i][j]}));

            for k in [1..Nrows(eigenvectors[i])] do
                u := CertainRows(eigenvectors[i], [k]);
                for l in [1..Nrows(eigenvectors[j])] do
                    v := CertainRows(eigenvectors[j], [l]);

                    x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);

                    if not ConvertSparseMatrixToMatrix(x)[1] in space then
                        Error("Fusion table is not obeyed");
                    fi;
                od;
            od;
        od;
    od;

    end );

InstallGlobalFunction(MAJORANA_TestEvecs,

    function(rep)

    local   u, v, i, j, k, ev, x, y;

    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], rep.field);

        for j in [1..3] do
            ev := MAJORANA_FusionTable[1, j + 1];
            for k in [1..Nrows(rep.evecs[i, j])] do
                v := CertainRows(rep.evecs[i, j], [k]);

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

    function(u, evecs, innerproducts, algebraproducts, setup)

    local   dim,            # size of setup.coords
            field,
            evals,
            new,
            a,              # first eigenvalue
            b,              # second eigenvalue
            ev_a,           # a - eigenvectors
            ev_b,           # b - eigenvectors
            ev,             # new eigenvalue
            i,              # loop over T
            j,
            k,
            x,              # product of eigenvectors
            y;              # product of x with u

    field := algebraproducts[1]!.ring;

    for i in setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], field);

        for evals in [[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]] do

            new := [0,0,0];

            for j in [1..3] do
                new[j] := SparseMatrix(0, dim, [], [], field);
            od;

            ev_a := evecs[i][evals[1]];
            ev_b := evecs[i][evals[2]];

            for j in [1..Nrows(ev_a)] do
                a := CertainRows(ev_a, [j]);
                for k in [1..Nrows(ev_b)] do
                    b := CertainRows(ev_b, [k]);
                    MAJORANA_FuseEigenvectorsNoForm(  a, b, i, evals, new,
                                                innerproducts,
                                                algebraproducts,
                                                setup );
                od;
            od;

            for j in [1..3] do
                ev := MAJORANA_FusionTable[1][j + 1];

                # new[j] := EchelonMatDestructive(new[j]).vectors;

                for k in [1..Nrows(new[j])] do
                    a := CertainRows(new[j], [k]);
                    x := MAJORANA_AlgebraProduct(u, a, algebraproducts, setup);
                    if x <> ev*a and x <> false then
                        Error("Algebra does not obey the fusion rules");
                    fi;
                od;
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

        MAJORANA_TestFusionAxis(u, rep.evecs[i], rep.innerproducts, rep.algebraproducts, rep.setup);
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

InstallGlobalFunction(MAJORANA_TestOrthogonality,

    function(rep)

    # Tests that eigenspaces are orthogonal with respect to the inner product

        local   errorortho, # list of indices which do not obey orthogonality of eigenvectors
                u,          # vector with 1 in j th position
                a,          # first eigenvalue
                b,          # second eigenvalue
                ev_a,       # list of a - eigenvectors
                ev_b,       # list of b - eigenvectors
                i,
                j,          # loop over T
                k,
                v,          # a - eigenvector
                w,          # b - eigenvector
                x;          # inner product

        errorortho := [];

        for i in rep.setup.orbitreps do

            u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], rep.field);

            for a in [1..3] do

                # orthogonality with 1-eigenvectors

                ev_a := rep.evecs[i][a];

                for j in [1..Nrows(ev_a)] do
                    v := CertainRows(ev_a, [j]);
                    x := MAJORANA_InnerProduct(u, v, rep.innerproducts, rep.setup);

                    if (x <> false) and (x <> 0) then

                        Add(errorortho, [i,0,a,u,v]);
                    fi;
                od;

                # orthogonality with all other eigenvectors

                for b in [a+1..3] do

                    ev_b := rep.evecs[i][b];

                    for j in [1..Nrows(ev_a)] do
                        v := CertainRows(ev_a, [j]);
                        for k in [1..Nrows(ev_b)] do
                            w := CertainRows(ev_b, [k]);

                            x := MAJORANA_InnerProduct(v, w, rep.innerproducts, rep.setup);

                            if (x <> false) and (x <> 0) then
                                Add(errorortho, [i,a,b,v,w]);
                            fi;
                        od;
                    od;
                od;
            od;
        od;

        if Size(errorortho) > 0 then Error("Orthog"); fi;

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

                    u := SparseMatrix(1, dim, [[l[1]]], [[1]], rep.field);
                    v := SparseMatrix(1, dim, [[l[2]]], [[1]], rep.field);
                    w := SparseMatrix(1, dim, [[k]], [[1]], rep.field);

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

    gram := MAJORANA_FillGramMatrix(Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0), rep);

    if MAJORANA_PositiveDefinite(ConvertSparseMatrixToMatrix(gram), rep.field) < 0 then
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

    end );

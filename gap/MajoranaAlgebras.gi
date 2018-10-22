#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

InstallGlobalFunction( MAJORANA_ChangeFieldOfRep,

    function(rep, field)

    local i, j;

    rep.field := field;

    for i in [1..Size(rep.algebraproducts)] do
        if not rep.algebraproducts[i] in [false, fail] then
            rep.algebraproducts[i] := rep.algebraproducts[i]*One(field);
            rep.algebraproducts[i]!.ring := field;
        fi;
    od;

    for i in [1..Size(rep.innerproducts)] do
        if not rep.innerproducts[i] in [false, fail] then
            rep.innerproducts[i] := rep.innerproducts[i]*One(field);
        fi;
    od;

    for i in rep.setup.orbitreps do
        rep.evecs[i] := rep.evecs[i]*One(field);
        for j in [1..3] do
            rep.evecs[i][j]!.ring := field;
        od;
    od;

    rep.setup.nullspace.vectors := rep.setup.nullspace.vectors*One(field);
    rep.setup.nullspace.vectors!.ring := field;

    rep.mat!.ring := field;
    rep.mat!.entries := List(rep.mat!.entries, x -> x*One(field));
    rep.vec!.ring := field;
    rep.vec!.entries := List(rep.vec!.entries, x -> x*One(field));

    end );

InstallGlobalFunction( MAJORANA_SetValue,

    function(rep, function_field, indets, vals)

    local i, j;

    for i in [1..Size(rep.algebraproducts)] do
        if not rep.algebraproducts[i] in [false, fail] then
            rep.algebraproducts[i]!.entries[1] := List(rep.algebraproducts[i]!.entries[1], x -> Value(x, indets, vals));
        fi;
    od;

    for i in [1..Size(rep.innerproducts)] do
        if not rep.innerproducts[i] in [false, fail] then
            rep.innerproducts[i] := Value(rep.innerproducts[i], indets, vals);
        fi;
    od;

    for i in rep.setup.orbitreps do
        for j in [1..3] do
            rep.evecs[i][j]!.entries := List(rep.evecs[i][j]!.entries, x -> List(x, y -> Value(y, indets, vals)));
        od;
    od;

    rep.mat!.entries := List(rep.mat!.entries, x -> List(x, y -> Value(y, indets, vals)) );
    rep.vec!.entries := List(rep.vec!.entries, x -> List(x, y -> Value(y, indets, vals)));

    rep.setup.nullspace.vectors!.entries := List(rep.setup.nullspace.vectors!.entries, x -> List(x, y -> Value(y, indets, vals)));

    MAJORANA_ChangeFieldOfRep(rep, CoefficientsRing(function_field));

    end );






InstallGlobalFunction(MajoranaRepresentation,

function(arg)

    local   rep, unknowns, main;

    if Size(arg) = 2 then arg[3] := "AllAxioms"; fi;

    rep :=  MAJORANA_SetUp(arg[1], arg[2],  arg[3]);

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

##
## The main part of the algorithm that is looped over in the main function.
##

InstallGlobalFunction(MAJORANA_MainLoop,

    function(rep)

    MAJORANA_AxiomM1(rep);

    MAJORANA_Fusion(rep);

    if MAJORANA_Dimension(rep) in [0,1] then return true; fi;

    MAJORANA_EigenvectorsAlgebraUnknowns(rep);

    MAJORANA_AxiomM1(rep);

    MAJORANA_Fusion(rep);

    if MAJORANA_Dimension(rep) in [0,1] then return true; fi;

    return MAJORANA_UnknownAlgebraProducts(rep);

    end);

##
## Uses axiom M1 to find inner products
##

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(rep)

    local   dim, mat, vec, i, j, k, u, v, w, x, y, z, eq, unknowns;

    if not false in rep.innerproducts then
        return;
    fi;

    Info(   InfoMajorana, 50, "Axiom M1");

    dim := Size(rep.setup.coords);
    unknowns := Positions(rep.innerproducts, false);

    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);

    for i in Filtered([1..Size(rep.algebraproducts)], x -> not rep.algebraproducts[x] in [false, fail]) do
        x := rep.algebraproducts[i];
        for j in [rep.setup.pairreps[i], Reversed(rep.setup.pairreps[i])] do

            u := SparseMatrix(1, dim, [[ j[1] ]], [[ 1 ]], Rationals);
            v := SparseMatrix(1, dim, [[ j[2] ]], [[ 1 ]], Rationals);

            for k in Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0) do

                w := SparseMatrix(1, dim, [[ k ]], [[ 1 ]], Rationals);

                y := MAJORANA_AlgebraProduct(v, w, rep.algebraproducts, rep.setup);

                if not y in [fail, false] then

                    eq := MAJORANA_SeparateInnerProduct(w, x, unknowns, rep);
                    eq := eq - MAJORANA_SeparateInnerProduct(y, u, unknowns, rep);

                    if Size(eq[1]!.indices[1]) = 1 then
                        z := MAJORANA_SingleInnerSolution(  eq, mat, vec,
                                                            unknowns,
                                                            rep.innerproducts);

                        mat := z.mat; vec := z.vec; unknowns := z.unknowns;

                        if unknowns = [] then
                            MAJORANA_CheckNullSpace(rep); return;
                        fi;
                    elif eq[1]!.indices[1] <> [] then
                        eq := eq*(1/eq[1]!.entries[1, 1]);
                        if not _IsRowOfSparseMatrix(mat, eq[1]) then
                            mat := UnionOfRows(mat, eq[1]);
                            vec := UnionOfRows(vec, eq[2]);
                        fi;
                    fi;
                fi;
            od;
        od;
    od;

    x := MAJORANA_SolutionInnerProducts(mat,vec,unknowns,rep.innerproducts);

    if x.unknowns = [] then
        MAJORANA_CheckNullSpace(rep);
    fi;

    return rec( mat := x.mat, vec := x.vec, unknowns := x.unknowns);

    end );

##
## Finds new eigenvectors using the fusion rules
##

InstallGlobalFunction( MAJORANA_Fusion,

function(arg)

    local   i, j, k, a, b, dim, new, evals, rep, FUSE, u;

    rep := arg[1];
    FUSE := MAJORANA_FuseEigenvectors;

    if Size(arg) = 2 and arg[2] = false then
        FUSE := MAJORANA_FuseEigenvectorsNoForm;
    fi;

    dim := Size(rep.setup.coords);

    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        while true do
            Info(   InfoMajorana, 50, STRINGIFY("Fusion of ", i, " evecs")) ;

            new := [0,0,0];

            for j in [1..3] do
                new[j] := CopyMat(rep.evecs[i, j]);
            od;

            for evals in [[1,1], [1,2], [1,3], [2,3], [2,2], [3,3]] do
                for j in [1..Nrows(rep.evecs[i, evals[1]])] do

                    a := CertainRows(rep.evecs[i, evals[1]], [j]);

                    for k in [1..Nrows(rep.evecs[i, evals[2]])] do

                        b := CertainRows(rep.evecs[i, evals[2]], [k]);

                        FUSE(a, b, u, evals, new, rep.innerproducts, rep.algebraproducts, rep.setup);
                    od;
                od;
            od;

            for j in [1..3] do
                new[j] := MAJORANA_BasisOfEvecs(new[j]);
            od;

            if ForAll([1..3], j -> Nrows(new[j]) = Nrows(rep.evecs[i, j])) then
                break;
            fi;

            rep.evecs[i] := new;

            MAJORANA_IntersectEigenspaces(rep);
            MAJORANA_IntersectEigenspaces(rep);
        od;
    od;

    end );

##
## Uses eigenvectors, the resurrection principle and nullspace vectors to find
## algebra product values
##

InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(arg)

    local   dim, x, y, i, j, k, l, evals, mat, vec, unknowns, u, a, b, c, bad, list, evecs_a, evecs_b, index, n, rep, evals_list;

    rep := arg[1];
    evals_list := [[1,2],[2,1],[1,3],[2,3]];

    if Size(arg) = 2 and arg[2] = false then evals_list := [[1,2], [1,3]]; fi;

    dim := Size(rep.setup.coords);

    n := Size(Positions(rep.algebraproducts, false));

    x := MAJORANA_EigenvectorsAlgebraUnknowns(rep);

    if not false in rep.algebraproducts then return true; fi;

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    x := MAJORANA_NullspaceUnknowns(mat, vec, unknowns, rep);

    if not false in rep.algebraproducts then return true; fi;

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    Info(   InfoMajorana, 50, "Building resurrection");

    for evals in evals_list do
        for i in rep.setup.orbitreps do

            evecs_a := rep.evecs[i, evals[1]];
            evecs_b := rep.evecs[i, evals[2]];

            u := SparseMatrix(1, dim, [[i]], [[1]], rep.field);

            list := [,];

            list[1] := List([1..Nrows(evecs_a)], i -> []);
            list[2] := List([1..Nrows(evecs_a)], i -> []);

            for j in [1..Nrows(evecs_a)] do
                bad := MAJORANA_FindBadIndices(CertainRows(evecs_a,[j]), rep.algebraproducts, rep.setup);
                for k in [1..Nrows(evecs_b)] do
                    x := Size(Intersection(bad, evecs_b!.indices[k]));
                    if x = 1 then
                        Add(list[1, j], k);
                    elif x > 1 then
                        Add(list[2, j], k);
                    fi;
                od;
            od;

            for index in [1,2] do
                for j in [1..Nrows(evecs_a)] do
                    c := CertainRows(evecs_a, [j]);
                    for k in list[index, j] do
                        b := CertainRows(evecs_b, [k]);

                        bad := MAJORANA_FindBadIndices(c, rep.algebraproducts, rep.setup);

                        for l in [1..Nrows(evecs_a)] do

                            a := CertainRows(evecs_a, [l]);

                            if CertainColumns(a, bad) = CertainColumns(b, bad) then

                                x := MAJORANA_Resurrection(  u, a, b, c, evals,
                                                        unknowns, rep);

                                if x <> false and x[1]!.indices[1] <> [] then
                                    if Size(x[1]!.indices[1]) = 1 then
                                        y := MAJORANA_SolveSingleSolution( x,
                                                            mat, vec, unknowns,
                                                            rep.algebraproducts,
                                                            rep.setup);

                                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;

                                        if not false in rep.algebraproducts then return true; fi;
                                    elif not _IsRowOfSparseMatrix(mat, x[1]) then
                                        mat := UnionOfRows(mat, x[1]);
                                        vec := UnionOfRows(vec, x[2]);
                                    fi;
                                fi;
                            fi;
                        od;

                        if Nrows(mat) > Ncols(mat) or Nrows(mat) > 8000 then
                            x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);

                            mat := x.mat; vec := x.vec; unknowns := x.unknowns;

                            if not false in rep.algebraproducts then return true; fi;
                        fi;

                    od;
                od;
            od;
        od;
    od;

    x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    if n = Size(Positions(rep.algebraproducts, false)) then
        return MAJORANA_AllConjugates(mat, vec, unknowns, rep);
    else
        return rec(mat := mat, vec := vec, unknowns := unknowns);
    fi;

    end );

##
## Finds the indices i such that v_i*v is not known
##

InstallGlobalFunction(MAJORANA_FindBadIndices,

    function(v, algebraproducts, setup)

    local   i,
            j,
            k,
            dim,
            list,
            bad;

    bad := [];
    dim := Size(setup.coords);
    list := [1..dim];

    for i in v!.indices[1] do
        for j in list do
            k :=  setup.pairorbit[i, j];

            if k < 0 then k := -k; fi;

            if algebraproducts[k] = false then
                Add(bad,j);
                list := Difference(list,[j]);
            fi;
        od;
    od;

    Sort(bad);

    return bad;

    end );

##
## Add a new eigenvector to a list of existing eigenvectors if not already a member
##

InstallGlobalFunction(MAJORANA_AddEvec,

    function(mat, x)

    if x!.indices[1] = [] then return mat; fi;

    x!.entries[1] := x!.entries[1]/x!.entries[1, Size(x!.entries[1])];

    if _IsRowOfSparseMatrix(mat, x) then
        return mat;
    else
        return UnionOfRows(mat, x);
    fi;

    end);

##
## Given two eigenvectors, if possible, finds product and adds it to appropriate set of evecs.
##

InstallGlobalFunction( MAJORANA_FuseEigenvectors,

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

    if x = false then return; fi;

    if x = fail then return; fi;

    if evals = [2,2] then
        y := MAJORANA_InnerProduct(a,b,innerproducts,setup);

        if y = false then return; fi;

        new[1] := MAJORANA_AddEvec(new[1], x - (1/4)*u*y);
    elif evals = [3,3] then
        y := MAJORANA_InnerProduct(a,b,innerproducts,setup);

        if y = false then return; fi;

        z := MAJORANA_AlgebraProduct(u,x,algebraproducts, setup);

        if z in [false,  fail] then return; fi;

        new[2] := MAJORANA_AddEvec(new[2], z - (1/32)*u*y);
        new[1] := MAJORANA_AddEvec(new[1], x + (3/32)*u*y - 4*z);
    else
        new[pos] := MAJORANA_AddEvec(new[pos],x);
    fi;

    end );

##
## Returns true if we have full eigenspace decomposition, returns false otherwise
##

InstallGlobalFunction( MAJORANA_CheckBasis,

    function(dim, evecs, rep)

    local i, basis;

    if rep.innerproducts = false then return false; fi;

    if Sum(List(evecs, Nrows)) + Nrows(rep.setup.nullspace.vectors) < dim - 1 then
        return false;
    fi;

    return true;

    end );

##
## Takes a matrix <mat> and conjugates it with respect to the signed perm <g>
##

InstallGlobalFunction( MAJORANA_ConjugateVec,

    function(mat,g)

    local   i,
            k,
            sign,
            res,
            pos;

    if ForAll(mat!.indices[1], i -> g[i] = i) then return mat; fi;

    res := SparseMatrix(1, Ncols(mat), [[]], [[]], mat!.ring);

    for i in [1..Size(mat!.indices[1])] do

        k := g[mat!.indices[1, i]];

        sign := 1;

        if k < 0 then k := -k; sign := -1; fi;

        pos := PositionSorted(res!.indices[1], k);

        Add(res!.indices[1], k, pos);
        Add(res!.entries[1], sign*mat!.entries[1, i], pos);
    od;

    return res;

    end );

##
## Takes two vectors and returns their algebra product if known and false or fail if not
##

InstallGlobalFunction(  MAJORANA_AlgebraProduct,

        function(u,v,algebraproducts,setup) # If all the relevant products are known, returns the algebra product of u and v. If not, returns 0

        local   i,      # loop over u
                j,      # loop over v
                k,      # pair orbit index
                x,      # algebra product
                g,      # conjugating element
                sign,   # correct sign of 5A axes
                vec,    # output vec
                vecs,
                elts,
                dim,
                pos;

        dim := Ncols(u);;

        vec := SparseMatrix(1, dim, [[]], [[]], algebraproducts[1]!.ring);

        elts := [];
        vecs := [];

        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do

                k := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

                if k > 0 then
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;

                x := algebraproducts[k];

                if not x in [fail, false] then

                    g := setup.pairconj[u!.indices[1, i], v!.indices[1, j]];

                    pos := Position(elts,g);

                    if pos <> fail then
                        AddRow(x!.indices[1],  sign*u!.entries[1, i]*v!.entries[1, j]*x!.entries[1],  vecs[pos]!.indices, vecs[pos]!.entries, 1);
                    else
                        Add(elts, g);
                        Add(vecs, CopyMat(sign*u!.entries[1, i]*v!.entries[1, j]*x));
                    fi;
                elif x = false then
                    # cannot calculate product
                    return false;
                else # product with a vector in the nullspace
                    return fail;
                fi;
            od;
        od;

        for i in [1..Size(elts)] do
            x := MAJORANA_ConjugateVec(vecs[i], setup.pairconjelts[elts[i]]);
            AddRow(x!.indices[1], x!.entries[1], vec!.indices, vec!.entries, 1);
        od;

        return RemoveMatWithHeads(vec, setup.nullspace);

        end );

##
## Takes two vectors and returns their inner product if known and false or fail if not
##

InstallGlobalFunction(  MAJORANA_InnerProduct,

    function(u, v, innerproducts, setup)

        local   i,              # loop over u
                j,              # loop over v
                k,              # pair orbit index
                sign,           # correct for 5A axes
                sum;            # output value

        sum := Zero(u!.ring);

        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do
                k := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

                if k > 0 then
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;

                if innerproducts[k] <> false then
                    sum := sum + sign*u!.entries[1, i]*v!.entries[1, j]*innerproducts[k];
                else
                    return false;
                fi;
            od;
        od;

        return sum;

        end );

##
## For a given range of basis vectors, fill in the Gram matrix wrt those vectors
##

InstallGlobalFunction(MAJORANA_FillGramMatrix,

function(range, rep)

    local   i, j, k, mat, l;

    l := Length(range);

    mat := SparseZeroMatrix(l, l, rep.field);

    for i in [1..l] do
        for j in [i..l] do

            k := rep.setup.pairorbit[range[i], range[j]];

            if k > 0 then
                SetEntry(mat, i, j, rep.innerproducts[k]);
                SetEntry(mat, j, i, rep.innerproducts[k]);
            else
                SetEntry(mat, i, j, -rep.innerproducts[-k]);
                SetEntry(mat, j, i, -rep.innerproducts[-k]);
            fi;
        od;
    od;

    return mat;

    end );

InstallGlobalFunction(MAJORANA_SeparateInnerProduct,

    function(u,v,unknowns,rep)

    local   row,            # record values of unknowns
            sum,            # record values of knowns
            dim,            # size of coordinates
            i,              # index for dim of u
            j,              # index for dim of v
            m,              # orbit of i,j
            pos,            # position of m in unknowns
            sign;           # correct sign of 5A axes

    dim := Size(rep.setup.coords);

    sum := SparseZeroMatrix(1, 1, rep.field);
    row := SparseZeroMatrix(1, Size(unknowns), rep.field);

    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            m := rep.setup.pairorbit[u!.indices[1][i]][v!.indices[1][j]];

            if m > 0 then
                sign := 1;
            else
                sign := -1;
                m := -m;
            fi;

            if rep.innerproducts[m] <> false then
                AddToEntry(sum, 1, 1, - sign*u!.entries[1][i]*v!.entries[1][j]*rep.innerproducts[m]);
            else
                pos := Position(unknowns,m);
                AddToEntry(row, 1, pos, sign*u!.entries[1][i]*v!.entries[1][j]);
            fi;
        od;
    od;

    return [row,sum];

    end );

InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(rep)

    local   i,          # loop over representatives
            j,
            unknowns,
            ev,         # loop over eigenvalues
            mat,        # matrix of unknowns
            vec,        # vector of knowns
            u,          # vector with 1 in j th position
            v,          # eigenvector
            x,          # result of SeparateAlgebraProduct
            y,          # result of SolutionAlgProducts
            dim;        # size of setup.coords

    dim := Size(rep.setup.coords);

    unknowns := [];

    mat := SparseMatrix(0, Size(unknowns), [], [], rep.field);
    vec := SparseMatrix(0, dim, [], [], rep.field);

    Info( InfoMajorana, 50, "Building eigenvector unknowns");

    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], rep.field);

        for ev in [1..3] do
            for j in [1..Nrows(rep.evecs[i, ev])] do

                v := CertainRows(rep.evecs[i, ev], [j]);

                x := MAJORANA_SeparateAlgebraProduct(u, v, unknowns, rep);

                if x <> fail then

                    x[2] := x[2] + MAJORANA_FusionTable[1, ev + 1]*v;

                    if Size(x[1]!.indices[1]) = 1 then
                        y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns,
                                                        rep.algebraproducts,
                                                        rep.setup);

                        if not false in rep.algebraproducts then return true; fi;

                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;

                    elif x[1]!.indices[1] <> [] and not _IsRowOfSparseMatrix(mat, x[1]) then
                        mat := UnionOfRows(mat, x[1]);
                        vec := UnionOfRows(vec, x[2]);
                    fi;
                fi;
            od;
        od;
    od;

    y := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep.algebraproducts, rep.setup);

    return y;

    end);

InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,unknowns,rep)

    local   row,        # record values of unknowns
            sum,        # record values of knowns
            i,          # index for dim of u
            j,          # index for dim of v
            l,          # ordered version of [i,j]
            k,
            g,
            sign,
            elts,
            vecs,
            x,          # vector with 1 in the ith position
            y,
            pos,        # position of unknown product
            dim;        # dimension

    dim := Size(rep.setup.coords);

    row := SparseZeroMatrix(1, Size(unknowns), rep.field);
    sum := SparseZeroMatrix(1, dim, rep.field);

    elts := [];
    vecs := [];

    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            k := rep.setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

            if k > 0 then
                sign := 1;
            else
                sign := -1;
                k := -k;
            fi;

            x := rep.algebraproducts[k];

            if x = fail then return fail; fi;

            if x <> false then

                g := rep.setup.pairconj[u!.indices[1, i], v!.indices[1, j]];

                pos := Position(elts,g);

                if pos <> fail then
                    vecs[pos] := vecs[pos] - sign*u!.entries[1, i]*v!.entries[1, j]*x;
                else
                    Add(elts,g);
                    Add(vecs,- sign*u!.entries[1, i]*v!.entries[1, j]*x);
                fi;
            else

                l := [u!.indices[1, i], v!.indices[1, j]];
                Sort(l);

                pos := Position(unknowns,l);

                if pos = fail then
                    Add(unknowns, l); pos := Size(unknowns);
                fi;

                AddToEntry(row, 1, pos, u!.entries[1, i]*v!.entries[1, j]);
            fi;
        od;
    od;

    for i in [1..Size(elts)] do
        sum := sum + MAJORANA_ConjugateVec(vecs[i],rep.setup.pairconjelts[elts[i]]);
    od;

    return [row, RemoveMatWithHeads(sum, rep.setup.nullspace)];

    end);

InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns)

    local   output,     # output row
            len,        # length of row
            i,          # loop over length of row
            x,y,
            k,
            sign,       # corrects sign of 5A axis
            pos;        # position of new product

    if ForAll(g, i -> g[i] = i) then return row; fi;

    len     := Ncols(row);
    output  := SparseZeroMatrix(1, len, Rationals);

    for i in [1..Size(row!.indices[1])] do

        x := unknowns[row!.indices[1, i]];
        y := g{x};

        sign := 1;

        if y[1] < 0 then sign := -sign; y[1] := -y[1]; fi;
        if y[2] < 0 then sign := -sign; y[2] := -y[2]; fi;

        Sort(y);

        k := Position(unknowns,y);

        if k = fail then Add(unknowns, y); k := Size(unknowns); fi;

        pos := PositionSorted(output!.indices[1], k);

        Add(output!.indices[1], k, pos);
        Add(output!.entries[1], sign*row!.entries[1, i], pos);
    od;

    return output;

    end);

InstallGlobalFunction(MAJORANA_BasisOfEvecs,

    function(mat)

    local ech, dim;

    dim := Ncols(mat);

    ech := EchelonMatDestructive(CertainColumns(mat, [dim, dim - 1..1]));

    return CertainColumns(ech.vectors, [dim, dim - 1..1]);

    end);

InstallGlobalFunction( MAJORANA_AllConjugates,

    function(mat, vec, unknowns, rep)

    local x, i, new_mat, new_vec, y, nonzero, g, conj;

    Info(   InfoMajorana, 50, "All conjugates") ;

    x := EchelonMatTransformationDestructive(CertainColumns(mat, [Size(unknowns), Size(unknowns) - 1..1]));

    mat := CertainColumns(x.vectors, [Size(unknowns), Size(unknowns) - 1..1]);
    vec := x.coeffs*vec;

    new_mat := CopyMat(mat);
    new_vec := CopyMat(vec);

    for g in rep.setup.conjelts do
        for i in [1 .. Nrows(mat)] do
            if mat!.indices[i] <> [] then
                conj := [,];

                conj[1] := MAJORANA_ConjugateRow(CertainRows(mat, [i]), g, unknowns );
                conj[2] := MAJORANA_ConjugateVec(CertainRows(vec, [i]), g);

                new_mat := UnionOfRows(new_mat, conj[1]);
                new_vec := UnionOfRows(new_vec, conj[2]);
            fi;
        od;

        if Nrows(new_mat) > Ncols(new_mat) then

            x := MAJORANA_SolutionAlgProducts(new_mat, new_vec, unknowns, rep.algebraproducts, rep.setup);

            if not false in rep.algebraproducts then return true; fi;

            nonzero := Filtered([1..Nrows(x.mat)], j -> x.mat!.indices[j] <> []);

            new_mat := CertainRows(x.mat, nonzero);
            new_vec := CertainRows(x.vec, nonzero);

            y := MAJORANA_RemoveKnownAlgProducts(mat, vec, unknowns, rep.algebraproducts, rep.setup);

            mat := y.mat; vec := y.vec; unknowns := y.unknowns;

            if not false in rep.algebraproducts then return true; fi;

        fi;
    od;

    x := MAJORANA_SolutionAlgProducts(new_mat,new_vec,unknowns, rep.algebraproducts, rep.setup);

    return rec(mat := x.mat, vec := x.vec, unknowns := x.unknowns);

    end );

InstallGlobalFunction(MAJORANA_Resurrection,

    function(u, a, b, c, evals, unknowns, rep)

    local   x, y, z, ev, res,  n, i;

    res := MAJORANA_SeparateAlgebraProduct(b, c, unknowns, rep);

    if res[1]!.indices[1] = [] or res = fail then return false; fi;

    ev := MAJORANA_FusionTable[evals[1] + 1, evals[2] + 1];

    res := ev*res;

    x := MAJORANA_AlgebraProduct(c, a - b, rep.algebraproducts, rep.setup);

    if x = fail then return false; fi;

    z := MAJORANA_SeparateAlgebraProduct(u, x, unknowns, rep);

    if z = fail then return; fi;

    res := res + z;

    n := a - b;

    if evals[1] = 2 then
        i := u!.indices[1, 1];

        y := MAJORANA_InnerProduct(n - SparseMatrix(1, n!.ncols, [[i]], [[GetEntry(n, 1, i)]], Rationals), c, rep.innerproducts, rep.setup);

        if y <> false then
            res[2] := res[2] + (1/4)*y*u;
        else
            Error();
            return false;
        fi;
    fi;

    return res;

    end );

InstallGlobalFunction( MAJORANA_NullspaceUnknowns,

    function(mat, vec, unknowns, rep)

    local   i, j,
            u,
            v,
            x,
            y,
            dim;

    if Nrows(rep.setup.nullspace.vectors) = 0 then
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    Info( InfoMajorana, 50, "Building nullspace unknowns" );

    dim := Size(rep.setup.coords);

    x := MAJORANA_Orbits(rep.generators, dim, rep.setup);

    for i in x.orbitreps do
        u := SparseMatrix(1, dim, [[i]], [[1]], mat!.ring);

        for j in [1..Nrows(rep.setup.nullspace.vectors)] do

            v := CertainRows(rep.setup.nullspace.vectors, [j]);

            if  ForAny(rep.setup.pairorbit[i], k -> rep.algebraproducts[AbsInt(k)] = false) and
                ForAll(rep.setup.pairorbit[i], k -> rep.algebraproducts[AbsInt(k)] <> fail) then   # TODO think I can improve this, only an issue if a fail is actually hit during separation

                x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,rep);

                if x <> fail and Size(x[1]!.indices[1]) = 1 then

                    y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns,
                                                        rep.algebraproducts,
                                                        rep.setup);

                    if not false in rep.algebraproducts then return true; fi;

                    mat := y.mat; vec := y.vec; unknowns := y.unknowns;

                elif x <> fail and x[1]!.indices[1] <> [] then
                    if not _IsRowOfSparseMatrix(mat, x[1]) then
                        mat := UnionOfRows(mat, x[1]);
                        vec := UnionOfRows(vec, x[2]);
                    fi;
                fi;
            fi;
        od;

        if Nrows(mat) > 8000 then
            y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);

            if not false in rep.algebraproducts then return; fi;

            mat := y.mat; vec := y.vec; unknowns := y.unknowns;
        fi;
    od;

    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);

    return rec( mat := y.mat, vec := y.vec, unknowns := y.unknowns);

    end );

MAJORANA_SolutionMatVecs_Whatever := MAJORANA_SolutionMatVecs;
InstallGlobalFunction( MAJORANA_SolutionAlgProducts,

    function( mat, vec, unknowns, algebraproducts, setup)

    local   sol,        # solution of system
            sign,       # correct sign of 5A axes
            i,          # loop over <unknowns>
            x,
            nonzero;

    if ForAll(mat!.indices, x -> x = []) then
        return rec( mat := SparseMatrix(0, Ncols(mat), [], [], algebraproducts[1]!.ring),
                    vec := SparseMatrix(0, Ncols(vec), [], [], algebraproducts[1]!.ring),
                    unknowns := unknowns    );
    fi;

    mat!.ncols := Size(unknowns);

    Info(   InfoMajorana, 40,
            STRINGIFY("Solving a ", Nrows(mat), " x ", Ncols(mat), " matrix") );

    if algebraproducts[1]!.ring = Rationals then
        for i in [1..Nrows(mat)] do
            x := _FoldList2(mat!.entries[i], DenominatorRat, LcmInt);
            mat!.entries[i] := mat!.entries[i]*x;
            vec!.entries[i] := vec!.entries[i]*x;
        od;
    fi;

    sol := MAJORANA_SolutionMatVecs_Whatever(mat,vec);

    Info(   InfoMajorana, 40, "Solved it!" );

    if ForAll(sol.solutions, x -> x = fail) then
        return rec( mat := sol.mat, vec := sol.vec, unknowns := unknowns);
    fi;

    for i in [1..Size(unknowns)] do

        if sol.solutions[i] <> fail then

            MAJORANA_RecordSolution(    sol.solutions[i], unknowns[i],
                                        algebraproducts, setup );
        fi;
    od;

    Unbind(sol.solutions);

    x := MAJORANA_RemoveKnownAlgProducts(   sol.mat, sol.vec, unknowns,
                                            algebraproducts, setup    );

    nonzero := Filtered([1..Nrows(x.mat)], j -> x.mat!.indices[j] <> []);

    mat := CertainRows(x.mat, nonzero);
    vec := CertainRows(x.vec, nonzero);
    unknowns := x.unknowns;

    x := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, algebraproducts, setup);

    mat := x.mat;
    vec := x.vec;
    unknowns := x.unknowns;

    return rec( mat := mat, vec := vec, unknowns := unknowns );

    end );

InstallGlobalFunction( MAJORANA_SolveSingleSolution,

    function(x, mat, vec, unknowns, algebraproducts, setup)

    local   elm,
            y,
            switch,
            nonzero,
            i;

    Info( InfoMajorana, 60, "Solved a single solution");

    elm := x[1]!.entries[1][1];
    x := x*(1/elm);

    MAJORANA_RecordSolution(    x[2], unknowns[x[1]!.indices[1, 1]],
                                algebraproducts, setup );

    y := MAJORANA_RemoveKnownAlgProducts(   mat, vec, unknowns,
                                            algebraproducts, setup );

    if Nrows(y.mat) > 0 then
        nonzero := Filtered([1..Nrows(y.mat)], j -> y.mat!.indices[j] <> []);

        mat := CertainRows(y.mat, nonzero);
        vec := CertainRows(y.vec, nonzero);
        unknowns := y.unknowns;

        switch := true;

        while switch = true do

            if unknowns = [] then
                return rec( mat := mat, vec := vec, unknowns := unknowns );
            fi;

            switch := false;

            for i in [1..Nrows(mat)] do
                if Size(mat!.indices[i]) = 1 then
                    switch := true;
                    elm := mat!.entries[i, 1];
                    MAJORANA_RecordSolution(    CertainRows(vec, [i])*(1/elm),
                                                unknowns[mat!.indices[i, 1]],
                                                algebraproducts, setup);
                fi;;
            od;

            if switch = true then
                Info( InfoMajorana, 60, "Solved a new single solution");
            fi;

            x := MAJORANA_RemoveKnownAlgProducts(   mat, vec, unknowns,
                                                    algebraproducts, setup );

            nonzero := Filtered([1..Nrows(x.mat)], j -> x.mat!.indices[j] <> []);

            mat := CertainRows(x.mat, nonzero);
            vec := CertainRows(x.vec, nonzero);
            unknowns := x.unknowns;
        od;

    fi;

    x := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, algebraproducts, setup);

    return rec( mat := x.mat, vec := x.vec, unknowns := x.unknowns );

    end );

InstallGlobalFunction( MAJORANA_RecordSolution,

    function( v, x, algebraproducts, setup)

    local   y,
            g,
            sign;

    y := setup.pairorbit[x[1], x[2]];
    g := SP_Inverse(setup.pairconjelts[setup.pairconj[x[1], x[2]]]);

    sign := 1;

    if y < 0 then sign := -1; y := -y; fi;

    if algebraproducts[y] = false then
        algebraproducts[y] := sign*MAJORANA_ConjugateVec(v,g);
        algebraproducts[y] := RemoveMatWithHeads(algebraproducts[y], setup.nullspace);
    fi;

    end );

##
## Takes a system [mat, vec] of unknown algebra products and removes
## from the system any variables which have already been found
##

InstallGlobalFunction( MAJORANA_RemoveKnownAlgProducts,

    function( mat, vec, unknowns, algebraproducts, setup)

    local   unsolved,
            i,
            j,
            elm,
            x,
            y,
            sign,
            g,
            switch,
            pos,
            prod;

    if Nrows(mat) = 0 then

        unknowns := [];
        mat := SparseMatrix(0, Size(unknowns), [], [], mat!.ring);
        vec := SparseMatrix(0, Size(setup.coords), [], [], mat!.ring);

        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    unsolved := [];

    switch := false;

    for i in [1..Size(unknowns)] do

        x := unknowns[i];

        y := setup.pairorbit[x[1], x[2]];

        sign := 1;

        if y < 0 then sign := -1; y := -y; fi;

        prod := algebraproducts[y];

        if prod <> false then

            switch := true;

            g := setup.pairconjelts[setup.pairconj[x[1], x[2]]];

            prod := MAJORANA_ConjugateVec(prod,g);

            for j in [1..Nrows(vec)] do
                pos := Position(mat!.indices[j], i);
                if pos <> fail then
                    elm := mat!.entries[j, pos];
                    AddRow( prod!.indices[1],-sign*elm*prod!.entries[1],
                            vec!.indices, vec!.entries, j);
                fi;
            od;
        else
            Add(unsolved,i);
        fi;
    od;

    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};

    return rec( mat := mat, vec := vec, unknowns := unknowns);

    end );

##
## Takes a system of linear equations whose unknowns are inner product values and
## removes any unknowns from the system if there function has already been found
##

InstallGlobalFunction( MAJORANA_RemoveKnownInnProducts,

    function(mat, vec, unknowns, innerproducts)

    local   unsolved, i, j, elm, prod, nonzero;

    unsolved := [];

    if Nrows(mat) = 0 then
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    for i in [1..Size(unknowns)] do
        prod := innerproducts[unknowns[i]];

        if prod <> false then
            for j in [1..Nrows(vec)] do
                elm := GetEntry(mat, j, i);

                if elm <> 0 then
                    AddToEntry(vec, j, 1, -elm*prod);
                fi;
            od;
        else
            Add(unsolved, i);
        fi;
    od;

    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};

    nonzero := Filtered([1..Nrows(mat)], j -> mat!.indices[j] <> []);

    mat := CertainRows(mat, nonzero);
    vec := CertainRows(vec, nonzero);

    return rec( mat := mat, vec := vec, unknowns := unknowns);

    end );

##
## Takes a linear equation with just one unknown inner product value and
## records the new value.
##

InstallGlobalFunction( MAJORANA_SingleInnerSolution,

    function(eq, mat, vec, unknowns, innerproducts)

    local x;

    x := unknowns[eq[1]!.indices[1][1]];

    if eq[2]!.entries[1] = [] then
        innerproducts[x] := Zero(mat!.ring);
    else
        innerproducts[x] := eq[2]!.entries[1, 1]/eq[1]!.entries[1, 1];
    fi;

    return MAJORANA_RemoveKnownInnProducts(mat, vec, unknowns, innerproducts);

    end );

##
## Takes a system of linear equations where are the unknowns are inner products,
## solves the system and records any new values.
##

InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, unknowns, innerproducts)

    local   sol,    # solution of system
            i,      # loop over <unknowns>
            x,      # element of <unknowns>
            nonzero;

    if Nrows(mat) = 0 then
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    sol := MAJORANA_SolutionMatVecs(mat,vec);

    for i in [1..Size(sol.solutions)] do
        if sol.solutions[i] <> fail then
            x := unknowns[i];
            if sol.solutions[i]!.entries[1] = [] then
                innerproducts[x] := Zero(mat!.ring);
            else
                innerproducts[x] := sol.solutions[i]!.entries[1, 1];
            fi;
        fi;
    od;

    return MAJORANA_RemoveKnownInnProducts( sol.mat, sol.vec, unknowns, innerproducts);

    end );

##
## Calculate the nullspace from the Gram matrix of the form and removes the
## relevant vectors from algebra products and eigenvectors
##

InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(rep)

    local   dim, gram, null, unknowns, i, j, x;

    if fail in rep.innerproducts then return; fi;

    dim := Size(rep.setup.coords);

    gram := MAJORANA_FillGramMatrix([1..dim], rep.innerproducts, rep.setup);
    null := KernelEchelonMatDestructive(gram, [1..dim]).relations;;
    null := ReversedEchelonMatDestructive(null);

    rep.setup.nullspace := null;

    if null.heads = [] then return; fi;

    for i in [1..Size(rep.setup.pairreps)] do
        x := Filtered([1..dim], j -> i in rep.setup.pairorbit[j]);
        if ForAll(x, j -> rep.setup.nullspace.heads[j] <> 0) then
            rep.setup.pairreps[i] := fail;
            rep.algebraproducts[i] := fail;
            rep.innerproducts[i] := fail;
        fi;
    od;

    for i in [1..Size(rep.algebraproducts)] do
        if not rep.algebraproducts[i] in [false, fail] then
            rep.algebraproducts[i] := RemoveMatWithHeads(rep.algebraproducts[i], null);
        fi;
    od;

    for i in rep.setup.orbitreps do
        for j in [1..3] do
            rep.evecs[i, j] := RemoveMatWithHeads(rep.evecs[i, j], null);
            rep.evecs[i, j] := MAJORANA_BasisOfEvecs(rep.evecs[i, j]);
        od;
    od;

    end );

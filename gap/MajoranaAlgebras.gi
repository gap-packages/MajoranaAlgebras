#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

##
## The main Majorana representation function.
##

InstallGlobalFunction(MajoranaRepresentation,

function(arg)

    local   rep, unknowns, main;

    # The default axiomatic setting is to assume all axioms as given in Seress (2012)
    if Size(arg) = 2 then arg[3] := "AllAxioms"; fi;

    # Run the setup function
    rep :=  MAJORANA_SetUp(arg[1], arg[2],  arg[3]);

    # While there are still unknown algebra product loop over the main part of the algorithm
    while true do

        unknowns := Positions(rep.algebraproducts, false);

        main := MAJORANA_MainLoop(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );
        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.innerproducts, false)), " unknown inner products ") );

        if not false in rep.algebraproducts then
            # We have completely constructed the algebra
            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            # No more algebra products have been found - the algebra is incomplete
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

    MAJORANA_EigenvectorsAlgebraUnknowns(rep);

    MAJORANA_AxiomM1(rep);

    MAJORANA_Fusion(rep);

    return MAJORANA_UnknownAlgebraProducts(rep);

    end);

##
## Uses axiom M1 to find inner products
##

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(rep)

    local   dim, mat, vec, i, j, k, u, v, w, x, y, z, eq, unknowns;

    if not false in rep.innerproducts then return; fi;

    Info(   InfoMajorana, 50, "Axiom M1");

    dim := Size(rep.setup.coords);

    # Create an empty system of unknowns where the ideterminants are the indices of the unknown inner products
    unknowns := Positions(rep.innerproducts, false);
    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);

    # Loop over all known representative algebra products
    for i in Filtered([1..Size(rep.algebraproducts)], x -> not rep.algebraproducts[x] in [false, fail]) do
        x := rep.algebraproducts[i];
        # Do this both for the pair rep and also for the reverse of the pair rep
        for j in [rep.setup.pairreps[i], Reversed(rep.setup.pairreps[i])] do

            u := SparseMatrix(1, dim, [[ j[1] ]], [[ 1 ]], Rationals);
            v := SparseMatrix(1, dim, [[ j[2] ]], [[ 1 ]], Rationals);

            # Loop over all elements of the spanning set
            for k in Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0) do

                w := SparseMatrix(1, dim, [[ k ]], [[ 1 ]], Rationals);

                y := MAJORANA_AlgebraProduct(v, w, rep.algebraproducts, rep.setup);

                if not y in [fail, false] then

                    # Use axiom M1 (i.e. that (uv, w) = (u,vw)) to find a new equation
                    eq := MAJORANA_SeparateInnerProduct(w, x, unknowns, rep.innerproducts, rep.setup);
                    eq := eq - MAJORANA_SeparateInnerProduct(y, u, unknowns, rep.innerproducts, rep.setup);

                    if Size(eq[1]!.indices[1]) = 1 then
                        # If there is only one unknown value in the equation then immediately record this value
                        z := MAJORANA_SingleInnerSolution(  eq, mat, vec, unknowns, rep.innerproducts);

                        mat := z.mat; vec := z.vec; unknowns := z.unknowns;

                        # If all inner products have been found then calculate the nullspace of the alg
                        if unknowns = [] then
                            MAJORANA_CheckNullSpace(rep); return;
                        fi;
                    elif eq[1]!.indices[1] <> [] then
                        # Otherwise, add this equation to the system of unknowns
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

    # Solve this system of linear equations
    x := MAJORANA_SolutionInnerProducts(mat, vec, unknowns, rep.innerproducts);

    if x.unknowns = [] then MAJORANA_CheckNullSpace(rep); fi;

    return rec( mat := x.mat, vec := x.vec, unknowns := x.unknowns);

    end );

##
## Finds new eigenvectors using the fusion rules
##

InstallGlobalFunction( MAJORANA_Fusion,

function(arg)

    local   i, j, k, a, b, new, evals, rep, FUSE, u;

    rep := arg[1];
    FUSE := MAJORANA_FuseEigenvectors;

    # If specified by the optional second argument, do not assume that we have a Frobenius form
    if Size(arg) = 2 and arg[2] = false then FUSE := MAJORANA_FuseEigenvectorsNoForm; fi;

    # Loop over representatives of the orbits of G on T
    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

        # Continue to loop until no new eigenvectors have been found
        while true do
            Info(   InfoMajorana, 50, STRINGIFY("Fusion of ", i, " evecs")) ;

            # A list in which to store the new eigenvectors found from fusion
            new := List( [1, 2, 3], j -> CopyMat(rep.evecs[i, j]) );

            # For the following pairs of eigenvalues
            for evals in [[1,1], [1,2], [1,3], [2,3], [2,2], [3,3]] do
                for j in [1..Nrows(rep.evecs[i, evals[1]])] do
                    a := CertainRows(rep.evecs[i, evals[1]], [j]);

                    for k in [1..Nrows(rep.evecs[i, evals[2]])] do
                        b := CertainRows(rep.evecs[i, evals[2]], [k]);

                        FUSE(a, b, u, evals, new, rep);
                    od;
                od;
            od;

            # Find a basis of the new eigenspaces
            new := List( new, MAJORANA_BasisOfEvecs );

            # If no new eigenvectors have been found then break
            if ForAll([1..3], j -> Nrows(new[j]) = Nrows(rep.evecs[i, j])) then
                break;
            fi;

            rep.evecs[i] := new;

            # Calculate the intersection of the eigenspaces to find new nullspace vectors
            # Do not understand why we have to do this twice!
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

    # If specified by the optional second argument, do not assume that we have a Frobenius form
    if Size(arg) = 2 and arg[2] = false then evals_list := [[1,2], [1,3]]; fi;

    dim := Size(rep.setup.coords);
    n := Size(Positions(rep.algebraproducts, false));

    # Use the known eigenvectors to calculate more algebra products
    x := MAJORANA_EigenvectorsAlgebraUnknowns(rep);
    if not false in rep.algebraproducts then return true; fi;
    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    # Use the known nullspace vectors to calculate more algebra products
    x := MAJORANA_NullspaceUnknowns(mat, vec, unknowns, rep.algebraproducts, rep.setup, rep.group);
    if not false in rep.algebraproducts then return true; fi;
    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    # Use the known eigenvectors and the resurrection principle to calculate more algebra products
    Info(   InfoMajorana, 50, "Building resurrection");

    for evals in evals_list do
        for i in rep.setup.orbitreps do

            evecs_a := rep.evecs[i, evals[1]];
            evecs_b := rep.evecs[i, evals[2]];

            u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

            # For the each eigenvector a in <evecs_a> create a list of indices of the eigenvectors b in
            # <evecs_b> such that the product ab involves only one unknown product. We will calculate
            # with these eigenvectors first.

            list := [,];
            list[1] := List([1..Nrows(evecs_a)], i -> []);
            list[2] := List([1..Nrows(evecs_a)], i -> []);

            for j in [1..Nrows(evecs_a)] do
                bad := MAJORANA_FindBadIndices(CertainRows(evecs_a,[j]), rep);
                for k in [1..Nrows(evecs_b)] do
                    x := Size(Intersection(bad, evecs_b!.indices[k]));
                    if x = 1 then
                        Add(list[1, j], k);
                    elif x > 1 then
                        Add(list[2, j], k);
                    fi;
                od;
            od;

            # Loop over eigevectors a, b and c where a and c are <evals[1]>-eigenvectors
            # and b is an <evals[2]>-eigenvectors.
            for index in [1,2] do
                for j in [1..Nrows(evecs_a)] do
                    c := CertainRows(evecs_a, [j]);
                    for k in list[index, j] do
                        b := CertainRows(evecs_b, [k]);

                        bad := MAJORANA_FindBadIndices(c, rep);

                        for l in [1..Nrows(evecs_a)] do

                            a := CertainRows(evecs_a, [l]);

                            # If this is satisfied then the product (a - b)*c is known
                            if CertainColumns(a, bad) = CertainColumns(b, bad) then

                                # Perform resurrection on these eigenvectors
                                x := MAJORANA_Resurrection(  u, a, b, c, evals, unknowns,rep);

                                if x <> false and x[1]!.indices[1] <> [] then
                                    # If this equation has only one unknown then record this value immediately
                                    if Size(x[1]!.indices[1]) = 1 then
                                        y := MAJORANA_SolveSingleSolution( x, mat, vec, unknowns, rep);

                                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;

                                        if not false in rep.algebraproducts then return true; fi;
                                    elif not _IsRowOfSparseMatrix(mat, x[1]) then
                                    # Otherwise, add the equation to the system of linear equations

                                        mat := UnionOfRows(mat, x[1]);
                                        vec := UnionOfRows(vec, x[2]);
                                    fi;
                                fi;
                            fi;
                        od;

                        # If the matrix is over a certain size then for performance reasons, solve it already
                        if Nrows(mat) > Ncols(mat) or Nrows(mat) > 8000 then
                            x := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

                            mat := x.mat; vec := x.vec; unknowns := x.unknowns;

                            if not false in rep.algebraproducts then return true; fi;
                        fi;

                    od;
                od;
            od;
        od;
    od;

    x := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;

    # If no new algebra products have been found then also add the images of the currect equations under the group action
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

    function(v, rep)

    local i, j, k, list, bad;

    bad := [];
    list := [1 .. Size(rep.setup.coords)];

    # Loop over the non-zero indices of v
    for i in v!.indices[1] do
        # Loop over all vectors in coords
        for j in list do
            k := MAJORANA_UnorderedOrbitalRep(rep.setup.orbitalstruct, [i, j]);
            k := rep.setup.pairrepsmap[ k ];

            if k < 0 then k := -k; fi;

            if rep.algebraproducts[k] = false then
                # Add this index to the list of unknown indices and remove it
                # from <list> as we don't need to check it again.
                Add(bad,j);
                list := Difference(list, [j]);
            fi;
        od;
    od;

    return SortedList(bad);

    end );

##
## Add a new eigenvector to a list of existing eigenvectors if not already a member
##

InstallGlobalFunction(MAJORANA_AddEvec,

    function(mat, x)

    # Don't add a zero row
    if x!.indices[1] = [] then return mat; fi;

    # Scale the vector so that the last entry is equal to one
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

    function(a, b, u, evals, new, rep)

    local   new_ev, pos, x, y, z;

    # Find the eigenvalue of the new vector
    new_ev := MAJORANA_FusionTable[evals[1] + 1, evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;

    # Find the product of the eigenvectors a and b
    x := MAJORANA_AlgebraProduct(a,b,rep.algebraproducts,rep.setup);

    if x in [false, fail] then return; fi;

    # If the eigenvalues are both equal to 1/4 or to 1/32 then we need to
    # do more work to recover the eigenvectors
    if evals = [2,2] then
        y := MAJORANA_InnerProduct(a,b,rep.innerproducts,rep.setup);

        if y = false then return; fi;

        new[1] := MAJORANA_AddEvec(new[1], x - (1/4)*u*y);
    elif evals = [3,3] then
        y := MAJORANA_InnerProduct(a,b,rep.innerproducts,rep.setup);

        if y = false then return; fi;

        z := MAJORANA_AlgebraProduct(u,x,rep.algebraproducts, rep.setup);

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

    # If there is no inner product then we cannot tell if we have a full eigenspace decomposition
    if rep.innerproducts = false then return false; fi;

    if Sum(List(evecs, Nrows)) + Nrows(rep.setup.nullspace.vectors) < dim - 1 then
        return false;
    fi;

    return true;

    end );

##
## Takes a row vector <vec> and finds its image under the signed perm <g>
##

InstallGlobalFunction( MAJORANA_ConjugateVec,

    function(vec,g)

    local   i,
            k,
            sign,
            res,
            pos;

    # If g is the identity on vec then return
    if ForAll(vec!.indices[1], i -> g[i] = i) then return vec; fi;

    res := SparseMatrix(1, Ncols(vec), [[]], [[]], Rationals);

    # Loop over the non-zero indices of vec and add their image to res
    for i in [1..Size(vec!.indices[1])] do

        k := g[vec!.indices[1, i]];

        sign := 1;

        if k < 0 then k := -k; sign := -1; fi;

        pos := PositionSorted(res!.indices[1], k);

        Add(res!.indices[1], k, pos);
        Add(res!.entries[1], sign*vec!.entries[1, i], pos);
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
                pos;

        vec := SparseMatrix(1, Ncols(u), [[]], [[]],  Rationals);

        # <elts> will record the permutations that we have to conjugate by
        # <vecs> will record the corresponding vector that we are permuting
        elts := [];
        vecs := [];

        # Loop over the non-zero coefficients of u and v
        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do

                # Find the representative of the orbital containing (i,j)
                k := MAJORANA_UnorderedOrbitalRep(setup.orbitalstruct, [u!.indices[1, i], v!.indices[1, j]]);
                k := setup.pairrepsmap[ k ];

                # Adjust the sign
                if k > 0 then
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;

                x := algebraproducts[k];

                if not x in [fail, false] then

                    g := MAJORANA_OrbitalCanonizingElementInverse(setup.orbitalstruct,  [u!.indices[1, i], v!.indices[1, j]]);
                    g := ListSignedPerm(g, Ncols(u));

                    pos := Position(elts,g);

                    # If we have already seen this elt then add the product <x> to the corresponding
                    # vector in vecs, else add the elt and the product to these lists.
                    if pos <> fail then
                        AddRow(x!.indices[1],  sign*u!.entries[1, i]*v!.entries[1, j]*x!.entries[1],  vecs[pos]!.indices, vecs[pos]!.entries, 1);
                    else
                        Add(elts, g);
                        Add(vecs, CopyMat(sign*u!.entries[1, i]*v!.entries[1, j]*x));
                    fi;
                elif x = false then # cannot calculate product
                    return false;
                else # product with a vector in the nullspace
                    return fail;
                fi;
            od;
        od;

        # Now go over all group elements are permute their corresponding vectors, add the
        # result to the output vec as we go along

        for g in elts do
            x := MAJORANA_ConjugateVec(vecs[i], g);
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

        sum := 0;

        # Loop over the non-zero coefficients of u and v
        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do
                # Find the representative of the orbital containing the pair (i,j)
                k := MAJORANA_UnorderedOrbitalRep( setup.orbitalstruct, [u!.indices[1, i], v!.indices[1, j]] );
                k := setup.pairrepsmap[ k ];

                # Adjust the sign
                if k > 0 then
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;

                if not innerproducts[k] in [fail, false] then
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

function(range, innerproducts, setup)

    local   i, j, k, mat, l;

    l := Length(range);

    mat := SparseZeroMatrix(l, l, Rationals);

    for i in [1..l] do
        for j in [i..l] do

            # Find the representative of the orbital containing the pair <range{[i,j]}>
            k := MAJORANA_UnorderedOrbitalRep( setup.orbitalstruct, range{[i,j]} );
            k := setup.pairrepsmap[k];

            # Adjust for the sign
            if k > 0 then
                SetEntry(mat, i, j, innerproducts[k]);
                SetEntry(mat, j, i, innerproducts[k]);
            else
                SetEntry(mat, i, j, -innerproducts[-k]);
                SetEntry(mat, j, i, -innerproducts[-k]);
            fi;
        od;
    od;

    return mat;

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
            y;          # result of SolutionAlgProducts

    unknowns := [];

    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, Size(rep.setup.coords), [], [], Rationals);

    Info( InfoMajorana, 50, "Building eigenvector unknowns");

    # Loop over the representatives of the action of G on T
    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

        # Loop over the three eigenvalues 0, 1/4 and 1/32
        for ev in [1..3] do
            for j in [1..Nrows(rep.evecs[i, ev])] do

                v := CertainRows(rep.evecs[i, ev], [j]);

                # Create the equation u*evec = 0
                x := MAJORANA_SeparateAlgebraProduct(u, v, unknowns, rep.algebraproducts, rep.setup);

                if x <> fail then
                    # Adjust to make the equation u*evec = ev*evec
                    x[2] := x[2] + MAJORANA_FusionTable[1, ev + 1]*v;

                    # If this equation has only one unknown, record its value immediately
                    if Size(x[1]!.indices[1]) = 1 then
                        y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns, rep);

                        if not false in rep.algebraproducts then return true; fi;

                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;

                    elif x[1]!.indices[1] <> [] and not _IsRowOfSparseMatrix(mat, x[1]) then
                    # Otherwise add the equation to the linear system
                        mat := UnionOfRows(mat, x[1]);
                        vec := UnionOfRows(vec, x[2]);
                    fi;
                fi;
            od;
        od;
    od;

    y := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

    return y;

    end);

InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,unknowns,algebraproducts,setup)

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
            pos;       # position of unknown product

    row := SparseZeroMatrix(1, Size(unknowns), Rationals);
    sum := SparseZeroMatrix(1, Size(setup.coords), Rationals);

    # <elts> will record the permutations that we have to conjugate by
    # <vecs> will record the corresponding vector that we are permuting
    elts := [];
    vecs := [];

    # Loop over the non-zero coefficients of u and v
    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            # Find the representative of the orbital containing the pair (i,j)
            k := MAJORANA_UnorderedOrbitalRep( setup.orbitalstruct, [u!.indices[1, i], v!.indices[1, j]] );
            k := setup.pairrepsmap[k];

            # Adjust the sign
            if k > 0 then
                sign := 1;
            else
                sign := -1;
                k := -k;
            fi;

            x := algebraproducts[k];

            if x = fail then return fail; fi;

            if x <> false then

                g := MAJORANA_OrbitalCanonizingElementInverse(setup.orbitalstruct,  [u!.indices[1, i], v!.indices[1, j]]);
                g := ListSignedPerm(g, Size(setup.coords));

                pos := Position(elts,g);

                if pos <> fail then
                    vecs[pos] := vecs[pos] - sign*u!.entries[1, i]*v!.entries[1, j]*x;
                else
                    Add(elts,g);
                    Add(vecs,- sign*u!.entries[1, i]*v!.entries[1, j]*x);
                fi;
            else
                # Otherwise, record the coefficient of the unknown pair in the matrix
                # of the system of linear equations
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

    for g in elts do
        sum := sum + MAJORANA_ConjugateVec(vecs[i], g);
    od;

    return [row, RemoveMatWithHeads(sum, setup.nullspace)];

    end);

InstallGlobalFunction(MAJORANA_SeparateInnerProduct,

    function(u,v,unknowns,innerproducts,setup)

    local   row,            # record values of unknowns
            sum,            # record values of knowns
            dim,            # size of coordinates
            i,              # index for dim of u
            j,              # index for dim of v
            k,              # orbit of i,j
            pos,            # position of k in unknowns
            sign;           # correct sign of 5A axes

    dim := Size(setup.coords);

    sum := SparseZeroMatrix(1, 1, Rationals);
    row := SparseZeroMatrix(1, Size(unknowns), Rationals);

    # Loop over the nonzero coefficients of u and v
    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            k := MAJORANA_UnorderedOrbitalRep( setup.orbitalstruct, [u!.indices[1, i], v!.indices[1, j]] );

            if setup.pairrepsmap[k] = fail then Error(); fi;

            k := setup.pairrepsmap[k];

            # Adjust the sign
            if k > 0 then
                sign := 1;
            else
                sign := -1;
                k := -k;
            fi;

            # If the product is known the calculate as usual. Otherwise, add the coefficient
            # to the matrix of the system of linear equations.
            if innerproducts[k] <> false then
                AddToEntry(sum, 1, 1, - sign*u!.entries[1, i]*v!.entries[1, j]*innerproducts[k]);
            else
                pos := Position(unknowns,k);
                AddToEntry(row, 1, pos, sign*u!.entries[1, i]*v!.entries[1, j]);
            fi;
        od;
    od;

    return [row,sum];

    end );

InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns)

    local   output,     # output row
            len,        # length of row
            i,          # loop over length of row
            y,
            k,
            sign,       # corrects sign of 5A axis
            pos;        # position of new product

    # If the elt is the identity then return
    if ForAll(g, i -> g[i] = i) then return row; fi;

    len     := Ncols(row);
    output  := SparseZeroMatrix(1, len, Rationals);

    # Loop over the non-zero coefficients of row
    for i in [1..Size(row!.indices[1])] do

        # Find the image of the corresponding element of unknowns
        y := g{ unknowns[row!.indices[1, i]] };

        # Adjust the sign
        sign := 1;

        if y[1] < 0 then sign := -sign; y[1] := -y[1]; fi;
        if y[2] < 0 then sign := -sign; y[2] := -y[2]; fi;

        Sort(y);

        # find the position of the image in unknowns
        k := Position(unknowns,y);
        if k = fail then Add(unknowns, y); k := Size(unknowns); fi;

        # Add the value to the ouput row
        pos := PositionSorted(output!.indices[1], k);

        Add(output!.indices[1], k, pos);
        Add(output!.entries[1], sign*row!.entries[1, i], pos);
    od;

    return output;

    end);

#
# Calculates the reversed echelon form of a matrix
#

InstallGlobalFunction(MAJORANA_BasisOfEvecs,

    function(mat)

    local ech, dim;

    dim := Ncols(mat);

    ech := EchelonMatDestructive(CertainColumns(mat, [dim, dim - 1..1]));

    return CertainColumns(ech.vectors, [dim, dim - 1..1]);

    end);

##
## Takes a system of linear equations and calculates the image of each row under
## the action of each of the elements in conjelts.
##

InstallGlobalFunction( MAJORANA_AllConjugates,

    function(mat, vec, unknowns, rep)

    local x, i, new_mat, new_vec, y, nonzero, g, conj;

    Info(   InfoMajorana, 50, "All conjugates") ;

    # TODO Why do we do this?
    x := EchelonMatTransformationDestructive(CertainColumns(mat, [Size(unknowns), Size(unknowns) - 1..1]));

    mat := CertainColumns(x.vectors, [Size(unknowns), Size(unknowns) - 1..1]);
    vec := x.coeffs*vec;

    # These matrices will be what we add the conjugated rows to
    new_mat := SparseMatrix( 0, Ncols(mat), [], [], Rationals );
    new_vec := SparseMatrix( 0, Ncols(vec), [], [], Rationals );

    # Loop over group elements and matrix rows
    for g in rep.setup.conjelts do
        for i in [1 .. Nrows(mat)] do
            if mat!.indices[i] <> [] then
                conj := [,];

                # Calculate the images under g
                conj[1] := MAJORANA_ConjugateRow(CertainRows(mat, [i]), g, unknowns );
                conj[2] := MAJORANA_ConjugateVec(CertainRows(vec, [i]), g);

                new_mat := UnionOfRows(new_mat, conj[1]);
                new_vec := UnionOfRows(new_vec, conj[2]);
            fi;
        od;

        # If the matrix is sufficiently big then for performance reasons, solve already

        if Nrows(new_mat) > Ncols(new_mat) then

            x := MAJORANA_SolutionAlgProducts( new_mat, new_vec, unknowns, rep );

            if not false in rep.algebraproducts then return true; fi;

            # Take only the nonzero rows of new_mat and new_vec
            nonzero := Filtered([1..Nrows(x.mat)], j -> x.mat!.indices[j] <> []);

            new_mat := CertainRows(x.mat, nonzero);
            new_vec := CertainRows(x.vec, nonzero);

            # And also remove the known products from the original system
            y := MAJORANA_RemoveKnownAlgProducts(mat, vec, unknowns, rep);

            mat := y.mat; vec := y.vec; unknowns := y.unknowns;

            if not false in rep.algebraproducts then return true; fi;
        fi;
    od;

    x := MAJORANA_SolutionAlgProducts(new_mat, new_vec, unknowns, rep);

    return rec(mat := x.mat, vec := x.vec, unknowns := x.unknowns);

    end );

##
## Performs the resurrection priciple on the eigenvectors a, b and c where
## a, c are <evals[1]>-eigenvectors and b in a <evals[2]>-eigenvector.
##

InstallGlobalFunction(MAJORANA_Resurrection,

    function(u, a, b, c, evals, unknowns, rep)

    local   x, y, z, ev, res, diff, i, p;

    ev := MAJORANA_FusionTable[evals[1] + 1, evals[2] + 1];

    # Calculate the product b*c
    res := ev*MAJORANA_SeparateAlgebraProduct(b, c, unknowns, rep.algebraproducts, rep.setup);

    # If this product has known unknowns then it is not useful to us
    if res[1]!.indices[1] = [] or res = fail then return false; fi;

    # Calculate the product (a - b)*c
    x := MAJORANA_AlgebraProduct(c, a - b, rep.algebraproducts, rep.setup);
    if x = fail then return false; fi;

    # Calculate the product u*((a - b)*c)
    z := MAJORANA_SeparateAlgebraProduct(u, x, unknowns, rep.algebraproducts, rep.setup);
    if z = fail then return; fi;

    res := res + z;
    diff := a - b;

    # If a and c are 1/4-eigenvectors then we need to calculate the projection
    # of (a - b)*c onto the 1-eigenspace of u.
    if evals[1] = 2 then
        i := u!.indices[1, 1];

        p := SparseMatrix(1, diff!.ncols, [[i]], [[GetEntry(diff, 1, i)]], Rationals);
        y := MAJORANA_InnerProduct(diff - p , c, rep.innerproducts, rep.setup);

        if y <> false then
            res[2] := res[2] + (1/4)*y*u;
        else
            return false;
        fi;
    fi;

    return res;

    end );

InstallGlobalFunction( MAJORANA_NullspaceUnknowns,

    function(mat, vec, unknowns, rep)

    local   i, j, gens, u, v, x, y, dim;

    if Nrows(rep.setup.nullspace.vectors) = 0 then
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    Info( InfoMajorana, 50, "Building nullspace unknowns" );

    dim := Size(rep.setup.coords);

    # Calculate the orbits of G on the spanning set coords
    gens := GeneratorsOfGroup(rep.group);
    gens := List( gens, g -> MAJORANA_FindPerm(g, rep, rep) );
    x := MAJORANA_Orbits(gens, dim, rep.setup);

    # Loop over the representatives of these orbits
    for i in x.orbitreps do
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

        # Loop over the rows of the nullspace
        for j in [1..Nrows(rep.setup.nullspace.vectors)] do
            v := CertainRows(rep.setup.nullspace.vectors, [j]);

            # Calculate the equation u*v = 0
            x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,rep.algebraproducts,rep.setup);

            # If the equation has only one unknown then immediately record this value
            if x <> fail and Size(x[1]!.indices[1]) = 1 then

                y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns,rep);

                if not false in rep.algebraproducts then return true; fi;

                mat := y.mat; vec := y.vec; unknowns := y.unknowns;

            elif x <> fail and x[1]!.indices[1] <> [] then
                # Otherwise add this equation to the system of linear equations
                if not _IsRowOfSparseMatrix(mat, x[1]) then
                    mat := UnionOfRows(mat, x[1]);
                    vec := UnionOfRows(vec, x[2]);
                fi;
            fi;
        od;

        # If the matrix is too big then for performance reasons, solve already
        if Nrows(mat) > 8000 then
            y := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

            if not false in rep.algebraproducts then return; fi;

            mat := y.mat; vec := y.vec; unknowns := y.unknowns;
        fi;
    od;

    y := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

    return rec( mat := y.mat, vec := y.vec, unknowns := y.unknowns);

    end );

MAJORANA_SolutionMatVecs_Whatever := MAJORANA_SolutionMatVecs;
InstallGlobalFunction( MAJORANA_SolutionAlgProducts,

    function( mat, vec, unknowns, rep)

    local   sol,        # solution of system
            sign,       # correct sign of 5A axes
            i,          # loop over <unknowns>
            x,
            nonzero;

    # If the matrix is zero then return
    if ForAll(mat!.indices, x -> x = []) then
        return rec( mat := SparseMatrix(0, Ncols(mat), [], [], Rationals),
                    vec := SparseMatrix(0, Ncols(vec), [], [], Rationals),
                    unknowns := unknowns    );
    fi;

    mat!.ncols := Size(unknowns);

    Info(   InfoMajorana, 40,
            STRINGIFY("Solving a ", Nrows(mat), " x ", Ncols(mat), " matrix") );

    # Turn the matrix into an integer matrix
    for i in [1..Nrows(mat)] do
        x := _FoldList2(mat!.entries[i], DenominatorRat, LcmInt);
        mat!.entries[i] := mat!.entries[i]*x;
        vec!.entries[i] := vec!.entries[i]*x;
    od;

    # Solve the system of linear equations
    sol := MAJORANA_SolutionMatVecs_Whatever(mat,vec);

    Info(   InfoMajorana, 40, "Solved it!" );

    # If no new solutions have been found then return
    if ForAll(sol.solutions, x -> x = fail) then
        return rec( mat := sol.mat, vec := sol.vec, unknowns := unknowns);
    fi;

    # Otherwise, record the new solutions
    for i in [1..Size(unknowns)] do
        if sol.solutions[i] <> fail then
            MAJORANA_RecordSolution( sol.solutions[i], unknowns[i], rep );
        fi;
    od;

    Unbind(sol.solutions);

    # Adjust the system of linear equations to take into account the new known products
    x := MAJORANA_RemoveKnownAlgProducts( sol.mat, sol.vec, unknowns, rep );

    # Take out any zero rows
    nonzero := Filtered([1..Nrows(x.mat)], j -> x.mat!.indices[j] <> []);

    mat := CertainRows(x.mat, nonzero);
    vec := CertainRows(x.vec, nonzero);
    unknowns := x.unknowns;

    return MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

    end );

InstallGlobalFunction( MAJORANA_SolveSingleSolution,

    function(eqn, mat, vec, unknowns, rep)

    local   elm, y, switch, nonzero, i, ind;

    Info( InfoMajorana, 60, "Solved a single solution");

    # Divide through by the coefficient of the indeterminant
    elm := eqn[1]!.entries[1, 1];
    eqn := eqn/elm;

    # Record the new algebra product
    MAJORANA_RecordSolution( eqn[2], unknowns[eqn[1]!.indices[1, 1]], rep );

    # Reduce the system of linear equations using this new product
    y := MAJORANA_RemoveKnownAlgProducts(   mat, vec, unknowns, rep );

    # Remove any nonzero rows
    nonzero := Filtered([1..Nrows(y.mat)], j -> y.mat!.indices[j] <> []);

    mat := CertainRows(y.mat, nonzero);
    vec := CertainRows(y.vec, nonzero);
    unknowns := y.unknowns;

    # If there still remain non-zero rows of the matrix then look for more single solutions
    if Nrows(y.mat) > 0 then
        switch := true;

        # While we continue to find new products
        while switch = true do

            if unknowns = [] then
                return rec( mat := mat, vec := vec, unknowns := unknowns );
            fi;

            switch := false;

            for i in [1..Nrows(mat)] do
                if Size(mat!.indices[i]) = 1 then
                    # We have found a new single solution
                    switch := true;
                    elm := mat!.entries[i, 1];
                    ind := mat!.indices[i, 1];
                    MAJORANA_RecordSolution( CertainRows(vec, [i])*(1/elm), unknowns[ind], rep);
                fi;;
            od;

            if switch = true then
                Info( InfoMajorana, 60, "Solved a new single solution");
            fi;

            # Reduce the system of linear equations using this new product
            y := MAJORANA_RemoveKnownAlgProducts( mat, vec, unknowns, rep );

            # Remove any nonzero rows
            nonzero := Filtered([1..Nrows(y.mat)], j -> y.mat!.indices[j] <> []);

            mat := CertainRows(y.mat, nonzero);
            vec := CertainRows(y.vec, nonzero);
            unknowns := y.unknowns;
        od;
    fi;

    y := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, rep);

    return rec( mat := y.mat, vec := y.vec, unknowns := y.unknowns );

    end );

##
## Takes a vector that is the product of a pair of vectors given by x and records
## it in the correct place in the list rep.algebraproducts.
##

InstallGlobalFunction( MAJORANA_RecordSolution,

    function( vec, x, rep)

    local   y, g, sign;

    y := rep.setup.pairorbit[x[1], x[2]];
    g := MAJORANA_OrbitalCanonizingElement( rep.setup.orbitalstruct, x);
    g := ListSignedPerm(g, Size(rep.setup.coords));
    y := rep.setup.pairrepsmap[g{x}];

    # Adjust the sign
    sign := 1;
    if y < 0 then sign := -1; y := -y; fi;

    # Record the new product
    if rep.algebraproducts[y] = false then
        rep.algebraproducts[y] := sign*MAJORANA_ConjugateVec(vec,g);
        rep.algebraproducts[y] := RemoveMatWithHeads(rep.algebraproducts[y], rep.setup.nullspace);
    fi;

    end );

##
## Takes a system [mat, vec] of unknown algebra products and removes
## from the system any variables which have already been found
##

InstallGlobalFunction( MAJORANA_RemoveKnownAlgProducts,

    function( mat, vec, unknowns, rep)

    local   unsolved,
            i,
            j,
            elm,
            x,
            y,
            sign,
            g,
            pos,
            prod;

    if Nrows(mat) = 0 then
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    unsolved := [];

    # Loop over the unknown algebra products
    for i in [1..Size(unknowns)] do

        # Find the representative of the orbital containing the unknown value
        x := unknowns[i];
        y := MAJORANA_UnorderedOrbitalRep(rep.setup.orbitalstruct, x);

        # Adjust the sign
        sign := 1;
        if y < 0 then sign := -1; y := -y; fi;

        prod := rep.algebraproducts[y];

        # If the product is now known the remove its value from the rhs
        if prod <> false then

            g := MAJORANA_OrbitalCanonizingElementInverse(rep.setup.orbitalstruct, x);
            g := ListSignedPerm(g, Size(rep.setup.coords));

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
            # Otherwise, we want to keep this column of the matrix
            Add(unsolved,i);
        fi;
    od;

    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};

    # We don't want to remove the non-zero rows because we want to preserve the
    # row indexing when we use this inside MAJORANA_AllConjugates.

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

    # Loop over the unknown inner product values
    for i in [1..Size(unknowns)] do
        prod := innerproducts[unknowns[i]];

        # If the product is now known then loop over the rows and remove this value from the rhs
        if prod <> false then
            for j in [1..Nrows(vec)] do
                elm := GetEntry(mat, j, i);

                if elm <> 0 then
                    AddToEntry(vec, j, 1, -elm*prod);
                fi;
            od;
        else
            # Otherwise, we want to keep this column of the matrix
            Add(unsolved, i);
        fi;
    od;

    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};

    # Remove any zero rows
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

    x := unknowns[eq[1]!.indices[1, 1]];

    if eq[2]!.entries[1] = [] then
        innerproducts[x] := 0;
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

    # Solve the system of linear equations
    sol := MAJORANA_SolutionMatVecs(mat,vec);

    # Record any new solutions that have been found
    for i in [1..Size(sol.solutions)] do
        if sol.solutions[i] <> fail then
            x := unknowns[i];
            if sol.solutions[i]!.entries[1] = [] then
                innerproducts[x] := 0;
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

    # Calculate the gram matrix wrt to the spanning set coords and find its kernel
    gram := MAJORANA_FillGramMatrix([1..dim], rep.innerproducts, rep.setup);
    null := KernelEchelonMatDestructive(gram, [1..dim]).relations;;
    null := ReversedEchelonMatDestructive(null);

    rep.setup.nullspace := null;

    # If the kernel is trivial the return
    if null.heads = [] then return; fi;

    # Find the orbits that only involve nullspace vectors and set these to fail
    # because we don't need to find their values
    for i in [1..Size(rep.setup.pairreps)] do
        # TODO wow this is ugly
        x := Filtered([1..dim], j -> ForAny( [1..dim], k -> rep.setup.pairrepsmap[ MAJORANA_UnorderedOrbitalRep( [j,k] ) ] = i ) );

        if ForAll(x, j -> rep.setup.nullspace.heads[j] <> 0) then
            rep.setup.pairreps[i] := fail;
            rep.algebraproducts[i] := fail;
            rep.innerproducts[i] := fail;
        fi;
    od;

    # Quotient out the algebra products and evecs by the nullspace vectors
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

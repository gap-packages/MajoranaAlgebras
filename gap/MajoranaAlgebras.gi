#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

################################################################################
##
## The main Majorana representation function.
##
################################################################################

InstallGlobalFunction(MajoranaRepresentation,

function(arg)

    local   input, index, options, rep, unknowns, main;

    input := arg[1]; index := arg[2];

    if Size(arg) = 2 then
        options := rec( );
    else
        options := arg[3];
    fi;

    # Run the setup function
    rep :=  MAJORANA_SetUp(input, index, options);

    # Find representations of the maximal subgroups of G and embed them
    if IsBound(options.embeddings) and options.embeddings = true and Size(rep.group) > 120 then
        MAJORANA_MaximalSubgps(rep, options);
    fi;

    # While there are still unknown algebra product loop over the main part of the algorithm
    while true do

        unknowns := Positions(rep.algebraproducts, false);

        main := MAJORANA_MainLoop(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );
        if IsBound(rep.innerproducts) then
            Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.innerproducts, false)), " unknown inner products ") );
        fi;

        if not false in rep.algebraproducts then
            # We have completely constructed the algebra
            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            # No more algebra products have been found - the algebra is incomplete
            Info( InfoMajorana, 10, "Fail" );
            rep.system := main;
            return rep;
        fi;
    od;

    end );

################################################################################
##
## The main loop functions
##
################################################################################

##
## The main part of the algorithm that is looped over in the main function.
##

InstallGlobalFunction(MAJORANA_MainLoop,

    function(rep)

    MAJORANA_AxiomM1(rep);

    MAJORANA_Fusion(rep);

    return MAJORANA_UnknownAlgebraProducts(rep);

    end);

##
## Uses axiom M1 to find inner products
##

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(rep)

    local   dim, system, i, j, k, u, v, w, x, y, z, eq;

    # If the algebra does not have a form
    if not IsBound(rep.innerproducts) then return; fi;

    # Or if all inner product values have been found
    if not false in rep.innerproducts then return; fi;

    Info(   InfoMajorana, 50, "Axiom M1");

    dim := Size(rep.setup.coords);

    # Create an empty system of unknowns where the indeterminants are the indices of the unknown inner products
    system := rec(  unknowns := Positions(rep.innerproducts, false) );
    system.mat := SparseMatrix(0, Size(system.unknowns), [], [], Rationals);
    system.vec := SparseMatrix(0, 1, [], [], Rationals);

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
                    eq := MAJORANA_SeparateInnerProduct(w, x, system.unknowns, rep.innerproducts, rep.setup);
                    eq := eq - MAJORANA_SeparateInnerProduct(y, u, system.unknowns, rep.innerproducts, rep.setup);

                    if Size(eq[1]!.indices[1]) = 1 then
                        # If there is only one unknown value in the equation then immediately record this value
                        MAJORANA_SingleInnerSolution( eq, system, rep.innerproducts );

                        # If all inner products have been found then calculate the nullspace of the alg
                        if system.unknowns = [] then
                            MAJORANA_CheckNullSpace(rep); return;
                        fi;
                    elif eq[1]!.indices[1] <> [] then
                        # Otherwise, add this equation to the system of unknowns
                        eq := eq*(1/eq[1]!.entries[1, 1]);
                        if not _IsRowOfSparseMatrix(system.mat, eq[1]) then
                            system.mat := UnionOfRows(system.mat, eq[1]);
                            system.vec := UnionOfRows(system.vec, eq[2]);
                        fi;
                    fi;
                fi;
            od;
        od;
    od;

    # Solve this system of linear equations
    MAJORANA_SolutionInnerProducts( system, rep.innerproducts );

    if system.unknowns = [] then
        MAJORANA_CheckNullSpace(rep);
    fi;

    return system;

    end );

##
## Finds new eigenvectors using the fusion rules
##

InstallGlobalFunction( MAJORANA_Fusion,

function(rep)

    local   i, ev, a, b, dim, new, FUSE, u;

    FUSE := MAJORANA_FuseEigenvectors;

    # If specified by the optional second argument, do not assume that we have a Frobenius form
    if not IsBound(rep.innerproducts) then FUSE := MAJORANA_FuseEigenvectorsNoForm; fi;

    # Loop over representatives of the orbits of G on T
    for i in rep.setup.orbitreps do
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

        # Continue to loop until no new eigenvectors have been found
        while true do
            Info(   InfoMajorana, 50, STRINGIFY("Fusion of ", i, " evecs")) ;

            # A record in which to store the new eigenvectors found from fusion
            new := ShallowCopy(rep.evecs[i]);

            # For pairs of eigenvalues
            for ev in UnorderedTuples(RecNames(rep.evecs[i]), 2) do
                for a in Iterator( rep.evecs[i].(ev[1]) ) do
                    for b in Iterator( rep.evecs[i].(ev[2]) ) do
                        # Fuse eigenvectors a and b
                        FUSE(a, b, u, ev, new, rep);
                    od;
                od;
            od;

            # Find a basis of the new eigenspaces
            for ev in RecNames(new) do
                new.(ev) := MAJORANA_BasisOfEvecs(new.(ev));
            od;

            # If no new eigenvectors have been found then break
            if ForAll(RecNames(new), ev -> Nrows(new.(ev)) = Nrows(rep.evecs[i].(ev))) then
                break;
            fi;

            rep.evecs[i] := new;

            # Calculate the intersection of the eigenspaces to find new nullspace vectors
            MAJORANA_IntersectEigenspaces(rep);
        od;
    od;

    end );

##
## Uses eigenvectors, the resurrection principle and nullspace vectors to find
## algebra product values
##

InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(rep)

    local   x, i, j, k, evals, system, u, a, b, c, bad, list, evecs_a, evecs_b, index, n, evals_list;

    evals_list := [["0", "1/4"], ["1/4", "0"], ["0", "1/32"], ["1/4", "1/32"]];

    # If specified by the optional second argument, do not assume that we have a Frobenius form
    if not IsBound(rep.innerproducts) then evals_list := [["0", "1/4"], ["0", "1/32"]]; fi;

    # Keep track of the number of unknown algebra products
    n := Size(Positions(rep.algebraproducts, false));

    # Setup the system of linear equations
    system := rec(  unknowns := [],
                    mat := SparseMatrix(0, 0, [], [], Rationals),
                    vec := SparseMatrix(0, Size(rep.setup.coords), [], [], Rationals) );

    # Use the known eigenvectors to calculate more algebra products
    MAJORANA_EigenvectorsAlgebraUnknowns(system, rep);
    if not false in rep.algebraproducts then return true; fi;

    # Use the known nullspace vectors to calculate more algebra products
    MAJORANA_NullspaceUnknowns(system, rep);
    if not false in rep.algebraproducts then return true; fi;

    # Use the known eigenvectors and the resurrection principle to calculate more algebra products
    Info(   InfoMajorana, 50, "Building resurrection");

    # For certain pairs of eigenvalues
    for evals in evals_list do
        for i in rep.setup.orbitreps do

            evecs_a := rep.evecs[i].(evals[1]);
            evecs_b := rep.evecs[i].(evals[2]);

            u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

            list := MAJORANA_ListOfBadIndicesForResurrection(evecs_a, evecs_b, rep);

            # Loop over eigevectors a, b and c where a and c are <evals[1]>-eigenvectors
            # and b is an <evals[2]>-eigenvectors.
            for index in [1,2] do
                for j in [1..Nrows(evecs_a)] do
                    c := CertainRows(evecs_a, [j]);
                    for k in list[index, j] do
                        b := CertainRows(evecs_b, [k]);

                        bad := MAJORANA_FindBadIndices( c, rep );

                        for a in Iterator( evecs_a ) do

                            # If this is satisfied then the product (a - b)*c is known
                            if CertainColumns(a, bad) = CertainColumns(b, bad) then

                                # Perform resurrection on these eigenvectors
                                x := MAJORANA_Resurrection(  u, a, b, c, evals, system.unknowns, rep);

                                if x <> fail and x[1]!.indices[1] <> [] then
                                    # If this equation has only one unknown then record this value immediately
                                    if Size(x[1]!.indices[1]) = 1 then
                                        MAJORANA_SolveSingleSolution( x, system, rep );

                                        if not false in rep.algebraproducts then return true; fi;
                                    elif not _IsRowOfSparseMatrix(system.mat, x[1]) then
                                        # Otherwise, add the equation to the system of linear equations
                                        system.mat := UnionOfRows(system.mat, x[1]);
                                        system.vec := UnionOfRows(system.vec, x[2]);
                                    fi;
                                fi;
                            fi;
                        od;
                        # If the matrix is over a certain size then for performance reasons, solve it already
                        if Nrows(system.mat) > Ncols(system.mat) or Nrows(system.mat) > 8000 then
                            MAJORANA_SolutionAlgProducts(system, rep);

                            if not false in rep.algebraproducts then return true; fi;
                        fi;

                    od;
                od;
            od;
        od;
    od;

    MAJORANA_SolutionAlgProducts(system, rep);

    # If no new algebra products have been found then also add the images of the current equations under the group action
    if n = Size(Positions(rep.algebraproducts, false)) then
        return MAJORANA_AllConjugates(system, rep);
    else
        return system;
    fi;

    end );

################################################################################
##
## Functions used in fusion
##
################################################################################

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

    local   ev, x, y, z;

    # Find the eigenvalue of the new vector
    ev := MAJORANA_FusionTable[ evals ][1];

    # Find the product of the eigenvectors a and b
    x := MAJORANA_AlgebraProduct(a, b, rep.algebraproducts, rep.setup);

    if x in [fail, false] then return; fi;

    # If the eigenvalues are both equal to 1/4 or to 1/32 then we need to
    # do more work to recover the eigenvectors
    if evals = ["1/4", "1/4"] then
        y := MAJORANA_InnerProduct(a, b, rep.innerproducts, rep.setup);

        if y = false then
            return;
        fi;

        new.("0") := MAJORANA_AddEvec(new.("0"), x - (1/4)*u*y);
    elif evals = ["1/32", "1/32"] then
        y := MAJORANA_InnerProduct(a, b, rep.innerproducts, rep.setup);

        if y = false then
            return;
        fi;

        z := MAJORANA_AlgebraProduct(u, x, rep.algebraproducts, rep.setup);

        if z in [false,  fail] then
            return;
        fi;

        new.("1/4") := MAJORANA_AddEvec(new.("1/4"), z - (1/32)*u*y);
        new.("0") := MAJORANA_AddEvec(new.("0"), x + (3/32)*u*y - 4*z);
    else
        new.(String(ev)) := MAJORANA_AddEvec(new.(String(ev)), x);
    fi;

    end );

##
## Performs the same function as <MAJORANA_FuseEigenvectors> but does not assume
## the existence of a Frobenius form
##

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
## Calculates the reversed echelon form of a matrix
##

InstallGlobalFunction(MAJORANA_BasisOfEvecs,

    function(mat)

    local ech, dim;

    dim := Ncols(mat);

    ech := EchelonMatDestructive(CertainColumns(mat, [dim, dim - 1..1]));

    return CertainColumns(ech.vectors, [dim, dim - 1..1]);

    end);

##
## Adds any nonzero intersection of eigenspaces to the nullspace
##

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

################################################################################
##
## Functions used in MAJORANA_UnknownAlgebraProducts
##
################################################################################

##
##  Uses known eigenvectors to find more algebra products
##

InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(system, rep)

    local   i,          # loop over representatives
            ev,         # loop over eigenvalues
            u,          # vector with 1 in j th position
            v,          # eigenvector
            x;          # result of SeparateAlgebraProduct

    Info( InfoMajorana, 50, "Building eigenvector unknowns");

    # Loop over the representatives of the action of G on T
    for i in rep.setup.orbitreps do

        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);

        # Loop over the three eigenvalues 0, 1/4 and 1/32
        for ev in rep.eigenvalues do
            for v in Iterator(rep.evecs[i].(String(ev))) do

                # Create the equation u*evec = 0
                x := MAJORANA_SeparateAlgebraProduct(u, v, system.unknowns, rep.algebraproducts, rep.setup);

                if x <> fail then

                    # Adjust to make the equation u*evec = ev*evec
                    x[2] := x[2] + ev*v;

                    # If this equation has only one unknown, record its value immediately
                    if Size(x[1]!.indices[1]) = 1 then
                        MAJORANA_SolveSingleSolution( x, system, rep );

                        if not false in rep.algebraproducts then return true; fi;

                    elif x[1]!.indices[1] <> [] and not _IsRowOfSparseMatrix(system.mat, x[1]) then
                        # Otherwise add the equation to the linear system
                        system.mat := UnionOfRows(system.mat, x[1]);
                        system.vec := UnionOfRows(system.vec, x[2]);
                    fi;
                fi;
            od;
        od;
    od;

    MAJORANA_SolutionAlgProducts(system, rep);

    end);

##
## Uses known nullspace vectors to find new algebra products
##

InstallGlobalFunction( MAJORANA_NullspaceUnknowns,

    function(system, rep)

    local   i, j, u, v, x, dim;

    if (rep.setup.nullspace.vectors) = 0 then return; fi;

    Info( InfoMajorana, 50, "Building nullspace unknowns" );

    dim := Size(rep.setup.coords);

    # Calculate the orbits of G on the spanning set coords
    x := MAJORANA_Orbits(rep.generators, rep.setup);

    # TODO these should be in place but makes it slower :(
    # Append( rep.setup.conjelts, x.conjelts );
    # rep.setup.conjelts := DuplicateFreeList(rep.setup.conjelts);

    # Loop over the representatives of these orbits
    for i in x.orbitreps do
        if  ForAny(rep.setup.pairorbit[i], k -> rep.algebraproducts[AbsInt(k)] = false) then

            u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);

            for v in Iterator(rep.setup.nullspace.vectors) do
                # Calculate the equation u*v = 0
                x := MAJORANA_SeparateAlgebraProduct(u,v,system.unknowns,rep.algebraproducts,rep.setup);

                # If the equation has only one unknown then immediately record this value
                if x <> fail and Size(x[1]!.indices[1]) = 1 then

                    MAJORANA_SolveSingleSolution( x, system, rep);
                    if not false in rep.algebraproducts then return true; fi;

                elif x <> fail and x[1]!.indices[1] <> [] then
                    # Otherwise add this equation to the system of linear equations
                    if not _IsRowOfSparseMatrix(system.mat, x[1]) then
                        system.mat := UnionOfRows(system.mat, x[1]);
                        system.vec := UnionOfRows(system.vec, x[2]);
                    fi;
                fi;
            od;
        fi;

        # If the matrix is too big then for performance reasons, solve already
        if Nrows(system.mat) > 8000 then
            MAJORANA_SolutionAlgProducts(system, rep);

            if not false in rep.algebraproducts then return; fi;
        fi;
    od;

    MAJORANA_SolutionAlgProducts(system, rep);

    end );

##
## Performs the resurrection priciple on the eigenvectors a, b and c where
## a, c are <evals[1]>-eigenvectors and b in a <evals[2]>-eigenvector.
##

InstallGlobalFunction(MAJORANA_Resurrection,

    function(u, a, b, c, evals, unknowns, rep)

    local   x, y, z, eqn, diff, proj, i, ev;

    # Calculate the product ev*b*c
    ev := MAJORANA_FusionTable[ evals ][1];
    eqn := ev*MAJORANA_SeparateAlgebraProduct(b, c, unknowns, rep.algebraproducts, rep.setup);

    # If this product has known unknowns then it is not useful to us
    if eqn[1]!.indices[1] = [] or eqn = fail then return fail; fi;

    # Calculate the product (a - b)*c
    x := MAJORANA_AlgebraProduct(c, a - b, rep.algebraproducts, rep.setup);
    if x = fail then return fail; fi;

    # Calculate the product u*((a - b)*c)
    z := MAJORANA_SeparateAlgebraProduct(u, x, unknowns, rep.algebraproducts, rep.setup);
    if z = fail then return fail; fi;

    eqn := eqn + z;
    diff := a - b;

    # If a and c are 1/4-eigenvectors then we need to calculate the projection
    # of (a - b)*c onto the 1-eigenspace of u.
    if evals[1] = "1/4" then
        i := u!.indices[1, 1];

        proj := SparseMatrix(1, diff!.ncols, [[i]], [[GetEntry(diff, 1, i)]], Rationals);
        y := MAJORANA_InnerProduct(diff - proj, c, rep.innerproducts, rep.setup);

        if y <> false then
            eqn[2] := eqn[2] + (1/4)*y*u;
        else
            return fail;
        fi;
    fi;

    return eqn;

    end );

##
## Takes a system of linear equations and calculates the image of each row under
## the action of each of the elements in conjelts.
##

InstallGlobalFunction( MAJORANA_AllConjugates,

    function(system, rep)

    local i, new, g, conj, x;

    Info(   InfoMajorana, 50, "All conjugates") ;

    new := rec( mat := SparseMatrix( 0, 0, [], [], Rationals ),
                vec := SparseMatrix( 0, Ncols(system.vec), [], [], Rationals ),
                unknowns := system.unknowns );

    # Loop over group elements and matrix rows
    for g in DuplicateFreeList(rep.setup.conjelts) do
        for i in [1 .. Nrows(system.mat)] do
            if system.mat!.indices[i] <> [] then
                conj := [,];

                # Calculate the images under g
                conj[1] := MAJORANA_ConjugateRow(CertainRows(system.mat, [i]), g, new.unknowns );
                conj[2] := MAJORANA_ConjugateVec(CertainRows(system.vec, [i]), g);

                new.mat := UnionOfRows(new.mat, conj[1]);
                new.vec := UnionOfRows(new.vec, conj[2]);
            fi;
        od;

        # If the matrix is sufficiently big then for performance reasons, solve already
        if Nrows(new.mat) > Ncols(new.mat) then

            MAJORANA_SolutionAlgProducts(new, rep);

            if not false in rep.algebraproducts then return true; fi;

            # And also remove the known products from the original system
            MAJORANA_RemoveKnownAlgProducts(system, rep, false);
        fi;
    od;

    MAJORANA_SolutionAlgProducts(new, rep);

    return new;

    end );

################################################################################
##
## The product functions
##
################################################################################

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

    # Loop over the non-zero coefficients of u and v in reverse order
    for i in Reversed([1..Size(u!.indices[1])]) do
        for j in Reversed([1..Size(v!.indices[1])]) do

            # Find the representative of the orbital containing (i,j)
            k := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

            # Adjust the sign
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

                # If we have already seen this elt then add the product <x> to the corresponding
                # vector in vecs, else add the elt and the product to these lists.
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

    # Now go over all group elements are permute their corresponding vectors, add the
    # result to the output vec as we go along
    for i in [1..Size(elts)] do
        x := MAJORANA_ConjugateVec(vecs[i], setup.pairconjelts[elts[i]]);
        AddRow(x!.indices[1], x!.entries[1], vec!.indices, vec!.entries, 1);
    od;

    if Nrows(setup.nullspace.vectors) > 0 then
        return RemoveMatWithHeads(vec, setup.nullspace);
    else
        return vec;
    fi;

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
            k := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

            # Adjust the sign
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

################################################################################
##
## Functions for finding indices that give unknown algebra products
##
################################################################################

##
## Finds the indices i such that v_i*v is not known
##

InstallGlobalFunction(MAJORANA_FindBadIndices,

    function(v, rep)

    local   i, j, k, list, bad;

    bad := [];
    list := [1..Size(rep.setup.coords)];

    # Loop over the non-zero coefficients of v
    for i in v!.indices[1] do
        # Loop over all vectors in coords
        for j in list do
            k := rep.setup.pairorbit[i, j];

            if k < 0 then k := -k; fi;

            if rep.algebraproducts[k] = false then
                # Add this index to the list of unknown indices and remove it
                # from <list> as we don't need to check it again.
                Add(bad,j);
                list := Difference(list,[j]);
            fi;
        od;
    od;

    Sort(bad);

    return bad;

    end );

##
## For the each eigenvector a in <evecs_a> create a list of indices of the eigenvectors b in
## <evecs_b> such that the product ab involves only one unknown product. We will calculate
## with these eigenvectors first.
##

InstallGlobalFunction(MAJORANA_ListOfBadIndicesForResurrection,

    function(evecs_a, evecs_b, rep)

    local list, j, k, x, bad;

    list := [,];

    list[1] := List([1..Nrows(evecs_a)], i -> []);
    list[2] := List([1..Nrows(evecs_a)], i -> []);

    for j in [1..Nrows(evecs_a)] do
        bad := MAJORANA_FindBadIndices( CertainRows(evecs_a,[j]), rep );
        for k in [1..Nrows(evecs_b)] do
            x := Size(Intersection(bad, evecs_b!.indices[k]));
            if x = 1 then
                Add(list[1, j], k);
            elif x > 1 then
                Add(list[2, j], k);
            fi;
        od;
    od;

    return list;

end );

################################################################################
##
## The conjugating functions
##
################################################################################

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

    # If g is the identity on vec then return
    if ForAll(mat!.indices[1], i -> g[i] = i) then return mat; fi;

    res := SparseMatrix(1, Ncols(mat), [[]], [[]], Rationals);

    # Loop over the non-zero indices of vec and add their image to res
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

InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns)

    local   output,     # output row
            i,          # loop over length of row
            y,
            k,
            sign,       # corrects sign of 5A axis
            pos;        # position of new product

    # If the elt is the identity then return
    if ForAll(g, i -> g[i] = i) then return row; fi;

    output  := SparseMatrix(1, Ncols(row), [[]], [[]], Rationals);

    # Loop over the non-zero coefficients of row
    for i in [1..Size(row!.indices[1])] do

        # Find the image of the corresponding element of unknowns
        y := g{ unknowns[row!.indices[1, i]] };
        Sort(y);

        # Adjust the sign
        sign := 1;
        if y[1] < 0 then sign := -sign; y[1] := -y[1]; fi;
        if y[2] < 0 then sign := -sign; y[2] := -y[2]; fi;

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

################################################################################
##
## Ancilliary functions for finding unknown algebra products
##
################################################################################

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
            pos;        # position of unknown product

    row := SparseMatrix( 1, Size(unknowns), [[]], [[]], Rationals);
    sum := SparseMatrix( 1, Size(setup.coords), [[]], [[]], Rationals);

    # <elts> will record the permutations that we have to conjugate by
    # <vecs> will record the corresponding vector that we are permuting
    elts := [];
    vecs := [];

    # Loop over the non-zero coefficients of u and v
    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            # Find the representative of the orbital containing the pair (i,j)
            k := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

            # Adjust the sign
            if k > 0 then
                sign := 1;
            else
                sign := -1;
                k := -k;
            fi;

            x := algebraproducts[k];

            if x = fail then return fail; fi;

            # If product is known then calculate as usual
            if x <> false then

                g := setup.pairconj[u!.indices[1, i], v!.indices[1, j]];

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

    for i in [1..Size(elts)] do
        x := MAJORANA_ConjugateVec(vecs[i], setup.pairconjelts[elts[i]]);
        AddRow(x!.indices[1], x!.entries[1], sum!.indices, sum!.entries, 1);
    od;

    if Nrows(setup.nullspace.vectors) > 0 then
        sum := RemoveMatWithHeads(sum, setup.nullspace);
    fi;

    return [row, sum];

    end);

MAJORANA_SolveSystem_Whatever := MAJORANA_SolveSystem;
InstallGlobalFunction( MAJORANA_SolutionAlgProducts,

    function( system, rep)

    local   sol,        # solution of system
            sign,       # correct sign of 5A axes
            i,          # loop over <unknowns>
            d,
            nonzero;

    # If the matrix is zero then return
    if ForAll(system.mat!.indices, x -> x = []) then return; fi;

    system.mat!.ncols := Size(system.unknowns);

    Info( InfoMajorana, 40, STRINGIFY("Solving a ", Nrows(system.mat), " x ", Ncols(system.mat), " matrix") );

    # Turn the matrix into an integer matrix
    for i in [1..Nrows(system.mat)] do
        d := _FoldList2(system.mat!.entries[i], DenominatorRat, LcmInt);
        system.mat!.entries[i] := system.mat!.entries[i]*d;
        system.vec!.entries[i] := system.vec!.entries[i]*d;
    od;

    # Solve the system of linear equations
    MAJORANA_SolveSystem_Whatever(system);

    Info(   InfoMajorana, 40, "Solved it!" );

    # If no new solutions have been found then return
    if ForAll(system.solutions, x -> x = fail) then
        Unbind(system.solutions);
        return;
    fi;

    # Otherwise, record the new solutions
    for i in [1..Size(system.unknowns)] do
        if system.solutions[i] <> fail then
            MAJORANA_RecordSolution( system.solutions[i], system.unknowns[i], rep );
        fi;
    od;

    if not false in rep.algebraproducts then return system; fi;

    Unbind(system.solutions);

    # Adjust the system of linear equations to take into account the new known products
    MAJORANA_RemoveKnownAlgProducts(system, rep);

    MAJORANA_SolutionAlgProducts(system, rep);

    end );

InstallGlobalFunction( MAJORANA_SolveSingleSolution,

    function(x, system, rep)

    local   elm, ind, y, switch, i;

    Info( InfoMajorana, 60, "Solved a single solution");

    # Divide through by the coefficient of the indeterminant
    elm := x[1]!.entries[1, 1];
    x := x/elm;

    # Record the new algebra product
    MAJORANA_RecordSolution( x[2], system.unknowns[x[1]!.indices[1, 1]], rep );

    # Reduce the system of linear equations using this new product
    MAJORANA_RemoveKnownAlgProducts( system, rep );

    # While we continue to find new products
    switch := true;

    while switch = true do

        if Nrows(system.mat) = 0 or system.unknowns = [] then return; fi;

        switch := false;

        for i in [1..Nrows(system.mat)] do
            if Size(system.mat!.indices[i]) = 1 then
            # We have found a new single solution
                switch := true;
                ind := system.unknowns[system.mat!.indices[i, 1]];
                elm := system.mat!.entries[i, 1];
                MAJORANA_RecordSolution( CertainRows(system.vec, [i])*(1/elm), ind, rep);
            fi;
        od;

        if switch = true then
            Info( InfoMajorana, 60, "Solved a new single solution");
            # Reduce the system of linear equations using this new product
            MAJORANA_RemoveKnownAlgProducts( system, rep );
        fi;
    od;

    MAJORANA_SolutionAlgProducts(system, rep);

    end );

##
## Takes a vector that is the product of a pair of vectors given by x and records
## it in the correct place in the list rep.algebraproducts.
##

InstallGlobalFunction( MAJORANA_RecordSolution,

    function( v, x, rep)

    local   y,
            g,
            sign;

    y := rep.setup.pairorbit[x[1], x[2]];
    g := SP_Inverse(rep.setup.pairconjelts[rep.setup.pairconj[x[1], x[2]]]);

    # Adjust the sign
    sign := 1;
    if y < 0 then sign := -1; y := -y; fi;

    # Record the new product
    if rep.algebraproducts[y] = false then
        v := sign*MAJORANA_ConjugateVec(v,g);
        rep.algebraproducts[y] := RemoveMatWithHeads(v, rep.setup.nullspace);
    fi;

    end );

##
## Takes a system [mat, vec] of unknown algebra products and removes
## from the system any variables which have already been found
##

InstallGlobalFunction( MAJORANA_RemoveKnownAlgProducts,

    function( arg )

    local   system, rep, nonzero,
            unsolved,
            i,
            j,
            elm,
            x,
            y,
            sign,
            g,
            pos,
            prod,
            nonzero_rows;

    system := arg[1]; rep := arg[2];

    nonzero := true;

    if Length(arg) > 2 then nonzero := arg[3]; fi;

    if Nrows(system.mat) = 0 then return; fi;

    unsolved := [];

    # Loop over the unknown algebra products
    for i in [1..Size(system.unknowns)] do

        # Find the representative of the orbital containing the unknown value
        x := system.unknowns[i];
        y := rep.setup.pairorbit[x[1], x[2]];

        # Adjust the sign
        sign := 1;
        if y < 0 then sign := -1; y := -y; fi;

        prod := rep.algebraproducts[y];

        # If the product is now known the remove its value from the rhs
        if prod <> false then
            g := rep.setup.pairconjelts[rep.setup.pairconj[x[1], x[2]]];
            prod := MAJORANA_ConjugateVec(prod,g);

            for j in [1..Nrows(system.vec)] do
                pos := Position(system.mat!.indices[j], i);
                if pos <> fail then
                    elm := system.mat!.entries[j, pos];
                    AddRow( prod!.indices[1],-sign*elm*prod!.entries[1],
                            system.vec!.indices, system.vec!.entries, j);
                fi;
            od;
        else
            # Otherwise, we want to keep this column of the matrix
            Add(unsolved,i);
        fi;
    od;

    system.mat := CertainColumns(system.mat, unsolved);
    system.unknowns := system.unknowns{unsolved};

    if nonzero = true then
        # Take out any zero rows
        nonzero_rows := Filtered([1..Nrows(system.mat)], j -> system.mat!.indices[j] <> []);
        system.mat := CertainRows(system.mat, nonzero_rows);
        system.vec := CertainRows(system.vec, nonzero_rows);
    fi;

    end );

################################################################################
##
## Ancilliary functions for finding unknown inner products
##
################################################################################

InstallGlobalFunction(MAJORANA_SeparateInnerProduct,

    function(u,v,unknowns,innerproducts,setup)

    local   row,            # record values of unknowns
            sum,            # record values of knowns
            i,              # index for dim of u
            j,              # index for dim of v
            m,              # orbit of i,j
            pos,            # position of m in unknowns
            sign;           # correct sign of 5A axes

    sum := SparseZeroMatrix(1, 1, Rationals);
    row := SparseZeroMatrix(1, Size(unknowns), Rationals);

    # Loop over the nonzero coefficients of u and v
    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do

            # Find the orbital containing the pair (i,j)
            m := setup.pairorbit[u!.indices[1, i], v!.indices[1, j]];

            # Adjust the sign
            if m > 0 then
                sign := 1;
            else
                sign := -1;
                m := -m;
            fi;

            # If the product is known the calculate as usual. Otherwise, add the coefficient
            # to the matrix of the system of linear equations.
            if innerproducts[m] <> false then
                AddToEntry(sum, 1, 1, - sign*u!.entries[1, i]*v!.entries[1, j]*innerproducts[m]);
            else
                pos := Position(unknowns,m);
                AddToEntry(row, 1, pos, sign*u!.entries[1, i]*v!.entries[1, j]);
            fi;
        od;
    od;

    return [row,sum];

    end );

##
## Takes a system of linear equations whose unknowns are inner product values and
## removes any unknowns from the system if there function has already been found
##

InstallGlobalFunction( MAJORANA_RemoveKnownInnProducts,

    function(system, innerproducts)

    local   unsolved, i, j, elm, prod, nonzero_rows;

    unsolved := [];

    if Nrows(system.mat) = 0 then return; fi;

    # Loop over the unknown inner product values
    for i in [1..Size(system.unknowns)] do
        prod := innerproducts[system.unknowns[i]];

        # If the product is now known then loop over the rows and remove this value from the rhs
        if prod <> false then
            for j in [1..Nrows(system.vec)] do
                elm := GetEntry(system.mat, j, i);

                if elm <> 0 then
                    AddToEntry(system.vec, j, 1, -elm*prod);
                fi;
            od;
        else
            Add(unsolved, i);
        fi;
    od;

    # Otherwise, we want to keep this column of the matrix
    system.mat := CertainColumns(system.mat, unsolved);
    system.unknowns := system.unknowns{unsolved};

    # Remove any zero rows
    nonzero_rows := Filtered([1..Nrows(system.mat)], j -> system.mat!.indices[j] <> []);
    system.mat := CertainRows(system.mat, nonzero_rows);
    system.vec := CertainRows(system.vec, nonzero_rows);

    end );

##
## Takes a linear equation with just one unknown inner product value and
## records the new value.
##

InstallGlobalFunction( MAJORANA_SingleInnerSolution,

    function(eq, system, innerproducts)

    local x;

    x := system.unknowns[eq[1]!.indices[1, 1]];

    if eq[2]!.entries[1] = [] then
        innerproducts[x] := 0;
    else
        innerproducts[x] := eq[2]!.entries[1, 1]/eq[1]!.entries[1, 1];
    fi;

    MAJORANA_RemoveKnownInnProducts(system, innerproducts);

    end );

##
## Takes a system of linear equations where are the unknowns are inner products,
## solves the system and records any new values.
##

InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( system, innerproducts)

    local   i, x;

    if Nrows(system.mat) = 0 then return; fi;

    # Solve the system of linear equations
    MAJORANA_SolveSystem(system);

    # Record any new solutions that have been found
    for i in [1..Size(system.solutions)] do
        if system.solutions[i] <> fail then
            x := system.unknowns[i];
            if system.solutions[i]!.entries[1] = [] then
                innerproducts[x] := 0;
            else
                innerproducts[x] := system.solutions[i]!.entries[1, 1];
            fi;
        fi;
    od;

    Unbind(system.solutions);

    MAJORANA_RemoveKnownInnProducts( system, innerproducts);

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
            k := setup.pairorbit[range[i], range[j]];

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

##
## Calculate the nullspace from the Gram matrix of the form and removes the
## relevant vectors from algebra products and eigenvectors
##

InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(rep)

    local   dim, gram, null, unknowns, i, ev, x;

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
        x := Filtered([1..dim], j -> i in rep.setup.pairorbit[j]);
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
        for ev in RecNames(rep.evecs[i]) do
            rep.evecs[i].(ev) := RemoveMatWithHeads(rep.evecs[i].(ev), null);
            rep.evecs[i].(ev) := MAJORANA_BasisOfEvecs(rep.evecs[i].(ev));
        od;
    od;

    end );

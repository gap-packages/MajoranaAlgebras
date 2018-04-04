# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
# This is a GAP prototype which already works quite a bit faster than any code
# inside GAP, for reasonably large systems (thousands of equations)
#
# If you find any bugs, email markus.pfeiffer@morphism.de
#
# TODO:
# * Carry denominator forward
# * actually only solve the solvable variables.
# * Make a better implementation of the padics code. Its currently pretty brittle
#   and hacky
# * More tests
# * look at flint that has some of this functionality
# * Implement in C or Rust (or Julia)?
# * Parallelisation strategies
# * use meataxe64
#

# Just to make sure we're not shooting ourselves
# in the foot with inconsistent entries.
CheckSystem := function(system)
    local b, r, c;

    Info( InfoMajoranaLinearEq, 5,
          " testing system of equation structure" );
    if not IsPrime(system.p) then
        Error("p is not prime");
    fi;

    for r in system.int_mat do
        for c in r do
            if DenominatorRat(c) <> 1 then
                Error("Non-1-denominator in system.int_mat");
            fi;
        od;
    od;

    for r in system.int_vecs do
        for c in r do
            if DenominatorRat(c) <> 1 then
                Error("Non-1-denominator in system.int_vecs");
            fi;
        od;
    od;
    Info( InfoMajoranaLinearEq, 5,
          " success.");
end;

InstallGlobalFunction( MAJORANA_Padic_Presolve,
function(system)
    local v, r, i;

    system.mat_mod_p := system.int_mat * Z(system.p)^0;
    ConvertToMatrixRep(system.mat_mod_p);

    system.echelon := EchelonMatTransformation(system.mat_mod_p);

    # FIXME: refactor
    v := List( system.echelon.vectors
             , x -> PositionsProperty(x, y -> not IsZero(y) ) );
    system.solvable_rows := PositionsProperty( v, z -> Length(z) = 1);

    system.interesting_rows := Set(Concatenation( List( system.echelon.coeffs,
                                               r -> PositionsProperty( r, x -> not IsZero(x)) ) ));
    system.solvable_variables := Concatenation( Filtered( v
                                                        , z -> Length(z) = 1) );
    system.unsolvable_variables := Difference([1..system.number_variables], system.solvable_variables);
    system.unsolvable_rows := Difference([1..system.number_equations], system.solvable_rows);
    system.uninteresting_rows := Difference([1..system.number_equations], system.interesting_rows);

    return system;
end );

InstallGlobalFunction(MAJORANA_SetupMatVecsSystem_Padic,
function(mat, vecs, p, precision, max_iter)
    local system, mmults, vmults, lcm;
    system := rec( mat := mat
                 , vecs := vecs
                 , number_variables := Length(mat[1])
                 , number_equations := Length(mat) );

    #  MakeIntSystem(system);
    Info(InfoMajoranaLinearEq, 5,
         "MakeIntSystem2: computing denominator lcms" );

    mmults := List(system.mat, x -> _FoldList2(x, DenominatorRat, LcmInt));
    vmults := List(system.vecs, x -> _FoldList2(x, DenominatorRat, LcmInt));
    lcm := _FoldList2(Concatenation(mmults, vmults), IdFunc, LcmInt);

    Info(InfoMajoranaLinearEq, 5,
         "MakeIntSystem2: lcm: ", lcm);

    system.lcm := lcm;
    system.int_mat := system.mat * lcm;
    system.int_vecs := system.vecs * lcm;

    system.p := p;
    system.precision := precision;
    system.padic_family := PurePadicNumberFamily(p, precision);
    system.padic_iterations := max_iter;

    MAJORANA_Padic_Presolve(system);

    # Transposingpalooza
    # FIXME: cleanup and only keep the ones we need
    system.transposed_int_mat := TransposedMat(system.int_mat);
    system.transposed_coeffs := TransposedMat(system.echelon.coeffs);
    system.transposed_vecs := TransposedMat(system.int_vecs);
    system.lifted_coeffs := List(system.transposed_coeffs, y -> List(y, IntFFESymm));

    system.solution_denominator := 1;

    return system;
end);

InstallGlobalFunction( MAJORANA_SolutionIntMatVec_Padic,
function(system, vi)
    local
        p,

        # These are *integer* vectors
        tmp_soln, soln,
        residue,

        # These are vectors in GF(p)
        vec_p,
        soln_p, soln_pp,

        done, iterations,
        soln_padic,
        ppower, sol, x, y, i,
        k, old_denom, denom, vecd, iter;

    p := system.p;

    # Accumulator for integer solution
    # FIXME: only solved variables? Should cut down on memory use
    #        and how many p-adic expansions we have to actually turn
    #        into denominators
    soln := ListWithIdenticalEntries(system.number_variables, 0);

    # These are the *integer* residuals of the RHS
    # initially this is the RHS we're solving for
    # FIXME: for testing, only do one RHS right now
    #        this will trivially generalise fortunately
    #        but there might be a point in not solving all
    #        RHS at the same time, in case we discover enough 
    #        of the denominator to not have to approximate?
    residue := MutableCopyMat(system.transposed_vecs[vi])
               * system.solution_denominator;

    done := false;
    iterations := 0;
    ppower := 1; # this is p^iterations, FIXME: can we avoid this?

    # digits in the p-adic approximation to the solution
    soln_padic := List([1..system.number_variables], x -> PadicNumber(system.padic_family, 0));

    while (not done) do
        iterations := iterations + 1;

        if iterations mod 100 = 0 then
            Info(InfoMajoranaLinearEq, 5, STRINGIFY(iterations, " iterations"));
        fi;

        # solve the system mod p
        vec_p := Z(p)^0 * residue;
        soln_p := vec_p * system.transposed_coeffs;
        soln_p := List( system.echelon.heads,
                      function(x)
                          if x > 0 then
                              return soln_p[x];
                          else
                              return Zero(soln_p[1]);
                          fi;
                      end);

        # Convert the solution from GF(p) to integers -p/2..p/2-1
        y := List(soln_p, IntFFESymm);

        # they are the coefficients of the p-adic expansion of the denominator
        # the below is slow, and hence replaced by the hack below that.
        # soln_padic := soln_padic + List(y, c -> PadicNumber(fam, ppower * -c));
        soln_padic := soln_padic + List(y, c -> PadicNumber(system.padic_family, [iterations, c mod system.padic_family!.modulus ] ) );
        AddRowVector(soln, y, ppower);

        residue := (residue - y * system.transposed_int_mat);
        residue{ system.interesting_rows } := residue{ system.interesting_rows } / p;

        ppower := ppower * p;

        Info(InfoMajoranaLinearEq, 10, "soln:    ", soln);
        Info(InfoMajoranaLinearEq, 10, "y:       ", y);
        Info(InfoMajoranaLinearEq, 10, "residue: ", residue);

        # Solution found?
        if IsZero( residue ) then
            Info(InfoMajoranaLinearEq, 5,
                 "found an integer solution");

            # FIXME: I don't like this state struct design at the moment
            system.int_solution := soln;
            return true;
        else
            if iterations > system.precision then
                Info(InfoMajoranaLinearEq, 5,
                     "reached iteration limit, trying to compute denominator");
                # Compute the least common denominator of them all

                denom := PadicDenominatorList( soln_padic{ system.solvable_variables }, system.padic_iterations);
                Info(InfoMajoranaLinearEq, 5,
                     "found denominator: ", denom);

                if denom = fail then
                    Info( InfoMajoranaLinearEq, 10,
                          "failed to find denominator trying to increase p-adic precision");
                    Info( InfoMajoranaLinearEq, 10,
                          "failed to solve rhs ", vi);
                    Info( InfoMajoranaLinearEq, 10,
                          "rhs not zero: ", system.rhns, " ", residue{ system.rhns } );

                    if IsBound(system.solution_denominator) then
                        Info( InfoMajoranaLinearEq, 10,
                              "got a denominator, so running with that");
                        system.int_solution := soln;
                        return true;
                    fi;

                    Error("Failed to compute denominator");
                    # FIXME: adjust system, i.e. we could increase precision?
                    # MAJORANA_SolutionIntMatVec_Padic(system);
                    return false;
                elif denom = 1 then
                    Error("Denominator 1 occurred. This should not happen.");
                    return false;
                else
                    Info(InfoMajoranaLinearEq, 5,
                         "solving system after multiplying rhs by denominator.");

                    system.solution_denominator := system.solution_denominator * denom;
                    system.rhns := Intersection( system.solvable_rows,
                                                 PositionsProperty(residue, x -> not IsZero(x)));

                    Info(InfoMajoranaLinearEq, 5,
                         "failed to solve rhs ", vi);
                    Info(InfoMajoranaLinearEq, 5,
                         "rhs not zero: ", system.rhns, " ", residue{ system.rhns });

                    # try again with new denominator
                    MAJORANA_SolutionIntMatVec_Padic(system, vi);
                    return true;
                fi;
            fi;
        fi;
    od;
end );

# FIXME: This has to be compatible with MAJORANA_SolutionMatVecs
# FIXME: Make compatible with sparse matrices.
InstallGlobalFunction( MAJORANA_SolutionMatVecs_Padic,
function(mat, vecs)
    local vi, system, res;

    res := rec();
    res.solutions := [];
    res.mat := mat;
    res.vec := vecs;

    mat := ConvertSparseMatrixToMatrix(mat);
    vecs := ConvertSparseMatrixToMatrix(vecs);
    if Length(mat) = 0 or Length(mat[1]) = 0 then
        return res;
    fi;


    system := MAJORANA_SetupMatVecsSystem_Padic( mat, vecs
                                                 , MAJORANA_Padic_Prime
                                                 , MAJORANA_Padic_Precision
                                                 , MAJORANA_Padic_Iterations );
    if Length(system.solvable_variables) > 0 then
    for vi in [1..Length(system.transposed_vecs)] do
        MAJORANA_SolutionIntMatVec_Padic(system, vi);
        Add(res.solutions, system.int_solution / system.solution_denominator );
    od;

    res.solutions := TransposedMatMutable(res.solutions);
    res.solutions := List(res.solutions, x -> SparseMatrix([x], Rationals));
    fi;
    res.solutions{ system.unsolvable_variables } := ListWithIdenticalEntries(Length(system.unsolvable_variables), fail);
   # res.mat := SparseMatrix(system.mat{ system.uninteresting_rows }, Rationals);
   # res.vec := SparseMatrix(system.vecs{ system.uninteresting_rows }, Rationals);

    # Debugging
    res.system := system;
    return res;
end);

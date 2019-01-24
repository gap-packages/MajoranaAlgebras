##
##
## A fix for the UnionOfRows bug in some versions of Gauss
##
##

Info(InfoMajorana, 10, "Warning, rebinding `UnionOfRows`, because it was broken sometime");
MakeReadWriteGVar("UnionOfRows");
UnbindGlobal("UnionOfRows");
BindGlobal("UnionOfRows",
    function( A, B )
      return SparseMatrix( A!.nrows + B!.nrows, A!.ncols, Concatenation( A!.indices, B!.indices ), Concatenation( A!.entries, B!.entries ), A!.ring );
    end );


##
## Takes as its input a record <system> with the components <mat>, <vec> and <unknowns>
## and attempts to solve the system of linear equations where each row of <mat> is
## a linear combinations of indeterminants indexed by <unknowns> and the corresponding
## row of <vec> gives the value of this linear combination.
##
## Changes the argument <system> in situ. Adds the component <solutions> which is a
## list where the ith position is the value of the ith indeterminant (given by
## <system.unknowns>) if found, and fail if not. The components <mat> and <vec>
## are changed to give only the rows that still contain non-zero coefficients
## of unsolved indeterminants.
##
## The argument <vec> may be an n x m matrix, in which case the solutions of the system
## will be row vectors of length m. Here n is the number of rows of <mat>.
##

InstallGlobalFunction(MAJORANA_SolveSystem,

function(system)

    local   n,
            res,
            sol,
            i,
            pos,
            elm,
            rowlist;

    n := Ncols(system.mat);

    res := EchelonMatTransformationDestructive(system.mat);
    system.mat := res.vectors;
    system.vec := res.coeffs*system.vec;

    sol := [1..n]*0;
    rowlist := [];

    for i in Reversed([1 .. n]) do

        pos := res.heads[i];

        if pos = 0 then
            sol[i] := fail;
        elif Size(system.mat!.indices[pos]) > 1 then
            Add(rowlist, pos);
            sol[i] := fail;
        else
            sol[i] := CertainRows(system.vec, [pos]);
        fi;
    od;

    system.solutions := sol;

    end );

##
## Takes as input a (list of lists) matrix A. If A is positive semidefinite
## then will return [L,D] such that A= LDL^T. Otherwise returns false.
##
## Note: does not test if matrix is square or symmetric.
##

InstallGlobalFunction(MAJORANA_LDLTDecomposition,

    function(A)

    local   B,      # input matrix
            n,      # size of matrix
            L,      # output lower triangular matrix
            D,      # output diagonal matrix
            i,      # loop over rows of matrix
            j,      # loop over columns of matrix
            k,      # loop over diagonals
            sum;    # sum used in values of L and D

    B := ShallowCopy(A); n := Size(B); L := NullMat(n,n); D := NullMat(n,n);

    for i in [1..n] do
        sum := [];
        for j in [1..i-1] do
            Add(sum, L[i, j]*L[i, j]*D[j, j]);
        od;

        D[i, i] := B[i, i] - Sum(sum);

        if D[i, i] = 0 then
                for j in [i+1..n] do
                    sum := [];
                    for k in [1..i-1] do
                        Add(sum, L[j, k]*L[i, k]*D[k, k]);
                    od;
                    if B[j, i] - Sum(sum) = 0 then
                        L[j, i]:=0;
                    else
                        return false;
                    fi;
                od;
                L[i, i]:=1;
        else
            for j in [i+1..n] do
                sum := [];
                for k in [1..i-1] do
                    Add(sum, L[j, k]*L[i, k]*D[k, k]);
                od;
                L[j, i] := (B[j, i] - Sum(sum))/D[i, i];
            od;
            L[i, i] := 1;
        fi;
    od;

    return Concatenation([L],[D]);

    end );

##
## Takes a matrix <mat> and returns 1 or 0 is <mat> is positive definite or
## positive semidefinite and returns -1 otherwise.
##

InstallGlobalFunction(MAJORANA_PositiveDefinite,

    function(mat)

    local   L,          # decomposition of matrix
            Diagonals,  # list of diagonals from decomposition
            i;          # loop over sze of matrix

    L := MAJORANA_LDLTDecomposition(mat);

    if L = false then
        return -1;
    fi;

    Diagonals := [];

    for i in [1..Size(mat)] do
        Append(Diagonals,[L[2][i, i]]);
    od;

    if ForAny(Diagonals, x->x<0) then
        return -1;
    elif ForAny(Diagonals, x->x=0) then
        return 0;
    else
        return 1;
    fi;

    end );

InstallGlobalFunction(_FoldList2,
    function(list, func, op)
    local k, s, old_s, r, i, len, n, nh, res, r1, r2;


    len := Length(list);
    # FIXME: We don't know a default value to return here.
    if len = 0 then
        Error("Lists of length 0 are not supported by this function");
    elif len = 1 then
        return func(list[1]);
    fi;

    res := List(list, func);
    k := len;
    s := 1;
    while k > 1 do
        r := k mod 2;
        old_s := s;
        k := QuoInt(k, 2);
        s := s * 2;
        i := s;
        while i <= k * s do
            if IsBound(res[i-old_s]) then
                r1 := res[i-old_s];
            else
                r1 := 1;
            fi;
            if IsBound(res[i]) then
                r2 := res[i];
            else
                r2 := 1;
            fi;
            res[i] := op(r1, r2);
            res[i-old_s] := 0;
            i := i + s;
        od;
        if r = 1 then
            k := k + 1;
            res[i] := res[i-old_s];
        fi;
    od;
    return res[ k * s ];
end );

##
## Takes a matrix <mat> and a row vector <row>, both in sparse matrix format.
## If <row> is already a row of <mat>, return true. Otherwise return false.
##

InstallGlobalFunction(_IsRowOfSparseMatrix,

    function(mat, row)

    local pos;

    pos := Positions(mat!.indices, row!.indices[1]);

    if ForAny(pos, i -> mat!.entries[i] = row!.entries[1]) then
        return true;
    else
        return false;
    fi;

    end);

##
## Performs the same matrix echelon reduction as <EchelonMatDestructive>
## but outputs the RREF of <mat> where the diagonal submatrix is on the right.
##

InstallGlobalFunction(ReversedEchelonMatDestructive,

    function(mat)

    local ncols, ech;

    ncols := Ncols(mat);;

    ech := EchelonMatDestructive(CertainColumns(mat, [ncols, ncols - 1 .. 1]));

    ech.vectors := CertainColumns(ech.vectors, [ncols, ncols - 1 .. 1]);
    ech.heads   := Reversed(ech.heads);

    return ech;

    end );

InstallGlobalFunction( SumIntersectionSparseMat,

    function(mat1, mat2)

    local mat, n, v, sum, row, int, i, M1, M2;


        M1 := CopyMat(mat1);
        M2 := CopyMat(mat2);

    # Basic checks on input
    if Nrows(M1) = 0 then
        return [ M2, M1 ];
    elif Nrows(M2) = 0 then
        return [ M1, M2 ];
    elif Ncols(M1) <> Ncols(M2) then
        Error("Matrices must have the same number of columns");
    elif RingOfDefinition(M1) <> RingOfDefinition(M2) then
        Error("Matrices must be defined over the same ring");
    fi;

    n := Ncols(M1);

    # Set up the matrix for Zassenhaus' algorithm
    mat := SparseMatrix(0, 2*n, [], [], M1!.ring);
    for v in Iterator(M1) do
        mat := UnionOfRows( mat, UnionOfColumns(v, v) );
    od;
    for v in Iterator(M2) do
        mat := UnionOfRows( mat, v );
    od;

    # Transform <mat> into echelon form
    mat := EchelonMatDestructive(mat);

    # Extract the basis for the sum
    sum := SparseMatrix(0, n, [], [], M1!.ring);
    for i in [1 .. n] do
        if mat.heads[i] <> 0 then
            row := CertainRows( mat.vectors, [mat.heads[i]] );
            sum := UnionOfRows( sum, CertainColumns(row, [1..n]) );
        fi;
    od;

    # Extract the basis for the intersection
    int := SparseMatrix(0, n, [], [], M1!.ring);
    for i in [n+1 .. 2*n] do
        if mat.heads[i] <> 0 then
            row := CertainRows(mat.vectors, [mat.heads[i]]);
            row := CertainColumns(row, [n+1 .. 2*n]);
            int := UnionOfRows(int, row );
        fi;
    od;

    return [ sum, int ];

end );

##
## Takes <mat>, a generic sparse matrix, and <null>, the output from the <EchelonMat> or
## <ReversedEchelonMat> functions such that the <null.vectors> has the same number of
## columns as <mat>. Returns that matrix <mat> that has been reduced wrt <null.vectors>.
##

InstallGlobalFunction(RemoveMatWithHeads,

    function(mat, null)

    local i, j, k, x;

    for i in [1 .. Nrows(mat)] do
        for j in mat!.indices[i] do
            k := null.heads[j];
            if k <> 0 then
                x := -GetEntry(mat, i, j);
                AddRow(null.vectors!.indices[k], x*null.vectors!.entries[k],
                        mat!.indices, mat!.entries, i);
            fi;
        od;
    od;

    return mat;

    end );

##
## Implement an iterator that runs over the rows of a sparse matrix
##

BindGlobal( "NextIterator_SparseMatrix", function( iter )

    iter!.pos := iter!.pos + 1;
    return CertainRows( iter!.matrix, [iter!.pos] );

    end );

InstallOtherMethod( Iterator, "for a sparse matrix",
    [ IsSparseMatrix ],
    function( M )
        local iter;
        iter := rec(    matrix := M,
                        pos := 0 );

        iter.NextIterator   := NextIterator_SparseMatrix;
        iter.IsDoneIterator := iter -> ( iter!.pos = Nrows(iter!.matrix) );
        iter.ShallowCopy    := iter -> rec( pos := iter!.pos,
                                            matrix := iter!.matrix ) ;

        return IteratorByFunctions(iter);

    end );

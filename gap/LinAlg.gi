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

    function(A, ring) # Takes as input a matrix A. If A is positive semidefinite then will return [L,D] such that A= LDL^T. Else returns 0. Note: does not test if matrix is square or symmetric.

    local   B,      # input matrix
            n,      # size of matrix
            L,      # output lower triangular matrix
            D,      # output diagonal matrix
            i,      # loop over rows of matrix
            j,      # loop over columns of matrix
            k,      # loop over diagonals
            sum;    # sum used in values of L and D

    B := ShallowCopy(A); n := Size(B); L := NullMat(n,n)*One(ring); D := NullMat(n,n)*One(ring);

    for i in [1..n] do
        sum := [];
        for j in [1..i-1] do
            Add(sum, L[i, j]*L[i, j]*D[j, j]);
        od;

        D[i][i] := B[i][i] - Sum(sum)*One(ring);

        if D[i][i] = Zero(ring) then
                for j in [i+1..n] do
                    sum := [];
                    for k in [1..i-1] do
                        Add(sum, L[j, k]*L[i, k]*D[k, k]);
                    od;
                    if B[j][i] - Sum(sum) = Zero(ring) then
                        L[j][i]:= Zero(ring);
                    else
                        return false;
                    fi;
                od;
                L[i][i] := One(ring);
        else
            for j in [i+1..n] do
                sum := [];
                for k in [1..i-1] do
                    Add(sum, L[j, k]*L[i, k]*D[k, k]);
                od;
                L[j][i] := (B[j][i] - Sum(sum)*One(ring))/D[i][i];
            od;
            L[i][i] := One(ring);
        fi;
    od;

    return Concatenation([L],[D]);

    end );

##
## Takes a matrix <mat> and returns 1 or 0 is <mat> is positive definite or
## positive semidefinite and returns -1 otherwise.
##

InstallGlobalFunction(MAJORANA_PositiveDefinite,

    function(GramMatrix, field) # Check returns 1, 0, -1 if Gram matrix is positive definite, positive semidefinite or neither respectively

    local   L,          # decomposition of matrix
            Diagonals,  # list of diagonals from decomposition
            i;          # loop over sze of matrix

    L := MAJORANA_LDLTDecomposition(GramMatrix, field);

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

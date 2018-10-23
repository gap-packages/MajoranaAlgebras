
InstallGlobalFunction(MAJORANA_SolutionMatVecs,

function(mat,vec) # Takes as input two matrices, the second being interpreted as a vector of vectors. Returns record where solutions[i] gives the value of unknown variable i if found, and fail otherwise

    local   n,
            res,
            sol,
            i,
            pos,
            elm,
            rowlist;

    n := Ncols(mat);

    res := EchelonMatTransformationDestructive(mat);
    mat := res.vectors;
    vec := res.coeffs*vec;

    sol := [1..n]*0;
    rowlist := [];

    for i in Reversed([1 .. n]) do

        pos := res.heads[i];

        if pos = 0 then
            sol[i] := fail;
        elif Size(mat!.indices[pos]) > 1 then
            Add(rowlist, pos);
            sol[i] := fail;
        else
            sol[i] := CertainRows(vec, [pos]);
        fi;
    od;

    return rec( solutions := sol,
                mat := CertainRows(mat, rowlist),
                vec := CertainRows(vec, rowlist)  );

    end );

InstallGlobalFunction(MAJORANA_LDLTDecomposition,

    function(A) # Takes as input a matrix A. If A is positive semidefinite then will return [L,D] such that A= LDL^T. Else returns 0. Note: does not test if matrix is square or symmetric.

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

InstallGlobalFunction(MAJORANA_PositiveDefinite,

    function(GramMatrix) # Check returns 1, 0, -1 if Gram matrix is positive definite, positive semidefinite or neither respectively

    local   L,          # decomposition of matrix
            Diagonals,  # list of diagonals from decomposition
            i;          # loop over sze of matrix

    L := MAJORANA_LDLTDecomposition(GramMatrix);

    if L = false then
        return -1;
    fi;

    Diagonals := [];

    for i in [1..Size(GramMatrix)] do
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

InstallGlobalFunction(_IsRowOfSparseMatrix,

    function(mat, row)

    local pos;

    pos := Position(mat!.indices, row!.indices[1]);

    if pos = fail then
        return false;
    fi;

    if mat!.entries[pos] = row!.entries[1] then
        return true;
    else
        return false;
    fi;

    end);

InstallGlobalFunction(ReversedEchelonMatDestructive,

    function(mat)

    local ncols, ech;

    ncols := Ncols(mat);;

    ech := EchelonMatDestructive(CertainColumns(mat, [ncols, ncols - 1 .. 1]));

    ech.vectors := CertainColumns(ech.vectors, [ncols, ncols - 1 .. 1]);
    ech.heads   := Reversed(ech.heads);

    return ech;

    end );

InstallGlobalFunction(RemoveMatWithHeads,

    function(mat, null)

    local v, i, j, k, x;

    v := null.vectors;

    for j in PositionsProperty(null.heads, x -> x <> 0) do
        k := null.heads[j];

        for i in [1..Nrows(mat)] do
            x := -GetEntry(mat, i, j);

            if x <> 0 then
                AddRow(v!.indices[k], x*v!.entries[k], mat!.indices, mat!.entries, i);
            fi;
        od;
    od;

    return mat;

    end );

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

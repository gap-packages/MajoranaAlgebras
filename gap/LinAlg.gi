 
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
                Add(sum, L[i][j]*L[i][j]*D[j][j]);
            od;

            D[i][i] := B[i][i] - Sum(sum);

            if D[i][i] = 0 then
                    for j in [i+1..n] do
                        sum := [];
                        for k in [1..i-1] do
                            Add(sum, L[j][k]*L[i][k]*D[k][k]);
                        od;
                        if B[j][i] - Sum(sum) = 0 then
                            L[j][i]:=0;
                        else
                            return false;
                        fi;
                    od;
                    L[i][i]:=1;
            else
                for j in [i+1..n] do
                    sum := [];
                    for k in [1..i-1] do
                        Add(sum, L[j][k]*L[i][k]*D[k][k]);
                    od;
                    L[j][i] := (B[j][i] - Sum(sum))/D[i][i];
                od;
                L[i][i] := 1;
            fi;
        od;

        return Concatenation([L],[D]);
        end
    );
    
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
        Append(Diagonals,[L[2][i][i]]);
    od;

    if ForAny(Diagonals, x->x<0) then
        return -1;
    elif ForAny(Diagonals, x->x=0) then
        return 0;
    else
        return 1;
    fi;
    
    end

    );

InstallGlobalFunction(MAJORANA_ReversedEchelonForm,
function( mat )
    local nrows,     # number of rows in <mat>
          ncols,     # number of columns in <mat>
          vectors,   # list of basis vectors
          i,         # loop over rows
          j,         # loop over columns
          x,         # a current element
          nzheads,   # list of non-zero heads
          row,       # the row of current interest
          inv,       # inverse of a matrix entry
          temp;

    nrows:= Length( mat );
    ncols:= Length( mat[1] );
    nzheads := [];
    vectors := [];
    
    i := 1;

    while i <= nrows do
    
        if ForAll(mat[i], x -> x = 0) then 
        
            i := i + 1;
        
        else
    
            if not i in nzheads then 
            
                # Reduce the row with the known basis vectors.
                
                for j in [ 1 .. Length(nzheads) ] do
                    x := mat[i][ncols + 1 - nzheads[j]];
                    if x <> 0 then
                      mat[i] := mat[i] - mat[ nzheads[j] ]*x;
                    fi;
                od;

                j := PositionNot( Reversed(mat[i]), 0 );
                
                if j <= nrows and j < ncols + 1 then

                    # We found a new basis vector.

                    mat[i] := mat[i]/mat[i][ncols + 1 - j] ;
                    
                    if j = i then 
                        
                        temp := ShallowCopy(mat[i]);
                        
                        i := i + 1;
                        
                    elif j > i then 
                    
                        # Swap rows i and j 
                    
                        temp := ShallowCopy(mat[i]);
                        mat[i] := ShallowCopy(mat[j]);
                        mat[j] := ShallowCopy(temp);
                        
                    elif j < i then 
                    
                        # Swap rows i and n - j 
                        
                        temp := ShallowCopy(mat[i]);
                        mat[i] := ShallowCopy(mat[j]);
                        mat[j] := ShallowCopy(temp);
                        
                        i := i + 1;
                    fi;

                    Add( nzheads, j );
                    Add( vectors, temp );
                else;
                    i := i + 1;
                fi;
            else
                i := i + 1;
            fi;
        fi;
    od;
    
    for i in [1..nrows] do 
        for j in [i + 1..nrows] do
            mat[i] := mat[i] - mat[j]*mat[i][ncols + 1 - j];
        od; 
    od;
    
    end );
    
InstallGlobalFunction(_FoldList2,
    function(list, func, op)
    local k, s, old_s, r, i, len, n, nh, res, r1, r2;


    len := Length(list);
    if len = 0 then
        return 1;
    elif len = 1 then
        return list[1];
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

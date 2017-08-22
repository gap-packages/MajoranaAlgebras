# Takes input matrix, returns a matrix whose rows form a basis of the nullspace of mat

InstallGlobalFunction(MAJORANA_NullSpace,

    function(mat) 

        local   A,      # input matrix
                B;      # output matrix

        A := ShallowCopy(mat);
        
        A := NullspaceMat(TransposedMat(A));
        
        if A <> [] then 
        
            B := List( A, ShallowCopy );
        
            MAJORANA_ReversedEchelonForm(B);
            
            return [[],B];
            
        else
            return [[],[]];            
        fi;

        end

        );
        
InstallGlobalFunction(MAJORANA_SolutionMatVecs1,
    
    function(mat,vec)
    
    local   A,
            B,
            C,
            id,
            error,
            m,
            n,
            dim,
            i,
            j,
            list,
            pos_1,
            pos_2,
            sol,
            unsolved;
    
    m := Size(mat);
    n := Size(mat[1]);
    
    dim := Size(vec[1]);
    
    id := IdentityMat(m);
    
    A := NullMat(n,n);
    B := NullMat(n,dim);
    C := NullMat(n,m);
    
    # Put matrix A into row echelon form
    
    for i in [1..n] do
    
        list := [];
        
        for j in [1..m] do
            Add(list,mat[j][i]);
        od;
        
        pos_1 := Position(list,1);
        
        if pos_1 <> fail then 
        
            A[i] := ShallowCopy(mat[pos_1]);
            B[i] := ShallowCopy(vec[pos_1]);
            C[i] := ShallowCopy(id[pos_1]);
                
            for j in [pos_1 + 1 .. m] do 
                if mat[j][i] <> 0 then 
                
                    mat[j] := mat[j] - mat[pos_1];
                    vec[j] := vec[j] - vec[pos_1];
                    id[j] := id[j] - id[pos_1];
                    
                    if ForAny(mat[j], x -> x <> 0) then 
                        pos_2 := PositionNonZero(mat[j]);
                        
                        vec[j] := vec[j]/mat[j][pos_2];
                        id[j]  := id[j]/mat[j][pos_2]; 
                        mat[j] := mat[j]/mat[j][pos_2]; # change leading elt to 1
                    fi;
                fi;
            od;   
            
            mat[pos_1] := [1..n]*0;
            vec[pos_1] := [1..dim]*0; 
            id[pos_1] := [1..m]*0;           
        fi;
    od;
    
    # Check if we can solve system
    
    error := [];
    
    for i in [1..m] do
        if ForAll(mat[i], x -> x = 0 ) and ForAny(vec[i], x -> x <> 0) then 
            Add(error, i);
        fi;
    od;
    
    if Size(error) >0 then
        # no solutions
        return [error,C,mat,vec];
    fi;
    
    Unbind(mat);
    Unbind(vec);
    Unbind(C);
    
    # solve system
     
    sol := [1..n]*0;
    unsolved := [];
    
    Error("Pause 3");

    for i in Reversed([1..n]) do
        if i in unsolved then
            
            for j in [1..i - 1] do
                if A[j][i] <> 0 then 
                    Add(unsolved,j);
                    sol[j] := [];
                fi;
            od;
            
        elif A[i][i] = 0 then 
        
            sol[i] := [];
            Append(unsolved,[i]);
            
            for j in [1..i - 1] do
                if A[j][i] <> 0 then 
                    Add(unsolved,j);
                    sol[j] := [];
                fi;
            od;
            
        else
            list:=[];
            
            j:= i + 1;
            
            sol[i] := B[i];
            
            for j in [1 .. i - 1] do
                if A[j][i] <> 0 then 
                    B[j] := B[j] - B[i]*A[j][i];
                    A[j] := A[j] - A[i]*A[j][i];
                fi;
            od;
        fi;
    od;

    return [sol,unsolved];
    
    end );
    
    
InstallGlobalFunction(MAJORANA_SolutionMatVecs,

    function(mat,vec) # Takes as input two matrices, the second being interpreted as a vector of vectors. Returns a list of size four if system is inconsistent, otherwise returns a list of size 4

        local   m,
                n,
                res,
                sol,
                unsolved,
                heads,
                i,
                j,
                pos,
                relations;
        
        m := Size(mat);
        n := Size(mat[1]);

        res := SemiEchelonMatTransformationDestructive(mat);
        mat := List(res.vectors,ShallowCopy);
        vec := res.coeffs*vec;
        
        sol := [1..n]*0;
        unsolved := [];
        relations := [[],[]];
        
        heads := res.heads;
        
        for i in Reversed([1 .. n]) do 
        
            pos := heads[i];
            
            if pos = 0 then 
                Add(unsolved,i);   
                sol[i] := [];             
            else
                for j in [i + 1 .. n] do
                    if mat[pos][j] <> 0 then 
                        if j in unsolved then
                            Add(relations[1],mat[pos]);
                            Add(relations[2],vec[pos]);
                            Add(unsolved,i);
                            break;
                        else
                            vec[pos] := vec[pos] - mat[pos][j]*sol[j];
                            mat[pos][j] := 0;
                        fi;
                    fi;
                od; 
                
                if not i in unsolved then 
                    sol[i] := vec[pos];
                else
                    sol[i] := [];
                fi;
            fi;
        od;
        
        Display("Solved it!");
        
        return [sol,unsolved,relations];

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

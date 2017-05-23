
#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

# Creates list of indexes of ProductList[1] whose algebra products are not known

BindGlobal( "MAJORANA_ExtractUnknownAlgebraProducts",

function(AlgebraProducts, ProductList)

    local   unknowns,   # list of unknown algebra products
            i,          # loop over ProductList[1]
            j,          # loop over ProductList[1]
            k,          # pair orbit index
            dim;        # size of ProductList[1]
    
    unknowns := [];
    dim := Size(AlgebraProducts[1]);
    
    for i in [1..dim] do
        for j in [i..dim] do 
            
            k := AbsInt(ProductList[3][i][j]);
        
            if AlgebraProducts[k] = false then 
                Add(unknowns,[i,j]);
            fi;
        od;
    od;

    return AsSet(unknowns);
end);

# Table of fusion rules

BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4,     0, 1/32]
            ,[ 1/32, 1/32,  1/32, 1/4 ] ]);

# This creates new eigenvectors through fusion rules

InstallGlobalFunction( MAJORANA_Fusion,

function(a, b, j, AlgebraProducts, EigenVectors, GramMatrix, ProductList)

    local   u,                  # vector with 1 in the j th position
            x,                  # new eigenvector
            y,                  # inner product of u and x
            z,                  # algebra product of u and x
            NewEigenVectors,    # list of new eigenvectors
            k,                  # run through ev_a-eigenvectors
            l,                  # run through ev_b-eigenvectors
            ev_a,               # first eigenvector
            ev_b,               # second eigenvector
            ev,                 # new eigenvector
            pos,                # index of new eigenvector
            dim,                # size of coordinates              
            FusionError;        # list of indexes which do not obey fusion
            
    dim := Size(ProductList[1]);

    ev := MAJORANA_FusionTable[a+1][b+1];
    
    NewEigenVectors := [];
    FusionError := [];

    u := [1..dim]*0;; u[j] := 1;

    ev_a := EigenVectors[j][a];
    ev_b := EigenVectors[j][b];

    # the 1/4,1/4 case is special
    if (a=2) and (b=2) then
        for k in [1..Size(ev_a)] do
            for l in [1..Size(ev_b)] do

                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );

                if x <> false then
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                    if y <> false then
                        x := x - y*u;
                        
                        z := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList);

                        if (z <> false) and (z <> x*0) then
                            Add(FusionError,[j,k,l]);
                        else
                            Add(NewEigenVectors,x);
                        fi;
                    fi;
                fi;
            od;
        od;
    # the 1/32, 1/32 case is even more special    
    elif (a=3) and (b=3) then
        for k in [1..Size(ev_a)] do
            for l in [1..Size(ev_b)] do
                
                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );

                if x <> false then
                    
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                    
                    if y <> false then
                        
                        x := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList );
                        
                        if x <> false then 
                            x := x - y*u;
                            
                            z := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList);
                        
                            if (z <> false) and ( z <> x/4) then
                                Add(FusionError,[j,k,l]);
                            else
                                Add(NewEigenVectors,x);
                            fi;
                        fi;
                    fi;
                    
                fi;
                
            od;
        od;
    else
        for k in [1..Size(ev_a)] do
            for l in [1..Size(ev_b)] do

                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );
                
                if x <> false then 
                
                    z := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList );

                    if (z <> false) and (z <> ev * x) then
                        Add(FusionError,[j,k,l]);
                    else
                        Add(NewEigenVectors, x);
                    fi;
                fi;
            od;
        od;
    fi;

    if Size(FusionError) > 0 then
        return [false, FusionError];
    else
        pos := Position(MAJORANA_FusionTable[1],ev) - 1;
        return [true, NewEigenVectors, pos];
    fi;
end);        
            

InstallGlobalFunction(MAJORANA_FullFusion,

    function(AlgebraProducts,EigenVectors, GramMatrix, ProductList)
    
    local   j,                  # loop over T
            k,                  # loop over pairs of eigenvalues
            ev,                 # a pair of eigenvalues
            x,                  # result of fusion
            new,                # new eigenvectors
            switch;             # 1 if new vectors have been found
            
    switch := 0;

    for j in ProductList[10] do
        
        new := [ [], [], [] ];
        
        for k in [[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]] do
        
            ev := [,];
        
            ev[1] := MAJORANA_FusionTable[1][k[1] + 1];
            ev[2] := MAJORANA_FusionTable[1][k[2] + 1];

            x := MAJORANA_Fusion(k[1],k[2],j,AlgebraProducts,EigenVectors, GramMatrix, ProductList);
            
            if x[1] then
                Append(new[x[3]], x[2]);
            else
                return [false, 
                        STRINGIFY( "Fusion of ", 
                            ev[1], ",", ev[2], 
                            " eigenvectors does not hold" ),
                        x[2] ];
            fi;
        od;
        
        if ForAny(new, x -> Size(x) > 0) then 
            switch := 1;
        fi;


        for k in [1..3] do
            Append(EigenVectors[j][k],new[k]);
            EigenVectors[j][k]:=ShallowCopy(BaseMat(EigenVectors[j][k]));
        od;        
        
    od;
    
    return [true,switch];
    
    end );

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
        fi;
        
        return B;

        end

        );
        
InstallGlobalFunction(MAJORANA_SolutionMatVecs,
    
    function(mat,vec)
    
    local   A,
            B,
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
    
    A := NullMat(n,n);
    B := NullMat(n,dim);
    
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
                
            for j in [pos_1 + 1 .. m] do 
                if mat[j][i] <> 0 then 
                
                    mat[j] := mat[j] - mat[pos_1];
                    vec[j] := vec[j] - vec[pos_1];
                    
                    if ForAny(mat[j], x -> x <> 0) then 
                        pos_2 := PositionNonZero(mat[j]);
                        
                        vec[j] := vec[j]/mat[j][pos_2];
                        mat[j] := mat[j]/mat[j][pos_2]; # change leading elt to 1
                    fi;
                fi;
            od;   
            
            mat[pos_1] := [1..n]*0;
            vec[pos_1] := [1..dim]*0;            
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
        return [error,,A,B];
    fi;
    
    Unbind(mat);
    Unbind(vec);
    
    # solve system
     
    sol := [1..n]*0;
    unsolved := [];

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
        
InstallGlobalFunction(MAJORANA_Append,

    function(x,mat,vec)

    local   i,          # loop over size of x
            pos;        # position of first non zero elt of row
    
    for i in [1..Size(x[1])] do
        if not x[1][i] in mat then 
            pos := PositionNonZero(x[1][i]);
            
            Add(mat,x[1][i]/x[1][i][pos]);
            Add(vec,x[2][i]/x[1][i][pos]);
            
        fi;
    od;
    
    end);

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
    
InstallGlobalFunction( MAJORANA_ConjugateVector,

    function(v,g,ProductList)
    
    local   i,      # loop over vector
            dim,    # length of vector
            vec,    # output vector
            pos_1,  # position of conjugated element in longcoords
            pos_2;  # position of conjugated element in coords
            
    dim := Size(v);
    
    vec := [1..dim]*0;
    
    for i in [1..dim] do 
        
        pos_1 := Position(ProductList[2],ProductList[1][i]^g);
        pos_2 := ProductList[5][pos_1];
        
        if pos_2 > 0 then 
            vec[pos_2] := v[i];
        else
            vec[-pos_2] := -v[i];
        fi;
    od;
    
    vec := MAJORANA_RemoveNullSpace(vec, ProductList[6]);
    
    return vec;
    
    end );

InstallGlobalFunction(  MAJORANA_AlgebraProduct,

        function(u,v,AlgebraProducts,list) # If all the relevant products are known, returns the algebra product of u and v. If not, returns 0

        # list should be of the form [ProductList[1],ProductList[2],ProductList[3],ProductList[4],ProductList[5],ProductList[6],ProductList[7],ProductList[9]]

        local   i,      # loop over u 
                j,      # loop over v
                k,      # pair orbit index
                x,      # algebra product
                g,      # conjugating element
                sign,   # correct sign of 5A axes
                vec,    # output vec
                dim;    # size of vectors 

        vec:=[];
        dim:=Size(u);
        vec:=[1..dim]*0;

        if ForAll(u,x-> x= 0 ) or ForAll(v,x->x=0) then
            return u*0;
        fi;

        for i in [1..dim] do
            if u[i] <> 0 then 
                for j in [1..dim] do
                    if v[j] <> 0 then 
                    
                        k := list[3][i][j];
                        
                        if k > 0 then 
                            sign := 1;
                        else
                            sign ::= -1;
                        fi;
                        
                        x := AlgebraProducts[AbsInt(k)];
                        
                        if x <> false then
                            
                            g := list[4][i][j];
                            
                            vec := vec + sign*u[i]*v[j]*MAJORANA_ConjugateVector(x,g,list);
                        else
                            if u[i] <> 0 and v[j] <> 0 then
                                # cannot calculate product
                                return false;
                            fi;
                        fi;
                    fi;
                od;
            fi;
        od;
        
        vec := MAJORANA_RemoveNullSpace(vec, list[6]);
                
        return vec;
        
        end

        );

InstallGlobalFunction(  MAJORANA_InnerProduct,

    function(u, v, GramMatrix, ProductList) # If all the relevant products are known, returns the algebra product of u and v. If not, returns [0]

        local   i,      # loop over u 
                j,      # loop over v
                k,      # pair orbit index
                sign,   # correct for 5A axes
                sum;    # output value

        sum := 0;

        for i in [1..Size(u)] do
            if u[i] <> 0 then
                for j in [1..Size(v)] do
                    if v[j] <> 0 then
                    
                        k := ProductList[3][i][j];
                        
                        if k > 0 then 
                            sign := 1;
                        else
                            sign := -1;
                        fi;
                        
                        if GramMatrix[AbsInt(k)] <> false then
                            sum := sum + sign*u[i]*v[j]*GramMatrix[AbsInt(k)];
                        else
                            # cannot calculate product
                            return false;
                        fi;
                    fi;
                od;
            fi;
        od;
        
        return sum;
        
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
        
# Checks if bilinear and algebra products obey axiom M1, outputs a list which is empty if they do obey the axiom

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(GramMatrix,AlgebraProducts,list) 

    # list should be of the form [ProductList[1],ProductList[2],ProductList[3],ProductList[4],ProductList[5]]

        local   ErrorM1,    # list of indices which do not obey axiom M1
                j,          # loop over algebra products
                k,          # loop over ProductList[1]
                p,          # second product
                dim,        # size of ProductList[1]
                x,          # first inner product
                y,          # second inner product
                u,          # vectors
                w,          #
                v;          #

        dim:=Size(AlgebraProducts[1]);

        ErrorM1:=[];
        
        for j in [1..Size(AlgebraProducts)] do
            if AlgebraProducts[j] <> false then
                for k in [1..dim] do 
                    
                    u := NullMat(1,dim)[1];
                    u[list[7][j][1]] := 1;
                    
                    v := NullMat(1,dim)[1];
                    v[list[7][j][2]] := 1;
                    
                    w := NullMat(1,dim)[1];
                    w[k] := 1;
                    
                    p := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,list);
                    
                    if p <> false then
                        x := MAJORANA_InnerProduct(u,p,GramMatrix, list);
                        y := MAJORANA_InnerProduct(AlgebraProducts[j],w,GramMatrix, list);
                        
                        if x <> false and y <> false and x <> y then 
                            Add(ErrorM1,[j,k]);
                        fi;
                        
                    fi;
                od;
            fi;
        od;
                    
        return ErrorM1;

        end

        );
        
# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(GramMatrix,AlgebraProducts,EigenVectors,ProductList) 
        
    # list should be of the form [ProductList[1],ProductList[2],ProductList[3],ProductList[4],ProductList[5],ProductList[6]]
        
        local   errorfusion,    # list of indices which do not obey fusion rules
                dim,            # size of ProductList[1]
                a,              # first eigenvalue
                b,              # second eigenvalue
                ev_a,           # a - eigenvectors
                ev_b,           # b - eigenvectors
                ev,             # new eigenvalue
                u,              # vector with 1 in i th position
                j,              # loop over T 
                v,              # a - eigenvector
                w,              # b - eigenvector
                x,              # product of eigenvectors
                y,              # product of x with u
                z,              # inner product where needed
                x0;             # further product in 1/32 case

        errorfusion:=[];

        dim := Size(AlgebraProducts[1]);

        for j in ProductList[10] do

            u := [1..dim]*0; u[j]:=1;
            
            for a in [1..3] do 
                for b in [a..3] do 
                    if not [a,b] in [[2,2],[3,3]] then 
                        
                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := MAJORANA_FusionTable[a + 1][b + 1];

                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then
                                
                                    y:=MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                    
                                    if y <> false and y <> ev * x then 
                                        Add(errorfusion,[j,a,b,v,w]);
                                    fi;
                                fi;
                            od;
                        od;
                        
                    # 1/4,1/4 fusion is special
                        
                    elif [a,b] = [2,2] then 

                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := 0; 
                        
                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then
                                    
                                    z := MAJORANA_InnerProduct(u,x,GramMatrix, ProductList);
                                        
                                    if z <> false then 
                                        
                                        x := x - z*u;
                                        y := MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                        
                                        if y <> false and y <> ev * x then 
                                            Add(errorfusion,[j,a,b,v,w]);
                                        fi;
                                        
                                    fi;
                                fi;
                            od;
                        od;
                    
                    # 1/32,1/32 fusion is even more special
                        
                    else
                    
                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := 0; 
                        
                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then 
                                
                                    x0 := MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                    
                                    if x0 <> false then 
                                    
                                        y := MAJORANA_InnerProduct(u,x,GramMatrix,ProductList);
                    
                                        if y <> 0 then 
                                            
                                            x0 := x0 - y*u;
                                            
                                            z := MAJORANA_AlgebraProduct(u,x0,AlgebraProducts,ProductList);
                                            
                                            if (z <> false) and (z <> x0/4) then  
                                            
                                                Add(errorfusion,[j,a,b,v,w]);
                                            
                                            fi;
                                        fi;
                                    fi;
                                fi;
                            od;
                        od;
                    fi;
                od;
            od;
        od;
        
        return errorfusion;
        
        end
        
        );

InstallGlobalFunction(MAJORANA_TestOrthogonality,

    function(GramMatrix,AlgebraProducts,EigenVectors, ProductList) # Tests that eigenspaces are orthogonal with respect to the inner product

        local   errorortho, # list of indices which do not obey orthogonality of eigenvectors
                u,          # vector with 1 in j th position
                a,          # first eigenvalue
                b,          # second eigenvalue
                ev_a,       # list of a - eigenvectors
                ev_b,       # list of b - eigenvectors
                j,          # loop over T
                v,          # a - eigenvector
                w,          # b - eigenvector
                x;          # inner product
        
        errorortho := [];

        for j in ProductList[10] do

            u := [1..Size(AlgebraProducts[1])]*0; u[j]:=1;
            
            for a in [1..3] do 
            
                # orthogonality with 1-eigenvectors
                
                ev_a := EigenVectors[j][a];
                
                for v in ev_a do
                    x := MAJORANA_InnerProduct(u, v, GramMatrix, ProductList);
                    
                    if (x <> false) and (x <> 0) then 
                        Add(errorortho, [j,1,a,u,v]);
                    fi;
                od;
                
                # orthogonality with all other eigenvectors
                
                for b in [a+1..3] do 
                
                    ev_b := EigenVectors[j][b];
                    
                    for v in ev_a do
                        for w in ev_b do
                            x := MAJORANA_InnerProduct(v, w, GramMatrix, ProductList);
                            
                            if (x <> false) and (x <> 0) then 
                                Add(errorortho, [j,a,b,u,v]);
                            fi;
                        od;
                    od;
                od;
            od;
        od;
        
        return errorortho;
        
        end
        
        );        

InstallGlobalFunction(MAJORANA_AxiomM2,

        function(GramMatrix,AlgebraProducts,ProductList) # Tests that the algebra obeys axiom M2

        local   B,      # matrix of inner products
                dim,    # size of ProductList[1]
                j,      # loop through ProductList[1]
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
                x3;     #

        dim:=Size(AlgebraProducts[1]);

        B:=NullMat(dim^2,dim^2);

        for j in [1..dim] do
            for k in [1..dim] do
                for l in [1..dim] do
                    for m in [1..dim] do
                        
                        a := [1..dim]*0; a[j] := 1; 
                        b := [1..dim]*0; b[k] := 1;
                        c := [1..dim]*0; c[l] := 1;
                        d := [1..dim]*0; d[m] := 1;
                        
                        x0 := MAJORANA_AlgebraProduct(a,c,AlgebraProducts,ProductList);
                        x1 := MAJORANA_AlgebraProduct(b,d,AlgebraProducts,ProductList);
                        x2 := MAJORANA_AlgebraProduct(b,c,AlgebraProducts,ProductList);
                        x3 := MAJORANA_AlgebraProduct(a,d,AlgebraProducts,ProductList);
                    
                        B[dim*(j-1) + k][dim*(l-1) +m]:=
                              MAJORANA_InnerProduct(x0,x1,GramMatrix, ProductList)
                            - MAJORANA_InnerProduct(x2,x3,GramMatrix, ProductList);
                    od;
                od;
            od;
        od;
        
        return MAJORANA_PositiveDefinite(B);

        end

        );

InstallGlobalFunction(MAJORANA_FillGramMatrix,

function(GramMatrix, ProductList)

    local   i,                  # loop over ProductList[1]
            j,                  # loop over ProductList[1]
            k,                  # pair orbit index
            dim,                # size of ProductList[1]
            GramMatrixFull;     # output matrix

    dim := Size(ProductList[1]);
    
    GramMatrixFull := NullMat(dim,dim);
    
    for i in [1..dim] do 
        for j in [1..dim] do
            
            k := ProductList[3][i][j];
            
            if k > 0 then 
                GramMatrixFull[i][j] := GramMatrix[k];
            else
                GramMatrixFull[i][j] := -GramMatrix[-k];
            fi;
        od;
    od;

    return GramMatrixFull;

    end

    );
    
InstallGlobalFunction(MAJORANA_SeparateInnerProduct,

    function(u,v,UnknownInnerProducts,GramMatrix,ProductList)

    local   row,        # record values of unknowns 
            sum,        # record values of knowns
            dim,        # size of coordinates
            i,          # index for dim of u
            j,          # index for dim of v
            m,          # orbit of i,j
            pos,        # position of m in unknowns
            sign;       # correct sign of 5A axes
            
    dim := Size(ProductList[1]);
            
    sum := 0;
    row := [1..Size(UnknownInnerProducts)]*0;

    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then 
                
                    m := ProductList[3][i][j];
                    
                    if m > 0 then 
                        sign := 1;
                    else
                        sign := -1;
                    fi;
                    
                    m := AbsInt(m);

                    if GramMatrix[m] <> false then
                        sum := sum - sign*u[i]*v[j]*GramMatrix[m];
                    else
                        pos := Position(UnknownInnerProducts,m);
                        row[pos] := row[pos] + sign*u[i]*v[j];
                    fi;
                fi;
            od;
        fi;
    od;

    return [row,sum];

    end );


InstallGlobalFunction(MAJORANA_Orthogonality,

function(a,b,i,UnknownInnerProducts, EigenVectors, GramMatrix, ProductList)

    local   mat,                    # matrix of unknowns 
            vec,                    # vector of knowns
            ev_a,                   # a - eigenvectors
            ev_b,                   # b - eigenvectors
            u,                      # a - eigenvector
            v,                      # b - eigenvector
            x,                      # result of separate inner product 
            dim,                    # size of coordinates
            OrthogonalityError;     # list of vectors which do not obey orthogonality
            
    dim := Size(ProductList[1]);

    mat := [];
    vec := [];

    OrthogonalityError := [];

    if a = 0 then
    
        ev_b := EigenVectors[i][b];
        u := [1..dim]*0; u[i] := 1;

        for v in ev_b do

            x := MAJORANA_SeparateInnerProduct(u,v,UnknownInnerProducts,GramMatrix,ProductList);

            if ForAll(x[1], y -> y = 0) then
                if x[2] <> 0 then
                    Add(OrthogonalityError,[u,v]);
                fi;
            else
                Add(mat,x[1]);
                Add(vec,[x[2]]);
            fi;
        od;
    else
    
        ev_a := EigenVectors[i][a];
        ev_b := EigenVectors[i][b];
        
        for u in ev_a do
            for v in ev_b do

                x := MAJORANA_SeparateInnerProduct(u,v,UnknownInnerProducts,GramMatrix,ProductList);

                if ForAll(x[1], y -> y = 0) then
                    if x[2] <> 0 then
                        Add(OrthogonalityError,[u,v]);
                    fi;
                else
                    Add(mat,x[1]);
                    Add(vec,[x[2]]);
                fi;
                
            od;
        od;
    fi;
    
    if Size(OrthogonalityError) > 0 then 
        return [false, OrthogonalityError];
    else
        return [true,[mat,vec]];
    fi;
    
    end

    );
    
InstallGlobalFunction(MAJORANA_FullOrthogonality,

    function(unknowns,EigenVectors,GramMatrix, ProductList)
    
    local   i,          # loop over T
            j,          # loop over eigenvalues
            k,          # loop over eigenvalues
            x,          # res of orthogonality
            mat,        # matrix of unknown values
            vec,        # vector of known values
            ev;         # pair of eigenvalues          
            
    mat := [];
    vec := [];

    for i in ProductList[10] do
        for j in [0..3] do 
            for k in [j+1..3] do

                x := MAJORANA_Orthogonality(j,k,i,unknowns,EigenVectors,GramMatrix, ProductList);
                
                if x[1] then 
                    MAJORANA_Append(x[2],mat,vec);
                else
                    ev := [,];
                
                    ev[1] := MAJORANA_FusionTable[1][j + 1];
                    ev[2] := MAJORANA_FusionTable[1][k + 1];
                    
                    return [false, 
                            STRINGIFY( "Orthogonality of "
                                , ev[1], ",", ev[2] 
                                , " eigenvectors does not hold"),
                            x[2]];
                fi;
            od;
        od;
    od;
    
    return [true, mat, vec];
    
    end );

InstallGlobalFunction(MAJORANA_MoreEigenvectors,

function(EigenVectors, AlgebraProducts, ProductList)

    local   i,          # loop over representatives
            j,          # loop over coordinates
            u,          # vector with 1 in i th position
            v,          # vector with 1 in j th position
            x,          # algebra product of u,v
            dim,        # size of coordinates
            mat,        # matrix of algebra products
            switch;     # 1 if new eigenvectors have been found
            
    dim := Size(ProductList[1]);
    
    switch := 0;
            
    for i in ProductList[10] do
    
        if Size(EigenVectors[i][1])+Size(EigenVectors[i][2])+Size(EigenVectors[i][3]) + 1 <> dim then
        
            u := [1..dim]*0; u[i] := 1;
            
            mat := [];

            for j in [1..dim] do
            
                v := [1..dim]*0; v[j] := 1;
                
                x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
                
                if x <> false then
                    Add(mat,x);
                else
                    mat := [];
                    break;
                fi;
            od;
            
            if mat <> [] then 
            
                switch := 1;

                EigenVectors[i][1]:=ShallowCopy(NullspaceMat(mat));
                EigenVectors[i][2]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim)/4));
                EigenVectors[i][3]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim)/32));
                EigenVectors[i][4]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim) ));

                if ProductList[6] <> [] and ProductList[6] <> false then 
                    for j in [1..4] do 
                        for x in EigenVectors[i][j] do
                            x := MAJORANA_RemoveNullSpace(x,ProductList[6]);
                        od;
                    od;
                fi;
                                       
                if Size(EigenVectors[i][4]) <> 1 and ProductList[6] <> false then
                
                    return [false, 1, i];
                    
                elif Size(EigenVectors[i][1]) + Size(EigenVectors[i][2]) + Size(EigenVectors[i][3]) + Size(EigenVectors[i][4]) > dim 
                    and ProductList[6] <> false then
                    
                    return [false, 2, i];
                fi;
            fi;
        fi; 
    od;
    
    return [true,switch];
    
    end );
    
InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(EigenVectors, AlgebraProducts, ProductList)

    local   i,          # loop over representatives
            ev,         # loop over eigenvalues
            unknowns,   # unknown algebra products
            mat,        # matrix of unknowns
            vec,        # vector of knowns
            table,      # table of eigenvalues
            u,          # vector with 1 in j th position
            v,          # eigenvector
            x,          # result of SeparateAlgebraProduct
            dim;        # size of ProductList[1]
    
    table := [0, 1/4, 1/32];
    
    dim := Size(AlgebraProducts[1]);
    
    for i in ProductList[10] do
    
        unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
                
        unknowns := Filtered( unknowns, x -> x[1] = i);
        
        if unknowns <> [] then
        
            mat := [];
            vec := [];
        
            for ev in [1..3] do 
            
                u := [1..dim]*0; u[i] := 1;
                
                for v in EigenVectors[i][ev] do
                    
                    x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,AlgebraProducts,ProductList);
                    
                    if ForAll(x[1], y -> y = 0) then 
                        if ForAny((x[2] + table[ev]*v), y -> y <> 0 ) then 
                            return [false,1,v];
                        fi;
                    elif not x[1] in mat then         
                        Add(mat, x[1]);
                        Add(vec, x[2] + table[ev]*v);       
                    fi;
                od;
            od;
        
            if mat <> [] then 
                x := MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
                
                if not x[1] and ProductList[6] <> false then 
                    return [false,2,x];
                fi;
            fi;
        fi;
    od;
            
    return [true];
    
    end);
    
InstallGlobalFunction(MAJORANA_RemoveNullSpace,

function(v,NullSp) 

    local   i,      # loop over nullspace
            j,      # leading coefficient (from rhs)
            dim;    # size of coefficients
    
    dim := Size(v);


    if NullSp <> [] and NullSp <> false then 
        for i in [1..Size(NullSp)] do
            j := Position(Reversed(NullSp[i]),1);
            if v[dim - j + 1] <> 0 then 
                v := v - v[dim - j + 1]*NullSp[i];
            fi;
        od;
    fi;
    
    return v;
    
    end
    
    );

InstallGlobalFunction(MAJORANA_UnknownsAxiomM1,
    
    function(k,l, GramMatrix, AlgebraProducts, ProductList)
    
    local j,dim,res,x,y,z,u,v,w;

    j := 1;
    
    dim := Size(AlgebraProducts[1]);
    
    res := false;
    
    while j < Size(AlgebraProducts) + 1 do 
    
        if AlgebraProducts[j] <> false and AlgebraProducts[j][k] <> 0 then 
        
            u := [1..dim]*0; u[l] := 1;
            v := [1..dim]*0; v[ProductList[7][j][1]] := 1;
        
            x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
            if x <> false then
                
                u := [1..dim]*0; u[ProductList[7][j][2]] := 1;
            
                y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                
                if y <> false then 
                    
                    v := StructuralCopy(AlgebraProducts[j]); v[k] := 0;
                    
                    w := [1..dim]*0; w[l] := 1;
                    
                    z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList);
                    
                    if z <> false then 
                                           
                        res := (y - z)/AlgebraProducts[j][k]; 
                                          
                    fi;
                fi;
                
            fi;
            
            if res <> false then 
                
                j := Size(AlgebraProducts) +  1;
            
            else
            
                u := [1..dim]*0; u[l] := 1;
                v := [1..dim]*0; v[ProductList[7][j][2]] := 1;
            
                x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
                if x <> false then
                    
                    u := [1..dim]*0; u[ProductList[7][j][1]] := 1;
                
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                    
                    if y <> false then 
                        
                        v := StructuralCopy(AlgebraProducts[j]); v[k] := 0;
                        
                        w := [1..dim]*0; w[l] := 1;
                        
                        z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList);
                        
                        if z <> false then 
                                               
                            res := (y - z)/AlgebraProducts[j][k];  
                                              
                        fi;
                    fi;
                    
                fi;
            fi;
            
            if res <> false then 
                j := Size(AlgebraProducts) + 1;
            else
                j := j + 1;
            fi;
            
        else
        
            j := j + 1;
            
        fi;
    od;
    
    return res;
    
    end );
    
InstallGlobalFunction(MAJORANA_FullUnknownsAxiomM1,

function(GramMatrix,AlgebraProducts,ProductList)

    local   switch,         # switch for while loop
            count,          # counts any newly found products
            i,              # loop over orbitals
            j,              # loop over an orbital 
            k,              # position of first element
            l,              # position of second element
            res,            # result of unknowns axiom M1
            sign;           # correct sign of 5A axes
    
    switch := 1;
            
    while switch = 1 do
        
        count := 0;

        for i in [1..Size(ProductList[9])] do
            if GramMatrix[i] = false then
            
                j := 1;
    
                while j < Size(ProductList[9][i]) + 1 do 
                    
                    k := ProductList[5][Position(ProductList[2],ProductList[9][i][j][1])];
                    l := ProductList[5][Position(ProductList[2],ProductList[9][i][j][2])];
                    
                    if k*l > 0 then 
                        sign := 1;
                    else
                        sign := -1;
                    fi;
                    
                    res := MAJORANA_UnknownsAxiomM1(AbsInt(k),AbsInt(l),GramMatrix,AlgebraProducts,ProductList);
                    
                    if res = false then 
                        j := j + 1;
                    else
                        GramMatrix[i] := sign*res;
                        count := count + 1;
                        
                        j := Size(ProductList[9][i]) + 1;
                    fi;        
                od;

            fi;
        od;
        
        if count = 0 then 
            switch := 0;
            break;
        fi;
    od; 

    end );
    
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
    
InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,UnknownAlgebraProducts,AlgebraProducts,ProductList)
    
    local   row,        # record values of unknowns
            sum,        # record values of knowns
            i,          # index for dim of u
            j,          # index for dim of v
            l,          # ordered version of [i,j]
            x,          # vector with 1 in the ith position
            y,          # vector with 1 in the jth position
            pos,        # position of unknown product 
            dim;        # dimension
    
    row := [1..Size(UnknownAlgebraProducts)]*0;
    sum := 0;
    
    dim := Size(AlgebraProducts[1]);
    
    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then
                
                    l := [Minimum([i,j]),Maximum([i,j])];
                    
                    if not l in UnknownAlgebraProducts then 
                        
                        x := [1..dim]*0; x[i] := 1;
                        y := [1..dim]*0; y[j] := 1;
                        
                        sum := sum - u[i]*v[j]*MAJORANA_AlgebraProduct(x,y,AlgebraProducts,ProductList);
                    else
                        
                        pos := Position(UnknownAlgebraProducts,l);
                        row[pos] := row[pos] + u[i]*v[j]; 
                    fi;
                fi;
            od;
        fi;
    od;
       
    return [row,sum];
    
    end);
    
InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns, ProductList)
    
    local   output,     # output row
            len,        # length of row
            i,          # loop over length of row
            j,          # first elt of original product
            k,          # second elt of original product
            x,          # new product
            pos_1,      # position of first element conjugated
            pos_2,      # position of second element conjugated
            sign,       # corrects sign of 5A axis
            pos;        # position of new product
            
    
    len := Size(row);
    output := [1..len]*0;
    
    for i in [1..Size(row)] do
        if row[i] <> 0 then 
    
            j := unknowns[i][1];
            k := unknowns[i][2];
            
            x := [0,0];
            
            pos_1 := Position(ProductList[2],(ProductList[1][j])^g);
            pos_2 := Position(ProductList[2],(ProductList[1][k])^g);
            
            x[1] := ProductList[5][pos_1];
            x[2] := ProductList[5][pos_2];
            
            if x[1]*x[2] < 0 then 
                sign := -1;
            else
                sign := 1;
            fi;
            
            x[1] := AbsInt(x[1]);
            x[2] := AbsInt(x[2]);
            
            Sort(x);
            
            pos := Position(unknowns,x);
            output[pos] := sign*row[i];
            
        fi;
    od;
    
    return output;
    
    end);     
    
InstallGlobalFunction(MAJORANA_Resurrection,

    function(i,ev_a, ev_b, EigenVectors,UnknownAlgebraProducts,AlgebraProducts,ProductList,GramMatrix)
    
    local a, b, sum, row, alpha, beta, gamma, mat, vec, bad, list, j, m, k, l, dim, u, v, x, y, z, g, sign, record;
    
    dim := Size(AlgebraProducts[1]);
    
    mat := [];
    vec := [];
    record := [];
    
    for a in [1..Size(EigenVectors[i][ev_a])] do 
        for b in [1..Size(EigenVectors[i][ev_b])] do 
        
            if MAJORANA_AlgebraProduct(EigenVectors[i][ev_a][a],EigenVectors[i][ev_b][b],AlgebraProducts,ProductList) = false then 
        
                row := [1..Size(UnknownAlgebraProducts)]*0;
                sum := [];
                
                beta := EigenVectors[i][ev_a][a];
                gamma := EigenVectors[i][ev_b][b];
                
                list := [];
                bad := [];
                
                x := MAJORANA_SeparateAlgebraProduct(MAJORANA_FusionTable[ev_a + 1][ev_b + 1]*beta,gamma,UnknownAlgebraProducts,AlgebraProducts,ProductList);
                
                row := row + x[1];
                sum := sum + x[2];
                   
                # find a suitable alpha   
                   
                for m in [1..dim] do
                    for k in [1..dim] do 
                        if gamma[k] <> 0 and [Minimum([m,k]),Maximum([m,k])] in UnknownAlgebraProducts then 
                            Add(bad,m);
                        fi;
                    od;
                od;
                
                bad := DuplicateFreeList(bad);
                
                Sort(bad);

                for m in [1..Size(EigenVectors[i][ev_b])] do 
                    Add(list, StructuralCopy(EigenVectors[i][ev_b][m]{bad}));
                od;

                if ev_a = ev_b then 
                    list[a] := [1..Size(bad)]*0;
                fi;
                
                x := SolutionMat(list,beta{bad});
                
                if x <> fail then 
                    
                    alpha := [];
                    
                    for m in [1..Size(x)] do
                        if ev_a <> ev_b or m <> a then
                            alpha := alpha + x[m]*EigenVectors[i][ev_b][m];
                        fi;
                    od;
                
                    y := MAJORANA_AlgebraProduct((alpha - beta),gamma,AlgebraProducts,ProductList);
                    
                    u := [1..dim]*0; u[i] := 1;
                        
                    z := MAJORANA_SeparateAlgebraProduct(u,y,UnknownAlgebraProducts,AlgebraProducts,ProductList);  
                        
                    row := row + z[1];
                    sum := sum + z[2];
                
                    if ev_b = 2 then 
                    
                        u := [1..dim]*0; u[i] := 1; 
                        x := MAJORANA_InnerProduct(alpha,gamma,GramMatrix,ProductList);
                    
                        if x <> false then 
                            sum := sum + u*x/4;
                        else
                            row := [];
                        fi;                            
                    fi;
                
                    if row <> [] then 
                        if ForAll(row, x -> x = 0) then 
                            if ForAny( sum, x -> x <> 0) then
                                Error("Resurrection error");
                            fi;
                        else
                            
                            for g in DuplicateFreeList(ProductList[12]) do
                                if g <> false then 
                                
                                    Add(mat,MAJORANA_ConjugateRow(row,g,UnknownAlgebraProducts,ProductList));
                                    Add(vec,MAJORANA_ConjugateVector(sum,g,ProductList));
                                    
                                    Add(record,[i,ev_a,ev_b,alpha,beta,gamma,g]);
                                fi;
                            od;
                        fi;  
                    fi;
                fi;
            fi;
        od;
    od;
    
    return [mat,vec,record];
        
end );

InstallGlobalFunction(MAJORANA_FullResurrection,

    function(EigenVectors,UnknownAlgebraProducts,AlgebraProducts,ProductList,GramMatrix)
    
    local mat, vec, record, j, k, l, x, g, t;
    
    mat := [];
    vec := [];
    record := [];
    
    t := Size(EigenVectors);
    
    for j in ProductList[10] do 
    
        if EigenVectors[j][1] = false then 
        
            g := ProductList[12][j];
            x := ProductList[13][j];
            
            k := ProductList[10][x];
            
            for l in [1..3] do 
            
                EigenVectors[j][l] := [];
                
                for x in EigenVectors[k][l] do
                    Add(EigenVectors[j][l], MAJORANA_ConjugateVector(x,g,ProductList));
                od;
            od;
        fi;
    
        for k in [1..3] do
            for l in [1..2] do 
                if [k,l] <> [2,2] then 
            
                    x := MAJORANA_Resurrection(j,k,l,EigenVectors,UnknownAlgebraProducts,AlgebraProducts,ProductList,GramMatrix);
                    
                    if x[1] <> [] then 
                        MAJORANA_Append(x,mat,vec);
                        Append(record, x[3]);
                    fi;
                fi;
            od;
        od;
    od;
    
    return [mat,vec,record];
    
    end);
    
InstallGlobalFunction(MAJORANA_NullSpaceAlgebraProducts,

    function(UnknownAlgebraProducts, AlgebraProducts, ProductList)
    
    local i, m, j, k, row, sum, dim, y, mat, vec, a, x, record;
    
    dim := Size(ProductList[6][1]);
    
    mat := [];
    vec := [];
    record := [];
    
    for m in [1..Size(UnknownAlgebraProducts)] do
    
        j := UnknownAlgebraProducts[m][1];
    
        a := [1..dim]*0; a[j] := 1;        
    
        for k in [1..Size(ProductList[6])] do
        
            row := [1..Size(UnknownAlgebraProducts)]*0;
            sum := [];
            
            x := MAJORANA_SeparateAlgebraProduct(a,ProductList[6][k],UnknownAlgebraProducts,AlgebraProducts,ProductList);

            if ForAll(x[1], x -> x = 0) then 
                if ForAny( x[2] , y -> y <> 0) then 
                    Error("Nullspace"); 
                fi;
            else
                Add(mat,x[1]);
                Add(vec,x[2]);
                Add(record,[k,"n"]);
            fi;
        od;
    od;
    
    return [mat,vec,record];
    
    end );
    
InstallGlobalFunction( MAJORANA_PairConjElements,

    function(ProductList)
    
    local   x,      # input elements
            y,      # representative elements
            z,      # index of elements which will also have g
            list,   # list of indices z
            i,      # first basis element
            j,      # second basis element
            k,      # orbital of elements
            l,
            table,
            pos_1,
            pos_2,
            dim,    # size of coordinates
            g;      # conjugating element
            
    table := [[],[1],[1,2],[1,3],[1,2,3,4]];
    
    dim := Size(ProductList[1]);
    
    for i in [1..dim] do
        for j in [1..dim] do
            if ProductList[4][i][j] = 0 then 
            
                x := [ProductList[1][i],ProductList[1][j]];
                y := [0,0];
                
                # create list of other indices which are going to have this elt
                
                list := [];
                
                for k in table[Order(x[1])] do 
                    for l in table[Order(x[2])] do
                        Add(list,[x[1]^k,x[2]^l]);
                    od;
                od;
                
                k := AbsInt(ProductList[3][i][j]);
            
                y := ProductList[9][k][1];
                
                l := 1;
                
                while l  < Size(list) + 1 do 
            
                    g := RepresentativeAction(ProductList[8],y,list[l],OnPairs);
                
                    if g <> fail then
                    
                        for z in list do
                            
                            pos_1 := AbsInt(ProductList[5][Position(ProductList[2],z[1])]);
                            pos_2 := AbsInt(ProductList[5][Position(ProductList[2],z[2])]);
                            
                            ProductList[4][pos_1][pos_2] := g;
                            ProductList[4][pos_2][pos_1] := g;
                        od;
                        
                        l := Size(list) + 1;
                        
                    else 
                             
                        g := RepresentativeAction(ProductList[8],y,Reversed(list[l]),OnPairs);
                    
                        if g <> fail then 
                        
                            for z in list do
                            
                                pos_1 := AbsInt(ProductList[5][Position(ProductList[2],z[1])]);
                                pos_2 := AbsInt(ProductList[5][Position(ProductList[2],z[2])]);
                                
                                ProductList[4][pos_1][pos_2] := g;
                                ProductList[4][pos_2][pos_1] := g;
                            od;
                            
                            l := Size(list) + 1;
                            
                        else
                            l := l + 1;
                        fi;
                    fi;
                od;
            fi;
        od;
    od;
    
    end );
    
InstallGlobalFunction( MAJORANA_PairOrbits,

    function(ProductList)
    
    local   i,    # first basis element
            j,    # second basis element
            dim,  # size of coordinates  
            k;    # loop through orbitals
            
    dim := Size(ProductList[1]);
    
    for i in [1..dim] do
        for j in [1..dim] do
        
            if ProductList[3][i][j] = 0 then 
                
                k := 1;
                
                while k < Size(ProductList[9]) + 1 do
                
                    if [ProductList[1][i],ProductList[1][j]] in ProductList[9][k] then
                    
                        if [ProductList[1][i],ProductList[1][j]] in ProductList[14] then
                        
                            ProductList[3][i][j] := -k;
                            ProductList[3][j][i] := -k;
                        
                        else
                    
                            ProductList[3][i][j] := k;
                            ProductList[3][j][i] := k;
                        fi;
                        
                        k := Size(ProductList[9]) + 1;
                        
                    else
                        k := k + 1;
                    fi;
                od;
            fi;
        od;
    od;
    
        
    end );
    
InstallGlobalFunction(MAJORANA_PairRepresentatives,

    function(ProductList)
    
    local   i,          # loop through orbitals
            x,          # representative elements
            y,          # positions of representatives
            list,       # list of equivalent pairs of elts
            j,          # orders of first element
            k,          # orders of second element
            pos_1,      # position of first element
            pos_2,      # position of second element
            table,      # table of orders of basis elements
            z;          # an equivalent pair of elements       
            
    table := [[],[1],[1,2],[1,3],[1,2,3,4]];

    for i in [1..Size(ProductList[9])] do
            
        x := ProductList[9][i][1];
        y := [Position(ProductList[1],x[1]), Position(ProductList[1],x[2])];
        
        Add(ProductList[7], y);
        
        # create list of other indices which are going to have this elt
                
        list := [];
        
        for j in table[Order(x[1])] do 
            for k in table[Order(x[2])] do
                Add(list,[x[1]^j,x[2]^k]);
            od;
        od;
        
        # put conj elt and pair orbit for representatives 
        
        for z in list do
                            
            pos_1 := AbsInt(ProductList[5][Position(ProductList[2],z[1])]);
            pos_2 := AbsInt(ProductList[5][Position(ProductList[2],z[2])]);
            
            ProductList[4][pos_1][pos_2] := ();
            ProductList[4][pos_2][pos_1] := ();
            
            ProductList[3][pos_1][pos_2] := i;
            ProductList[3][pos_2][pos_1] := i;
        od;
    od;
    
    end);
    
InstallGlobalFunction( MAJORANA_MakeVector,

    function( pos, vals, dim)
    
    local   vec,    # output vector
            i;      # loop over input
            
    vec := [1..dim]*0;;
    
    for i in [1..Size(pos)] do
        vec[pos[i]] := vec[pos[i]] + vals[i];
    od;
    
    return vec;
    
    end );
    
InstallGlobalFunction( MAJORANA_OutputError,

    function(message, error, OutputList)
    
    local   output,         # output vector
            i;              # loop over output list
            
    output := ["Error", message, error];
    
    for i in [1..5] do 
        Add(output, OutputList[i]);
    od;
    
    return StructuralCopy(output);
    
    end );
    
InstallGlobalFunction( MAJORANA_SolutionAlgProducts,

    function( mat, vec, UnknownAlgebraProducts, AlgebraProducts, ProductList)
    
    local   Solution,   # solution of system
            i,          # loop over <UnknownAlgebraProducts>
            x,          # element of <UnknownAlgebraProducts>
            y,          # orbit of x
            g;          # conj element of x
    
    if mat <> [] then
        
        Solution := MAJORANA_SolutionMatVecs(mat,vec);

        if Size(Solution) = 2 then
            for i in [1..Size(UnknownAlgebraProducts)] do
                if not i in Solution[2] then 
                
                    x := UnknownAlgebraProducts[i]; 
                    
                    y := ProductList[3][x[1]][x[2]];
                    g := ProductList[4][x[1]][x[2]];
                                                    
                    if y > 0 then 
                        AlgebraProducts[y] := MAJORANA_ConjugateVector(Solution[1][i],Inverse(g),ProductList);
                    else
                        AlgebraProducts[-y] := -MAJORANA_ConjugateVector(Solution[1][i],Inverse(g),ProductList);
                    fi;
                fi;
            od;
            
            return [true];
            
        else
            return [false, Solution];
        fi;
    fi;
    
    end );
    
InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, UnknownInnerProducts, GramMatrix)
    
    local   Solution,   # solution of system
            i,          # loop over <UnknownInnerProducts>
            x;          # element of <UnknownInnerProducts>    
    
    Solution := MAJORANA_SolutionMatVecs(mat,vec);

    if Size(Solution) = 2 then                    
        
        for i in [1..Size(Solution[1])] do
            if not i in Solution[2] then
    
                x:=UnknownInnerProducts[i]; 

                GramMatrix[x]:=Solution[1][i][1];
            fi;
        od;
        
        if Size(Solution[2]) = Size(Solution[1]) then
            return [false, Solution];
        fi;
        
        return [true];
    else
    
        return [false, Solution];
    fi;
    
    end );
    
InstallGlobalFunction( MAJORANA_MergeOrbitalsAxes,

    function(ProductList, OrbitalsT)
    
    local   i,      # loop over orbitals
            x,      # representative of orbital
            y,      # orders of representatives
            table,  # table of orders of axes
            j,      # loop over orders of first axis
            k,      # loop over orders of second axis
            l;      # loop over orbitals
    
    i := Size(OrbitalsT) + 1;
            
    while i < Size(ProductList[9]) + 1 do
        
        x := ProductList[9][i][1];

        y := [Order(x[1]),Order(x[2])];
        
        table  := [[],[1],[1,2],[1,3],[1,2,3,4]];
        
        for j in table[y[1]] do
            for k in table[y[2]] do  
                      
                if not [x[1]^j, x[2]^k] in ProductList[9][i] then
                    
                    l := i + 1;
                    
                    while l < Size(ProductList[9]) + 1 do
                    
                        if [x[1]^j, x[2]^k] in ProductList[9][l] then 
                        
                            Append(ProductList[9][i],ProductList[9][l]);
                           
                            if y = [5,5] then
                                if [j,k] in [[1,2],[2,1],[1,3],[3,1],[2,4],[4,2],[3,4],[4,3]] then 
                                    Append(ProductList[14], ProductList[9][l]);
                                fi;
                            elif y[2] = 5 and k in [2,3] then 
                                Append(ProductList[14], ProductList[9][l]);
                            fi;
                           
                            Remove(ProductList[9],l);
                           
                            l := Size(ProductList[9]) + 1;
                           
                        else
                            l := l + 1;
                        fi;
                        
                    od;
                fi;
            od;
        od;
        
        i := i + 1;       
        
    od;
        
    end );    
    
InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(GramMatrix,AlgebraProducts,EigenVectors,ProductList)
    
    local   GramMatrixFull,     # full gram matrix
            x,                  # result of positive definite
            i,                  # loop over orbitals
            j,                  # loop over representatives
            k;                  # loop over eigenvalues
    
        if ProductList[6] = false then 
            if ForAll(GramMatrix, x -> x <> false) then 
                GramMatrixFull := MAJORANA_FillGramMatrix(GramMatrix, ProductList);

                x := MAJORANA_PositiveDefinite(GramMatrixFull);

                if x < 0 then
                    return false;
                elif x = 0 then
                    ProductList[6] := MAJORANA_NullSpace(GramMatrixFull);
                fi;
            fi;

            if ProductList[6] <> [] then

                # Change alg products to get rid of any axes not in the basis
                
                for i in [1..Size(ProductList[9])] do
                    if AlgebraProducts[i] <> false then
                        AlgebraProducts[i]:= MAJORANA_RemoveNullSpace(AlgebraProducts[i], ProductList[6]);
                    fi;
                od;

                # Change evecs to get rid of any axes not in the basis

                for j in ProductList[10] do
                    for k in [1..3] do                        
                        for x in [1..Size(EigenVectors[j][k])] do
                            EigenVectors[j][k][x] := MAJORANA_RemoveNullSpace(EigenVectors[j][k][x],ProductList[6]);
                        od;                            
                    od;
                od;
            fi;
        fi;
        
        return true;
    
    end );
    
InstallGlobalFunction(MAJORANA_MainSteps,

    function(i,GramMatrix,AlgebraProducts,EigenVectors,ProductList,Output,OutputList,maindimensions)

    local   x,          # result of calculations
            dim,        # size of coordinates
            error,      # record error
            switch,     # 1 if new values have been found
            falsecount,
            count,
            newfalsecount;
            
    dim := Size(ProductList[1]);
    
    falsecount := [0,0];
    
    falsecount[1] := Size(Positions(AlgebraProducts,false));
    falsecount[2] := Size(Positions(GramMatrix,false));

    count := 0;
    
    if ForAny(falsecount, x -> x > 0) then 
        switch := 1;
    else
        switch := 0;
    fi;
    
    while switch = 1 do 
    
        count := count + 1;

                                ## STEP 5: INNER PRODUCTS M1 ##
                                
        
                                    
        if falsecount[2] > 0 then 
            MAJORANA_FullUnknownsAxiomM1(GramMatrix,AlgebraProducts,ProductList);
        fi;
        
        x := MAJORANA_CheckNullSpace(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
        
        if x = false then
            Output[i] := MAJORANA_OutputError("The inner product is not positive definite"
                                , []
                                , OutputList);                                     
            break;
        fi;
        
                                    ## STEP 6: FUSION ##                                        
                                

        # Use these eigenvectors and the fusion rules to find more

        if ForAny(maindimensions, x -> x < dim - 1) then                
        
            x := MAJORANA_FullFusion(AlgebraProducts,EigenVectors, GramMatrix, ProductList);
            
            if not x[1] and ProductList[6] <> false then 
                Output[i] := MAJORANA_OutputError(x[2],
                                x[3],
                                OutputList);
                break;
            fi;
            
        fi;
        
                                    ## STEP 7: PRODUCTS FROM EIGENVECTORS ##

        # Check fusion and M1

        error := MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,ProductList);

        if Size(error) > 0 and ProductList[6] <> false then
            Output[i] := MAJORANA_OutputError("Algebra does not obey axiom M1"
                            , error
                            , OutputList);
            break;
        fi;

        error := MAJORANA_TestFusion(GramMatrix, AlgebraProducts, EigenVectors,ProductList);

        if Size(error) > 0 and ProductList[6] <> false then
            Output[i] := MAJORANA_OutputError("Algebra does not obey fusion rules"
                         , error
                         , OutputList);
            break;
        fi;

        # Use eigenvectors to find more products
        
        x := MAJORANA_EigenvectorsAlgebraUnknowns(EigenVectors,AlgebraProducts,ProductList);
        
        if not x[1] then 
            if x[2] = 1 then 
                Output[i] := MAJORANA_OutputError( "Error eigenvectors algebra unknowns"
                        , x[3]
                        , OutputList);
                break;
            elif x[2] = 2 then 
                Output[i] := MAJORANA_OutputError("Inconsistent system of unknown algebra products step 7"
                        , x[3]
                        , OutputList);
                break;
            fi;
        fi;
        
                                            ## STEP 8: MORE EVECS II ##

        # Check if we have full espace decomp, if not find it

        x := MAJORANA_MoreEigenvectors(EigenVectors, AlgebraProducts, ProductList);
        
        if not x[1] then 
            if x[2] = 1 then 
                Output[i] := MAJORANA_OutputError("Algebra does not obey axiom M5"
                    , []
                    , OutputList);
                break;
            elif x[2] = 2 then 
                Output[i] := MAJORANA_OutputError("Algebra does not obey axiom M4"
                    , []
                    , OutputList);
                break;
            fi;
        fi;
    
        if  Size(Positions(AlgebraProducts,false)) = falsecount[1] and
            Size(Positions(GramMatrix,false)) = falsecount[2] then
            
            switch := 0;
            
            break;
        else
            falsecount[1] := Size(Positions(AlgebraProducts,false));
            falsecount[2] := Size(Positions(GramMatrix,false)); 
        fi;    
    
    od;
        
    end );
        
InstallGlobalFunction(MajoranaRepresentation,

function(G,T)

    local   # Seress
            ProductList, error, OrbitsT, 

            # indexing and temporary variables
            i, j, k, x, y,

            # Step 0 - Set Up
            Output, t, SizeOrbitals, OrbitalsT, 
            
            # Step 1 - Shape
            Shape, RepsSquares6A, Unknowns3X,

            # Step 2 - Possible shapes
            Binaries, master,

            # Step 3 - Products and evecs I
            GramMatrix, GramMatrixFull, AlgebraProducts, EigenVectors, sign,

            # Step 4 - More products and evecs
            h, s, dim, 

            # Step 6 - More inner products
            unknowns, mat, vec, 
            
            vals, pos, OutputList, record, 

            
            falsecount, newfalsecount, maindimensions, newdimensions, switchmain, count;     
            
 

                                            ## STEP 0: SETUP ##

    Output:=[];

    t:=Size(T);

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;
    
    # Set up ProductList
    
    ProductList := [1..13]*0;
        
    ProductList[8]  := G;
    ProductList[10] := [];
    ProductList[12] := [1..t]*0;
    ProductList[13] := [1..t]*0;
    
    # Construct orbits of G on T
    
    ProductList[11] := OrbitsDomain(G,T);
    
    for x in ProductList[11] do
        Add(ProductList[10], Position( T, Representative(x)));
    od;
    
    for i in [1..t] do
        
        j := 1;
        
        while j < Size(ProductList[11]) + 1 do 
            if T[i] in ProductList[11][j] then 
            
                ProductList[13][i] := j;
                
                k := ProductList[10][j];
                
                ProductList[12][i] := RepresentativeAction(G,T[k],T[i]);
                
                j := Size(ProductList[11]) + 1;;
                
            else
                j := j + 1;
            fi;
        od;
    od;

    # Construct orbitals of G on T x T

    x := OrbitsDomain(G,Cartesian(T,T),OnPairs);
    
    OrbitalsT := [];
    
    for i in [1..Size(x)] do
        Add(OrbitalsT, ShallowCopy(x[i]));
    od;
    
    i := 1;
    
    while i < Size(OrbitalsT) do 
    
        if not [OrbitalsT[i][1][2],OrbitalsT[i][1][1]] in OrbitalsT[i] then
        
            j := i + 1;
            
            while j < Size(OrbitalsT) + 1 do
            
                if  [OrbitalsT[i][1][2],OrbitalsT[i][1][1]]  in OrbitalsT[j] then
                
                    Append(OrbitalsT[i],OrbitalsT[j]);
                    Remove(OrbitalsT,j);
                        
                    j := Size(OrbitalsT) + 1;
                    
                else
                    
                    j := j + 1;
                fi;
            od;
        fi;
        
        i := i + 1;
        
    od;

                                        ## STEP 1: SHAPE ##

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    Shape:=NullMat(1,Size(OrbitalsT))[1];

    RepsSquares6A:=[];

    for i in [1..Size(OrbitalsT)] do
    
        x:=Representative(OrbitalsT[i]);
        
        if Order(x[1]*x[2]) = 1 then
            Shape[i]:="1A";
        fi;
        if Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            Shape[i]:="2A";
        fi;
        if Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            Shape[i]:="2B";
        fi;
        if Order(x[1]*x[2]) = 3 then
            Shape[i]:="3X";
        fi;
        if Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            Shape[i]:="4A";
        fi;
        if Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            Shape[i]:="4B";
        fi;
        if Order(x[1]*x[2]) = 5 then
            Shape[i]:="5A";
        fi;
        if Order(x[1]*x[2])=6 then
            Shape[i]:="6A";
            Add(RepsSquares6A,(x[1]*x[2])^2);
        fi;
    od;

    # Check for inclusions of 3A in 6A

    for i in [1..Size(OrbitalsT)] do
        if Shape[i][1] = '3' then
            j := 0;
            while j < Size(OrbitalsT[i]) do
                j := j + 1;
                if OrbitalsT[i][j][1]*OrbitalsT[i][j][2] in RepsSquares6A then
                    Shape[i]:="3A";;
                    j:=Size(OrbitalsT[i])+1;;
                fi;
            od;
        fi;
    od;

    Unknowns3X:=[];

    for i in [1..Size(OrbitalsT)] do
        if Shape[i] = ['3','X'] then
            Add(Unknowns3X,i);
        fi;
    od;

                                            ## STEP 2: POSSIBLE SHAPES ##

    # Run through possibilities for unknowns

    Binaries:=AsList(FullRowSpace(GF(2),Size(Unknowns3X)));

    for i in [1..Size(Binaries)] do
        
        Output[i]:=[];

        master:=1;
        
        ProductList[6] := false;
        ProductList[1] := StructuralCopy(T);

        while master = 1 do

            # Add new values in the shape

            for j in [1..Size(Unknowns3X)] do
                k:=Unknowns3X[j];
                if Binaries[i][j] = 1*Z(2) then
                    Shape[k]:="3A";
                else
                    Shape[k]:="3C";
                fi;
            od;

            # Create lists of 3A, 4A and 5A axes

            for j in [1..Size(OrbitalsT)] do
                if Shape[j]=['3','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
                        x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(ProductList[1],Set([x,x^2]));
                    od;
                fi;
            od;
            
            for j in [1..Size(OrbitalsT)] do
                if Shape[j]=['4','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
                        x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(ProductList[1],Set([x,x^3]));
                    od;
                fi;
            od;
            
            for j in [1..Size(OrbitalsT)] do 
                if Shape[j]=['5','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
                        x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(ProductList[1],Set([x,x^2,x^3,x^4]));
                    od;
                fi;
            od;

            ProductList[1] := DuplicateFreeList(ProductList[1]);
            
            dim := Size(ProductList[1]);

            for j in [t+1..dim] do
                ProductList[1][j] := ProductList[1][j][1];
            od;
            
            ProductList[2]:=StructuralCopy(T);
            ProductList[5]:=[1..t];

            for j in [t+1..dim] do
                
                x := ProductList[1][j];
            
                if Order(x) = 3 then 
            
                    Append(ProductList[5],[j,j]);
                    Append(ProductList[2],[x,x^2]);
                    
                elif Order(x) = 4 then 

                    Append(ProductList[5],[j,j]);
                    Append(ProductList[2],[x,x^3]);
                    
                elif Order(x) = 5 then 
                    Append(ProductList[5],[j,-j,-j,j]);
                    Append(ProductList[2],[x,x^2,x^3,x^4]); 
                fi;
            od;
            
            ProductList[2] := Flat(ProductList[2]);

            x:=Orbits(G,Cartesian(ProductList[1],ProductList[1]),OnPairs);
            
            ProductList[9] := [];
    
            for j in [1..Size(x)] do
                Add(ProductList[9], ShallowCopy(x[j]));
            od;

            # This is a bit of a patch, ask Markus tomorrow

            j:=1;

            while j < Size(ProductList[9]) + 1 do
                if Order(ProductList[9][j][1][1]) = 2 and Order(ProductList[9][j][1][2]) = 2 then
                    Remove(ProductList[9],j);
                else
                    j := j+1;
                fi;
            od;

            ProductList[9] := Concatenation(OrbitalsT,ProductList[9]);
            
            j := Size(OrbitalsT) + 1;
            
            while j < Size(ProductList[9]) + 1 do 
                if not [ProductList[9][j][1][2],ProductList[9][j][1][1]] in ProductList[9][j] then
                    k := j + 1;
                    
                    while k < Size(ProductList[9]) +1 do
                    
                        if  [ProductList[9][j][1][2],ProductList[9][j][1][1]]  in ProductList[9][k] then
                        
                            if Order(ProductList[9][j][1][1]) < Order(ProductList[9][j][1][2]) then 
                        
                                Append(ProductList[9][j],ProductList[9][k]);
                                Remove(ProductList[9],k);
                            
                            else 
                            
                                Append(ProductList[9][k],ProductList[9][j]);
                                Remove(ProductList[9],j);
                                
                                j := j - 1;
                                
                            fi;
                            
                            k := Size(ProductList[9]) + 1;
                            
                        else
                            
                            k := k + 1;
                        fi;
                    od;
                                        
                fi;
                
                j := j + 1;
                
            od;
            
            ProductList[14] := [];
            
            MAJORANA_MergeOrbitalsAxes(ProductList, OrbitalsT);

            SizeOrbitals:=Size(ProductList[9]);

            ProductList[7] := [];
            ProductList[4] := NullMat(dim,dim);
            ProductList[3] := NullMat(dim,dim);
            
            MAJORANA_PairRepresentatives(ProductList);
            MAJORANA_PairOrbits(ProductList);
            MAJORANA_PairConjElements(ProductList);
            
                                        ## STEP 3: PRODUCTS AND EVECS I ##


            # Set up algebra product and gram matrices
            
            OutputList := [0,0,0,0,0];

            AlgebraProducts := NullMat(1,SizeOrbitals)[1];
            GramMatrix := NullMat(1,SizeOrbitals)[1];

            for j in [1..SizeOrbitals] do
                AlgebraProducts[j]:=false;
                GramMatrix[j]:=false;
            od;

            # Set up eigenvector matrix

            EigenVectors:=NullMat(t,3);

            for j in [1..t] do
                if j in ProductList[10] then
                    for k in [1..3] do
                        EigenVectors[j][k] := [];
                    od;
                else
                    for k in [1..3] do
                        EigenVectors[j][k] := false;
                    od;
                fi;
            od;
            
            OutputList[1] := Shape;
            OutputList[2] := GramMatrix;
            OutputList[3] := AlgebraProducts;
            OutputList[4] := EigenVectors;
            OutputList[5] := ProductList;
            
            # Start filling in values and products!

            # (2,2) products and eigenvectors from IPSS10

            # Add eigenvectors from IPSS10

            for j in ProductList[10] do
                for k in [1..t] do

                    x := ProductList[3][j][k];

                    if Shape[x] = ['2','A'] then
                    
                        pos := [j, k, 0];
                        pos[3] := Position(T,T[j]*T[k]);
                        
                        vals := [-1/4, 1, 1];
                        
                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));
                        
                        vals := [0, 1, -1];

                        Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['2','B'] then
                    
                        pos := [k];
                        vals := [1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['3','A'] then
                    
                        pos := [j, k, 0, 0];

                        pos[3] := Position(T, T[j]*T[k]*T[j]);
                        pos[4] := ProductList[5][Position(ProductList[2],T[j]*T[k])];

                        vals := [-10/27, 32/27, 32/27, 1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));
                        
                        vals := [-8/45, -32/45, -32/45, 1];

                        Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

                        vals := [0, 1, -1, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['3','C'] then
                    
                        pos := [j, k, 0];

                        pos[3] := Position(T, T[j]*T[k]*T[j]);

                        vals := [-1/32, 1, 1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                        vals := [0, 1, -1];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['4','A'] then
                        
                        pos := [j, k, 0, 0, 0];
                        pos[3] := Position(T, T[j]*T[k]*T[j]);
                        pos[4] := Position(T, T[k]*T[j]*T[k]);
                        pos[5] := ProductList[5][Position(ProductList[2],T[j]*T[k])];

                        vals := [-1/2, 2, 2, 1, 1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                        vals := [-1/3, -2/3, -2/3, -1/3, 1]; 

                        Add(EigenVectors[j][2], MAJORANA_MakeVector(pos,vals,dim));

                        vals := [0, 1, -1, 0, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['4','B'] then
                        
                        pos := [j, k, 0, 0, 0];
                        pos[3] := Position(T, T[j]*T[k]*T[j]);
                        pos[4] := Position(T, T[k]*T[j]*T[k]);
                        pos[5] := Position(T, (T[j]*T[k])^2);
                        
                        vals := [-1/32, 1, 1, 1/8, -1/8];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos,vals,dim));

                        vals := [0, 1, -1, 0, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos,vals,dim));

                    elif Shape[x] = ['5','A'] then
                    
                        pos := [j, k, 0, 0, 0, 0];
                        pos[3] := Position(T, T[j]*T[k]*T[j]);
                        pos[4] := Position(T, T[k]*T[j]*T[k]);
                        pos[5] := Position(T, T[j]*T[k]*T[j]*T[k]*T[j]);
                        pos[6] := ProductList[5][Position(ProductList[2],T[j]*T[k])]; 

                        if pos[6] < 0 then
                            pos[6] := -pos[6];
                            sign := -1;
                        else
                            sign := 1;
                        fi;

                        vals := [3/512, -15/128, -15/128, -1/128, -1/128, sign*1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));

                        vals := [-3/512, 1/128, 1/128, 15/128, 15/128, sign*1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));
                        
                        vals := [0, 1/128, 1/128, -1/128, -1/128, sign*1];

                        Add(EigenVectors[j][2], MAJORANA_MakeVector(pos, vals, dim));

                        vals := [0, 1, -1, 0, 0, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

                        vals := [0, 0, 0, 1, -1, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

                    elif Shape[x] = ['6','A'] then

                        pos := [j, k, 0, 0, 0, 0, 0, 0];
                        pos[3] := Position(T, T[j]*T[k]*T[j]);
                        pos[4] := Position(T, T[k]*T[j]*T[k]);
                        pos[5] := Position(T, T[j]*T[k]*T[j]*T[k]*T[j]);
                        pos[6] := Position(T, T[k]*T[j]*T[k]*T[j]*T[k]);
                        pos[7] := Position(T, (T[j]*T[k])^3);
                        pos[8] := ProductList[5][Position(ProductList[2],(T[j]*T[k])^2)];

                        vals := [2/45, -256/45, -256/45, -32/45, -32/45, -32/45, 32/45, 1];

                        Add(EigenVectors[j][1], MAJORANA_MakeVector(pos, vals, dim));

                        vals := [-8/45, 0, 0, -32/45, -32/45, -32/45, 32/45, 1];

                        Add(EigenVectors[j][2], MAJORANA_MakeVector(pos, vals, dim));
                        
                        vals := [0, 1, -1, 0, 0, 0, 0, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));

                        vals := [0, 0, 0, 1, -1, 0, 0, 0];

                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));
                        
                        # put in products of 2A and 3A axes
                        
                        x := ProductList[3][pos[7]][pos[8]];
                        
                        AlgebraProducts[x] := [1..dim]*0;
                        
                        GramMatrix[x] := 0;
                    fi;
                od;
                
                # 1/32 eigenvectors from conjugation
                
                for k in [t+1..dim] do 
                
                    h := ProductList[1][k];
                    
                    if not h^T[j] in [h,h^2] then 
                    
                        pos := [k, 0];
                        pos[2] := ProductList[5][Position(ProductList[2],h^T[j])];
                        
                        if pos[2] < 0 then
                            pos[2] := -pos[2];
                            sign := -1;
                        else
                            sign := 1;
                        fi;
                        
                        vals := [1,-sign*1];
                        
                        Add(EigenVectors[j][3], MAJORANA_MakeVector(pos, vals, dim));
                    fi;
                od;
            od;
                
            # Products from IPSS10

            for j in [1..SizeOrbitals] do

                x := ProductList[7][j][1];
                y := ProductList[7][j][2];

                if Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 2 then

                    if Shape[j] = ['1','A'] then

                        AlgebraProducts[j] := NullMat(1,dim)[1];

                        AlgebraProducts[j][x] := 1;

                        GramMatrix[j] := 1;

                    elif Shape[j] = ['2','A'] then

                        pos := [x, y, 0];
                        pos[3] := Position(T,T[x]*T[y]);

                        vals := [1/8, 1/8, -1/8];
                        
                        AlgebraProducts[j] := MAJORANA_MakeVector( pos, vals, dim);

                        GramMatrix[j] := 1/8;

                    elif Shape[j] = ['2','B'] then

                        AlgebraProducts[j] := NullMat(1,dim)[1];

                        GramMatrix[j] := 0;

                    elif Shape[j] = ['3','A'] then
                    
                        pos := [x, y, 0, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        pos[4] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                        vals := [1/16, 1/16, 1/32, -135/2048];
                        
                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j] := 13/256;

                    elif Shape[j] = ['3','C'] then
                    
                        pos := [x, y, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        
                        vals := [1/64, 1/64, -1/64];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j]:=1/64;

                    elif Shape[j] = ['4','A'] then

                        pos := [x, y, 0, 0, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        pos[4] := Position(T,T[y]*T[x]*T[y]);
                        pos[5] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                        vals := [3/64, 3/64, 1/64, 1/64, -3/64]; 
    
                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j] := 1/32;

                    elif Shape[j] = ['4','B'] then
                    
                        pos := [x, y, 0, 0, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        pos[4] := Position(T,T[y]*T[x]*T[y]);
                        pos[5] := Position(T,(T[x]*T[y])^2);

                        vals := [1/64, 1/64, -1/64, -1/64, 1/64];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j]:=1/64;

                    elif  Shape[j] = ['5','A'] then
                    
                        pos := [x, y, 0, 0, 0, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        pos[4] := Position(T,T[y]*T[x]*T[y]);
                        pos[5] := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                        pos[6] := ProductList[5][Position(ProductList[2],T[x]*T[y])];

                        if pos[6] < 0 then
                            pos[6] := -pos[6];
                            sign := -1;
                        else
                            sign := 1;
                        fi;
                        
                        vals := [3/128, 3/128, -1/128, -1/128, -1/128, sign];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j]:=3/128;

                    elif Shape[j] = ['6','A'] then
                    
                        pos := [x, y, 0, 0, 0, 0, 0, 0];
                        pos[3] := Position(T,T[x]*T[y]*T[x]);
                        pos[4] := Position(T,T[y]*T[x]*T[y]);
                        pos[5] := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                        pos[6] := Position(T,T[y]*T[x]*T[y]*T[x]*T[y]);
                        pos[7] := Position(T,(T[x]*T[y])^3);
                        pos[8] := ProductList[5][Position(ProductList[2],(T[x]*T[y])^2)];
                        
                        vals := [1/64, 1/64, -1/64, -1/64, -1/64, -1/64, 1/64, 45/2048];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j]:=5/256;
                    fi;

                # 2,3 products

                elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 3 then
                    if ProductList[1][x]*ProductList[1][y] in T then 

                        s := ProductList[1][x]; h := ProductList[1][y];

                        # Inside a 3A algebra
                        
                        pos := [x, 0, 0, y];
                        pos[2] := Position(T,s*h);
                        pos[3] := Position(T,s*h*h);
                        
                        vals := [2/9, -1/9, -1/9, 5/32];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j]:=1/4;
                    fi;

                # 2,4 products

                elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 4 then

                    s := ProductList[1][x];
                    h := ProductList[1][y];

                    if s*h in T then

                        # Inside a 4A algebra
                        
                        pos := [x, 0, 0, 0, y];
                        pos[2] := Position(T,s*h);
                        pos[3] := Position(T,s*h*h*h);
                        pos[4] := Position(T,s*h*h);
                        
                        vals := [5/16, -1/8, -1/8, -1/16, 3/16];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j] := 3/8;
                    fi;

                # (2,5) values

                elif Order(ProductList[1][x]) = 2 and Order(ProductList[1][y]) = 5 then

                    s := ProductList[1][x];
                    h := ProductList[1][y];

                    if s*h in T then

                        # Inside a 5A algebra
                        
                        pos := [0, 0, 0, 0, y];
                        pos[1] := Position(T,s*h);
                        pos[2] := Position(T,s*h^4);
                        pos[3] := Position(T,s*h^2);
                        pos[4] := Position(T,s*h^3);
                        
                        vals := [7/4096, 7/4096, -7/4096, -7/4096, 7/32];

                        AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                        GramMatrix[j] := 0;
                    fi;
                    
                elif x = y then 
                
                    h := ProductList[1][x];
                    
                    if Order(h) = 3 then    # (3,3) values
                        
                        AlgebraProducts[j] := NullMat(1,dim)[1];
                        AlgebraProducts[j][x] := 1;

                        GramMatrix[j] := 8/5;
                        
                    elif Order(h) = 4 then  # (4,4) values
                    
                        AlgebraProducts[j] := NullMat(1,dim)[1];
                        AlgebraProducts[j][x] := 1;

                        GramMatrix[j] := 2;
                        
                    elif Order(h) = 5 then  # (5,5) values
                    
                        k := 1;

                        while k < t+1 do

                            if T[k]*h in T then

                                s:=T[k]; 
                                
                                pos := [k, 0, 0, 0, 0];
                                pos[2] := Position(T,s*h); 
                                pos[3] := Position(T,s*h^2); 
                                pos[4] := Position(T,s*h^3); 
                                pos[5] := Position(T,s*h^4);
                                
                                vals := [1,1,1,1,1]*(175/524288);

                                AlgebraProducts[j] := MAJORANA_MakeVector(pos, vals, dim);

                                k:=t+1;
                            else
                                k:=k+1;
                            fi;
                        od;

                        GramMatrix[j] := 875/2^(19);
                    fi;
                fi;
            od;

                                        ## STEP 4: MAIN LOOP ##
            
            maindimensions:=[];

            for j in ProductList[10] do
                for k in [1..3] do
                    if Size(EigenVectors[j][k]) > 0 then
                        EigenVectors[j][k]:=ShallowCopy(BaseMat(EigenVectors[j][k]));
                    fi;
                od;
                Add(maindimensions,Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3])+1);
            od;
            
            falsecount := [0,0];
            
            if false in GramMatrix then
                falsecount[1] := Size(Positions(GramMatrix,false));
            fi;
            
            if false in AlgebraProducts then
                falsecount[2] := Size(Positions(AlgebraProducts,false));
            fi;
            
            if ForAll(maindimensions, x -> x = dim) and falsecount = [0,0] then 
                switchmain := 1;
            else
                switchmain := 0;
            fi;
            
            count := 1;
            
            while switchmain = 0 do 
                
                count := count + 1;
                
                MAJORANA_MainSteps(i,GramMatrix,AlgebraProducts,EigenVectors,ProductList,Output,OutputList,maindimensions);
                
                if false in AlgebraProducts then 
                
                                                ## STEP 9: RESURRECTION PRINCIPLE ##
                
                    # Add nullspace to eigenvectors
                    
                    if ProductList[6] <> [] and ProductList[6] <> false then
                        for j in ProductList[10] do 
                            for k in [1..3] do
                                Append(EigenVectors[j][1],ProductList[6]);
                            od;
                        od;
                    fi;
                    
                    # put eigenvectors into reversed echelon form 
                    
                    for j in ProductList[10] do 
                        for k in [1..3] do 
                            if EigenVectors[j][k] <> [] then
                                MAJORANA_ReversedEchelonForm(EigenVectors[j][k]);
                            fi;
                        od;
                    od;
                    
                    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
                            
                    x := MAJORANA_FullResurrection(EigenVectors,unknowns,AlgebraProducts,ProductList,GramMatrix);
                    
                    mat := x[1]; vec := x[2]; record := x[3];
                    
                    if ProductList[6] <> [] and ProductList[6] <> false then 
                                
                        x := MAJORANA_NullSpaceAlgebraProducts(unknowns, AlgebraProducts, ProductList);
                        
                        MAJORANA_Append(x,mat,vec);
                    fi;
                    
                    if mat <> [] then 
                        x := MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
                            
                        if not x[1] and ProductList[6] <> false then 
                            Output[i] := MAJORANA_OutputError("Inconsistent system of unknown algebra products step 8"
                                            , [x[2],mat,vec,record]
                                            , OutputList);
                            break;
                        fi;
                    fi;
                
                fi;

                newdimensions := [];
                
                for j in ProductList[10] do 
                    Add(newdimensions, Size(EigenVectors[j][1]) 
                                        + Size(EigenVectors[j][2]) 
                                        + Size(EigenVectors[j][3]) + 1);
                od;
                
                newfalsecount := [0,0];
                
                if false in GramMatrix then
                    newfalsecount[1] := Size(Positions(GramMatrix,false));
                fi;
                
                if false in AlgebraProducts then
                    newfalsecount[2] := Size(Positions(AlgebraProducts,false));
                fi;
                
                if ForAll(newdimensions, x -> x = dim) and newfalsecount = [0,0] then
                    break;
                elif newdimensions = maindimensions and newfalsecount = falsecount then 

                    Output[i] := StructuralCopy(["Fail"
                                , "Missing values"
                                ,
                                , OutputList[1]
                                , OutputList[2]
                                , OutputList[3]
                                , OutputList[4]
                                , OutputList[5]]);
                    break;
                else
                    maindimensions := StructuralCopy(newdimensions);
                    falsecount := StructuralCopy(newfalsecount);
                fi;
                
                if count > 7 then 
                    switchmain := 1;
                fi;
                
            od;
            
            if Output[i] <> [] then 
                break;
            fi;
            
                                    ## STEP 10: INNER PRODUCTS FROM ORTHOGONALITY ##
                    
                                            
                # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products
                
                unknowns:=[];

                for j in [1..SizeOrbitals] do
                    if GramMatrix[j] = false then
                        Add(unknowns,j);
                    fi;
                od;
            
                if Size(unknowns) > 0 then

                    x := MAJORANA_FullOrthogonality(unknowns,EigenVectors,GramMatrix, ProductList);
                    
                    if not x[1] then 
                        Output[i] := MAJORANA_OutputError( x[2]
                                        , x[3]
                                        , OutputList);
                        break;
                    elif x[2] <> [] then 
                    
                        Error("Inner");
                    
                        y := MAJORANA_SolutionInnerProducts(x[2], x[3], unknowns, GramMatrix);
                        
                        if not y[1] then 
                            if Size(y[2]) <> 2 then 
                                Output[i] := MAJORANA_OutputError("Inconsistent system of unknown inner products"
                                            , [mat,vec]
                                            , OutputList);
                            fi;
                            break;
                        fi;
                    fi;
                fi;
                
                x := MAJORANA_CheckNullSpace(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
                
                if x = false then
                    Output[i] := MAJORANA_OutputError("The inner product is not positive definite"
                                        , []
                                        , OutputList);    
                    break;
                fi;


                                        ## STEP 11: CHECK ALGEBRA ##
                                    
            # Check bilinear form is positive definite
            
            GramMatrixFull := MAJORANA_FillGramMatrix(GramMatrix, ProductList);

            if MAJORANA_PositiveDefinite(GramMatrixFull) <0 then
                Output[i] := MAJORANA_OutputError("Gram Matrix is not positive definite"
                                , []
                                , OutputList);
                            
            fi;

            # Check that all triples obey axiom M1

            error:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,ProductList);

            if Size(error)>0 then
                Output[i] := MAJORANA_OutputError("Algebra does not obey axiom M1"
                                    , error
                                    , OutputList);
            fi;

            # Check that eigenvectors obey the fusion rules

            error := MAJORANA_TestFusion(GramMatrix,AlgebraProducts,EigenVectors,ProductList);

            if ForAny(error,x->Size(x)>0) then
                Output[i] := MAJORANA_OutputError("Algebra does not obey fusion rules"
                                    , error
                                    , OutputList);
                break;
            fi;

            # Check that the eigenspaces are orthogonal

            error := MAJORANA_TestOrthogonality(GramMatrix,AlgebraProducts,EigenVectors,ProductList);

            if Size(error) > 0 then
                Output[i] := MAJORANA_OutputError("Eigenspaces are not orthogonal with respect to the inner product"
                            , error
                            , OutputList);
                break;
            fi;

            # Check M2

            # error:=MAJORANA_AxiomM2(GramMatrix,AlgebraProducts,ProductList);

            # if error = -1 then
            #    Output[i] := MAJORANA_OutputError("Algebra does not obey axiom M2"
            #                        , error
            #                        , OutputList);
            #   break;
            #fi;

            Output[i] := StructuralCopy(["Success"
                        ,
                        ,
                        , OutputList[1]
                        , OutputList[2]
                        , OutputList[3]
                        , OutputList[4]
                        , OutputList[5]]);
            master:=0;
        od;
        
    od;    

    return Output;

    end );

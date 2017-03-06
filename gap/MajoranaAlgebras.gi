#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

InstallGlobalFunction(MAJORANA_LDLTDecomposition,

function(A) # Takes as input a matrix A. If A is positive semidefinite then will return [L,D] such that A= LDL^T. Else returns 0. Note: does not test if matrix is square or symmetric.

        local B, n, L, D, i, j, k, temp;

        B:=ShallowCopy(A); n:=Size(B); L:=NullMat(n,n); D:=NullMat(n,n);

        for i in [1..n] do
            temp:=[];
            for j in [1..i-1] do
                Append(temp,[L[i][j]*L[i][j]*D[j][j]]);
            od;

            D[i][i] := B[i][i] - Sum(temp);

            if D[i][i] =0 then
                    for j in [i+1..n] do
                        temp:=[];
                        for k in [1..i-1] do
                            Append(temp,[L[j][k]*L[i][k]*D[k][k]]);
                        od;
                        if B[j][i] - Sum(temp)= 0 then
                            L[j][i]:=0;
                        else
                            return D;
                        fi;
                    od;
                    L[i][i]:=1;
            else
                for j in [i+1..n] do
                    temp:=[];
                    for k in [1..i-1] do
                        Append(temp,[L[j][k]*L[i][k]*D[k][k]]);
                    od;
                    L[j][i]:=(B[j][i] - Sum(temp))/D[i][i];
                od;
                L[i][i]:=1;
            fi;
        od;

        return Concatenation([L],[D]);
        end
    );

InstallGlobalFunction(MAJORANA_SolutionMatVecs,

function(mat,vecs) # Takes as input two matrices, the second being interpreted as a vector of vectors. Returns a list of size four if system is inconsistent, otherwise returns a list of size 4

        local A, C, n, m, d, absd, B, i, j, k, x, imax, temp, tempv, tempi, sol, list, newmat, newvec, pos, p, unsolved, zeros, error;

        A:=StructuralCopy(mat);
        B:=StructuralCopy(vecs);

        n:=Size(A);
        m:=Size(A[1]);

        if n<m then
            A{[n+1..m]} := NullMat(m-n,m);
        elif m<n then
            for i in [1..n] do
                A[i]{[m+1..n]} := NullMat(1,n-m)[1];
            od;
        fi;

        p:=Maximum(n,m);

        d:=NullMat(1,p)[1];

        C:=IdentityMat(p);

        # Put matrix in row echelon form

        i:=1;

        while i <= p do

            for j in [i..p] do
                d[j]:=A[j][i];
            od;

            absd:=List(d,x->AbsoluteValue(x));

            imax:=Position(absd,Maximum(absd));

            if d[imax] = 0 then

                k:=i+1;

                while k <= p do
                    if A[i][k] <> 0 then
                        C[i] := C[i]/A[i][k];
                        B[i] := B[i]/A[i][k];
                        A[i] := A[i]/A[i][k];


                        k:=p+1;
                    else
                        k:=k+1;
                    fi;
                od;

                i:=i+1;

            else

                # Swap rows i and imax

                temp:=ShallowCopy(A[imax]); tempv:=ShallowCopy(B[imax]); tempi:=ShallowCopy(C[imax]);
                A[imax]:=ShallowCopy(A[i]); A[i]:=ShallowCopy(temp);
                B[imax]:=ShallowCopy(B[i]); B[i]:=ShallowCopy(tempv);
                C[imax]:=ShallowCopy(C[i]); C[i]:=ShallowCopy(tempi);

                for k in [i+1..p] do

                    x:=A[k][i]/A[i][i];

                    B[k]:=B[k] - x*B[i];
                    C[k]:=C[k] - x*C[i];
                    A[k]:=A[k] - x*A[i];

                od;

                C[i]:=C[i]/A[i][i];
                B[i]:=B[i]/A[i][i];
                A[i]:=A[i]/A[i][i];


                d[i]:=0;

                i:=i+1;

            fi;

        od;

        # Check if we can solve the system of equations

        newmat:=NullMat(p,p);
        newvec:=NullMat((p),Size(vecs[1]));
        error:=[];

        for i in [1..p] do
            if ForAll(A[i],x->x=0) then
                if ForAny(B[i],x->x<> 0 ) then

                    Append(error,[i]);

                fi;
            else
                pos:=Position(A[i],1);
                newmat[pos] := StructuralCopy(A[i]);
                newvec[pos] := StructuralCopy(B[i]);
            fi;
        od;

        if Size(error) >0 then
            # no solutions
            return [error,C,A,B];
        fi;

        zeros:=NullMat(1,p)[1];
        sol:=NullMat(1,m)[1];
        unsolved:=[];

        if newmat[m] = zeros then
            # sol[m] is unknown

            sol[m]:=[];
            Append(unsolved,[m]);

        else
            sol[m]:=newvec[m];
        fi;

        for i in [1..m-1] do
            if newmat[m-i] = zeros then
                sol[m-i] := [];
                Append(unsolved,[m-i]);
            else
                list:=[];
                j:=m-i+1;
                while j<=m do
                    if newmat[m-i][j] <> 0 then
                        if not j in unsolved then
                            Append(list,[newmat[m-i][j]*sol[j]]);
                            j:=j+1;
                        else
                            sol[m-i] := [];
                            Append(unsolved,[m-i]);
                            j:=m+1;
                        fi;
                    else
                        j:=j+1;
                    fi;
                od;
                if not m-i in unsolved then
                    sol[m-i]:=(newvec[m-i] - Sum(list))/newmat[m-i][m-i];
                fi;
            fi;
        od;



        return [sol,unsolved];


        end

        );

InstallGlobalFunction(MAJORANA_NullSpace,

        function(mat) # Takes as input matrix, returns a matrix whose rows form a basis of the nullspace of mat

        local A, C, n, m, d, absd, i, j, k, x, imax, temp, tempi, basis, basic, free, vec;

        A:=ShallowCopy(mat);

        n:=Size(A);
        m:=Size(A[1]);

        d:=NullMat(1,n)[1];

        C:=IdentityMat(n);

        # Put matrix in row echelon form

        i:=1;

        while i <= n do

            for j in [i..n] do
                d[j]:=A[j][i];
            od;

            absd:=List(d,x->AbsoluteValue(x));

            imax:=Position(absd,Maximum(absd));

            if d[imax] = 0 then

                k:=i+1;

                while k <= m do
                    if A[i][k] <> 0 then

                        # Turn leading coefficient of row i (at pos k) into 1

                        C[i] := C[i]/A[i][k];
                        A[i] := A[i]/A[i][k];

                        k:=m+1;
                    else
                        k:=k+1;
                    fi;
                od;

                i:=i+1;

            else

                # Swap rows i and imax

                temp:=ShallowCopy(A[imax]); tempi:=ShallowCopy(C[imax]);
                A[imax]:=ShallowCopy(A[i]); A[i]:=ShallowCopy(temp);
                C[imax]:=ShallowCopy(C[i]); C[i]:=ShallowCopy(tempi);

                for k in [i+1..n] do

                    x:=A[k][i]/A[i][i];

                    C[k]:=C[k] - x*C[i];
                    A[k]:=A[k] - x*A[i];

                od;

                C[i]:=C[i]/A[i][i];
                A[i]:=A[i]/A[i][i];


                d[i]:=0;

                i:=i+1;

            fi;

        od;

        # Compute null space

        basis:=[];

        basic:=[];
        free:=[];


        for i in [1..n] do
            Append(basic,[Position(A[i],1)]);
        od;

        for i in [0..m-1] do
            if not m-i in basic then
                Append(free,[m-i]);
            fi;
        od;

        for i in free do

            vec:=NullMat(1,m)[1];

            for j in [0..m-1] do
                if i = m-j then
                    vec[m-j] :=1;
                elif m-j in basic then
                    for k in [m-j+1..m] do
                        vec[m-j] := vec[m-j] - A[m-j][k]*vec[k];
                    od;
                fi;
            od;

            Append(basis,[vec]);
        od;

        for j in [1..Size(basis)] do
            basis[j]:=basis[j]/basis[j][m-j+1];

            for k in [1..j-1] do
                basis[j]:=basis[j] - basis[j][m-k+1]*basis[k];
            od;
        od;

        for j in [1..Size(basis)] do
            for k in [1..(Size(basis)-j)] do
                basis[j]:=basis[j] - basis[j][m-k]*basis[m-k];
            od;
        od;

        return basis;

        end

        );

InstallGlobalFunction(  MAJORANA_AlgebraProduct,

        function(u,v,basis,done) # If all the relevant products are known, returns the algebra product of u and v. If not, returns 0

        local i, j, sum;

        sum:=[];

        if ForAll(u,x-> x= 0 ) or ForAll(v,x->x=0) then
            return u*0;
        fi;

        for i in [1..Size(u)] do
            if u[i] <> 0 then
                for j in [1..Size(v)] do
                    if v[j] <> 0 then
                        if [i,j] in done then
                            Append(sum,[u[i]*v[j]*basis[i][j]]);
                        else

                            # cannot calculate product

                            return 0;
                        fi;
                    fi;
                od;
            fi;
        od;
        return Sum(sum);
        end

        );

InstallGlobalFunction(  MAJORANA_InnerProduct,

    function(u,v,gram,done) # If all the relevant products are known, returns the algebra product of u and v. If not, returns [0]

        local i, j, sum;

        sum:=[];

        for i in [1..Size(u)] do
            if u[i] <> 0 then
                for j in [1..Size(v)] do
                    if v[j] <> 0 then
                        if [i,j] in done then
                            Append(sum,[u[i]*v[j]*gram[i][j]]);
                        else

                            # cannot calculate product

                            return [0];
                        fi;
                    fi;
                od;
            fi;
        od;
        return Sum(sum);
        end

    );

InstallGlobalFunction(MAJORANA_PositiveDefinite,

        function(GramMatrix) # Check returns 1, 0, -1 if Gram matrix is positive definite, positive semidefinite or neither respectively

        local L, Diagonals, i;

        L:=MAJORANA_LDLTDecomposition(GramMatrix);

        Diagonals:=[];

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

InstallGlobalFunction(MAJORANA_AxiomM1,

        function(GramMatrix,AlgebraProducts) # Checks if bilinear and algebra products obey axiom M1, outputs a list which is empty if they do obey the axiom

        local sum1, sum2, ErrorM1, j, k, l, m, dim;

        dim:=Size(AlgebraProducts);

        ErrorM1:=[];

        for j in [1..dim] do
            for k in [1..dim] do
                for l in [1..dim] do

                    sum1:=[];
                    sum2:=[];

                    for m in [1..dim] do
                        if AlgebraProducts[k][l][m] <> 0 then
                            Append(sum1,[AlgebraProducts[k][l][m]*GramMatrix[j][m]]);
                        fi;
                        if AlgebraProducts[j][k][m] <> 0 then
                            Append(sum2,[AlgebraProducts[j][k][m]*GramMatrix[m][l]]);
                        fi;
                    od;

                    if Sum(sum1) <> Sum(sum2) then
                        Add(ErrorM1,[j,k,l]);
                    fi;

                od;
            od;
        od;

        return ErrorM1;

        end

        );

InstallGlobalFunction(MAJORANA_Fusion,

        function(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors) # Checks if algebra obeys the fusion rules, outputs list of six lists which are empty if it does obey fusion rules

        local errorfusion00, errorfusion02, errorfusion04, errorfusion22, errorfusion24, errorfusion44, a, t, j, k, l, x, y, z, x0;

        errorfusion00:=[];
        errorfusion02:=[];
        errorfusion04:=[];
        errorfusion22:=[];
        errorfusion24:=[];
        errorfusion44:=[];

        t:=Size(T);

        for j in [1..t] do

            a:=NullMat(1,Size(AlgebraProducts))[1]; a[j]:=1;

            # 0,0 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][1])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][1][l],AlgebraProducts,KnownAlgebraProducts);
                    if x<>0 then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts);
                        if y <> 0 then
                            if ForAny(y, z-> z <> 0) then
                                Add(errorfusion00,[j,k,l]);
                            fi;
                        fi;
                    fi;
                od;
            od;

            # 0,1/4 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][2][l],AlgebraProducts,KnownAlgebraProducts);
                    if x<>0 then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts);
                        if y <> 0 and y <> x/4 then
                            Add(errorfusion02,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 0,1/32 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][3][l],AlgebraProducts,KnownAlgebraProducts);
                    if x <> 0 then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts);
                        if y <> 0 and y <> x/32 then
                            Add(errorfusion04,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 1/4,1/4 fusion

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][2][l],AlgebraProducts,KnownAlgebraProducts);
                    if x<>0 then
                        if MAJORANA_InnerProduct(a,x,GramMatrix,KnownInnerProducts) <> [0] then

                            y:= x - MAJORANA_InnerProduct(a,x,GramMatrix,KnownInnerProducts)*a;
                            z:=MAJORANA_AlgebraProduct(a,y,AlgebraProducts,KnownAlgebraProducts);

                            if z <> 0 and ForAny(z, x -> x <> 0)  then
                                Add(errorfusion22,[j,k,l]);
                            fi;
                        fi;
                    fi;
                od;
            od;

            # 1/4,1/32 fusion

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][3][l],AlgebraProducts,KnownAlgebraProducts);
                    if x <> 0 then
                        if MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> x/32 then
                            Add(errorfusion24,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 1/32,1/32 fusion

            for k in [1..Size(EigenVectors[j][3])] do
                for l in [1..Size(EigenVectors[j][3])] do

                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][3][k],EigenVectors[j][3][l],AlgebraProducts,KnownAlgebraProducts);

                    if x <> 0 then
                        y:= x - MAJORANA_InnerProduct(a,x,GramMatrix,KnownInnerProducts)*a;
                        z:= y - 4*MAJORANA_AlgebraProduct(a,y,AlgebraProducts,KnownAlgebraProducts);
                        x0:=MAJORANA_AlgebraProduct(a,z,AlgebraProducts,KnownAlgebraProducts);
                        
                        if x0 <> 0 and ForAny(x0, x -> x <> 0)  then
                            Add(errorfusion44,[j,k,l]);
                        fi;
                    fi;
                od;
            od;
        od;

        return [errorfusion00,errorfusion02,errorfusion04,errorfusion22,errorfusion24,errorfusion44];

        end

        );

InstallGlobalFunction(MAJORANA_Orthogonality,

        function(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors) # Tests that eigenspaces are orthogonal with respect to the inner product

        local a, b, t, errorortho01, errorortho02, errorortho04, errorortho12, errorortho14, errorortho24, j, k, l;

        t:=Size(T);

        for j in [1..t] do

            a:=NullMat(1,Size(AlgebraProducts))[1];
            a[j]:=1;

            errorortho01:=[];
            errorortho02:=[];
            errorortho04:=[];
            errorortho12:=[];
            errorortho14:=[];
            errorortho24:=[];

            for k in [1..Size(EigenVectors[j][1])] do
                if MAJORANA_InnerProduct(EigenVectors[j][1][k],a,GramMatrix,KnownInnerProducts) <> 0 then
                    Append(errorortho01,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    if MAJORANA_InnerProduct(EigenVectors[j][1][k],EigenVectors[j][2][l],GramMatrix,KnownInnerProducts) <> 0 then
                        Append(errorortho02,[[j,k,l]]);
                    fi;
                od;
            od;

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    if MAJORANA_InnerProduct(EigenVectors[j][1][k],EigenVectors[j][3][l],GramMatrix,KnownInnerProducts) <> 0 then
                        Append(errorortho02,[[j,k,l]]);
                    fi;
                od;
            od;

            for k in [1..Size(EigenVectors[j][2])] do
                if MAJORANA_InnerProduct(EigenVectors[j][2][k],a,GramMatrix,KnownInnerProducts) <> 0 then
                    Append(errorortho12,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][3])] do
                if MAJORANA_InnerProduct(EigenVectors[j][3][k],a,GramMatrix,KnownInnerProducts) <> 0 then
                    Append(errorortho14,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    if MAJORANA_InnerProduct(EigenVectors[j][2][k],EigenVectors[j][3][l],GramMatrix,KnownInnerProducts) <> 0 then
                        Append(errorortho24,[[j,k,l]]);
                    fi;
                od;
            od;
        od;

        return [errorortho01,errorortho02,errorortho04,errorortho12,errorortho14,errorortho24];

        end

        );

InstallGlobalFunction(MAJORANA_AxiomM2,

        function(KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts) # Tests that the algebra obeys axiom M2

        local B, dim, L, i, j , k , l, m, Diagonals;

        dim:=Size(AlgebraProducts);

        B:=NullMat(dim^2,dim^2);

        for j in [1..dim] do
            for k in [1..dim] do
                for l in [1..dim] do
                    for m in [1..dim] do
                        B[dim*(j-1) + k][dim*(l-1) +m]:= MAJORANA_InnerProduct(AlgebraProducts[j][l],AlgebraProducts[k][m],GramMatrix,KnownInnerProducts) - MAJORANA_InnerProduct(AlgebraProducts[k][l],AlgebraProducts[m][j],GramMatrix,KnownInnerProducts);
                    od;
                od;
            od;
        od;

        L:=MAJORANA_LDLTDecomposition(B);

        Diagonals:=[];

        for i in [1..Size(B)] do
            Append(Diagonals,[L[2][i][i]]);
        od;

        if ForAny(Diagonals, x->x<0) then
            return -1;
        else
            return 1;
        fi;

        end

        );

InstallGlobalFunction(MAJORANA_Form,

        function(str) # Outputs value of bilinear product depending on shape of orbital

        if str = ['1','A'] then
            return 1;
        elif str = ['2','A'] then
            return 1/8;
        elif str = ['2','B'] then
            return 0;
        elif str = ['3','A'] then
            return 13/256;
        elif str = ['3','C'] then
            return 1/64;
        elif str = ['4','A'] then
            return 1/32;
        elif str = ['4','B'] then
            return 1/64;
        elif str = ['5','A'] then
            return 3/128;
        elif str = ['6','A'] then
            return 5/256;
        fi;
        end

        );

InstallGlobalFunction(MajoranaRepresentation,

function(G,T)

    local   # error checking
            ErrorFusion, ErrorM1, ErrorM2, ErrorOrthogonality,

            # indexing and temporary variables
            i, j, k, l, m, n, x, y, z,

            # Step 0 - Set Up
            Output, t, Orbitals, SizeOrbitals,

            # Step 1 - Shape
            Shape, RepsSquares6A, Unknowns3X,

            # Step 2 - Possible shapes
            Binaries, master, 3Aaxes, 4Aaxes, 5Aaxes, 5AaxesFixed, u, v, w,

            # Step 3 - Products and evecs I
            GramMatrix, GramMatrixT, LI, NullSpT, AlgebraProducts, EigenVectors, KnownInnerProducts, KnownAlgebraProducts, EigenVector, sign, x0, x1, xm1, x2, xm2, x3, x4, x5, x6, x7, x8, x2A, x3A,

            # Step 4 - More products and evecs
            h, s, xj, xk, xl, xik, xil, xjk, xjl, xkl, xx, L, Diagonals, NullSp, dim, a,

            # Step 5 - More evecs
            switch, Dimensions, NewDimensions, NewEigenVectors,

            # Step 6 - More inner products
            UnknownInnerProducts, mat, vec, sum, row, Solution,

            # Step 7 - More algebra products
            Alpha, Alpha2, Beta, Beta2, walpha, wbeta, c, Form, str1, str2, str3, str4, str5, zeros,  UnknownAlgebraProducts, record;


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

    # Construct orbitals of G on T

    Orbitals:=OrbitsDomain(G,Cartesian(T,T),OnPairs);

    SizeOrbitals:=Size(Orbitals);

                                        ## STEP 1: SHAPE ##

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    Shape:=NullMat(1,SizeOrbitals)[1];

    RepsSquares6A:=[];

    for i in [1..SizeOrbitals] do
        x:=Representative(Orbitals[i]);
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

    for i in [1..SizeOrbitals] do
        if Shape[i][1] = '3' then
            j:=0;
            while j<Size(Orbitals[i]) do
                j:=j+1;
                if Orbitals[i][j][1]*Orbitals[i][j][2] in RepsSquares6A then
                    Shape[i]:="3A";;
                    j:=Size(Orbitals[i])+1;;
                fi;
            od;
        fi;
    od;

    Unknowns3X:=[];

    for i in [1..SizeOrbitals] do
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

            3Aaxes:=[];
            4Aaxes:=[];
            5Aaxes:=[];

            for j in [1..SizeOrbitals] do
                if Shape[j]=['3','A'] then
                    for k in [1..Size(Orbitals[j])] do
                        Add(3Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));
                    od;
                fi;
                if Shape[j]=['4','A'] then
                    for k in [1..Size(Orbitals[j])] do
                        Add(4Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));
                    od;
                fi;
                if Shape[j]=['5','A'] then
                    for k in [1..Size(Orbitals[j])] do
                        Add(5Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));
                    od;
                fi;
            od;

            3Aaxes:=DuplicateFreeList(3Aaxes); u:=Size(3Aaxes);
            4Aaxes:=DuplicateFreeList(4Aaxes); v:=Size(4Aaxes);
            5Aaxes:=DuplicateFreeList(5Aaxes); w:=Size(5Aaxes);

            5AaxesFixed:=NullMat(w,0);

            for j in [1..w] do
                x:=5Aaxes[j].1;
                5AaxesFixed[j]:=[x,x^4];
            od;


                                        ## STEP 3: PRODUCTS AND EVECS I ##

            # Set up Gram matrix

            dim:= t+u+v+w;

            GramMatrix:=NullMat(dim,dim);

            KnownInnerProducts:=[];

            # Set up algebra product matrix

            AlgebraProducts:=NullMat(dim,dim);

            for j in [1..dim] do
                for k in [1..dim] do
                    AlgebraProducts[j][k]:=NullMat(1,dim)[1];
                od;
            od;

            KnownAlgebraProducts:=[];

            # Set up eigenvector matrix

            EigenVectors:=NullMat(t,3);

            for j in [1..t] do
                for k in [1..3] do
                    EigenVectors[j][k]:=[];
                od;
            od;

            # Start filling in values and products!

            l:=1;

            for j in [1..t+u+v+w] do
                for k in [1..t+u+v+w] do
                    Add(KnownInnerProducts,[j,k]);
                od;
            od;

            for j in [1..t] do
                for k in [1..t] do
                    Add(KnownAlgebraProducts,[j,k]);
                od;
            od;

            # (2,2) products and eigenvectors from IPSS10

            for j in [1..(t+u+v)] do
                GramMatrix[j][j]:=1;
                AlgebraProducts[j][j][j]:=1;
            od;

            for j in [(t+1)..dim] do
                Add(KnownAlgebraProducts,[j,j]);
            od;

            # Add products of 5A axes with themselves

            for j in [1..w] do
                GramMatrix[t+u+v+j][t+u+v+j]:=875/524288;

                l:=1;

                while l < t+1 do
                    h:=5AaxesFixed[j][1];
                    if T[l]*h in T then
                        x:=T[l]; x1:=Position(T,x*h); x2:=Position(T,x*h*h); x3:=Position(T,x*h*h*h); x4:=Position(T,x*h*h*h*h);

                        AlgebraProducts[t+u+v+j][t+u+v+j][l]:=175/524288;
                        AlgebraProducts[t+u+v+j][t+u+v+j][x1]:=175/524288;
                        AlgebraProducts[t+u+v+j][t+u+v+j][x2]:=175/524288;
                        AlgebraProducts[t+u+v+j][t+u+v+j][x3]:=175/524288;
                        AlgebraProducts[t+u+v+j][t+u+v+j][x4]:=175/524288;

                        l:=t+1;
                    else
                        l:=l+1;
                    fi;
                od;

            od;



            # Remaining products from IPSS10

            for j in [1..SizeOrbitals] do
                if Shape[j] = ['2','A'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); x2:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]);

                        GramMatrix[x0][x1]:=1/8;

                        AlgebraProducts[x0][x1][x0]:= 1/8;
                        AlgebraProducts[x0][x1][x1]:= 1/8;
                        AlgebraProducts[x0][x1][x2]:=-1/8;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[x2]:=1;
                        EigenVector[x0]:=-1/4;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[x2]:=-1;

                        Add(EigenVectors[x0][2],StructuralCopy(EigenVector));

                    od;
                fi;
                if Shape[j] = ['2','B'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]);

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                    od;
                fi;
                if Shape[j] = ['3','A'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x3:=Position(3Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));

                        GramMatrix[x0][x1]:=13/256;

                        AlgebraProducts[x0][x1][x0]:=1/16;
                        AlgebraProducts[x0][x1][x1]:=1/16;
                        AlgebraProducts[x0][x1][xm1]:=1/32;
                        AlgebraProducts[x0][x1][t+x3]:= - 135/2048;

                        AlgebraProducts[x0][t+x3][x0]:=2/9;
                        AlgebraProducts[x0][t+x3][x1]:=-1/9;
                        AlgebraProducts[x0][t+x3][xm1]:=-1/9;
                        AlgebraProducts[x0][t+x3][t+x3]:=5/32;

                        AlgebraProducts[t+x3][x0][x0]:=2/9;
                        AlgebraProducts[t+x3][x0][x1]:=-1/9;
                        AlgebraProducts[t+x3][x0][xm1]:=-1/9;
                        AlgebraProducts[t+x3][x0][t+x3]:=5/32;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-10/27;
                        EigenVector[x1]:=32/27;
                        EigenVector[xm1]:=32/27;
                        EigenVector[t+x3]:=1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-8/45;
                        EigenVector[x1]:=-32/45;
                        EigenVector[xm1]:=-32/45;
                        EigenVector[t+x3]:=1;

                        Add(EigenVectors[x0][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));

                    od;
                fi;
                if Shape[j] = ['3','C'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]);

                        GramMatrix[x0][x1]:=1/64;

                        AlgebraProducts[x0][x1][x0]:=1/64;
                        AlgebraProducts[x0][x1][x1]:=1/64;
                        AlgebraProducts[x0][x1][xm1]:=-1/64;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/32;
                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));
                    od;
                fi;
                if Shape[j] = ['4','A'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x2:=Position(T,Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);
                        x4:=Position(4Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));

                        GramMatrix[x0][x1]:=1/32;

                        AlgebraProducts[x0][x1][x0]:=3/64;
                        AlgebraProducts[x0][x1][x1]:=3/64;
                        AlgebraProducts[x0][x1][xm1]:=1/64;
                        AlgebraProducts[x0][x1][x2]:=1/64;
                        AlgebraProducts[x0][x1][t+u+x4]:= -3/64;

                        AlgebraProducts[x0][t+u+x4][x0]:=5/16;
                        AlgebraProducts[x0][t+u+x4][x1]:=-1/8;
                        AlgebraProducts[x0][t+u+x4][xm1]:=-1/8;
                        AlgebraProducts[x0][t+u+x4][x2]:=-1/16;
                        AlgebraProducts[x0][t+u+x4][t+u+x4]:= 3/16;

                        AlgebraProducts[t+u+x4][x0][x0]:=5/16;
                        AlgebraProducts[t+u+x4][x0][x1]:=-1/8;
                        AlgebraProducts[t+u+x4][x0][xm1]:=-1/8;
                        AlgebraProducts[t+u+x4][x0][x2]:=-1/16;
                        AlgebraProducts[t+u+x4][x0][t+u+x4]:= 3/16;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/2;
                        EigenVector[x1]:=2;
                        EigenVector[xm1]:=2;
                        EigenVector[x2]:=1;
                        EigenVector[t+u+x4]:=1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/3;
                        EigenVector[x1]:=-2/3;
                        EigenVector[xm1]:=-2/3;
                        EigenVector[x2]:=-1/3;
                        EigenVector[t+u+x4]:=1;

                        Add(EigenVectors[x0][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));

                    od;
                fi;
                if Shape[j] = ['4','B'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x2:=Position(T,Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);
                        x4:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);

                        GramMatrix[x0][x1]:=1/64;

                        AlgebraProducts[x0][x1][x0]:=1/64;
                        AlgebraProducts[x0][x1][x1]:=1/64;
                        AlgebraProducts[x0][x1][xm1]:=-1/64;
                        AlgebraProducts[x0][x1][x2]:=-1/64;
                        AlgebraProducts[x0][x1][x4]:= 1/64;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/32;
                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=1;
                        EigenVector[x2]:=1/8;
                        EigenVector[x4]:=-1/8;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));

                    od;
                fi;
                if Shape[j] = ['5','A'] then
                    for k in [1..Size(Orbitals[j])] do
                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x2:=Position(T,Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);
                        xm2:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]);
                        x5:=Position(5Aaxes,Group(Orbitals[j][k][1]*Orbitals[j][k][2]));

                        if Orbitals[j][k][1]*Orbitals[j][k][2] in 5AaxesFixed[x5] then
                            sign:=1;
                        else
                            sign:=-1;
                        fi;

                        GramMatrix[x0][x1]:=3/128;

                        AlgebraProducts[x0][x1][x0]:=3/128;
                        AlgebraProducts[x0][x1][x1]:=3/128;
                        AlgebraProducts[x0][x1][xm1]:=-1/128;
                        AlgebraProducts[x0][x1][x2]:=-1/128;
                        AlgebraProducts[x0][x1][xm2]:=-1/128;
                        AlgebraProducts[x0][x1][t+u+v+x5] := sign*1;

                        AlgebraProducts[x0][t+u+v+x5][x1]:=sign*7/4096;
                        AlgebraProducts[x0][t+u+v+x5][xm1]:=sign*7/4096;
                        AlgebraProducts[x0][t+u+v+x5][x2]:=-sign*7/4096;
                        AlgebraProducts[x0][t+u+v+x5][xm2]:=-sign*7/4096;
                        AlgebraProducts[x0][t+u+v+x5][t+u+v+x5] := 7/32;

                        AlgebraProducts[t+u+v+x5][x0][x1]:=sign*7/4096;
                        AlgebraProducts[t+u+v+x5][x0][xm1]:=sign*7/4096;
                        AlgebraProducts[t+u+v+x5][x0][x2]:=-sign*7/4096;
                        AlgebraProducts[t+u+v+x5][x0][xm2]:=-sign*7/4096;
                        AlgebraProducts[t+u+v+x5][x0][t+u+v+x5] := 7/32;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=3/512;
                        EigenVector[x1]:=-15/128;
                        EigenVector[xm1]:=-15/128;
                        EigenVector[x2]:=-1/128;
                        EigenVector[xm2]:=-1/128;
                        EigenVector[t+u+v+x5]:=sign*1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-3/512;
                        EigenVector[x1]:=1/128;
                        EigenVector[xm1]:=1/128;
                        EigenVector[x2]:=15/128;
                        EigenVector[xm2]:=15/128;
                        EigenVector[t+u+v+x5]:=sign*1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1/128;
                        EigenVector[xm1]:=1/128;
                        EigenVector[x2]:=-1/128;
                        EigenVector[xm2]:=-1/128;
                        EigenVector[t+u+v+x5]:=sign*1;

                        Add(EigenVectors[x0][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x2]:=1;
                        EigenVector[xm2]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));
                    od;
                fi;
                if Shape[j] = ['6','A'] then
                    for k in [1..Size(Orbitals[j])] do

                        x0:=Position(T,Orbitals[j][k][1]); x1:=Position(T,Orbitals[j][k][2]); xm1:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x2:=Position(T,Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);
                        xm2:=Position(T,Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]); x3:=Position(T,Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]*Orbitals[j][k][1]*Orbitals[j][k][2]);
                        x2A:=Position(T,(Orbitals[j][k][1]*Orbitals[j][k][2])^3); x3A:=Position(3Aaxes,Group((Orbitals[j][k][1]*Orbitals[j][k][2])^2));

                        GramMatrix[x0][x1]:=5/256;

                        AlgebraProducts[x0][x1][x0]:=1/64;
                        AlgebraProducts[x0][x1][x1]:=1/64;
                        AlgebraProducts[x0][x1][xm1]:=-1/64;
                        AlgebraProducts[x0][x1][x2]:=-1/64;
                        AlgebraProducts[x0][x1][xm2]:=-1/64;
                        AlgebraProducts[x0][x1][x3] := -1/64;
                        AlgebraProducts[x0][x1][x2A] := 1/64;
                        AlgebraProducts[x0][x1][t+x3A] := 45/2048;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=2/45;
                        EigenVector[x1]:=-256/45;
                        EigenVector[xm1]:=-256/45;
                        EigenVector[x2]:=-32/45;
                        EigenVector[xm2]:=-32/45;
                        EigenVector[x3]:=-32/45;
                        EigenVector[x2A]:=32/45;
                        EigenVector[t+x3A]:=1;

                        Add(EigenVectors[x0][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-8/45;
                        EigenVector[x2]:=-32/45;
                        EigenVector[xm2]:=-32/45;
                        EigenVector[x3]:=-32/45;
                        EigenVector[x2A]:=32/45;
                        EigenVector[t+x3A]:=1;

                        Add(EigenVectors[x0][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x2]:=1;
                        EigenVector[xm2]:=-1;

                        Add(EigenVectors[x0][3],StructuralCopy(EigenVector));
                    od;
                fi;
            od;

            LI:=1;

            GramMatrixT:=[];

            for j in [1..t] do
                GramMatrixT[j]:=GramMatrix[j]{[1..t]};
            od;

            x:=MAJORANA_PositiveDefinite(GramMatrixT);

            if x = -1 then
                Output[i]:=[StructuralCopy(Shape),"Error","Inner product not positive definite on A", StructuralCopy(GramMatrixT)];
                break;
            elif x = 0 then
                NullSpT:=MAJORANA_NullSpace(GramMatrixT);
                LI:=0;
            fi;

            if LI=0 then

                n:=Size(NullSpT);

                Display(n);

                T:=T{[1..t-n]};

                # Change Gram matrix to get rid of any axes not in basis

                GramMatrix:=GramMatrix{Concatenation([1..t-n],[t+1..dim])};

                for j in [1..dim-n] do
                    GramMatrix[j]:=GramMatrix[j]{Concatenation([1..t-n],[t+1..dim])};
                od;

                # Change alg products to get rid of any axes not in the basis

                AlgebraProducts:=AlgebraProducts{Concatenation([1..t-n],[t+1..dim])};

                for j in [1..dim-n] do
                    AlgebraProducts[j]:=AlgebraProducts[j]{Concatenation([1..t-n],[t+1..dim])};
                od;

                for j in [1..n] do
                    for k in [1..dim-n] do
                        for l in [1..dim-n] do
                            AlgebraProducts[k][l]:=AlgebraProducts[k][l] - NullSpT[j]*AlgebraProducts[k][l][dim-j+1];
                        od;
                    od;
                od;

                for j in [1..dim] do
                    for k in [1..dim] do
                        AlgebraProducts[j][k]:=AlgebraProducts[j][k]{Concatenation([1..t-n],[t+1..dim])};
                    od;
                od;

                # Change evecs to get rid of any axes not in the basis

                for j in [1..t] do
                    for k in [1..3] do
                        for l in [1..Size(EigenVectors[j][k])] do
                            for m in [1..Size(NullSpT)] do
                                EigenVectors[j][k][l]:=EigenVectors[j][k][l] - NullSpT[m]*EigenVectors[j][k][l][dim - m+1];
                            od;
                            EigenVectors[j][k][l]:=EigenVectors[j][k][l]{Concatenation([1..t-n],[t+1..dim])};
                        od;
                    od;
                od;

                t:=t-n;

                dim:=dim-n;



            else
                dim:=t+u+v+w;
            fi;



                                        ## STEP 4: MORE PRODUCTS ##

            # (2,3) values

            for j in [1..t] do
                for k in [1 .. u] do
                    x:=T[j]; h:=3Aaxes[k].1;

                    if x*h in T then
                        # D6

                        GramMatrix[j][k+t] := 1/4;
                        GramMatrix[k+t][j] := 1/4;

                        # Products and eigenvalues were done in 2,2 case

                        Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);
                        Append(KnownAlgebraProducts,[[j,k+t],[k+t,j]]);

                    elif Order(x*h)=3 then

                        # Case (2A,3A) in IPSS10

                        GramMatrix[j][k+t] := 1/9;
                        GramMatrix[k+t][j] := 1/9;

                        xj:=Position(3Aaxes,Group(x*h*x));
                        xk:=Position(3Aaxes,Group(h*x*h*h*x*h*h));
                        xl:=Position(3Aaxes,Group(h*h*x*h*x*h));

                        AlgebraProducts[j][k+t][j]:=1/9;
                        AlgebraProducts[j][k+t][k+t]:=5/64;
                        AlgebraProducts[j][k+t][xj+t]:=3/64;
                        AlgebraProducts[j][k+t][xk+t]:=-1/16;
                        AlgebraProducts[j][k+t][xl+t]:=-1/16;
                        AlgebraProducts[k+t][j][j]:=1/9;
                        AlgebraProducts[k+t][j][k+t]:=5/64;
                        AlgebraProducts[k+t][j][xj+t]:=3/64;
                        AlgebraProducts[k+t][j][xk+t]:=-1/16;
                        AlgebraProducts[k+t][j][xl+t]:=-1/16;

                        Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);
                        Append(KnownAlgebraProducts,[[j,k+t],[k+t,j]]);

                    elif Order(x*h)=4 then
                        if (x*h)^2 in T then

                            # Case (2A,3A) in IPSS10

                            GramMatrix[j][k+t] := 1/36;
                            GramMatrix[k+t][j] := 1/36;

                            xik:=Position(T,h*h*x*h);
                            xil:=Position(T,h*x*h*h);
                            xjk:=Position(T,x*h*h*x*h*x);
                            xjl:=Position(T,x*h*x*h*h*x);
                            xkl:=Position(T,h*x*h*x*h*h*x*h*h);
                            xx:=Position(T,h*h*x*h*x*h*h);
                            xj:=Position(3Aaxes,Group(x*h*x));
                            xk:=Position(3Aaxes,Group(h*h*x*h*x*h));
                            xl:=Position(3Aaxes,Group(h*x*h*x*h*h));

                            AlgebraProducts[j][k+t][j]:=1/45;
                            AlgebraProducts[j][k+t][k+t]:=1/64;
                            AlgebraProducts[j][k+t][xik]:=-1/90;
                            AlgebraProducts[j][k+t][xil]:=-1/90;
                            AlgebraProducts[j][k+t][xjk]:=-1/90;
                            AlgebraProducts[j][k+t][xjl]:=-1/90;
                            AlgebraProducts[j][k+t][xkl]:=1/45;
                            AlgebraProducts[j][k+t][xx]:=-1/45;
                            AlgebraProducts[j][k+t][xj+t]:=-1/64;
                            AlgebraProducts[j][k+t][xk+t]:=1/64;
                            AlgebraProducts[j][k+t][xl+t]:=1/64;

                            AlgebraProducts[k+t][j][j]:=1/45;
                            AlgebraProducts[k+t][j][k+t]:=1/64;
                            AlgebraProducts[k+t][j][xik]:=-1/90;
                            AlgebraProducts[k+t][j][xil]:=-1/90;
                            AlgebraProducts[k+t][j][xjk]:=-1/90;
                            AlgebraProducts[k+t][j][xjl]:=-1/90;
                            AlgebraProducts[k+t][j][xkl]:=1/45;
                            AlgebraProducts[k+t][j][xx]:=-1/45;
                            AlgebraProducts[k+t][j][xj+t]:=-1/64;
                            AlgebraProducts[k+t][j][xk+t]:=1/64;
                            AlgebraProducts[k+t][j][xl+t]:=1/64;

                            Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);
                            Append(KnownAlgebraProducts,[[j,k+t],[k+t,j]]);
                        else
                            GramMatrix[j][k+t] := -13/8192;
                            GramMatrix[k+t][j] := -13/8192;

                            # Means we have to add extra vectors on, might need to do this at the end, this is case (2B,3A) in IPSS10 ******************************

                            Print("Not 2-closed, have a Maj rep of s4 of shape (2B,3A) for (2,3) value with [j,k] =",[j,k]);

                            Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);
                        fi;
                    else
                        l:=1;
                        while l<t+1 do
                            s:=T[l];
                            if s*h in T and Order(x*s*h) = 2 then
                                Append(KnownInnerProducts,[[j,k+t],[k+t,j]]); l:=t+1; m:=1;
                                if x*s*h in T then

                                    # Need to know values of (a_t,a_s) and (a_s, a_{tsh}) and (a_t, a_{sh^2})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [s,x*s*h] in Orbitals[m] then str2:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;

                                    # Use values to work out inner product

                                    GramMatrix[j][k+t] := (64/135)*(-2*MAJORANA_Form(str1) + 4*MAJORANA_Form(str2) + MAJORANA_Form(str3)) + 1/45;
                                    GramMatrix[k+t][j]:=GramMatrix[j][k+t];

                                    Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);

                                else

                                    # Need to know values of (a_t,a_s) and (a_t, a_{sh^2})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;

                                    # Use values to work out inner product

                                    GramMatrix[j][k+t] := (64/135)*(2*MAJORANA_Form(str1) + MAJORANA_Form(str3));
                                    GramMatrix[k+t][j]:=GramMatrix[j][k+t];

                                    Append(KnownInnerProducts,[[j,k+t],[k+t,j]]);
                                fi;
                            else
                                l:=l+1;
                            fi;
                        od;
                    fi;
                od;
            od;

            l:=1;

            # (2,4) values

            for j in [1..t] do
                for k in [1 .. v] do
                    x:=T[j]; h:=4Aaxes[k].1;
                    if x*h in T then

                        # D8

                        GramMatrix[j][k+t+u] := 3/8;
                        GramMatrix[k+t+u][j] := 3/8;

                        # Prod done in 2,2 case

                        Append(KnownInnerProducts,[[j,k+t+u],[k+t+u,j]]);
                        Append(KnownAlgebraProducts,[[j,k+t+u],[k+t+u,j]]);

                    else
                        l:=1;
                        while l<t+1 do
                            if T[l]*h in T and Order(x*T[l]*h) = 2 then

                                s:=T[l]; Append(KnownInnerProducts,[[j,k+t+u],[k+t+u,j]]); l:=t+1; m:=1;

                                if x*s*h in T then

                                    # Need to know values of (a_t,a_s) and (a_s, a_{tsh}), (a_t, a_{sh^2}) and (a_t, a_{sh^3})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str2:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                        od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [s,x*s*h] in Orbitals[m] then str4:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;

                                    # Use values to work out inner product

                                    GramMatrix[j][k+t+u] := ( - 5*MAJORANA_Form(str1) + MAJORANA_Form(str2) +MAJORANA_Form(str3) + 8*MAJORANA_Form(str4))/3 +1/24;
                                    GramMatrix[k+t+u][j]:=GramMatrix[j][k+t+u];

                                    Append(KnownInnerProducts,[[j,k+t+u],[k+t+u,j]]);

                                else

                                    # Need to know values of (a_t,a_s), (a_t, a_{sh^2}) and (a_t, a_{sh^3})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str2:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;

                                    # Use values to work out inner product

                                    GramMatrix[j][k+t+u] := (3*MAJORANA_Form(str1) + MAJORANA_Form(str2) + 2*MAJORANA_Form(str3))/3;
                                    GramMatrix[k+t+u][j]:=GramMatrix[j][k+t+u];

                                    Append(KnownInnerProducts,[[j,k+t+u],[k+t+u,j]]);
                                fi;
                            else
                                l:=l+1;
                            fi;
                        od;
                    fi;
                od;
            od;
            l:=1;

            # (2,5) values

            for j in [1..t] do
                for k in [1 .. w] do
                    x:=T[j]; h:=5AaxesFixed[k][1];



                    if x*h in T then

                        # D10

                        GramMatrix[j][k+t+u+v] := 0;
                        GramMatrix[k+t+u+v][j] := 0;

                        # Products done in 2,2 case

                        Append(KnownInnerProducts,[[j,k+t+u+v],[k+t+u+v,j]]);
                        Append(KnownAlgebraProducts,[[j,k+t+u+v],[k+t+u+v,j]]);

                    else
                        l:=1;
                        while l<t+1 do
                            if T[l]*h in T and Order(x*T[l]*h) = 2 then
                                s:=T[l]; Append(KnownInnerProducts,[[j,k+t+u+v],[k+t+u+v,j]]); l:=t+1; m:=1;
                                if x*s*h in T then

                                    # Need to know values of (a_t,a_s) and (a_s, a_{tsh}), (a_t, a_{sh^2}), (a_t, a_{sh^3}) and (a_t, a_{sh^4})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [s,x*s*h] in Orbitals[m] then str2:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h] in Orbitals[m] then str4:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h*h] in Orbitals[m] then str5:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;


                                    # Use values to work out inner product

                                    GramMatrix[j][k+t+u+v] := (13*MAJORANA_Form(str1) - MAJORANA_Form(str2))/8 + (MAJORANA_Form(str3) +MAJORANA_Form(str4) +MAJORANA_Form(str5))/128 -1/2048;
                                    GramMatrix[k+t+u+v][j]:=GramMatrix[j][k+t+u+v];

                                    Append(KnownInnerProducts,[[j,k+t+u+v],[k+t+u+v,j]]);

                                else

                                    # Need to know values of (a_t,a_s), (a_t, a_{sh^2}), (a_t, a_{sh^3}) and (a_t, a_{sh^4})

                                    while m < SizeOrbitals+1  do
                                        if [x,s] in Orbitals[m] then str1:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h] in Orbitals[m] then str3:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h] in Orbitals[m] then str4:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;
                                    while m < SizeOrbitals+1  do
                                        if [x,s*h*h*h*h] in Orbitals[m] then str5:=Shape[m]; m:=SizeOrbitals+1; else m:=m+1; fi;
                                    od; m:=1;

                                    # Use values to work out inner product

                                    GramMatrix[j][k+t+u+v] := (-3*MAJORANA_Form(str1) + MAJORANA_Form(str3) +MAJORANA_Form(str4) +MAJORANA_Form(str5))/128;
                                    GramMatrix[k+t+u+v][j]:=GramMatrix[j][k+t+u+v];

                                    Append(KnownInnerProducts,[[j,k+t+u+v],[k+t+u+v,j]]);

                                fi;
                            else
                                l:=l+1;
                            fi;
                        od;
                    fi;
                od;
            od;

            l:=1;

            # (3,3) values

            for j in [1..u] do
                for n in [j+1..u] do
                    h:=3Aaxes[j].1; k:=3Aaxes[n].1; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if x*h in T and x*k in T then
                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*k);
                            x4:=Position(T,x*k*k);

                            if [x1,t+n] in KnownInnerProducts and [x2,t+n] in KnownInnerProducts then

                                GramMatrix[t+j][t+n]:= 64*( -3*GramMatrix[x1][t+n] + GramMatrix[x2][t+n])/135 + 2048*(GramMatrix[x1][x3] + GramMatrix[x1][x4])/1215 + 16/243;
                                GramMatrix[t+n][t+j]:=GramMatrix[t+j][t+n];

                                Append(KnownInnerProducts,[[t+j,t+n],[t+n,t+j]]);

                                l:=t+1;
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            l:=1;

            # (3,4) values

            for j in [1..u] do
                for n in [1..v] do
                    h:=3Aaxes[j].1; k:=4Aaxes[n].1; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if x*h in T and x*k in T then

                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*k);
                            x4:=Position(T,x*k*k);
                            x5:=Position(T,x*k*k*k);

                            if ForAll([[x1,t+u+n],[x2,t+u+n],[x3,t+u+n]], x-> x in KnownInnerProducts) then

                                l:=t+1;

                                GramMatrix[t+j][t+u+n]:= 64*( -4*GramMatrix[x1][t+u+n] + GramMatrix[x2][t+u+n] +4*GramMatrix[x1][x3] +2*GramMatrix[x1][x4] + 4*GramMatrix[x1][x5])/135 + 127/270;
                                GramMatrix[t+u+n][t+j]:=GramMatrix[t+j][t+u+n];

                                Append(KnownInnerProducts,[[t+j,t+u+n],[t+u+n,t+j]]);
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            l:=1;

            # (3,5) values

            for j in [1..u] do
                for n in [1..w] do
                    h:=3Aaxes[j].1; k:=5AaxesFixed[n][1]; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if x*h in T and x*k in T then
                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*k);
                            x4:=Position(T,x*k*k);
                            x5:=Position(T,x*k*k*k);
                            x6:=Position(T,x*k*k*k*k);

                            if ForAll([[x1,t+u+v+n],[x2,t+u+v+n]], x-> x in KnownInnerProducts) then

                                l:=t+1;

                                GramMatrix[t+j][t+u+v+n]:= 64*( -5*GramMatrix[x1][t+u+v+n] + GramMatrix[x2][t+u+v+n])/135 + 7*(GramMatrix[x1][x3] - GramMatrix[x1][x4] - GramMatrix[x1][x5] + GramMatrix[x1][x6])/270;
                                GramMatrix[t+u+v+n][t+j]:=GramMatrix[t+j][t+u+v+n];

                                Append(KnownInnerProducts,[[t+j,t+u+v+n],[t+u+v+n,t+j]]);
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            l:=1;

            # (4,4) values

            for j in [1..v] do
                for n in [j+1..v] do
                    h:=4Aaxes[j].1; k:=4Aaxes[n].1; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if x*h in T and x*k in T then

                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*h*h*h);
                            x4:=Position(T,x*k);
                            x5:=Position(T,x*k*k);
                            x6:=Position(T,x*k*k*k);

                            if ForAll([[x1,t+u+n],[x2,t+u+n],[x3,t+u+n]],x -> x in KnownInnerProducts) then

                                l:=t+1;

                                GramMatrix[t+u+j][t+u+n]:= ( -9*GramMatrix[x1][t+u+n] + GramMatrix[x2][t+u+n] + GramMatrix[x3][t+u+n] + 8*GramMatrix[x1][x4] +4*GramMatrix[x1][x5] +8*GramMatrix[x1][x6])/3 + 1/6;
                                GramMatrix[t+u+n][t+u+j]:=GramMatrix[t+u+j][t+u+n];

                                Append(KnownInnerProducts,[[t+u+j,t+u+n],[t+u+n,t+u+j]]);
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            l:=1;

            # (4,5) values

            for j in [1..v] do
                for n in [1..w] do
                    h:=4Aaxes[j].1; k:=5AaxesFixed[n][1]; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if  x*h in T and x*k in T then
                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*h*h*h);
                            x4:=Position(T,x*k);
                            x5:=Position(T,x*k*k);
                            x6:=Position(T,x*k*k*k);
                            x7:=Position(T,x*k*k*k*k);

                            if ForAll([[x1,t+u+v+n],[x2,t+u+v+n],[x3,t+u+v+n]],x -> x in KnownInnerProducts) then

                                l:=t+1;

                                GramMatrix[t+u+j][t+u+v+n]:= ( -11*GramMatrix[x1][t+u+v+n] + GramMatrix[x2][t+u+v+n] + GramMatrix[x3][t+u+v+n])/3+ 7*(GramMatrix[x1][x4] - GramMatrix[x1][x5] - GramMatrix[x1][x6]+ GramMatrix[x1][x6])/192;
                                GramMatrix[t+u+v+n][t+u+j]:=GramMatrix[t+u+j][t+u+v+n];

                                Append(KnownInnerProducts,[[t+u+j,t+u+v+n],[t+u+v+n,t+u+j]]);
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            l:=1;

            # (5,5) values

            for j in [1..w] do
                for n in [j+1..w] do
                    h:=5AaxesFixed[j][1]; k:=5AaxesFixed[n][1]; l:=1;
                    while l < t+1 do
                        x:=T[l];
                        if  x*h in T and x*k in T then
                            x1:=Position(T,x*h);
                            x2:=Position(T,x*h*h);
                            x3:=Position(T,x*h*h*h);
                            x4:=Position(T,x*h*h*h*h);
                            x5:=Position(T,x*k);
                            x6:=Position(T,x*k*k);
                            x7:=Position(T,x*k*k*k);
                            x8:=Position(T,x*k*k*k*k);

                            if ForAll([[x1,t++v+n],[x2,t+u+v+n],[x3,t+u+v+n],[x4,t+u+v+n]],x -> x in KnownInnerProducts) then

                                l:=t+1;

                                GramMatrix[t+u+v+j][t+u+v+n]:= ( 25*GramMatrix[x1][t+u+v+n] + GramMatrix[x2][t+u+v+n] + GramMatrix[x3][t+u+v+n] + GramMatrix[x4][t+u+v+n])/4048+ 7*(GramMatrix[x1][x5] - GramMatrix[x1][x6] - GramMatrix[x1][x7]+ GramMatrix[x1][x8])/32;
                                GramMatrix[t+u+v+n][t+j+u+v]:=GramMatrix[t+u+v+j][t+u+v+n];

                                Append(KnownInnerProducts,[[t+u+v+j,t+u+v+n],[t+u+v+n,t+u+v+j]]);
                            else
                                l:=l+1;
                            fi;
                        else
                            l:=l+1;
                        fi;
                    od;
                od;
            od;

            # Redo 2,3 vals

            ### Not really sure what this is, need to work out the theory again

#           for j in [1..t] do
#               for k in [1..u] do
#                   if not [j,t+k] in KnownInnerProducts then
#                       x:=T[j]; h:=3Aaxes[k].1;
#                       l:=1;
#                       while l < t+1 do
#                           s:=T[l];
#                           if s*h in T then
#                               if s*x in T then
#
#                                   x1:=Position(T,s*x);;
#
#                                   if [x1,k] in KnownInnerProducts then
#
#                                       x2:=Position(T,s*h);
#                                       x3:=Position(T,s*h^2);
#
#                                       GramMatrix[j][k]:=2/45 + 32*(GramMatrix[j][x2] + GramMatrix[j][x3] + GramMatrix[x1][x2] + GramMatrix[x1][x3])/45 - GramMatrix[x1][k];
#                                       GramMatrix[k][j]:=GramMatrix[j][k];
#
#                                       Append(KnownInnerProducts,[[j,k],[k,j]]);
#
#                                       l:=t+1;
#                                   else
#                                       l:=l+1;
#                                   fi;
#                               elif Order(s*x) = 2 and not s*x in T then
#
#                                   x2:=Position(T,s*h);
#                                   x3:=Position(T,s*h^2);
#
#                                   GramMatrix[j][k]:= 32*(GramMatrix[x2][j] + GramMatrix[x3][j])/45;
#                                   GramMatrix[k][j]:=GramMatrix[j][k];
#
#                                   Append(KnownInnerProducts,[[j,k],[k,j]]);
#
#                                   l:=t+1;
#                               elif Order(s*x) = 3 then
#                                   if GramMatrix[UnknownInnerProducts[j][1]][l] = 13/256 then
#
#                                       x1:=Position(3Aaxes,Group(x*s));
#                                       x2:=Position(T,s*h);
#                                       x3:=Position(T,s*h^2);
#
#                                       if [x1,k] in KnownInnerProducts and [x2,k] in KnownInnerProducts and  [x3,k] in KnownInnerProducts then
#
#                                           GramMatrix[j][k]:=1/36 - 27*GramMatrix[x1][k]/64 +3*(GramMatrix[x1][x2]+GramMatrix[x1][x2])/10 +32*(GramMatrix[j][x2]+GramMatrix[j][x3])/45;
#                                           GramMatrix[k][j]:=GramMatrix[j][k];
#
#                                           Append(KnownInnerProducts,[[j,k],[k,j]]);
#
#                                           l:=t+1;
#                                       else
#                                           l:=l+1;
#                                       fi;
#                                   else
#                                       l:=l+1;
#                                   fi;
#                               else
#                                   l:=l+1;
#                               fi;
#                           else
#                               l:=l+1;
#                           fi;
#                       od;
#                   fi;
#               od;
#           od;


                                        ## STEP 5: MORE EVECS ##

            # Find linearly independent subsets of eigenvectors

            Dimensions:=[];

            for j in [1..t] do
                for k in [1..3] do
                    if Size(EigenVectors[j][k]) > 0 then
                        EigenVectors[j][k]:=ShallowCopy(BaseMat(EigenVectors[j][k]));
                    fi;
                od;
                Append(Dimensions,[Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3])]);
            od;

            # Use these eigenvectors and the fusion rules to find more

            switch:=0;

            if ForAll(Dimensions,x->x=dim-1) then
                switch:=1;
            fi;

            NewEigenVectors:=NullMat(t,3);

            for j in [1..t] do
                for k in [1..3] do
                    NewEigenVectors[j][k]:=[];
                od;
            od;

            while switch=0 do

                for j in [1..t] do

                    a:=NullMat(1,dim)[1]; a[j]:=1;
                    zeros:=NullMat(1,dim)[1];

                    # 1, x fusion is a waste of time because a_0 obviously just preserves the evectors!

                    # 0,0 fusion

                    for k in [1..Size(EigenVectors[j][1])] do
                        for l in [1..Size(EigenVectors[j][1])] do

                            x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][1][l],AlgebraProducts,KnownAlgebraProducts);

                            if x<> 0 then
                                if MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> zeros and MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> 0 then
                                    Output[i]:=[StructuralCopy(Shape),"Error","Fusion of 0,0 eigenvectors does not hold",StructuralCopy(j),StructuralCopy(EigenVectors[j][1][k]),StructuralCopy(EigenVectors[j][1][l]),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts)];
                                    break;
                                fi;
                                Add(NewEigenVectors[j][1],x);
                            fi;
                        od;
                    od;

                    if Output[i] <> [] then
                        break;
                    fi;

                    # 0,1/4 fusion

                    for k in [1..Size(EigenVectors[j][1])] do
                        for l in [1..Size(EigenVectors[j][2])] do

                            x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][2][l],AlgebraProducts,KnownAlgebraProducts);

                            if x<> 0 then
                                if MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> x/4 and MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> 0 then
                                    Output[i]:=[StructuralCopy(Shape),"Error","Fusion of 0,1/4 eigenvectors does not hold",StructuralCopy(j),StructuralCopy(EigenVectors[j][1][k]),StructuralCopy(EigenVectors[j][2][l]),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts)];
                                    break;
                                fi;
                                Add(NewEigenVectors[j][2],x);
                            fi;
                        od;
                    od;

                    if Output[i] <> [] then
                        break;
                    fi;

                    # 0,1/32 fusion

                    for k in [1..Size(EigenVectors[j][1])] do
                        for l in [1..Size(EigenVectors[j][3])] do

                            x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][3][l],AlgebraProducts,KnownAlgebraProducts);

                            if x<> 0 then
                                if MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> x/32 and MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> 0 then
                                    Output[i]:=[StructuralCopy(Shape),"Error","Fusion of 0,1/32 eigenvectors does not hold",StructuralCopy(j),StructuralCopy(EigenVectors[j][1][k]),StructuralCopy(EigenVectors[j][3][l]),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts)];
                                    break;
                                fi;
                                Add(NewEigenVectors[j][3],x);
                            fi;
                        od;
                    od;

                    if Output[i] <> [] then
                        break;
                    fi;

                    # 1/4,1/32 fusion

                    for k in [1..Size(EigenVectors[j][2])] do
                        for l in [1..Size(EigenVectors[j][3])] do

                            x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][3][l],AlgebraProducts,KnownAlgebraProducts);

                            if x<> 0 then
                                if MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> x/32 and MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> 0 then
                                    Output[i]:=[StructuralCopy(Shape),"Error","Fusion of 1/4,1/32 eigenvectors does not hold",StructuralCopy(j),StructuralCopy(EigenVectors[j][2][k]),StructuralCopy(EigenVectors[j][3][l]),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts)];
                                    break;;
                                fi;
                                Add(NewEigenVectors[j][3],x);
                            fi;
                        od;
                    od;

                    if Output[i] <> [] then
                        break;
                    fi;

                    # 1/4,1/4 Fusion

                    for k in [1..Size(EigenVectors[j][2])] do
                        for l in [1..Size(EigenVectors[j][2])] do

                            x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][2][l],AlgebraProducts,KnownAlgebraProducts);

                            if x <> 0 then
                                y:=MAJORANA_InnerProduct(a,x,GramMatrix,KnownInnerProducts);
                                if y <> 0 then
                                    z:=x-y*a;
                                    if MAJORANA_AlgebraProduct(a,z,AlgebraProducts,KnownAlgebraProducts) <> zeros and MAJORANA_AlgebraProduct(a,x,AlgebraProducts,KnownAlgebraProducts) <> 0 then
                                        Output[i]:=[StructuralCopy(Shape),"Error","Fusion of 1/4,1/4 eigenvectors does not hold",StructuralCopy(j),StructuralCopy(EigenVectors[j][1][k]),StructuralCopy(EigenVectors[j][1][l]),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts),StructuralCopy(KnownInnerProducts),StructuralCopy(GramMatrix)];
                                    break;
                                    fi;
                                    Add(NewEigenVectors[j][1],z);
                                fi;
                            fi;
                        od;
                    od;

                    if Output[i] <> [] then
                        break;
                    fi;

                    Append(EigenVectors[j][1],NewEigenVectors[j][1]);
                    Append(EigenVectors[j][2],NewEigenVectors[j][2]);
                    Append(EigenVectors[j][3],NewEigenVectors[j][3]);

                od;

                if Output[i] <> [] then
                    break;
                fi;

                NewDimensions:=[];

                for j in [1..t] do
                    for k in [1..3] do
                        if Size(EigenVectors[j][k]) > 0 then
                            EigenVectors[j][k]:=ShallowCopy(BaseMat(EigenVectors[j][k]));
                        fi;
                    od;
                    Append(NewDimensions,[Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3])]);
                od;

                if NewDimensions = Dimensions then
                    switch := 1;
                elif ForAll(NewDimensions,x->x=dim-1) then
                    switch := 1;
                else
                    Dimensions:=StructuralCopy(NewDimensions);
                fi;

            od;

            if Output[i] <> [] then
                break;
            fi;

                                        ## STEP 6: MORE INNER PRODUCTS ##

            UnknownInnerProducts:=[];

            for j in [1..dim] do
                for k in [j..dim] do
                    if not [j,k] in KnownInnerProducts then
                        Add(UnknownInnerProducts,[j,k]);
                    fi;
                od;
            od;

            # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products

            if Size(UnknownInnerProducts) > 0 then

                mat:=[];
                vec:=[];

                for j in [1..t] do

                    # 1- eigenvectors and 0-eigenvectors

                    for k in [1..Size(EigenVectors[j][1])] do

                        sum:=[];

                        Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                        for m in [1..j] do

                            if [m,j] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][1][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[m,j])]:=EigenVectors[j][1][k][m];
                            fi;

                        od;

                        for m in [j+1..dim] do

                            if [j,m] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][1][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[j,m])]:=EigenVectors[j][1][k][m];
                            fi;

                        od;

                        Add(vec,Sum(sum));

                    od;

                    # 1- eigenvectors and 1/4-eigenvectors

                    for k in [1..Size(EigenVectors[j][2])] do

                        sum:=[];

                        Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                        for m in [1..j] do

                            if [m,j] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][2][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[m,j])]:=EigenVectors[j][2][k][m];
                            fi;

                        od;

                        for m in [j+1..dim] do

                            if [j,m] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][2][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[j,m])]:=EigenVectors[j][2][k][m];
                            fi;

                        od;

                        Add(vec,Sum(sum));

                    od;

                    # 1- eigenvectors and 1/32-eigenvectors

                    for k in [1..Size(EigenVectors[j][3])] do

                        sum:=[];

                        Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                        for m in [1..j] do

                            if [m,j] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][3][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[m,j])]:=EigenVectors[j][3][k][m];
                            fi;

                        od;

                        for m in [j+1..dim] do

                            if [j,m] in KnownInnerProducts then
                                Add(sum,-EigenVectors[j][3][k][m]*GramMatrix[j][m]);
                            else
                                mat[Size(mat)][Position(UnknownInnerProducts,[j,m])]:=EigenVectors[j][3][k][m];
                            fi;

                        od;

                        Add(vec,Sum(sum));

                    od;

                    # 0-eigenvectors and 1/4-eigenvectors

                    for k in [1..Size(EigenVectors[j][1])] do
                        for l in [1..Size(EigenVectors[j][2])] do

                            sum:=[];
                            Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                            for m in [1..dim] do
                                for n in [1..m] do

                                    if [n,m] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][1][k][m]*EigenVectors[j][2][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[n,m])]:=EigenVectors[j][1][k][m]*EigenVectors[j][2][l][n];
                                    fi;

                                od;

                                for n in [m+1..dim] do

                                    if [m,n] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][1][k][m]*EigenVectors[j][2][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[m,n])]:=EigenVectors[j][1][k][m]*EigenVectors[j][2][l][n];
                                    fi;
                                od;
                            od;

                            Add(vec,Sum(sum));

                        od;
                    od;

                    # 0-eigenvectors and 1/32-eigenvectors

                    for k in [1..Size(EigenVectors[j][1])] do
                        for l in [1..Size(EigenVectors[j][3])] do

                            sum:=[];
                            Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                            for m in [1..dim] do
                                for n in [1..m] do

                                    if [n,m] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][1][k][m]*EigenVectors[j][3][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[n,m])]:=EigenVectors[j][1][k][m]*EigenVectors[j][3][l][n];
                                    fi;

                                od;

                                for n in [m+1..dim] do

                                    if [m,n] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][1][k][m]*EigenVectors[j][3][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[m,n])]:=EigenVectors[j][1][k][m]*EigenVectors[j][3][l][n];
                                    fi;
                                od;
                            od;

                            Add(vec,Sum(sum));

                        od;
                    od;

                    # 1/4-eigenvectors and 1/32-eigenvectors

                    for k in [1..Size(EigenVectors[j][2])] do
                        for l in [1..Size(EigenVectors[j][3])] do

                            sum:=[];
                            Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                            for m in [1..dim] do
                                for n in [1..m] do

                                    if [n,m] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][2][k][m]*EigenVectors[j][3][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[n,m])]:=EigenVectors[j][2][k][m]*EigenVectors[j][3][l][n];
                                    fi;

                                od;

                                for n in [m+1..dim] do

                                    if [m,n] in KnownInnerProducts then
                                        Add(sum,-EigenVectors[j][2][k][m]*EigenVectors[j][3][l][n]*GramMatrix[m][n]);
                                    else
                                        mat[Size(mat)][Position(UnknownInnerProducts,[m,n])]:=EigenVectors[j][2][k][m]*EigenVectors[j][3][l][n];
                                    fi;
                                od;
                            od;

                            Add(vec,Sum(sum));

                        od;
                    od;
                od;

                Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                if Size(Solution) = 2 then
                    if Size(Solution[2])>0 then
                        Output[i]:=[StructuralCopy(Shape),"Fail","Missing inner product values",StructuralCopy(KnownInnerProducts),StructuralCopy(GramMatrix)];
                    else
                        for k in [1..Size(Solution[1])] do
                            x:=UnknownInnerProducts[k][1]; y:=UnknownInnerProducts[k][2];
                            GramMatrix[x][y]:=Solution[1][k];
                            GramMatrix[y][x]:=Solution[1][k];
                            Append(KnownInnerProducts,[[x,y],[y,x]]);
                        od;
                    fi;
                else
                    Output[i]:=[StructuralCopy(Shape),"Error","Inconsistent system of unknown inner products"];
                fi;
            fi;

            if Size(Output[i]) > 0 then
                break;
            fi;

            # Check that GramMatrix matrix is pd

            L:=MAJORANA_LDLTDecomposition(GramMatrix);

            Diagonals:=[];

            for j in [1..Size(GramMatrix)] do
                Append(Diagonals,[L[2][j][j]]);
            od;

            if ForAny(Diagonals, x->x<0) then
                Output[i]:=[StructuralCopy(Shape),"Error","The inner product is not positive definite",StructuralCopy(3Aaxes), StructuralCopy(4Aaxes), StructuralCopy(5Aaxes), StructuralCopy(5AaxesFixed), StructuralCopy(GramMatrix)];
                break;
            elif ForAny(Diagonals, x->x=0) then
                NullSp:=MAJORANA_NullSpace(GramMatrix);
                LI:=0;
            else
                LI:=1;
            fi;

            if LI=0 then

                dim:=t+u+v+w-Size(NullSp);

                for j in [1..Size(NullSp)] do
                    NullSp[j]:=NullSp[j]/NullSp[j][dim+Size(NullSp)-j+1];

                    for k in [1..j-1] do
                        NullSp[j]:=NullSp[j] - NullSp[j][n-k+1]*NullSp[k];
                    od;
                od;

                # Change alg products to get rid of any axes not in the basis

                AlgebraProducts:=AlgebraProducts{[1..dim]};

                for j in [1..dim] do
                    AlgebraProducts[j]:=AlgebraProducts[j]{[1..dim]};
                od;

                for j in [1..Size(NullSp)] do
                    for k in [1..dim] do
                        for l in [1..dim] do
                            AlgebraProducts[k][l]:=AlgebraProducts[k][l] - NullSp[j]*AlgebraProducts[k][l][dim+Size(NullSp)-j+1];
                        od;
                    od;
                od;

                for j in [1..dim] do
                    for k in [1..dim] do
                        AlgebraProducts[j][k]:=AlgebraProducts[j][k]{[1..dim]};
                    od;
                od;

                # Change evecs to get rid of any axes not in the basis

                for j in [1..t] do
                    for k in [1..3] do
                        for l in [1..Size(EigenVectors[j][k])] do
                            for m in [1..Size(NullSp)] do
                                EigenVectors[j][k][l]:=EigenVectors[j][k][l] - NullSp[m]*EigenVectors[j][k][l][dim+ Size(NullSp) - m+1];
                            od;
                            EigenVectors[j][k][l]:=EigenVectors[j][k][l]{[1..dim]};
                        od;
                    od;
                od;
            else
                dim:=t+u+v+w;
            fi;

                                        ## STEP 7: MORE PRODUCTS II ##

            # Check fusion and M1

            ErrorM1:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts);

            if Size(ErrorM1)>0 then
                Output[i]:=[StructuralCopy(Shape),"Error","Algebra does not obey axiom M1 step 7",StructuralCopy(GramMatrix),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts),StructuralCopy(ErrorM1)];
            fi;

            ErrorFusion:=MAJORANA_Fusion(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors);

            if ForAny(ErrorFusion, x->Size(x) > 0) then
                Output[i]:=[StructuralCopy(Shape),"Error","Algebra does not obey fusion rules step 7",StructuralCopy(GramMatrix),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts),StructuralCopy(EigenVectors),StructuralCopy(ErrorFusion)];
                break;
            fi;


            # Use eigenvectors to find more products

            for j in [1..t] do

                UnknownAlgebraProducts:=[];

                for k in [t+1..dim] do
                    if not [j,k] in KnownAlgebraProducts then
                        Add(UnknownAlgebraProducts,[j,k]);
                    fi;
                od;

                if Size(UnknownAlgebraProducts) > 0 then

                    mat:=NullMat(Size(Union(EigenVectors[j][1],EigenVectors[j][2],EigenVectors[j][3])),Size(UnknownAlgebraProducts));
                    vec:=NullMat(Size(mat),dim);

                    for k in [1..Size(EigenVectors[j][1])] do

                        sum:=[];

                        for l in [1..dim] do
                            if EigenVectors[j][1][k][l] <> 0 then
                                if [j,l] in KnownAlgebraProducts then
                                        Add(sum, EigenVectors[j][1][k][l]*AlgebraProducts[j][l]);
                                else
                                    mat[k][Position(UnknownAlgebraProducts,[j,l])] := EigenVectors[j][1][k][l];
                                fi;
                            fi;
                        od;

                        x:=Sum(sum);

                        if x <> 0 then
                            vec[k]:= -x;
                        fi;
                    od;

                    for k in [1..Size(EigenVectors[j][2])] do

                        sum:=[];

                        for l in [1..dim] do
                            if EigenVectors[j][2][k][l] <> 0 then
                                if [j,l] in KnownAlgebraProducts then
                                    Add(sum, EigenVectors[j][2][k][l]*AlgebraProducts[j][l]);
                                else
                                    mat[k+Size(EigenVectors[j][1])][Position(UnknownAlgebraProducts,[j,l])] := EigenVectors[j][2][k][l];
                                fi;
                            fi;
                        od;

                        vec[k+Size(EigenVectors[j][1])]:= (-1)*Sum(sum) + EigenVectors[j][2][k]/4;

                    od;

                    for k in [1..Size(EigenVectors[j][3])] do

                        sum:=[];

                        for l in [1..dim] do
                            if EigenVectors[j][3][k][l] <> 0 then
                                if [j,l] in KnownAlgebraProducts then
                                    Add(sum, EigenVectors[j][3][k][l]*AlgebraProducts[j][l]);
                                else
                                    mat[k+Size(EigenVectors[j][1])+Size(EigenVectors[j][2])][Position(UnknownAlgebraProducts,[j,l])] := EigenVectors[j][3][k][l];
                                fi;
                            fi;
                        od;

                        vec[k+Size(EigenVectors[j][1])+Size(EigenVectors[j][2])]:= (-1)*Sum(sum) + EigenVectors[j][3][k]/32;
                    od;


                    Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                    if Size(Solution) = 2 then
                            for k in [1..Size(Solution[1])] do
                                if not k in Solution[2] then

                                    x:=UnknownAlgebraProducts[k][1]; y:=UnknownAlgebraProducts[k][2];

                                    AlgebraProducts[x][y]:=Solution[1][k];
                                    AlgebraProducts[y][x]:=Solution[1][k];

                                    Append(KnownAlgebraProducts,[[x,y],[y,x]]);
                                fi;
                            od;
                    else
                        Output[i]:=[StructuralCopy(Shape),"Error","Inconsistent system of unknown algebra products step 7",StructuralCopy(GramMatrix),StructuralCopy(KnownAlgebraProducts),StructuralCopy(AlgebraProducts),StructuralCopy(EigenVectors),StructuralCopy(mat),StructuralCopy(vec),StructuralCopy(Solution),StructuralCopy(UnknownAlgebraProducts)];
                        break;
                    fi;
                fi;
            od;

            if Size(Output[i])>0 then
                break;
            fi;

                                        ## STEP 8: RESURRECTION PRINCIPLE I ##

            # Use resurrection principle (Step 7 Seress)

            UnknownAlgebraProducts:=[];

            for j in [1..t] do
                for k in [j..dim] do
                    if not [j,k] in KnownAlgebraProducts then
                        Append(UnknownAlgebraProducts,[[j,k]]);
                    fi;
                od;
            od;

            if Size(UnknownAlgebraProducts) > 0 then

                mat:=[];
                vec:=[];
                record:=[];

                for j in [1..t] do

                    a:=NullMat(1,dim)[1];
                    a[j]:=1;

                    Alpha:=[];
                    Alpha2:=[];
                    Beta2:=[];

                    for k in [1..Size(EigenVectors[j][1])] do
                        if Size(Positions(EigenVectors[j][1][k]{[t+1..dim]},0)) = dim-t-1 then
                            Append(Alpha,[k]);
                        elif Size(Positions(EigenVectors[j][1][k]{[t+1..dim]},0)) = dim-t then
                            Append(Alpha2,[k]);
                        fi;
                    od;

                    for k in [1..Size(EigenVectors[j][2])] do
                        if Size(Positions(EigenVectors[j][2][k]{[t+1..dim]},0)) = dim-t then
                            Append(Beta2,[k]);
                        fi;
                    od;

                    for k in Alpha do

                        Beta:=[];

                        l:=t+1;

                        while l <= dim do
                            if EigenVectors[j][1][k][l] <> 0 then
                                c:=l;
                                l:=dim+1;
                            else
                                l:=l+1;
                            fi;
                        od;

                        for l in [1..Size(EigenVectors[j][2])] do
                            if Size(Positions(EigenVectors[j][2][l]{[t+1..dim]},0)) = dim-t-1 and EigenVectors[j][2][l][c] <>0 then
                                Append(Beta,[l]);
                            fi;
                        od;

                        EigenVectors[j][1][k]:=EigenVectors[j][1][k]/EigenVectors[j][1][k][c];

                        for m in Beta do

                            EigenVectors[j][2][m]:=EigenVectors[j][2][m]/EigenVectors[j][2][m][c];

                            walpha:=Concatenation(EigenVectors[j][1][k]{[1..t]},NullMat(1,dim-t)[1]);
                            wbeta:=Concatenation(EigenVectors[j][2][m]{[1..t]},NullMat(1,dim-t)[1]);

                            ## first part

                            for l in Alpha2 do

                                sum:=[];
                                row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                # calculate lhs

                                x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][l],(walpha - wbeta),AlgebraProducts,KnownAlgebraProducts);

                                for n in [1..dim] do
                                    if x[n] <> 0 then
                                        if [j,n] in KnownAlgebraProducts then
                                            Append(sum,-[AlgebraProducts[j][n]*x[n]]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[j,n])]:=x[n];
                                        fi;
                                    fi;
                                od;

                                # calculate rhs

                                for n in [1..t] do
                                    if EigenVectors[j][1][l][n] <> 0 then
                                        if [n,c] in KnownAlgebraProducts then
                                            Append(sum,[-EigenVectors[j][1][l][n]*AlgebraProducts[n][c]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[n,c])]:=row[Position(UnknownAlgebraProducts,[n,c])] + EigenVectors[j][1][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                Append(mat,[row]);
                                Append(vec,[Sum(sum) - MAJORANA_AlgebraProduct(wbeta,EigenVectors[j][1][l],AlgebraProducts,KnownAlgebraProducts)/4]);
                                Append(record,[[1,j,k,l,m]]);

                            od;

                            ## Second part

                            for l in Beta2 do

                                sum:=[];
                                row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                # calculate lhs

                                x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][l],(walpha - wbeta),AlgebraProducts,KnownAlgebraProducts);

                                for n in [1..dim] do
                                    if x[n] <> 0 then
                                        if [j,n] in KnownAlgebraProducts then
                                            Append(sum,-[AlgebraProducts[j][n]*x[n]]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[j,n])]:=x[n];
                                        fi;
                                    fi;
                                od;

                                # calculate rhs

                                for n in [1..t] do
                                    if EigenVectors[j][2][l][n] <> 0 then
                                        if [n,c] in KnownAlgebraProducts then
                                            Append(sum,[EigenVectors[j][2][l][n]*AlgebraProducts[n][c]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[n,c])]:=row[Position(UnknownAlgebraProducts,[n,c])] + EigenVectors[j][2][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                x:= MAJORANA_InnerProduct(EigenVectors[j][2][m],EigenVectors[j][2][l],GramMatrix,KnownInnerProducts);

                                if x<>0 then
                                    Append(mat,[row]);
                                    Append(vec,[Sum(sum)+MAJORANA_AlgebraProduct(EigenVectors[j][2][l],walpha,AlgebraProducts,KnownAlgebraProducts)/4 - x*a/4]);
                                    Append(record,[[2,j,k,l,m]]);
                                fi;
                            od;
                        od;
                    od;
                od;

                Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                if Size(Solution)  = 2 then
                    if Size(Solution[2]) = 0 then
                        for k in [1..Size(UnknownAlgebraProducts)] do
                            x:=UnknownAlgebraProducts[k][1]; y:=UnknownAlgebraProducts[k][2];
                            AlgebraProducts[x][y]:=Solution[1][k];
                            AlgebraProducts[y][x]:=Solution[1][k];
                            Append(KnownAlgebraProducts,[[x,y],[y,x]]);
                        od;
                    else
                        Output[i]:=[Shape,"Fail","Missing algebra product values",GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors];
                        Output[i]:=StructuralCopy(Output[i]);
                        break;
                    fi;
                else
                    Output[i]:=[Shape,"Error","Inconsistent system of unknown algebra products"];
                    Output[i]:=StructuralCopy(Output[i]);
                fi;
            fi;

                                        ## STEP 9: MORE EVECS II ##

            # Step 8 - check if we have full espace decomp, if not find it

            for j in [1..t] do
                if Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3]) + 1 <> dim then
                    mat:=[];

                    for k in [1..dim] do
                        Append(mat,[AlgebraProducts[j][k]]);
                    od;

                    mat:=TransposedMat(mat);

                    EigenVectors[j][1]:=MAJORANA_NullSpace(mat);
                    EigenVectors[j][2]:=MAJORANA_NullSpace(mat - IdentityMat(dim)/4);
                    EigenVectors[j][3]:=MAJORANA_NullSpace(mat - IdentityMat(dim)/32);
                    EigenVectors[j][4]:=MAJORANA_NullSpace(mat - IdentityMat(dim) );

                    if Size(EigenVectors[j][4]) <> 1 then
                        Output[i]:=[Shape,"Error","Algebra does not obey axiom M5",GramMatrix,AlgebraProducts,EigenVectors];
                        Output[i]:=StructuralCopy(Output[i]);
                        break;
                    elif Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3]) + Size(EigenVectors[j][4]) <> dim then
                        Output[i]:=[Shape,"Error","Algebra does not obey axiom M4",GramMatrix,AlgebraProducts,EigenVectors];
                        Output[i]:=StructuralCopy(Output[i]);
                        break;
                    fi;
                fi;
            od;

            if Size(Output[i]) > 0 then
                break;
            fi;

                                        ## STEP 10: RESURRECTION PRINCIPLE II ##

            # Check that eigenvectors obey the fusion rules

            ErrorFusion:=MAJORANA_Fusion(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors);

            if ForAny(ErrorFusion,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Algebra does not obey fusion rules",GramMatrix,AlgebraProducts,EigenVectors,ErrorFusion];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Use resurrection principle

            UnknownAlgebraProducts:=[];

            for j in [t+1..dim] do
                for k in [j+1..dim] do
                    if not [j,k] in KnownAlgebraProducts then
                        Append(UnknownAlgebraProducts,[[j,k]]);
                    fi;
                od;
            od;

            if Size(UnknownAlgebraProducts) > 0 then
                mat:=[];
                vec:=[];
                record:=[];

                for j in [1..t] do

                    a:=NullMat(1,dim)[1];
                    a[j]:=1;

                    Alpha:=[];
                    Alpha2:=[];
                    Beta2:=[];

                    for k in [1..Size(EigenVectors[j][1])] do
                        if Size(Positions(EigenVectors[j][1][k]{[t+1..dim]},0)) = dim-t-1 then
                            Append(Alpha,[k]);
                        fi;
                        Append(Alpha2,[k]);
                    od;

                    for k in [1..Size(EigenVectors[j][2])] do
                        Append(Beta2,[k]);
                    od;

                    for k in Alpha do

                        Beta:=[];

                        l:=t+1;

                        while l <= dim do
                            if EigenVectors[j][1][k][l] <> 0 then
                                c:= l;
                                l:=dim+1;
                            else
                                l:=l+1;
                            fi;
                        od;

                        for l in [1..Size(EigenVectors[j][2])] do
                            if Size(Positions(EigenVectors[j][2][l]{[t+1..dim]},0)) = dim-t-1 and EigenVectors[j][2][l][c] <>0 then
                                Append(Beta,[l]);
                            fi;
                        od;

                        EigenVectors[j][1][k]:=EigenVectors[j][1][k]/EigenVectors[j][1][k][c];

                        for m in Beta do

                            EigenVectors[j][2][m]:=EigenVectors[j][2][m]/EigenVectors[j][2][m][c];

                            walpha:=Concatenation(EigenVectors[j][1][k]{[1..t]},NullMat(1,dim-t)[1]);
                            wbeta:=Concatenation(EigenVectors[j][2][m]{[1..t]},NullMat(1,dim-t)[1]);

                            ## first part

                            for l in Alpha2 do

                                sum:=[];
                                row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                # calculate unknowns

                                for n in [1..c] do
                                    if EigenVectors[j][1][l][n] <> 0 then
                                        if [n,c] in KnownAlgebraProducts then
                                            Append(sum,-[AlgebraProducts[c][n]*EigenVectors[j][1][l][n]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[n,c])]:=EigenVectors[j][1][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                for n in [c+1..dim] do
                                    if EigenVectors[j][1][l][n] <> 0 then
                                        if [c,n] in KnownAlgebraProducts then
                                            Append(sum,-[AlgebraProducts[c][n]*EigenVectors[j][1][l][n]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[c,n])]:=EigenVectors[j][1][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                # calculate knowns

                                Append(sum,-[MAJORANA_AlgebraProduct(a,MAJORANA_AlgebraProduct(EigenVectors[j][1][l],(walpha - wbeta),AlgebraProducts,KnownAlgebraProducts),AlgebraProducts,KnownAlgebraProducts)]);

                                Append(sum,-[MAJORANA_AlgebraProduct(EigenVectors[j][1][l],wbeta,AlgebraProducts,KnownAlgebraProducts)]/4);

                                Append(mat,[row]);
                                Append(vec,[Sum(sum)]);
                                Append(record,[[1,j,k,l,m]]);

                            od;

                            ## Second part

                            for l in Beta2 do

                                sum:=[];
                                row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                # calculate unknowns

                                for n in [1..c] do
                                    if EigenVectors[j][2][l][n] <> 0 then
                                        if [n,c] in KnownAlgebraProducts then
                                            Append(sum,[AlgebraProducts[c][n]*EigenVectors[j][2][l][n]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[n,c])]:=-EigenVectors[j][2][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                for n in [c+1..dim] do
                                    if EigenVectors[j][2][l][n] <> 0 then
                                        if [c,n] in KnownAlgebraProducts then
                                            Append(sum,[AlgebraProducts[c][n]*EigenVectors[j][2][l][n]/4]);
                                        else
                                            row[Position(UnknownAlgebraProducts,[c,n])]:=-EigenVectors[j][2][l][n]/4;
                                        fi;
                                    fi;
                                od;

                                # calculate knowns

                                Append(sum,-[MAJORANA_AlgebraProduct(a,MAJORANA_AlgebraProduct(EigenVectors[j][2][l],(walpha - wbeta),AlgebraProducts,KnownAlgebraProducts),AlgebraProducts,KnownAlgebraProducts)]);

                                Append(sum,[MAJORANA_AlgebraProduct(EigenVectors[j][2][l],walpha,AlgebraProducts,KnownAlgebraProducts)]/4);

                                x:= MAJORANA_InnerProduct(EigenVectors[j][2][m],EigenVectors[j][2][l],GramMatrix,KnownInnerProducts);

                                Append(sum,-[a*x/4]);

                                if x<>0 then
                                    Append(mat,[row]);
                                    Append(vec,[Sum(sum)]);
                                    Append(record,[[2,j,k,l,m]]);
                                fi;
                            od;
                        od;
                    od;
                od;

                Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                if Size(Solution)  = 2 then
                    if Size(Solution[2]) = 0 then
                        for k in [1..Size(UnknownAlgebraProducts)] do

                            x:=UnknownAlgebraProducts[k][1]; y:=UnknownAlgebraProducts[k][2];

                            AlgebraProducts[x][y]:=Solution[1][k];
                            AlgebraProducts[y][x]:=Solution[1][k];

                            Append(KnownAlgebraProducts,[[x,y],[y,x]]);
                        od;
                    else
                        Output[i]:=[Shape,"Fail","Missing algebra products",GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors];
                        Output[i]:=StructuralCopy(Output[i]);
                        break;
                    fi;
                else
                    Output[i]:=[Shape,"Error","Inconsistent system of unknown algebra products",GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors];
                    Output[i]:=StructuralCopy(Output[i]);
                    break;
                fi;
            fi;

                                        ## STEP 11: CHECK ALGEBRA ##

            # Check bilinear form is positive definite

            if MAJORANA_PositiveDefinite(GramMatrix) <0 then
                Output[i]:=[Shape,"Error","Gram Matrix is not positive definite",GramMatrix, AlgebraProducts, EigenVectors];
                Output[i]:=StructuralCopy(Output[i]);
            fi;

            # Check that all triples obey axiom M1

            ErrorM1:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts);

            if Size(ErrorM1)>0 then
                Output[i]:=[Shape,"Error","Algebra does not obey axiom M1",GramMatrix,AlgebraProducts,ErrorM1];
                Output[i]:=StructuralCopy(Output[i]);
            fi;

            # Check that eigenvectors obey the fusion rules

            ErrorFusion:=MAJORANA_Fusion(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors);

            if ForAny(ErrorFusion,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Algebra does not obey fusion rules",GramMatrix,AlgebraProducts,EigenVectors,ErrorFusion];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Check that the eigenspaces are orthogonal

            ErrorOrthogonality:=MAJORANA_Orthogonality(T,KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts,EigenVectors);

            if ForAny(ErrorOrthogonality,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Eigenspaces are not orthogonal with respect to the inner product",GramMatrix,AlgebraProducts,EigenVector,ErrorOrthogonality];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Check M2

            ErrorM2:=MAJORANA_AxiomM2(KnownInnerProducts,GramMatrix,KnownAlgebraProducts,AlgebraProducts);

            if ErrorM2 = -1 then
                Output[i]:=[Shape,"Error","Algebra does not obey axiom M2",GramMatrix,AlgebraProducts,ErrorM2];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            Output[i]:=[Shape,"Success",GramMatrix,AlgebraProducts,EigenVectors];
            Output[i]:=StructuralCopy(Output[i]);

            master:=0;
        od;
    od;

    return Output;

    end

    );






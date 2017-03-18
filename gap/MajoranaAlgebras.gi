#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

# Temporary replacement for KnownAlgebraProducts. To be removed
# after all functions have been cleaned up.
BindGlobal( "MAJORANA_ExtractUnknownAlgebraProducts",
function(AlgebraProducts)
    local table, i, j;
    
    table := [];

    for i in [1..Size(AlgebraProducts)] do
		if AlgebraProducts[i] = false then
			Add(table, i);
		fi;
    od;
    return AsSet(table);
end);


BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4, [1,0], 1/32]
            ,[ 1/32, 1/32,  1/32, []] ]);

# This is the test function for fusion

InstallGlobalFunction( MAJORANA_TestFusion,
function(a, b, j, Shape, AlgebraProducts, EigenVectors, GramMatrix, ProductList, dim)
    local u, zeros, x, y, z, k, l, NewEigenVectors, ev_a, ev_b, ev, table;

    table := [0, 1/4, 1/32];
    ev := MAJORANA_FusionTable[a+1][b+1];

    zeros := [1..dim] * 0;
    u := ShallowCopy(zeros);
    u[j] := 1;
    NewEigenVectors := [];

    ev_a := EigenVectors[j][a];
    ev_b := EigenVectors[j][b];

    # the 1/4,1/4 case is special
    if (a=2) and (b=2) then
        for k in [1..Size(ev_a)] do
            for l in [1..Size(ev_b)] do

                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );

                if x <> false then
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList[3]);
                    if y <> false then
                        z := x - y*u;

                        if MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList) <> false and
                           MAJORANA_AlgebraProduct( u, z, AlgebraProducts, ProductList ) <> zeros then
                            return [false, StructuralCopy([ Shape
                                           , "Error"
                                           , STRINGIFY( "Fusion of ",
                                                        table[a], ",", table[b],
                                                        " eigenvectors does not hold" )
                                           , j
                                           , ev_a[k]
                                           , ev_b[l]
                                           , AlgebraProducts
                                           , GramMatrix ])];
                        else
                            Add(NewEigenVectors, z);
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
                
                    y := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList );

                    if y <> false and y <> ev * x then

                        return [false, StructuralCopy([ Shape
                                       , "Error"
                                       , STRINGIFY( "Fusion of ",
                                                    table[a], ",", table[b],
                                                    " eigenvectors does not hold" )
                                       , j
                                       , ev_a[k]
                                       , ev_b[l]
                                       , AlgebraProducts
                                       , ProductList ])];
                    else
                        Add(NewEigenVectors, x);
                    fi;
                fi;
            od;
        od;
    fi;

    return [true, NewEigenVectors];
end);

InstallGlobalFunction(MAJORANA_NullSpace,

        function(mat) # Takes as input matrix, returns a matrix whose rows form a basis of the nullspace of mat

        local A,B;

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

function(mat,vecs) # Takes as input two matrices, the second being interpreted as a vector of vectors. Returns a list of size four if system is inconsistent, otherwise returns a list of size 4

        local A, C, n, m, d, absd, B, i, j, k, x, imax, temp, tempv, tempi, sol, list, newmat, newvec, pos, p, unsolved, zeros, error;

        A:=StructuralCopy(mat);
        B:=StructuralCopy(vecs);

        i:=1;

        while i < Size(A) do
            if ForAll(A[i], x-> x =0) and ForAll(B[i], x-> x = 0) then
                Remove(A,i);
                Remove(B,i);
            else
                i:=i+1;
            fi;
        od;

        n:=Size(A);
        m:=Size(A[1]);

        if n<m then
            A{[n+1..m]} := NullMat(m-n,m);
            B{[n+1..m]} := NullMat(m-n,Size(vecs[1]));
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


InstallGlobalFunction(  MAJORANA_AlgebraProduct,

        function(u,v,AlgebraProducts,list) # If all the relevant products are known, returns the algebra product of u and v. If not, returns 0

        # list should be of the form [coordinates,longcoordinates,pairorbitlist,pairconjelements,positionlist,NullSp]

        local i, j, k, s, sum, dim, vec, x;

        sum:=[];
        dim:=Size(u);
        vec:=[1..dim]*0;

        if ForAll(u,x-> x= 0 ) or ForAll(v,x->x=0) then
            return u*0;
        fi;

        for i in [1..dim] do
			if u[i] <> 0 then 
				for j in [1..dim] do
					if v[j] <> 0 then 
						x:=AlgebraProducts[list[3][i][j]];
						if x <> false then
							for k in [1..dim] do
								s:=list[5][Position(list[2],list[1][k]^list[4][i][j])];
								if s > 0 then
									vec[s]:=u[i]*v[j]*x[k];
								else
									vec[-s]:=-u[i]*v[j]*x[k];
								fi;
							od;
							sum:=sum + vec;
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
        
        sum := MAJORANA_RemoveNullSpace(sum, list[6]);
                
        return sum;
        
        end

        );

InstallGlobalFunction(  MAJORANA_InnerProduct,

    function(u,v,GramMatrix, pairorbitlist) # If all the relevant products are known, returns the algebra product of u and v. If not, returns [0]

        local i, j, sum;

        sum:=[];

        for i in [1..Size(u)] do
            if u[i] <> 0 then
                for j in [1..Size(v)] do
                    if v[j] <> 0 then
                        if GramMatrix[pairorbitlist[i][j]] <> false then
                            Append(sum,[u[i]*v[j]*GramMatrix[pairorbitlist[i][j]]]);
                        else

                            # cannot calculate product

                            return false;
                        fi;
                    fi;
                od;
            fi;
        od;
        return Sum(sum);
        end );

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

        function(GramMatrix,AlgebraProducts,Orbitals,list,pairrepresentatives) # Checks if bilinear and algebra products obey axiom M1, outputs a list which is empty if they do obey the axiom

		# list should be of the form [coordinates,longcoordinates,pairorbitlist,pairconjelements,positionlist]

        local ErrorM1, j, k, p, dim, x, y, u, w,v;

        dim:=Size(AlgebraProducts[1]);

        ErrorM1:=[];
        
        for j in [1..Size(AlgebraProducts)] do
			if AlgebraProducts[j] <> false then
				for k in [1..dim] do 
					
					u := NullMat(1,dim)[1];
					u[pairrepresentatives[j][1]] := 1;
					
					v := NullMat(1,dim)[1];
					v[pairrepresentatives[j][2]] := 1;
					
					w := NullMat(1,dim)[1];
					w[k] := 1;
					
					p := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,list);
					
					if p <> false then
						x := MAJORANA_InnerProduct(u,p,GramMatrix, list[3]);
						y := MAJORANA_InnerProduct(AlgebraProducts[j],w,GramMatrix, list[3]);
						
						if x <> false and y <> false and x <> y then 
							Add(ErrorM1,[j,k]);
							Error("Axiom M1");
						fi;
						
					fi;
				od;
			fi;
		od;
					
        return ErrorM1;

        end

        );

InstallGlobalFunction(MAJORANA_Fusion,

        function(T,GramMatrix,AlgebraProducts,EigenVectors,ProductList) # Checks if algebra obeys the fusion rules, outputs list of six lists which are empty if it does obey fusion rules
		
		# list should be of the form [coordinates,longcoordinates,pairorbitlist,pairconjelements,positionlist,NullSp]
		
        local errorfusion00, errorfusion02, errorfusion04, errorfusion22, errorfusion24, errorfusion44, a, t, j, k, l, x, y, z, x0;

        errorfusion00:=[];
        errorfusion02:=[];
        errorfusion04:=[];
        errorfusion22:=[];
        errorfusion24:=[];
        errorfusion44:=[];

        t:=Size(T);

        for j in [1..t] do

            a:=NullMat(1,Size(AlgebraProducts[1]))[1]; a[j]:=1;

            # 0,0 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][1])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][1][l],AlgebraProducts,ProductList);
                    if x <> false then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,ProductList);
                        if y <> false then
                            if ForAny(y, z-> z <> 0) then
                                Add(errorfusion00,[j,k,l]);
                                Error("0,0 fusion");
                            fi;
                        fi;
                    fi;
                od;
            od;

            # 0,1/4 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][2][l],AlgebraProducts,ProductList);
                    if x<> false then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,ProductList);
                        if y <> false and y <> x/4 then
                            Add(errorfusion02,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 0,1/32 fusion

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][k],EigenVectors[j][3][l],AlgebraProducts,ProductList);
                    if x <> false then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,ProductList);
                        if y <> false and y <> x/32 then
                            Add(errorfusion04,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 1/4,1/4 fusion

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][2][l],AlgebraProducts,ProductList);
                    if x<> false then
                        if MAJORANA_InnerProduct(a,x,GramMatrix, ProductList[3]) <> false then

                            y:= x - MAJORANA_InnerProduct(a,x,GramMatrix, ProductList[3])*a;
                            z:=MAJORANA_AlgebraProduct(a,y,AlgebraProducts,ProductList);

                            if z <> false and ForAny(z, x -> x <> 0)  then
                                Add(errorfusion22,[j,k,l]);
                            fi;
                        fi;
                    fi;
                od;
            od;

            # 1/4,1/32 fusion

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][k],EigenVectors[j][3][l],AlgebraProducts,ProductList);
                    if x <> false then
                        y:=MAJORANA_AlgebraProduct(a,x,AlgebraProducts,ProductList);
                        if y <> false and y <> x/32 then
                            Add(errorfusion24,[j,k,l]);
                        fi;
                    fi;
                od;
            od;

            # 1/32,1/32 fusion ## I think this is mathematically wrong, need to look at again

#            for k in [1..Size(EigenVectors[j][3])] do
#                for l in [1..Size(EigenVectors[j][3])] do
#
#                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][3][k],EigenVectors[j][3][l],AlgebraProducts,ProductList);
#
#                    if x <> false then
#                        y:= x - MAJORANA_InnerProduct(a,x,GramMatrix, positionlist)*a;
#
#                        z:= y - 4*MAJORANA_AlgebraProduct(a,y,AlgebraProducts,ProductList);
#                        x0:=MAJORANA_AlgebraProduct(a,z,AlgebraProducts,ProductList);
#
#                        if x0 <> false and ForAny(x0, x -> x <> 0)  then
#                            Add(errorfusion44,[j,k,l]);
#                        fi;
#                    fi;
#                od;
#            od;
        od;

        return [errorfusion00,errorfusion02,errorfusion04,errorfusion22,errorfusion24,errorfusion44];

        end

        );

InstallGlobalFunction(MAJORANA_TestOrthogonality,

        function(T,GramMatrix,AlgebraProducts,EigenVectors,pairorbitlist) # Tests that eigenspaces are orthogonal with respect to the inner product

        local a, b, t, errorortho01, errorortho02, errorortho04, errorortho12, errorortho14, errorortho24, j, k, l, x;

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
                x := MAJORANA_InnerProduct(EigenVectors[j][1][k],a,GramMatrix, pairorbitlist);
                if x <> false and x <> 0 then
                    Append(errorortho01,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][2])] do
                    x := MAJORANA_InnerProduct(EigenVectors[j][1][k],EigenVectors[j][2][l],GramMatrix, pairorbitlist);
                    if x <> false and x <> 0 then
                        Append(errorortho02,[[j,k,l]]);
                    fi;
                od;
            od;

            for k in [1..Size(EigenVectors[j][1])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x := MAJORANA_InnerProduct(EigenVectors[j][1][k],EigenVectors[j][3][l],GramMatrix, pairorbitlist);
                    if x <> false and x <> 0 then
                        Append(errorortho02,[[j,k,l]]);
                    fi;
                od;
            od;

            for k in [1..Size(EigenVectors[j][2])] do
                x := MAJORANA_InnerProduct(EigenVectors[j][2][k],a,GramMatrix, pairorbitlist);
                if x <> false and x <> 0 then
                    Append(errorortho12,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][3])] do
                x := MAJORANA_InnerProduct(EigenVectors[j][3][k],a,GramMatrix, pairorbitlist);
                if x <> false and x <> 0 then
                    Append(errorortho14,[[j,k]]);
                fi;
            od;

            for k in [1..Size(EigenVectors[j][2])] do
                for l in [1..Size(EigenVectors[j][3])] do
                    x := MAJORANA_InnerProduct(EigenVectors[j][2][k],EigenVectors[j][3][l],GramMatrix, pairorbitlist);
                    if x <> false and x <> 0 then
                        Append(errorortho24,[[j,k,l]]);
                    fi;
                od;
            od;
        od;

        return [errorortho01,errorortho02,errorortho04,errorortho12,errorortho14,errorortho24];

        end

        );

InstallGlobalFunction(MAJORANA_AxiomM2,

        function(GramMatrix,AlgebraProducts,ProductList,pairorbitlist) # Tests that the algebra obeys axiom M2

        local B, dim, L, i, j , k , l, m, Diagonals, a,b,c,d,x0,x1,x2,x3;

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
                              MAJORANA_InnerProduct(x0,x1,GramMatrix, pairorbitlist)
                            - MAJORANA_InnerProduct(x2,x3,GramMatrix, pairorbitlist);
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

InstallGlobalFunction(MAJORANA_FillGramMatrix,

function(GramMatrix, Orbitals, longcoordinates, pairorbitlist, dim)

    local i, j, x, y, GramMatrixFull;

    GramMatrixFull := NullMat(dim,dim);
    
    for i in [1..dim] do 
		for j in [1..dim] do
			GramMatrixFull[i][j] := GramMatrix[pairorbitlist[i][j]];
		od;
	od;

    return GramMatrixFull;

    end

    );

InstallGlobalFunction(MAJORANA_Orthogonality,

function(a,b,n,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim)

    local mat, vec, record, sum, i, j, k, l, m, p;

    mat := [];
    vec := [];
    record := [];

    if a = 0 then

        for i in [1..Size(EigenVectors[n][b])] do

            sum:=[];

            Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

            for j in [1..dim] do
				if EigenVectors[n][b][i][j] <> 0 then 

					m := pairorbitlist[j][n];

					if GramMatrix[m] <> false then
						Add(sum,-EigenVectors[n][b][i][j]*GramMatrix[m]);
					else
						p := Position(UnknownInnerProducts,m);
						mat[Size(mat)][p]:=mat[Size(mat)][p] + EigenVectors[n][b][i][j];
					fi;
				fi;
            od;

            Add(vec, [Sum(sum)]);

            if ForAll(mat[Size(mat)], x -> x = 0) then
                if vec[Size(vec)] <> [0] then

                    Error("Orthogonality error");
                else
					Remove(vec);
					Remove(mat);
                fi;
            fi;

        od;

    else

        for i in [1..Size(EigenVectors[n][a])] do
            for j in [1..Size(EigenVectors[n][b])] do

                sum:=[];

                Add(mat,NullMat(1,Size(UnknownInnerProducts))[1]);

                for k in [1..dim] do
					if EigenVectors[n][a][i][k] <> 0 then
						for l in [1..dim] do
							if EigenVectors[n][b][j][l] <> 0 then 
							
								m := pairorbitlist[k][l];

								if GramMatrix[m] <> false then
									Add(sum,-EigenVectors[n][a][i][k]*EigenVectors[n][b][j][l]*GramMatrix[m]);
								else
									p := Position(UnknownInnerProducts,m);
									mat[Size(mat)][p]:=mat[Size(mat)][p] + EigenVectors[n][a][i][k]*EigenVectors[n][b][j][l];
								fi;
							fi;
						od;
					fi;
				od;


                Add(vec,[Sum(sum)]);

                if ForAll(mat[Size(mat)], x -> x = 0) then
                    if vec[Size(vec)] <> [0] then

                        Error("Orthogonality error");
					else
						Remove(vec);
						Remove(mat);
                    fi;
                fi;

            od;
        od;

    fi;

    return [mat,vec];

    end

    );
    
InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(j, ev, EigenVectors, UnknownAlgebraProducts, AlgebraProducts, pairrepresentatives, ProductList)

	local mat, vec, row, sum, table, l, m, x, a, b, dim;
    
    mat := []; vec :=[];
	
	table := [0, 1/4, 1/32];
	
	dim := Size(AlgebraProducts[1]);
	
	a := [1..dim]*0; a[j] := 1;
	
	for l in [1..Size(EigenVectors[j][ev])] do
					
		row := [1..Size(UnknownAlgebraProducts)]*0;
		sum := [];
		
		for m in [1..dim] do 
			if EigenVectors[j][ev][l][m] <> 0 then
				
				x := ProductList[3][j][m];
                				
                if AlgebraProducts[x] = false then 
                    if pairrepresentatives[x] = [j,m] then 
                        row[Position(UnknownAlgebraProducts,x)] := EigenVectors[j][ev][l][m];
                    else
                        row := [];
                        break;
                    fi;
                else 
                    b := [1..dim]*0; b[m] := 1;
                    sum := sum + EigenVectors[j][ev][l][m]*MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                fi;
			fi;
		od;
        
        if row <> [] then 
		
            x := -sum + EigenVectors[j][ev][l]*table[ev];

            if sum <> [] then
                if ForAll(row, x -> x = 0) then 
                    if ForAny( x , y -> y <> 0) then 
                        Error("Step 7 system 1");
                    fi;
                else
                    Add(mat,row);
                    Add(vec,x);
                fi;
            fi;
        fi;
		
	od;
    
    return [mat,vec];
	
	end);
    
InstallGlobalFunction(MAJORANA_RemoveNullSpace,

function(v,NullSp) 

    local i, dim;
    
    dim := Size(v);

    if Size(NullSp) > 0 then 
        for i in [1..Size(NullSp)] do
            if v[dim - i + 1] <> 0 then 
                v := v - v[dim - i + 1]*NullSp[i];
            fi;
        od;
    fi;
    
    return v;
    
    end
    
    );


InstallGlobalFunction(MAJORANA_UnknownsAxiomM1,
    
    function(k,l, GramMatrix, AlgebraProducts, ProductList, pairrepresentatives)
    
    local j,dim,res,x,y,z,u,v,w;

    j := 1;
    
    dim := Size(AlgebraProducts[1]);
    
    res := false;
    
    while j < Size(AlgebraProducts) + 1 do 
    
        if AlgebraProducts[j] <> false and AlgebraProducts[j][k] <> 0 then 
        
            u := [1..dim]*0; u[l] := 1;
            v := [1..dim]*0; v[pairrepresentatives[j][1]] := 1;
        
            x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
            if x <> false then
                
                u := [1..dim]*0; u[pairrepresentatives[j][2]] := 1;
            
                y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList[3]);
                
                if y <> false then 
                    
                    v := StructuralCopy(AlgebraProducts[j]); v[k] := 0;
                    
                    w := [1..dim]*0; w[l] := 1;
                    
                    z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList[3]);
                    
                    if z <> false then 
                                           
                        res := (y - z)/AlgebraProducts[j][k]; 
                                          
                    fi;
                fi;
                
            fi;
            
            if res <> false then 
                
                j := Size(AlgebraProducts) +  1;
            
            else
            
                u := [1..dim]*0; u[l] := 1;
                v := [1..dim]*0; v[pairrepresentatives[j][2]] := 1;
            
                x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
                if x <> false then
                    
                    u := [1..dim]*0; u[pairrepresentatives[j][1]] := 1;
                
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList[3]);
                    
                    if y <> false then 
                        
                        v := StructuralCopy(AlgebraProducts[j]); v[k] := 0;
                        
                        w := [1..dim]*0; w[l] := 1;
                        
                        z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList[3]);
                        
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
    
InstallGlobalFunction(MAJORANA_UnknownOrbital,

function(i,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives)

    local j, k, l, x, res;
    
    res := false;
    
    j := 1;
    
    while j < Size(Orbitals[i]) + 1 do 
        
        k := ProductList[5][Position(ProductList[2],Orbitals[i][j][1])];
        l := ProductList[5][Position(ProductList[2],Orbitals[i][j][2])];
        
        if k > 0 and l > 0 then 

            res := MAJORANA_UnknownsAxiomM1(k,l,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
        
        elif k > 0 and l < 0 then 
        
            x := MAJORANA_UnknownsAxiomM1(k,-l,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
        
            if x <> false then 
                res := - x;
            fi;
            
        elif k < 0 and l > 0 then 
        
            x := MAJORANA_UnknownsAxiomM1(-k,l,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
        
            if x <> false then 
                res := - x;
            fi; 
            
        else
        
            res := MAJORANA_UnknownsAxiomM1(-k,-l,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
        fi;
        
        if res = false then 
            j := j + 1;
        else
            j := Size(Orbitals[i]) + 1;
        fi;
        
    od;
    
    return res;
    
    end );
    
InstallGlobalFunction(MAJORANA_ReversedEchelonForm,
    function(mat)
    
    local i, j, k, l, m, n, list, temp;   
     
    n := Size(mat);
    m := Size(mat[1]);
    
    i := 1;
    j := 1;
    
    while i < n + 1 do
        while j < m do
        
            list := [1..n]*0;
        
            for k in [i..n] do 
                if mat[k][m - j + 1] <> 0 then
                    mat[k] := mat[k]/mat[k][m - j + 1];
                    list[k] := 1;
                fi;
            od;
            
            if ForAll(list, x -> x = 0) then 
                j := j + 1;
            else
            
                l := Position(list,1);
                list[l] := 0;
            
                # swap rows
            
                temp := ShallowCopy(mat[i]);
                mat[i] := ShallowCopy(mat[l]);
                mat[l] := ShallowCopy(temp);
                
                for k in Positions(list,1) do
                    mat[k] := mat[k] - mat[l];
                od;
                
                j := j + 1;
                
                break;
                
            fi;
        od;
        
        i := i + 1;
        
    od;
    
    for i in [1..n] do
        for j in [i..n-1] do 
            mat[i] := mat[i] - mat[j+1]*mat[i][m-j];
        od;
    od;
            
    end );
        
InstallGlobalFunction(MajoranaRepresentation,

function(G,T)

    local   # Seress
            coordinates, representatives, conjelements, orbitlist, pairrepresentatives, pairconjelements, pairorbitlist, longcoordinates, positionlist, ProductList,

            long3Aaxes, long4Aaxes, long5Aaxes,

            # error checking
            ErrorFusion, ErrorM1, ErrorM2, ErrorOrthogonality,

            # indexing and temporary variables
            i, j, k, l, m, n, x, y, z, b,

            # Step 0 - Set Up
            Output, t, Orbitals, SizeOrbitals, OrbitalsT, SizeOrbitalsT, orbits, SizeOrbits, OrbitsT, SizeOrbitsT,

            # Step 1 - Shape
            Shape, RepsSquares6A, Unknowns3X,

            # Step 2 - Possible shapes
            Binaries, master, 3Aaxes, 4Aaxes, 5Aaxes, 5AaxesFixed, u, v, w,

            # Step 3 - Products and evecs I
            GramMatrix, GramMatrixT, GramMatrixFull, LI, NullSpT, AlgebraProducts, EigenVectors,
            # KnownInnerProducts,
            # KnownAlgebraProducts,
            EigenVector, sign, x0, x1, xm1, x2, xm2, x3, x4, x5, x6, x7, x8, x2A, x3A,

            # Step 4 - More products and evecs
            h, s, xj, xk, xl, xik, xil, xjk, xjl, xkl, xx, L, Diagonals, NullSp, dim, a,

            # Step 5 - More evecs
            switch, Dimensions, NewDimensions, NewEigenVectors,

            # Step 6 - More inner products
            UnknownInnerProducts, mat, vec, sum, row, Solution,

            # Step 7 - More algebra products
            Alpha, Alpha2, Beta, Beta2, walpha, wbeta, c, Form, zeros,  UnknownAlgebraProducts, record,
            
            falsecount, newfalsecount, maindimensions, newdimensions, switchmain, count, fres; 

            # str1, str2, str3, str4, str5,
            
            
            
 

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

    x:=OrbitsDomain(G,Cartesian(T,T),OnPairs);
    
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
	
	 SizeOrbitalsT:=Size(OrbitalsT);

                                        ## STEP 1: SHAPE ##

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    Shape:=NullMat(1,SizeOrbitalsT)[1];

    RepsSquares6A:=[];

    for i in [1..SizeOrbitalsT] do
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

    for i in [1..SizeOrbitalsT] do
        if Shape[i][1] = '3' then
            j:=0;
            while j<Size(OrbitalsT[i]) do
                j:=j+1;
                if OrbitalsT[i][j][1]*OrbitalsT[i][j][2] in RepsSquares6A then
                    Shape[i]:="3A";;
                    j:=Size(OrbitalsT[i])+1;;
                fi;
            od;
        fi;
    od;

    Unknowns3X:=[];

    for i in [1..SizeOrbitalsT] do
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
                    Shape[k]:="3C";
                else
                    Shape[k]:="3A";
                fi;
            od;

            # Create lists of 3A, 4A and 5A axes

            3Aaxes:=[];
            4Aaxes:=[];
            5Aaxes:=[];

            for j in [1..SizeOrbitalsT] do
                if Shape[j]=['3','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
						x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(3Aaxes,Set([x,x^2]));
                    od;
                fi;
                if Shape[j]=['4','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
						x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(4Aaxes,Set([x,x^3]));
                    od;
                fi;
                if Shape[j]=['5','A'] then
                    for k in [1..Size(OrbitalsT[j])] do
						x := OrbitalsT[j][k][1]*OrbitalsT[j][k][2];
                        Add(5Aaxes,Set([x,x^2,x^3,x^4]));
                    od;
                fi;
            od;

            3Aaxes:=DuplicateFreeList(3Aaxes); u:=Size(3Aaxes);
            4Aaxes:=DuplicateFreeList(4Aaxes); v:=Size(4Aaxes);
            5Aaxes:=DuplicateFreeList(5Aaxes); w:=Size(5Aaxes);

			for j in [1..u] do
				3Aaxes[j] := 3Aaxes[j][1];
			od;
			
			for j in [1..v] do
				4Aaxes[j] := 4Aaxes[j][1];
			od;
			
			for j in [1..w] do
				5Aaxes[j] := 5Aaxes[j][1];
			od;

            dim:=t+u+v+w;

            coordinates:=[];

            Append(coordinates,T);
            Append(coordinates,3Aaxes);
            Append(coordinates,4Aaxes);
            Append(coordinates,5Aaxes);

            long3Aaxes:=[];
            long4Aaxes:=[];
            long5Aaxes:=[];

            positionlist:=[];

            for j in [1..t] do
                Append(positionlist,[j]);
            od;

            for j in [t+1..t+u] do
                Append(long3Aaxes, [coordinates[j], coordinates[j]^2]);
                Append(positionlist,[j,j]);
            od;

            for j in [t+u+1..t+u+v] do
                Append(long4Aaxes, [coordinates[j], coordinates[j]^3]);
                Append(positionlist,[j,j]);
            od;

            for j in [t+u+v+1..dim] do
                Append(long5Aaxes, [coordinates[j], coordinates[j]^2, coordinates[j]^3, coordinates[j]^4]);
                Append(positionlist,[j,-j,-j,j]);
            od;

            longcoordinates:=StructuralCopy(T);

            Append(longcoordinates,long3Aaxes);
            Append(longcoordinates,long4Aaxes);
            Append(longcoordinates,long5Aaxes);

            orbits:=Orbits(G,coordinates);

            SizeOrbits:=Size(orbits);

            OrbitsT:=Filtered(orbits,x->Order(Representative(x)) = 2);

            SizeOrbitsT := Size(OrbitsT);

            representatives:=[];

            for j in [1..SizeOrbitsT] do
                x := Representative(OrbitsT[j]);
                Add(representatives,Position(coordinates,x));
            od;

            for j in [1..SizeOrbits] do
                x := Representative(orbits[j]);
                if Order(x) <> 2 then
                    Add(representatives,Position(coordinates,x));
                fi;
            od;


            conjelements:=[];
            orbitlist:=[];

            for j in [1..dim] do
                k:=1;
                while k < SizeOrbits + 1 do
                    if coordinates[j] in orbits[k] then
                        Add(orbitlist,k);
                        Add(conjelements,RepresentativeAction(G,representatives[k],coordinates[j]));
                        k := SizeOrbits + 1;
                    else
                        k := k + 1;
                    fi;
                od;
            od;

            x:=Orbits(G,Cartesian(coordinates,coordinates),OnPairs);
            
            Orbitals := [];
    
			for j in [1..Size(x)] do
				Add(Orbitals, ShallowCopy(x[j]));
			od;

            # This is a bit of a patch, ask Markus tomorrow

            j:=1;

            while j < Size(Orbitals) + 1 do
                if Order(Orbitals[j][1][1]) = 2 and Order(Orbitals[j][1][2]) = 2 then
                    Remove(Orbitals,j);
                else
                    j := j+1;
                fi;
            od;

            Orbitals := Concatenation(OrbitalsT,Orbitals);
            
            j := SizeOrbitalsT + 1;
            
            while j < Size(Orbitals) + 1 do 
    
				if not [Orbitals[j][1][2],Orbitals[j][1][1]] in Orbitals[j] then
				
					k := j + 1;
					
					while k < Size(Orbitals) +1 do
					
						if  [Orbitals[j][1][2],Orbitals[j][1][1]]  in Orbitals[k] then
                        
                            if Order(Orbitals[j][1][1]) < Order(Orbitals[j][1][2]) then 
						
                                Append(Orbitals[j],Orbitals[k]);
                                Remove(Orbitals,k);
                            
                            else 
                            
                                Append(Orbitals[k],Orbitals[j]);
                                Remove(Orbitals,j);
                                
                                j := j - 1;
                                
                            fi;
							
							k := Size(Orbitals) + 1;
							
						else
							
							k := k + 1;
						fi;
					od;
										
				fi;
				
				j := j + 1;
				
			od;

            SizeOrbitals:=Size(Orbitals);

            pairrepresentatives:=[];

            for j in [1..SizeOrbitals] do
                x := Orbitals[j][1];
                Add(pairrepresentatives, [Position(coordinates,x[1]), Position(coordinates,x[2])]);
            od;

            pairconjelements:=NullMat(dim,dim);
            pairorbitlist := NullMat(dim,dim);

            for j in [1..dim] do
                for k in [1..dim] do
                    l:=1;
                    while l < SizeOrbitals + 1 do
                        if [coordinates[j],coordinates[k]] in Orbitals[l] then
                            pairorbitlist[j][k] := l;
                            
                            x := RepresentativeAction(G,[coordinates[pairrepresentatives[l][1]],coordinates[pairrepresentatives[l][2]]],[coordinates[j],coordinates[k]],OnPairs);
                            
                            if x <> fail then         
								pairconjelements[j][k] := x;
							else
								x := RepresentativeAction(G,[coordinates[pairrepresentatives[l][1]],coordinates[pairrepresentatives[l][2]]],[coordinates[k],coordinates[j]],OnPairs);
								pairconjelements[j][k] := x;
							fi;
								
                            l := SizeOrbitals + 1;
                        else
                            l := l+1;
                        fi;
                    od;
                od;
            od;

            ProductList:=[coordinates,longcoordinates,pairorbitlist,pairconjelements,positionlist,[]];


                                        ## STEP 3: PRODUCTS AND EVECS I ##

            # Set up Gram matrix

            GramMatrix :=

            NullMat(dim,dim);
            for j in [1..dim] do
                for k in [1..dim] do
                    GramMatrix[j][k] := false;
                od;
            od;

            # Set up algebra product and gram matrices

            AlgebraProducts := NullMat(1,SizeOrbitals)[1];
            GramMatrix := NullMat(1,SizeOrbitals)[1];

            for j in [1..SizeOrbitals] do
                AlgebraProducts[j]:=false;
                GramMatrix[j]:=false;
            od;

            # Set up eigenvector matrix

            EigenVectors:=NullMat(t,3);

            for j in [1..t] do
                for k in [1..3] do
                    EigenVectors[j][k]:=[];
                od;
            od;

            # Start filling in values and products!

            l:=1;

            # (2,2) products and eigenvectors from IPSS10


            # Add eigenvectors from IPSS10

            for j in [1..t] do

                x0:=j;

                for l in [1..t] do

                    k:=pairorbitlist[x0][l];

                    x1 := Position(pairorbitlist[x0],k);

                    if Shape[k] = ['2','A'] then

                        x2 := Position(T,T[x0]*T[x1]);

                        EigenVector:=[1..dim]*0;

                        EigenVector[x0] :=-1/4;
                        EigenVector[x1] :=1;
                        EigenVector[x2] :=1;

                        Add(EigenVectors[j][1],EigenVector);

                        EigenVector:=[1..dim]*0;

                        EigenVector[x1]:=1;
                        EigenVector[x2]:=-1;

                        Add(EigenVectors[j][2],EigenVector);

                    elif Shape[k] = ['2','B'] then

                        EigenVector := [1..dim]*0;

                        EigenVector[x1] := 1;

                        Add(EigenVectors[j][1],EigenVector);

                    elif Shape[k] = ['3','A'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);
                        x3 := positionlist[t + Position(long3Aaxes,T[x0]*T[x1])];

                        EigenVector := [1..dim]*0;;

                        EigenVector[x0]:=-10/27;
                        EigenVector[x1]:=32/27;
                        EigenVector[xm1]:=32/27;
                        EigenVector[x3]:=1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-8/45;
                        EigenVector[x1]:=-32/45;
                        EigenVector[xm1]:=-32/45;
                        EigenVector[x3]:=1;

                        Add(EigenVectors[j][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    elif Shape[k] = ['3','C'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/32;
                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    elif Shape[k] = ['4','A'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);
                        x2 := Position(T, T[x1]*T[x0]*T[x1]);
                        x4 := positionlist[t + 2*u + Position(long4Aaxes,T[x0]*T[x1])];

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/2;
                        EigenVector[x1]:=2;
                        EigenVector[xm1]:=2;
                        EigenVector[x2]:=1;
                        EigenVector[x4]:=1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/3;
                        EigenVector[x1]:=-2/3;
                        EigenVector[xm1]:=-2/3;
                        EigenVector[x2]:=-1/3;
                        EigenVector[x4]:=1;

                        Add(EigenVectors[j][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    elif Shape[k] = ['4','B'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);
                        x2 := Position(T, T[x1]*T[x0]*T[x1]);
                        x4 := Position(T, (T[x0]*T[x1])^2);

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-1/32;
                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=1;
                        EigenVector[x2]:=1/8;
                        EigenVector[x4]:=-1/8;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    elif Shape[k] = ['5','A'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);
                        x2 := Position(T, T[x1]*T[x0]*T[x1]);
                        xm2 := Position(T, T[x0]*T[x1]*T[x0]*T[x1]*T[x0]);
                        x5 := positionlist[t + 2*u + 2*v + Position(long5Aaxes,T[x0]*T[x1])];

                        if x5 < 0 then
                            x5 := -x5;
                            sign := -1;
                        else
                            sign := 1;
                        fi;

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=3/512;
                        EigenVector[x1]:=-15/128;
                        EigenVector[xm1]:=-15/128;
                        EigenVector[x2]:=-1/128;
                        EigenVector[xm2]:=-1/128;
                        EigenVector[x5]:=sign*1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-3/512;
                        EigenVector[x1]:=1/128;
                        EigenVector[xm1]:=1/128;
                        EigenVector[x2]:=15/128;
                        EigenVector[xm2]:=15/128;
                        EigenVector[x5]:=sign*1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1/128;
                        EigenVector[xm1]:=1/128;
                        EigenVector[x2]:=-1/128;
                        EigenVector[xm2]:=-1/128;
                        EigenVector[x5]:=sign*1;

                        Add(EigenVectors[j][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x2]:=1;
                        EigenVector[xm2]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    elif Shape[k] = ['6','A'] then

                        xm1 := Position(T, T[x0]*T[x1]*T[x0]);
                        x2 := Position(T, T[x1]*T[x0]*T[x1]);
                        xm2 := Position(T, T[x0]*T[x1]*T[x0]*T[x1]*T[x0]);
                        x3 := Position(T, T[x1]*T[x0]*T[x1]*T[x0]*T[x1]);
                        x2A := Position(T, (T[x0]*T[x1])^3);
                        x3A := positionlist[t + Position(long3Aaxes,(T[x0]*T[x1])^2)];

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=2/45;
                        EigenVector[x1]:=-256/45;
                        EigenVector[xm1]:=-256/45;
                        EigenVector[x2]:=-32/45;
                        EigenVector[xm2]:=-32/45;
                        EigenVector[x3]:=-32/45;
                        EigenVector[x2A]:=32/45;
                        EigenVector[x3A]:=1;

                        Add(EigenVectors[j][1],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x0]:=-8/45;
                        EigenVector[x2]:=-32/45;
                        EigenVector[xm2]:=-32/45;
                        EigenVector[x3]:=-32/45;
                        EigenVector[x2A]:=32/45;
                        EigenVector[x3A]:=1;

                        Add(EigenVectors[j][2],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x1]:=1;
                        EigenVector[xm1]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                        EigenVector:=NullMat(1,dim)[1];

                        EigenVector[x2]:=1;
                        EigenVector[xm2]:=-1;

                        Add(EigenVectors[j][3],StructuralCopy(EigenVector));

                    fi;
                od;
                
                # 1/32 eigenvectors from conjugation
                
                for k in [1..u] do 
                
                    h := 3Aaxes[k];
                    
                    if not h^T[j] in [h,h^2] then 
                    
                        l := positionlist[Position(longcoordinates,h^T[j])];
                        
                        EigenVector := [1..dim]*0;
                        
                        EigenVector[t+k] := 1;
                        EigenVector[l] := -1;
                        
                        Add(EigenVectors[j][3],EigenVector);
                    fi;
                od;
                
                for k in [1..v] do 
                
                    h := 4Aaxes[k];
                    
                    if not h^T[j] in [h,h^3] then 
                    
                        l := positionlist[Position(longcoordinates,h^T[j])];
                        
                        EigenVector := [1..dim]*0;
                        
                        EigenVector[t+u+k] := 1;
                        EigenVector[l] := -1;
                        
                        Add(EigenVectors[j][3],EigenVector);
                    fi;
                od;
                
                for k in [1..w] do 
                
                    h := 5Aaxes[k];
                    
                    if h^T[j] in [h^2,h^3] then 
                    
                        Display("Warning - a 5A axis is conj to its square or cube");
                        
                        EigenVector := [1..dim]*0;
                        
                        EigenVector[t+u+v+k] := 2;
                        
                        Add(EigenVectors[j][3],EigenVector);
                    
                    elif not h^T[j] in [h,h^4] then 
                    
                        l := positionlist[Position(longcoordinates,h^T[j])];
                        
                        EigenVector := [1..dim]*0;
                        
                        EigenVector[t+u+v+k] := 1;
                        
                        if l > 0 then 
                            EigenVector[l] := -1;
                        else
                            EigenVector[-l] := 1;
                        fi;
                        
                        Add(EigenVectors[j][3],EigenVector);
                    fi;
                od;
                
                
            od;

            # Products from IPSS10

            for j in [1..SizeOrbitals] do

                x := pairrepresentatives[j][1];
                y := pairrepresentatives[j][2];

                if Order(coordinates[x]) = 2 and Order(coordinates[y]) = 2 then

                    if Shape[j] = ['1','A'] then

                        AlgebraProducts[j] := NullMat(1,dim)[1];

                        AlgebraProducts[j][x] := 1;

                        GramMatrix[j] := 1;

                    elif Shape[j] = ['2','A'] then

                        x2:=Position(T,T[x]*T[y]);

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:= 1/8;
                        AlgebraProducts[j][y]:= 1/8;
                        AlgebraProducts[j][x2]:=-1/8;

                        GramMatrix[j] := 1/8;

                    elif Shape[j] = ['2','B'] then

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        GramMatrix[j]:=0;

                    elif Shape[j] = ['3','A'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);
                        x3 := positionlist[t + Position(long3Aaxes,T[x]*T[y])];

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/16;
                        AlgebraProducts[j][y]:=1/16;
                        AlgebraProducts[j][xm1]:=1/32;
                        AlgebraProducts[j][x3]:= - 135/2048;

                        GramMatrix[j]:=13/256;

                    elif Shape[j] = ['3','C'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/64;
                        AlgebraProducts[j][y]:=1/64;
                        AlgebraProducts[j][xm1]:=-1/64;

                        GramMatrix[j]:=1/64;

                    elif Shape[j] = ['4','A'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);
                        x2 := Position(T,T[y]*T[x]*T[y]);
                        x4 := positionlist[t + 2*u + Position(long4Aaxes,T[x]*T[y])];

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=3/64;
                        AlgebraProducts[j][y]:=3/64;
                        AlgebraProducts[j][xm1]:=1/64;
                        AlgebraProducts[j][x2]:=1/64;
                        AlgebraProducts[j][x4]:= -3/64;

                        GramMatrix[j]:=1/32;

                    elif Shape[j] = ['4','B'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);
                        x2 := Position(T,T[y]*T[x]*T[y]);
                        x2A := Position(T,(T[x]*T[y])^2);

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/64;
                        AlgebraProducts[j][y]:=1/64;
                        AlgebraProducts[j][xm1]:=-1/64;
                        AlgebraProducts[j][x2]:=-1/64;
                        AlgebraProducts[j][x4]:= 1/64;

                        GramMatrix[j]:=1/64;

                    elif  Shape[j] = ['5','A'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);
                        x2 := Position(T,T[y]*T[x]*T[y]);
                        xm2 := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                        x5 := positionlist[t + 2*u + 2*v + Position(long5Aaxes,T[x]*T[y])];

                        if x5 < 0 then
                            x5 := -x5;
                            sign := -1;
                        else
                            sign:=1;
                        fi;

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=3/128;
                        AlgebraProducts[j][y]:=3/128;
                        AlgebraProducts[j][xm1]:=-1/128;
                        AlgebraProducts[j][x2]:=-1/128;
                        AlgebraProducts[j][xm2]:=-1/128;
                        AlgebraProducts[j][x5] := sign*1;

                        GramMatrix[j]:=3/128;

                    elif Shape[j] = ['6','A'] then

                        xm1 := Position(T,T[x]*T[y]*T[x]);
                        x2 := Position(T,T[y]*T[x]*T[y]);
                        xm2 := Position(T,T[x]*T[y]*T[x]*T[y]*T[x]);
                        x3 := Position(T,T[y]*T[x]*T[y]*T[x]*T[y]);
                        x2A := Position(T,(T[x]*T[y])^3);
                        x3A := positionlist[t + Position(long3Aaxes,(T[x]*T[y])^2)];

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/64;
                        AlgebraProducts[j][y]:=1/64;
                        AlgebraProducts[j][xm1]:=-1/64;
                        AlgebraProducts[j][x2]:=-1/64;
                        AlgebraProducts[j][xm2]:=-1/64;
                        AlgebraProducts[j][x3] := -1/64;
                        AlgebraProducts[j][x2A] := 1/64;
                        AlgebraProducts[j][x3A] := 45/2048;

                        GramMatrix[j]:=5/256;
                    fi;

                # 2,3 products

                elif Order(coordinates[x]) = 2 and Order(coordinates[y]) = 3 then
                    if Order(coordinates[x]*coordinates[y]) in T then

                        s := coordinates[x]; h := coordinates[y];

                        # Inside a 3A algebra

                        x1 := Position(T,s*h);
                        xm1 := Position(T,s*h*h);

                        AlgebraProducts[j] := NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=2/9;
                        AlgebraProducts[j][x1]:=-1/9;
                        AlgebraProducts[j][xm1]:=-1/9;
                        AlgebraProducts[j][y]:=5/32;

                        GramMatrix[x][y]:=1/4;

                    elif Order(coordinates[x]*coordinates[y]) = 3 then

                        s := coordinates[x]; h := coordinates[y];

                        # Case (2A,3A) in IPSS10

                        xj := positionlist[t + Position(long3Aaxes,s*h*s)];
                        xk := positionlist[t + Position(long3Aaxes,h*s*h*h*s*h*h)];
                        xl := positionlist[t + Position(long3Aaxes,h*h*s*h*s*h)];

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/9;
                        AlgebraProducts[j][y]:=5/64;
                        AlgebraProducts[j][xj]:=3/64;
                        AlgebraProducts[j][xk]:=-1/16;
                        AlgebraProducts[j][xl]:=-1/16;

                        GramMatrix[j]:=1/9;

                    elif Order(coordinates[x]*coordinates[y]) =4 and (coordinates[x]*coordinates[y])^2 in T then

                        s := coordinates[x]; h := coordinates[y];

                        # Case (2A,3A) in IPSS10

                        xik:=Position(T,h*h*s*h);
                        xil:=Position(T,h*s*h*h);
                        xjk:=Position(T,s*h*h*s*h*s);
                        xjl:=Position(T,s*h*s*h*h*s);
                        xkl:=Position(T,h*s*h*s*h*h*s*h*h);
                        xx:=Position(T,h*h*s*h*s*h*h);
                        xj := positionlist[t + Position(long3Aaxes,s*h*s)];
                        xk := positionlist[t + Position(long3Aaxes,h*s*h*h*s*h*h)];
                        xl := positionlist[t + Position(long3Aaxes,h*h*s*h*s*h)];

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=1/45;
                        AlgebraProducts[j][y]:=1/64;
                        AlgebraProducts[j][xik]:=-1/90;
                        AlgebraProducts[j][xil]:=-1/90;
                        AlgebraProducts[j][xjk]:=-1/90;
                        AlgebraProducts[j][xjl]:=-1/90;
                        AlgebraProducts[j][xkl]:=1/45;
                        AlgebraProducts[j][xx]:=-1/45;
                        AlgebraProducts[j][xj]:=-1/64;
                        AlgebraProducts[j][xk]:=1/64;
                        AlgebraProducts[j][xl]:=1/64;

                        GramMatrix[j]:=1/36;

                    else
 
                        GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);

                    fi;

                # 2,4 products

                elif Order(coordinates[x]) = 2 and Order(coordinates[y]) = 4 then

                    s := coordinates[x];
                    h := coordinates[y];

                    if s*h in T then

                        # Inside a 4A algebra

                        x1 := Position(T,s*h);
                        xm1 := Position(T,s*h*h*h);
                        x2 := Position(T,s*h*h);

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x]:=5/16;
                        AlgebraProducts[j][x1]:=-1/8;
                        AlgebraProducts[j][xm1]:=-1/8;
                        AlgebraProducts[j][x2]:=-1/16;
                        AlgebraProducts[j][y]:= 3/16;

                        GramMatrix[j]:=3/8;

                    else
                        GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                    fi;

                # (2,5) values

                 elif Order(coordinates[x]) = 2 and Order(coordinates[y]) = 5 then

                    s := coordinates[x];
                    h := coordinates[y];

                    if s*h in T then

                        # Inside a 5A algebra

                        x1 := Position(T,s*h);
                        xm1 := Position(T,s*h^4);
                        x2 := Position(T,s*h^2);
                        xm2 := Position(T,s*h^3);

                        AlgebraProducts[j]:=NullMat(1,dim)[1];

                        AlgebraProducts[j][x1] := 7/4096;
                        AlgebraProducts[j][xm1] := 7/4096;
                        AlgebraProducts[j][x2] := -7/4096;
                        AlgebraProducts[j][xm2] := -7/4096;
                        AlgebraProducts[j][y] := 7/32;

                        GramMatrix[j]:=0;
                    else

                        GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);

                    fi;

                # 3,3 values

                elif Order(coordinates[x]) = 3 and Order(coordinates[y]) = 3 then

                    h := coordinates[x];
                    k := coordinates[y];

                    if x = y then

                        AlgebraProducts[j] := NullMat(1,dim)[1];

                        AlgebraProducts[j][x] := 1;

                        GramMatrix[j] := 8/5;

                    else

                        GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                    fi;

                # (3,4) values

                elif Order(coordinates[x]) = 3 and Order(coordinates[y]) = 4 then
                
                    GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);

                # (3,5) values

                elif Order(coordinates[x]) = 3 and Order(coordinates[y]) = 5 then

                    GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                    
                # (4,4) values

                elif Order(coordinates[x]) = 4 and Order(coordinates[y]) = 4 then

                if x = y then

                    AlgebraProducts[j] := NullMat(1,dim)[1];

                    AlgebraProducts[j][x] := 1;

                    GramMatrix[j] := 2;

                else
                
                    GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                
                fi;

                    
                # (4,5) values

                elif Order(coordinates[x]) = 4 and Order(coordinates[y]) = 5 then
                
                GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);

                # (5,5) values

                elif Order(coordinates[x]) = 5 and Order(coordinates[y]) = 5 then
       
                    h := coordinates[x];
                    k := coordinates[y];

                    if x = y then

                        l:=1;

                        while l < t+1 do

                            if T[l]*h in T then

                                s:=T[l]; x1:=Position(T,s*h); x2:=Position(T,s*h*h); x3:=Position(T,s*h*h*h); x4:=Position(T,s*h*h*h*h);

                                AlgebraProducts[j]:=NullMat(1,dim)[1];

                                AlgebraProducts[j][l]:=175/524288;
                                AlgebraProducts[j][x1]:=175/524288;
                                AlgebraProducts[j][x2]:=175/524288;
                                AlgebraProducts[j][x3]:=175/524288;
                                AlgebraProducts[j][x4]:=175/524288;

                                l:=t+1;
                            else
                                l:=l+1;
                            fi;
                        od;

                        GramMatrix[j] := 875/2^(19);

                    else
                    
						GramMatrix[j] := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                    fi;

                fi;
            od;
            
            while 1 = 1 do
                
                count := 0;
            
                for j in [1..SizeOrbitals] do
                    if GramMatrix[j] = false then 
                        x := MAJORANA_UnknownOrbital(j,Orbitals,GramMatrix,AlgebraProducts,ProductList,pairrepresentatives);
                        
                        if x <> false then 
                            count := count + 1;
                            GramMatrix[j] := x;
                        fi;
                    fi;
                od;
                
                if count = 0 then 
                    break;
                fi;
            od;            

            LI:=1;

            GramMatrixT:=MAJORANA_FillGramMatrix(GramMatrix,OrbitalsT,longcoordinates,pairorbitlist,t);

            x:=MAJORANA_PositiveDefinite(GramMatrixT);

            if x = -1 then
                Output[i]:=[StructuralCopy(Shape),"Error","Inner product not positive definite on A", StructuralCopy(GramMatrixT)];
                break;
            elif x = 0 then
                NullSpT:=MAJORANA_NullSpace(GramMatrixT);
                LI:=0;
            fi;

                                        ## STEP 4: MORE PRODUCTS ##



            
            maindimensions := [];
            
            maindimensions:=[];

            for j in [1..t] do
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


                                            ## STEP 5: MORE EVECS ##
                                            
                # Find linearly independent subsets of eigenvectors

                # Use these eigenvectors and the fusion rules to find more

                switch:=0;

                Dimensions := StructuralCopy(maindimensions);

                if ForAll(Dimensions,x->x=dim) then
                    switch:=1;
                fi;

                NewEigenVectors:=NullMat(t,3);

                for j in [1..t] do
                    for k in [1..3] do
                        NewEigenVectors[j][k]:=[];
                    od;
                od;

                while switch = 0 do
                    for j in [1..t] do
                        # 1, x fusion is a waste of time because a_0 obviously just preserves the evectors!
                        Output[i] := [];

                        # 0,0 fusion
                        fres := MAJORANA_TestFusion(1,1,j,Shape,AlgebraProducts,EigenVectors, GramMatrix, ProductList, dim);
                        if fres[1] then
                            Append(NewEigenVectors[j][1], fres[2]);
                        else
                            Output[i] := fres[2];
                            break;
                        fi;

                        # 0,1/4 fusion
                        fres := MAJORANA_TestFusion(1,2,j,Shape,AlgebraProducts,EigenVectors, GramMatrix, ProductList, dim);
                        if fres[1] then
                            Append(NewEigenVectors[j][2], fres[2]);
                        else
                            Output[i] := fres[2];
                            break;
                        fi;

                        # 0,1/32 fusion
                        fres := MAJORANA_TestFusion(1,3,j,Shape,AlgebraProducts,EigenVectors, GramMatrix, ProductList, dim);
                        if fres[1] then
                            Append(NewEigenVectors[j][3], fres[2]);
                        else
                            Output[i] := fres[2];
                            Error("0, 1/32 fusion");
                            break;
                        fi;

                        # 1/4,1/32 fusion
                        fres := MAJORANA_TestFusion(2,3,j,Shape,AlgebraProducts,EigenVectors, GramMatrix, ProductList, dim);
                        if fres[1] then
                            Append(NewEigenVectors[j][3], fres[2]);
                        else
                            Output[i] := fres[2];
                            break;
                        fi;

                        # 1/4,1/4 Fusion
                        fres := MAJORANA_TestFusion(2,2,j,Shape,AlgebraProducts,EigenVectors, GramMatrix, ProductList, dim);
                        if fres[1] then
                            Append(NewEigenVectors[j][1], fres[2]);
                        else
                            Output[i] := fres[2];
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
                # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products

                switch := 0;

                while switch = 0 do 
                
                    UnknownInnerProducts:=[];

                    for j in [1..SizeOrbitals] do
                        if GramMatrix[j] = false then
                            Add(UnknownInnerProducts,j);
                        fi;
                    od;
                
                    if Size(UnknownInnerProducts) = 0 then
                        
                        break;
                        
                    else

                        mat:=[];
                        vec:=[];

                        for j in [1..t] do

                            # 1- eigenvectors and 0-eigenvectors

                            x := MAJORANA_Orthogonality(0,1,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);

                            # 1- eigenvectors and 1/4-eigenvectors

                            x := MAJORANA_Orthogonality(0,2,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);

                            # 1- eigenvectors and 1/32-eigenvectors

                            x := MAJORANA_Orthogonality(0,3,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);

                            # 0-eigenvectors and 1/4-eigenvectors

                            x := MAJORANA_Orthogonality(1,2,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);

                            # 0-eigenvectors and 1/32-eigenvectors

                            x := MAJORANA_Orthogonality(1,3,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);

                            # 1/4-eigenvectors and 1/32-eigenvectors

                            x := MAJORANA_Orthogonality(2,3,j,UnknownInnerProducts,EigenVectors,GramMatrix, pairorbitlist,representatives, dim);

                            Append(mat,x[1]);
                            Append(vec,x[2]);
                        od;

                        Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                        if Size(Solution) = 2 then                    
                            
                                for k in [1..Size(Solution[1])] do
                                    if not k in Solution[2] then

                                        x:=UnknownInnerProducts[k]; 

                                        GramMatrix[x]:=Solution[1][k][1];
                                    fi;
                                od;
                                
                                if Size(Solution[2]) = Size(Solution[1]) then
                                    break;
                                fi;
                                
                             #   if Size(Solution[2]) = Size(Solution[1]) then 
                             #       Output[i] := [ StructuralCopy(Shape)
                             #                , "Fail"
                             #                , "Missing inner product values"
                             #                , StructuralCopy(GramMatrix)];
                             #       Error("Missing inner products");
                             #   
                             #       break;
                             #   fi;
                                
                        else
                            Output[i] := [ Shape
                                         , "Error"
                                         , "Inconsistent system of unknown inner products"
                                         , mat
                                         , vec
                                         , EigenVectors
                                         , AlgebraProducts
                                         , GramMatrix
                                         , UnknownInnerProducts ];
                            Output[i]:=StructuralCopy(Output[i]);
                            
                            break;
                        fi;
                    fi;
                od;
                
                if Size(Output[i]) > 0 then
                    break;
                fi;

                # Check that GramMatrix matrix is pd
                
                if ForAll(GramMatrix, x -> x <> false) then 
                    GramMatrixFull := MAJORANA_FillGramMatrix(GramMatrix, Orbitals, longcoordinates, pairorbitlist, dim);

                    L:=MAJORANA_LDLTDecomposition(GramMatrixFull);

                    Diagonals:=[];

                    for j in [1..Size(GramMatrixFull)] do
                        Append(Diagonals,[L[2][j][j]]);
                    od;

                    if ForAny(Diagonals, x->x<0) then
                        Output[i] := [ StructuralCopy(Shape)
                                     , "Error"
                                     , "The inner product is not positive definite"
                                     , StructuralCopy(3Aaxes)
                                     , StructuralCopy(4Aaxes)
                                     , StructuralCopy(5Aaxes)
                                     , StructuralCopy(5AaxesFixed)
                                     , StructuralCopy(GramMatrix) ];
                        break;
                    elif ForAny(Diagonals, x->x=0) then
                    
                        NullSp:=MAJORANA_NullSpace(GramMatrixFull);
                        
                        ProductList[6] := NullSp;
                        
                        LI:=0;
                    else
                        LI:=1;
                    fi;
                else
                    LI := 1;
                fi;

                if LI=0 then

                    dim:=t+u+v+w;
                    
                    n := Size(NullSp);

                    # Change alg products to get rid of any axes not in the basis

                    for j in [1..n] do
                        for k in [1..SizeOrbitals] do
                            if AlgebraProducts[k] <> false then
                                AlgebraProducts[k]:=AlgebraProducts[k] - NullSp[j]*AlgebraProducts[k][dim-j+1];
                            fi;
                        od;
                    od;

                    # Change evecs to get rid of any axes not in the basis

                    for j in [1..t] do
                        for k in [1..3] do
                            for l in [1..Size(EigenVectors[j][k])] do
                                for m in [1..n] do
                                    EigenVectors[j][k][l]:=EigenVectors[j][k][l] - NullSp[m]*EigenVectors[j][k][l][dim - m + 1];
                                od;
                            od;
                        od;
                    od;
                else
                    dim:=t+u+v+w;
                fi;

                                            ## STEP 7: MORE PRODUCTS II ##

                # Check fusion and M1

                ErrorM1:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,Orbitals,ProductList,pairrepresentatives);

                if Size(ErrorM1)>0 then
                    Output[i] := [ StructuralCopy(Shape)
                                 , "Error"
                                 , "Algebra does not obey axiom M1 step 7"
                                 , StructuralCopy(GramMatrix)
                                 , StructuralCopy(AlgebraProducts)
                                 , StructuralCopy(ErrorM1)];
                    Error("axiom M1");
                    break;
                fi;

                ErrorFusion:=MAJORANA_Fusion(T, GramMatrix, AlgebraProducts, EigenVectors,ProductList);

                if ForAny(ErrorFusion, x->Size(x) > 0) then
                
                    Output[i] := [ StructuralCopy(Shape)
                                 , "Error"
                                 , "Algebra does not obey fusion rules step 7"
                                 , StructuralCopy(GramMatrix)
                                 , StructuralCopy(AlgebraProducts)
                                 , StructuralCopy(EigenVectors)
                                 , StructuralCopy(ErrorFusion)];
                    Error("fusion");
                    break;
                fi;


                # Use eigenvectors to find more products
                            
                for j in [1..t] do
                
                    a := [1..dim]*0; a[j]:=1;
                
                    UnknownAlgebraProducts := [];
                
                    for k in [1..SizeOrbitals] do
                    
                        if AlgebraProducts[k] = false and pairrepresentatives[k][1] = j then 
                            Add(UnknownAlgebraProducts,k);
                        fi;
                            
                    od;
                    
                    if UnknownAlgebraProducts <> [] then
                    
                        mat := [];
                        vec := [];
                        
                        x := MAJORANA_EigenvectorsAlgebraUnknowns(j, 1, EigenVectors, UnknownAlgebraProducts, AlgebraProducts, pairrepresentatives, ProductList);
                        
                        if x <> 0 then 
                            Append(mat,x[1]);
                            Append(vec,x[2]);
                        fi;
                        
                        x := MAJORANA_EigenvectorsAlgebraUnknowns(j, 2, EigenVectors, UnknownAlgebraProducts, AlgebraProducts, pairrepresentatives, ProductList);
                        
                        if x <> 0 then 
                            Append(mat,x[1]);
                            Append(vec,x[2]);
                        fi;
                        
                        x := MAJORANA_EigenvectorsAlgebraUnknowns(j, 3, EigenVectors, UnknownAlgebraProducts, AlgebraProducts, pairrepresentatives, ProductList);
                        
                        if x <> 0 then 
                            Append(mat,x[1]);
                            Append(vec,x[2]);
                        fi;
                        
                        # use fact that if v in null space then a \cdot v = 0
                        
                        if LI = 0 then 
                        
                            for k in [1..Size(NullSp)] do
                                
                                row := [1..Size(UnknownAlgebraProducts)]*0;
                                sum := [];
                                
                                for l in [1..dim] do
                                    if NullSp[k][l] <> 0 then
                                    
                                        y := pairorbitlist[j][l]; 
                                        
                                        if AlgebraProducts[y] <> false then 
                                        
                                            a := [1..dim]*0; a[j] := 1;
                                            b := [1..dim]*0; b[l] := 1;
                                            
                                            sum := sum - NullSp[k][l]*MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                                            
                                        else
                                            
                                            if pairrepresentatives[y] = [j,l] or pairrepresentatives[y] = [l,j] then
                                                row[Position(UnknownAlgebraProducts,y)]:=NullSp[k][l];
                                            else
                                                row := [];
                                                break; ##### need to move on to next null space vector
                                            fi;
                                        fi;
                                    fi;
                                od;
                                
                                if sum <> [] and row <> [] then
                                    if ForAll(row, x -> x = 0) then 
                                        if ForAny( sum , y -> y <> 0) then 
                                            Error("Step 7 nullspace"); 
                                        fi;
                                    else
                                        Add(mat,row);
                                        Add(vec,sum);
                                    fi;
                                fi;
                            
                            od;
                        fi;
                                     
                        if mat <> [] then 
                            Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                            if Size(Solution) = 2 then
                                    for k in [1..Size(Solution[1])] do
                                        if not k in Solution[2] then

                                            x:=UnknownAlgebraProducts[k]; 

                                            AlgebraProducts[x]:=Solution[1][k];
                                        fi;
                                    od;
                            else
                                Output[i] := [ StructuralCopy(Shape)
                                             , "Error"
                                             , "Inconsistent system of unknown algebra products step 7"
                                             , StructuralCopy(GramMatrix)
                                             , [] # StructuralCopy(KnownAlgebraProducts)
                                             , StructuralCopy(AlgebraProducts)
                                             , StructuralCopy(EigenVectors)
                                             , StructuralCopy(mat)
                                             , StructuralCopy(vec)
                                             , StructuralCopy(Solution)
                                             , StructuralCopy(UnknownAlgebraProducts)];
                                Error("Step 7");
                                break;
                            fi;
                        fi;
                    fi;                        
                od;

                if Size(Output[i])>0 then
                    break;
                fi;
                
                                                   ## STEP 8: RESURRECTION PRINCIPLE I ##

                # Check fusion and M1

                ErrorM1:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,Orbitals,ProductList,pairrepresentatives);

                if Size(ErrorM1)>0 then
                    Output[i] := [ StructuralCopy(Shape)
                                 , "Error"
                                 , "Algebra does not obey axiom M1 step 8"
                                 , StructuralCopy(GramMatrix)
                                 , [] # StructuralCopy(KnownAlgebraProducts)
                                 , StructuralCopy(AlgebraProducts)
                                 , StructuralCopy(ErrorM1)];
                fi;

                ErrorFusion:=MAJORANA_Fusion(T, GramMatrix, AlgebraProducts,EigenVectors,ProductList);

                if ForAny(ErrorFusion, x->Size(x) > 0) then
                    Output[i] := [ StructuralCopy(Shape)
                                 , "Error"
                                 , "Algebra does not obey fusion rules step 8"
                                 , StructuralCopy(GramMatrix)
                                 , []
                                 , StructuralCopy(AlgebraProducts)
                                 , StructuralCopy(EigenVectors)
                                 , StructuralCopy(ErrorFusion)];
                    break;
                fi;

                # Use resurrection principle (Step 7 Seress)

                UnknownAlgebraProducts := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts);
                
                if ForAny(UnknownAlgebraProducts, x -> pairrepresentatives[x][1] < t+1) then

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

                                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][1][l],(walpha - wbeta),AlgebraProducts,ProductList);

                                    if x <> false then
                                        for n in [1..dim] do
                                            if x[n] <> 0 then
                                            
                                                y := pairorbitlist[j][n];
                                            
                                                if AlgebraProducts[y] <> false then
                                                
                                                    a := [1..dim]*0; a[j] := 1;
                                                    b := [1..dim]*0; b[n] := 1;
                                                    
                                                    Add(sum,-MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)*x[n]);
                                                else
                                                    if pairrepresentatives[y] = [j,n] or pairrepresentatives[y] = [n,j] then
                                                        row[Position(UnknownAlgebraProducts,y)]:=x[n];
                                                    else
                                                        row := [];
                                                        break; ##### need to move on to next alpha2
                                                    fi;
                                                fi;
                                            fi;
                                        od;
                                    fi;

                                    # calculate rhs

                                    if row <> [] then 
                                        for n in [1..t] do
                                            if EigenVectors[j][1][l][n] <> 0 then
                                            
                                                y := pairorbitlist[n][c];
                                                
                                                if AlgebraProducts[y] <> false then
                                                    
                                                    a := [1..dim]*0; a[c] := 1;
                                                    b := [1..dim]*0; b[n] := 1;
                                                    
                                                    Add(sum,-EigenVectors[j][1][l][n]*MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)/4);
                                                else
                                                    if pairrepresentatives[y] = [c,n] or pairrepresentatives[y] = [n,c] then
                                                        row[Position(UnknownAlgebraProducts,y)]:=row[Position(UnknownAlgebraProducts,y)] + EigenVectors[j][1][l][n]/4;
                                                    else
                                                        row := [];
                                                        break; ##### need to move on to next alpha2
                                                    fi;
                                                fi;
                                            fi;
                                        od;
                                    fi;
                                    
                                    x := Sum(sum) - MAJORANA_AlgebraProduct(wbeta,EigenVectors[j][1][l],AlgebraProducts,ProductList)/4;
                                    
                                    if sum <> [] and row <> [] then
                                        if ForAll(row, x -> x = 0) then 
                                            if ForAny( x , y -> y <> 0) then 
                                                Error("Step 8 system 1");
                                            fi;
                                        else
                                            Add(mat,row);
                                            Add(vec,x);
                                        fi;
                                    fi;       
                                    
                                od;

                                ## Second part

                                for l in Beta2 do

                                    sum:=[];
                                    row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                    # calculate lhs

                                    x:=MAJORANA_AlgebraProduct(EigenVectors[j][2][l],(walpha - wbeta),AlgebraProducts,ProductList);

                                    if x <> false then
                                        for n in [1..dim] do
                                            if x[n] <> 0 then
                                                
                                                y := pairorbitlist[j][n];
                                            
                                                if AlgebraProducts[y] <> false then
                                                
                                                    a := [1..dim]*0; a[j] := 1;
                                                    b := [1..dim]*0; b[n] := 1;
                                                
                                                    Add(sum,-MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)*x[n]);
                                                else
                                                    if pairrepresentatives[y] = [j,n] or pairrepresentatives[y] = [n,j] then
                                                        row[Position(UnknownAlgebraProducts,y)]:=row[Position(UnknownAlgebraProducts,y)] + x[n];
                                                    else
                                                        row := [];
                                                        break; # need to move on to next beta2
                                                    fi;
                                                fi;
                                            fi;
                                        od;
                                    fi;

                                    # calculate rhs

                                    if row <> [] then 

                                        for n in [1..t] do
                                            if EigenVectors[j][2][l][n] <> 0 then
                                            
                                                y := pairorbitlist[n][c];
                                            
                                                if AlgebraProducts[y] <> false then
                                                
                                                    a := [1..dim]*0; a[c] := 1;
                                                    b := [1..dim]*0; b[n] := 1;
                                                
                                                    Add(sum,[EigenVectors[j][2][l][n]*MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)/4]);
                                                else
                                                    if pairrepresentatives[y] = [c,n] or pairrepresentatives[y] = [n,c] then
                                                        row[Position(UnknownAlgebraProducts,[n,c])]:=row[Position(UnknownAlgebraProducts,[n,c])] + EigenVectors[j][2][l][n]/4;
                                                    else 
                                                        row := [];
                                                        break; # need to move on to next beta2
                                                    fi;
                                                fi;
                                            fi;
                                        od;
                                    fi;

                                    y:= MAJORANA_InnerProduct(EigenVectors[j][2][m],EigenVectors[j][2][l],GramMatrix, pairorbitlist);
                                    
                                    if y <> false then 
                                    
                                        x := Sum(sum) + MAJORANA_AlgebraProduct(EigenVectors[j][2][l],walpha,AlgebraProducts,ProductList)/4 - y*a/4;
                                        
                                        if sum <> [] and row <> [] then
                                            if ForAll(row, x -> x = 0) then 
                                                if ForAny( x , y -> y <> 0) then 
                                                    Error("Step 8 system 2"); 
                                                fi;
                                            else
                                                Add(mat,row);
                                                Add(vec,x);
                                            fi;
                                        fi;
                                    fi;
                                od;
                            od;
                        od;
                    od;

                    if mat <> [] then 
                        Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                        if Size(Solution)  = 2 then
                            if Size(Solution[2]) = 0 then
                                for k in [1..Size(UnknownAlgebraProducts)] do
                                    x:=UnknownAlgebraProducts[k]; 
                                    AlgebraProducts[x]:=Solution[1][k];
                                od;
                        #    else
                        #        Output[i] := [Shape,"Fail","Missing algebra product values",GramMatrix, [], AlgebraProducts,EigenVectors];
                        #        Error("Missing algebra product values");
                        #        Output[i] := StructuralCopy(Output[i]);
                        #        break;
                            fi;
                        else
                            Output[i] := [Shape,"Error","Inconsistent system of unknown algebra products",mat,vec,record,AlgebraProducts,EigenVectors];
                            Output[i] := StructuralCopy(Output[i]);
                            Error("Inconsistent system of unknown algebra products");
                            break;
                        fi;
                    fi;
                fi;

                                            ## STEP 9: MORE EVECS II ##

                # Step 8 - check if we have full espace decomp, if not find it

                for j in [1..t] do
                    
                    a := [1..dim]*0; a[j] := 1;
                
                    if Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3]) + 1 <> dim then
                        mat:=[];

                        for k in [1..dim] do
                            b := [1..dim]*0; b[k] := 1;
                            x := MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                            if x <> false then
                                Add(mat,x);
                            else
                                mat := [];
                                break;
                            fi;
                        od;

                        if mat <> [] then 

                            EigenVectors[j][1]:=ShallowCopy(NullspaceMat(mat));
                            EigenVectors[j][2]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim)/4));
                            EigenVectors[j][3]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim)/32));
                            EigenVectors[j][4]:=ShallowCopy(NullspaceMat(mat - IdentityMat(dim) ));

                            
                            if LI = 0 then 
                                for k in [1..4] do 
                                    for l in [1..Size(EigenVectors[j][k])] do
                                        MAJORANA_RemoveNullSpace(EigenVectors[j][k][l],NullSp);
                                    od;
                                od;
                            fi;
                                                   
                            if Size(EigenVectors[j][4]) <> 1 then
                                Output[i]:=[Shape,"Error","Algebra does not obey axiom M5",GramMatrix,AlgebraProducts,EigenVectors];
                                Output[i]:=StructuralCopy(Output[i]);
                                break;
                            elif Size(EigenVectors[j][1])+Size(EigenVectors[j][2])+Size(EigenVectors[j][3]) + Size(EigenVectors[j][4]) > dim then
                                Output[i]:=[Shape,"Error","Algebra does not obey axiom M4",GramMatrix,AlgebraProducts,EigenVectors];
                                Error("M4");
                                Output[i]:=StructuralCopy(Output[i]);
                                break;
                            fi;
                        fi;
                    fi;
                    
                od;

                if Size(Output[i]) > 0 then
                    break;
                fi;
                
                newdimensions := [];
                
                for j in [1..t] do 
                    Add(newdimensions, Size(EigenVectors[j][1]) + Size(EigenVectors[j][2]) + Size(EigenVectors[j][3]) + 1);
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

                    Output[i] := [Shape,"Fail","Missing values",GramMatrix, [], AlgebraProducts,EigenVectors];
                    Error("Missing values");
                    Output[i] := StructuralCopy(Output[i]);
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
            
             

                                        ## STEP 10: RESURRECTION PRINCIPLE II ##

            # Check that eigenvectors obey the fusion rules

            ErrorFusion:=MAJORANA_Fusion(T,GramMatrix,AlgebraProducts,EigenVectors,ProductList);

            if ForAny(ErrorFusion,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Algebra does not obey fusion rules",GramMatrix,AlgebraProducts,EigenVectors,ErrorFusion];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Use resurrection principle

            UnknownAlgebraProducts := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts);

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

                                for n in [1..dim] do
                                    if EigenVectors[j][1][l][n] <> 0 then
                                    
                                        y := pairorbitlist[c][n];
                                        
                                        if AlgebraProducts[y] <> false then
                                            
                                            a := [1..dim]*0; a[c] := 1;
                                            b := [1..dim]*0; b[n] := 1;
                                    
                                            Add(sum, -MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)*EigenVectors[j][1][l][n]/4);
                                        else
                                            if pairrepresentatives[y] = [c,n] or pairrepresentatives[y] = [n,c] then
                                                row[Position(UnknownAlgebraProducts,y)]:=EigenVectors[j][1][l][n]/4;
                                            else
                                                row := [];
                                                break; ##### need to move on to next alpha2
                                            fi;
                                        fi;
                                    fi;
                                od;
                                
                                if row <> [] then 
                                    # calculate knowns
                                    
                                    a := [1..dim]*0; a[j] := 1;

                                    Append(sum,-[MAJORANA_AlgebraProduct(a,MAJORANA_AlgebraProduct(EigenVectors[j][1][l],(walpha - wbeta),AlgebraProducts,ProductList),AlgebraProducts,ProductList)]);

                                    Append(sum,-[MAJORANA_AlgebraProduct(EigenVectors[j][1][l],wbeta,AlgebraProducts,ProductList)]/4);
                                    
                                    x := Sum(sum);
                                    
                                    if ForAll(row, x -> x = 0) then 
                                        if ForAny( x , y -> y <> 0) then 
                                            Error("Step 10 system 1"); 
                                        fi;
                                    else
                                        Add(mat,row);
                                        Add(vec,x);
                                    fi;
                                fi; 

                            od;

                            ## Second part

                            for l in Beta2 do

                                sum:=[];
                                row:=NullMat(1,Size(UnknownAlgebraProducts))[1];

                                # calculate unknowns

                                for n in [1..dim] do
                                    if EigenVectors[j][2][l][n] <> 0 then
                                    
                                        y := pairorbitlist[c][n];
                                    
                                        if AlgebraProducts[y] <> false then
                                        
                                            a := [1..dim]*0; a[c] := 1;
                                            b := [1..dim]*0; b[n] := 1;
                                        
                                            Add(sum,MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList)*EigenVectors[j][2][l][n]/4);
                                        else
                                            if pairrepresentatives[y] = [c,n] or pairrepresentatives[y] = [n,c] then
                                                row[Position(UnknownAlgebraProducts,y)]:=row[Position(UnknownAlgebraProducts,y)] - EigenVectors[j][2][l][n]/4;;
                                            else
                                                row := [];
                                                break; # need to move on to next beta2
                                            fi;
                                        fi;
                                    fi;
                                od;
                                
                                if row <> [] then 

                                    # calculate knowns
                                    
                                    a := [1..dim]*0; a[j] := 1;

                                    Append(sum,-[MAJORANA_AlgebraProduct(a,MAJORANA_AlgebraProduct(EigenVectors[j][2][l],(walpha - wbeta),AlgebraProducts,ProductList),AlgebraProducts,ProductList)]);

                                    Append(sum,[MAJORANA_AlgebraProduct(EigenVectors[j][2][l],walpha,AlgebraProducts,ProductList)]/4);

                                    x:= MAJORANA_InnerProduct(EigenVectors[j][2][m],EigenVectors[j][2][l],GramMatrix, pairorbitlist);

                                    Append(sum,-[a*x/4]);

                                    x := Sum(sum);
                                    
                                    if ForAll(row, x -> x = 0) then 
                                        if ForAny( x , y -> y <> 0) then 
                                            Error("Step 10 system 2"); 
                                        fi;
                                    else
                                        Add(mat,row);
                                        Add(vec,x);
                                    fi;
                                fi;
                            od;
                        od;
                    od;
                od;

                Solution:=MAJORANA_SolutionMatVecs(mat,vec);

                if Size(Solution)  = 2 then
                    if Size(Solution[2]) = 0 then
                        for k in [1..Size(UnknownAlgebraProducts)] do

                            x:=UnknownAlgebraProducts[k]; 

                            AlgebraProducts[x]:=Solution[1][k];
                        od;
                    else
                        Output[i]:=[Shape,"Fail","Missing algebra products",GramMatrix,[],AlgebraProducts,EigenVectors];
                        Output[i]:=StructuralCopy(Output[i]);
                        break;
                    fi;
                else
                    Output[i]:=[Shape,"Error","Inconsistent system of unknown algebra products",GramMatrix,[],AlgebraProducts,EigenVectors];
                    Output[i]:=StructuralCopy(Output[i]);
                    break;
                fi;
            fi;

                                        ## STEP 11: CHECK ALGEBRA ##

            # Check bilinear form is positive definite

            if MAJORANA_PositiveDefinite(GramMatrixFull) <0 then
                Output[i]:=[Shape,"Error","Gram Matrix is not positive definite",GramMatrix, AlgebraProducts, EigenVectors];
                Output[i]:=StructuralCopy(Output[i]);
            fi;

            # Check that all triples obey axiom M1

            ErrorM1:=MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,Orbitals,ProductList,pairrepresentatives);

            if Size(ErrorM1)>0 then
                Output[i]:=[Shape,"Error","Algebra does not obey axiom M1",GramMatrix,AlgebraProducts,ErrorM1];
                Output[i]:=StructuralCopy(Output[i]);
            fi;

            # Check that eigenvectors obey the fusion rules

            ErrorFusion:=MAJORANA_Fusion(T,GramMatrix,AlgebraProducts,EigenVectors,ProductList);

            if ForAny(ErrorFusion,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Algebra does not obey fusion rules",GramMatrix,AlgebraProducts,EigenVectors,ErrorFusion];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Check that the eigenspaces are orthogonal

            ErrorOrthogonality:=MAJORANA_TestOrthogonality(T,GramMatrix,AlgebraProducts,EigenVectors,pairorbitlist);

            if ForAny(ErrorOrthogonality,x->Size(x)>0) then
                Output[i]:=[Shape,"Error","Eigenspaces are not orthogonal with respect to the inner product",GramMatrix,AlgebraProducts,EigenVector,ErrorOrthogonality];
                Output[i]:=StructuralCopy(Output[i]);
                break;
            fi;

            # Check M2

            ErrorM2:=MAJORANA_AxiomM2(GramMatrix,AlgebraProducts,ProductList,pairorbitlist);

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

    end );

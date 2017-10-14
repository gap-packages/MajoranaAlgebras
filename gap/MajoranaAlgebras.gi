#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

# Creates list of indexes [i,j] where product of i th and j th coordinate vectors is not known

BindGlobal( "MAJORANA_ExtractUnknownAlgebraProducts",

function(AlgebraProducts, ProductList)

    local   unknowns,       # list of unknown algebra products
            i,              # loop over coordinates
            j,              # loop over coordinates
            k,              # pair orbit index
            dim;            # size of coordinates
    
    unknowns := [];
    dim := Size(ProductList[1]);
    
    for i in [1..dim] do
        for j in [i + 1..dim] do 
            
            k := ProductList[3][i][j];
        
            if k < 0 then 
                k := -k;
            fi;
        
            if AlgebraProducts[k] = false then 
                Add(unknowns,[i,j]);
            fi;
        od;
    od;

    return AsSet(unknowns);
end);

# Finds the indices i such that v_i*v is not known

InstallGlobalFunction(MAJORANA_FindBadIndices,
    
    function(v, unknowns)
    
    local   i,
            j,
            bad;            
    
    bad := [];
                    
    for i in [1..Size(v)] do
        if v[i] <> 0 then 
            for j in [1 .. i - 1] do 
                if [j,i] in unknowns then 
                    Add(bad,j);
                fi;
            od;
            for j in [i + 1 .. Size(v)] do 
                if [i,j] in unknowns then 
                    Add(bad,j);
                fi;
            od;
        fi;
    od;

    bad := DuplicateFreeList(bad);

    Sort(bad);
    
    return bad;
    
    end );        
    
# given two eigenvectors, if possible, finds product and adds it to appropriate set of evecs

InstallGlobalFunction( MAJORANA_FuseEigenvectors,

    function(a,b,i,evals,mat,other_mat, new, GramMatrix, AlgebraProducts, ProductList)
    
    local   dim,
            u,
            new_ev,
            pos,
            x,
            y,
            z;
         
    dim := Size(ProductList[1]);
    u := [1..dim]*0; u[i] := 1;
    
    new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
    
    x := MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                    
    if x <> false then 
        if evals = [2,2] then 
            y := MAJORANA_InnerProduct(a,b,GramMatrix,ProductList);
            
            if y <> false then 
                Add(new[1], x - (1/4)*u*y);
            fi;
        elif evals = [3,3] then 
            y := MAJORANA_InnerProduct(a,b,GramMatrix,ProductList);
            z := MAJORANA_AlgebraProduct(u,x,AlgebraProducts, ProductList);
            
            if y <> false and z <> false then 
                Add(new[2], z - (1/32)*u*y);
            elif y <> false then 
                Add(other_mat, x - (1/32)*u*y);
            fi;
                
        else
            Add(new[pos],x);
        fi;
    else
        Add(mat, b);
    fi;

    return [new, mat, other_mat];
    
    end );

# finds new eigenvectors using the fusion rules 

InstallGlobalFunction( MAJORANA_Fusion,

function(GramMatrix, AlgebraProducts, EigenVectors, ProductList)

    local   i,
            j,
            k,
            a,
            b,
            dim,
            unknowns,
            new,
            u,
            evals,
            new_ev,
            pos,
            mat,
            x,
            y,
            z,
            null,
            bad,
            other_mat;
            
    for j in ProductList[10][1] do 
        for k in [1..3] do 
            if EigenVectors[j][k] <> [] then
                MAJORANA_ReversedEchelonForm(EigenVectors[j][k]);
            fi;
        od;
    od;
    
    dim := Size(ProductList[1]);
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    for i in ProductList[10][1] do 
        
        new := [ [], [], [] ];
        other_mat := [];
    
        for evals in [[1,1],[1,2],[2,1],[1,3],[3,1],[2,3],[2,2],[3,3]] do
            for a in EigenVectors[i][evals[1]] do
                mat := [];
                for b in EigenVectors[i][evals[2]] do 
                    MAJORANA_FuseEigenvectors(a, b, i, evals, mat, other_mat, new, GramMatrix, AlgebraProducts, ProductList);
                od;
                
                if mat <> [] then 
                
                    bad := MAJORANA_FindBadIndices(a,unknowns);
                    null := NullspaceMat(List(mat, x -> x{bad}));
                    
                    for j in [1..Size(null)] do 
                        
                        b := null[j]*mat;
                        
                        MAJORANA_FuseEigenvectors(a, b, i, evals, mat, other_mat, new, GramMatrix, AlgebraProducts, ProductList);
                        
                    od;
                fi; 
            od;
        od;
        
        if other_mat <> [] then 
            u := [1..dim]*0; u[i] := 1;
        
            bad := MAJORANA_FindBadIndices(u,unknowns);
            null := NullspaceMat(List(other_mat, x -> x{bad}));
            
            for j in [1..Size(null)] do 
                z := MAJORANA_AlgebraProduct(u,null[j]*other_mat,AlgebraProducts, ProductList);
                Add(new[2], z);
            od;
        fi;
        
        for j in [1..3] do 
            Append(EigenVectors[i][j], new[j]);
            EigenVectors[i][j] := ShallowCopy(BaseMat(EigenVectors[i][j]));
        od;
        
    od;
    
    return [true];    
    
    end );     
         
InstallGlobalFunction(MAJORANA_Append,

    function(x,mat,vec)

    local   i,          # loop over size of x
            pos;        # position of first non zero elt of row
    
    for i in [1..Size(x[1])] do
    
        pos := PositionNonZero(x[1][i]);
        x[2][i] := x[2][i]/x[1][i][pos];
        x[1][i] := x[1][i]/x[1][i][pos];
    
        if not x[1][i] in mat then             
            Add(mat,x[1][i]);
            Add(vec,x[2][i]);
        fi;
    od;
    
    end);

InstallGlobalFunction( MAJORANA_FindVectorPermutation, 
    
    function(g,ProductList)
    
    local   dim,        # size of coordinates
            i,          # loop over group elements
            j,          # loop over coordinates
            list,       # list to build permutation
            perm,       # the permutation
            signlist,   # corrects signs of 5A axes
            pos_1,      # position of conjugated element in longcoordinates
            pos_2;      # corresponding position in coordinates
    
    dim := Size(ProductList[1]);
    
    signlist := ListWithIdenticalEntries(dim,1);
    
    if g = () then 
        return [(),signlist];
    else
        list := [1..dim]*0;
        for j in [1..dim] do 
        
            pos_1 := Position(ProductList[2],ProductList[1][j]^g);
            pos_2 := ProductList[5][pos_1];
            
            if pos_2 > 0 then 
                list[j] := pos_2;
            else
                list[j] := -pos_2;
                signlist[-pos_2] := -1;
            fi;
        od;
        
        perm := PermList(list);
    
        return [perm,signlist];
    fi;
    
    end);
 
    
InstallGlobalFunction( MAJORANA_ConjugateVector, 

    function(v,g,ProductList)
    
    local   i,              # loop over vector
            dim,            # length of vector
            vec,            # output vector
            perm,
            sign,           # corrects sign of 5A axes 
            pos_1,          # position of conjugated element in longcoords
            pos_2;          # position of conjugated element in coords
    
    if g[1] <> () then 
        
        perm := g[1];
        sign := g[2];
        
        dim := Size(v);
        
        vec := [1..dim]*0;
        
        for i in [1..dim] do 
        
            if v[i] <> 0 then 
            
                vec[i^perm] := sign[i^perm]*v[i];
                
            fi;
        od;
        
        return vec;
    else
        return v;
    fi;
    
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

        dim:=Size(u);
        vec:=[1..dim]*0;

        if ForAll(u, x -> x = 0) or ForAll(v, x -> x = 0) then
            return u*0;
        fi;

        for i in [1..dim] do
            if u[dim - i + 1] <> 0 then 
                for j in [1..dim] do
                    if v[dim - j + 1] <> 0 then 
                    
                        k := list[3][dim - i + 1][dim - j + 1];
                        
                        if k > 0 then 
                            sign := 1;
                        else
                            sign := -1;
                            k := -k;
                        fi;

                        x := AlgebraProducts[k];
                        
                        if x <> false then
                            
                            g := list[4][dim - i + 1][dim - j + 1][1];
                            
                            vec := vec + sign*u[dim - i + 1]*v[dim - j + 1]*MAJORANA_ConjugateVector(x,g,list);
                        else
                            # cannot calculate product
                            return false;
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

        local   i,              # loop over u 
                j,              # loop over v
                k,              # pair orbit index
                sign,           # correct for 5A axes
                sum,            # output value
                unknowns,
                x;
        
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
                            k := -k;
                        fi;
                        
                        if GramMatrix[k] <> false then
                            sum := sum + sign*u[i]*v[j]*GramMatrix[k];
                        else
                            # cannot calculate product
                            
                            unknowns := Positions(GramMatrix, false);
                            
                            x := MAJORANA_SeparateInnerProduct(u,v,unknowns,GramMatrix,ProductList);
                            
                            if ForAll(x[1], y -> y = 0) then 
                                return x[2]; 
                            else
                                return false;
                            fi;
                            
                        fi;
                    fi;
                od;
            fi;
        od;
        
        return sum;
        
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

    local   row,            # record values of unknowns 
            sum,            # record values of knowns
            dim,            # size of coordinates
            i,              # index for dim of u
            j,              # index for dim of v
            m,              # orbit of i,j
            pos,            # position of m in unknowns
            sign;           # correct sign of 5A axes
            
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
                        m := -m;
                    fi;

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
        
        u := [1..dim]*0; u[i] := 1;
        ev_a := [u];
    
     else
        ev_a := EigenVectors[i][a];
    fi;
    
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
    
    if Size(OrthogonalityError) > 0 then 
        return [false, OrthogonalityError];
    else
        return [true,[mat,vec]];
    fi;
    
    end

    );
    
InstallGlobalFunction(MAJORANA_FullOrthogonality,

    function(EigenVectors,GramMatrix, AlgebraProducts, ProductList)
    
    local   i,          # loop over T
            j,          # loop over eigenvalues
            k,          # loop over eigenvalues
            x,          # res of orthogonality
            mat,        # matrix of unknown values
            vec,        # vector of known values
            ev,         # pair of eigenvalues     
            unknowns;     
            
    mat := [];
    vec := [];
    
    unknowns:=Positions(GramMatrix,false);
    
    if Size(unknowns) > 0 then 
    
        for i in ProductList[10][1] do        
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
    fi;
    

    if mat <> [] then 
        MAJORANA_SolutionInnerProducts(mat,vec, unknowns, GramMatrix);
    fi;
    
    if not false in GramMatrix then 
        x := MAJORANA_CheckNullSpace(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
        
        if x = false then
            return [false, "The inner product is not positive definite", []];  
        fi;
    fi;

    return [true, mat, vec];
    
    end );

    
InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(AlgebraProducts, EigenVectors, ProductList)

    local   i,          # loop over representatives
            ev,         # loop over eigenvalues
            unknowns,   # unknown algebra products
            mat,        # matrix of unknowns
            vec,        # vector of knowns
            table,      # table of eigenvalues
            u,          # vector with 1 in j th position
            v,          # eigenvector
            x,          # result of SeparateAlgebraProduct
            y,          # result of SolutionAlgProducts
            dim;        # size of ProductList[1]
    
    dim := Size(ProductList[1]);
    
    table := [0, 1/4, 1/32];
    
    mat := [];
    vec := [];
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
    
    for i in ProductList[10][1] do 
        for ev in [1..3] do 
        
            if EigenVectors[i][ev] <> [] then
                MAJORANA_ReversedEchelonForm(EigenVectors[i][ev]);
                EigenVectors[i][ev] := List(EigenVectors[i][ev], x -> MAJORANA_RemoveNullSpace(x,ProductList[6]));
            fi;
            
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
    od;
    
    Display("EigenVectors unknowns");
    
    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
            
    return y;
    
    end);
    
InstallGlobalFunction(MAJORANA_PositionLastOne,

function(list)

    local   len,    # length of list
            i;      # loop over list

    len := Size(list);

    for i in [1..len] do 
        if list[len - i + 1] = 1 then 
            return len - i + 1;
        fi;
    od;
    
    return false;
    
    end );
    
InstallGlobalFunction(MAJORANA_RemoveNullSpace,

function(v,NullSp) 

    local   i,      # loop over nullspace
            j;      # leading coefficient (from rhs)
    
    if NullSp <> false then 
        for i in [1..Size(NullSp[2])] do
            
            j := NullSp[1][i];
            
            if v[j] <> 0 then 
                v := v - v[j]*NullSp[2][i];
            fi;
        od;
    fi;
    
    return v;
    
    end
    
    );

InstallGlobalFunction(MAJORANA_UnknownsAxiomM1,

    function(GramMatrix, AlgebraProducts, ProductList)
    
    local   dim,
            mat,
            vec,
            i,
            j,
            k,
            u,
            v,
            w,
            x,
            y,
            z,
            row,
            sum,
            pos,
            unknowns;
    
    dim := Size(ProductList[1]);
    
    mat := [];
    vec := [];
    
    unknowns := Positions(GramMatrix, false);
    
    for i in [1..dim] do 
        
        u := [1..dim]*0; u[i] := 1;
    
        for j in [1..Size(AlgebraProducts)] do 
            
            if AlgebraProducts[j] <> false then 
            
                pos := ProductList[7][j];
                
                for k in [pos,Reversed(pos)] do
                
                    v := [1..dim]*0; v[k[1]] := 1;
                    w := [1..dim]*0; w[k[2]] := 1;
                
                    row := [];
                    sum := 0;
                
                    x := MAJORANA_SeparateInnerProduct(u, AlgebraProducts[j], unknowns, GramMatrix, ProductList);
                    
                    row := row + x[1];
                    sum := sum + x[2];
                
                    y := MAJORANA_AlgebraProduct(u, v, AlgebraProducts, ProductList);
                    
                    if y <> false then 
                        z := MAJORANA_SeparateInnerProduct(y, w, unknowns, GramMatrix, ProductList);
                        
                        row := row - z[1];
                        sum := sum - z[2];
                        
                        if ForAny(row, x -> x <> 0) then 
                            Add(mat,row);
                            Add(vec,[sum]);
                        fi;
                    fi;     
                od;
            fi;
        od;
    od;
    
    if mat <> [] then 
        MAJORANA_SolutionInnerProducts(mat,vec,unknowns,GramMatrix);#
    fi;
    
    end );
    
    
InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,UnknownAlgebraProducts,AlgebraProducts,ProductList)
    
    local   row,        # record values of unknowns
            sum,        # record values of knowns
            i,          # index for dim of u
            j,          # index for dim of v
            l,          # ordered version of [i,j]
            k,
            g,
            sign,
            x,          # vector with 1 in the ith position
            y,
            pos,        # position of unknown product 
            dim;        # dimension
    
    dim := Size(ProductList[1]);
    
    row := [1..Size(UnknownAlgebraProducts)]*0;
    sum := [1..dim]*0;
    
    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then
                
                    if i < j then 
                        l := [i,j];
                    else
                        l := [j,i];
                    fi;
                    
                    if not l in UnknownAlgebraProducts then 
                        
                        k := ProductList[3][i][j];
                        g := ProductList[4][i][j][1];
                        
                        if k > 0 then 
                            sign := 1;
                        else
                            sign := -1;
                            k := -k;
                        fi;
                        
                        x := AlgebraProducts[k];

                        sum := sum - sign*u[i]*v[j]*MAJORANA_ConjugateVector(x,g,ProductList);
                    
                    else
                        
                        pos := Position(UnknownAlgebraProducts,l);
                        row[pos] := row[pos] + u[i]*v[j]; 
                    fi;
                fi;
            od;
        fi;
    od;
    
    sum := MAJORANA_RemoveNullSpace(sum,ProductList[6]);
       
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
    
    if g <> () then 
    
        len     := Size(row);
        output  := [1..len]*0;
        
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
                    if x[1] < 0 then 
                        x[1] := -x[1];
                    else
                        x[2] := -x[2];
                    fi;
                else
                    sign := 1;
                    if x[1] < 0 then 
                        x := -x;
                    fi;
                fi;
                
                Sort(x);
                
                pos := Position(unknowns,x);
                output[pos] := sign*row[i];
                
            fi;
        od;
    
        return output;
    
    else
    
        return row;
    fi;
    
    end);     
    
InstallGlobalFunction(MAJORANA_Resurrection,

    function(alpha,beta,gamma,i,evals,mat,vec,bad_mat,GramMatrix, AlgebraProducts, EigenVectors, ProductList)   
    
    local   dim,
            u,
            unknowns,
            row,
            sum,
            x,
            y,
            z,
            ev,
            g;
            
    dim := Size(ProductList[1]);
    u := [1..dim]*0; u[i] := 1;
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);

    row := [];
    sum := [];
    
    x := MAJORANA_AlgebraProduct(alpha - beta, gamma, AlgebraProducts, ProductList);
    
    if x = false then 
        Add(bad_mat, alpha - beta);
    else                                
        y := MAJORANA_SeparateAlgebraProduct(u, x, unknowns, AlgebraProducts, ProductList);
        
        row := row + y[1];
        sum := sum + y[2];
        
        ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
        
        z := ev*MAJORANA_SeparateAlgebraProduct(beta,gamma,unknowns,AlgebraProducts,ProductList);
            
        row := row + z[1];
        sum := sum + z[2];
        
        if evals[1] = 2 then 
            z := MAJORANA_InnerProduct(alpha,gamma,GramMatrix,ProductList);
            
            if z <> false then 
                sum := sum + (1/4)*u*z;
            else 
                row := [];
            fi;
        fi;                               
    fi;
    
    if row <> [] and ForAny(row, x -> x <> 0) then 
        for g in DuplicateFreeList(ProductList[12]) do

            if g <> false then 
            
                if sum = [] then 
                    Error("empty");
                fi;
                
                y := MAJORANA_ConjugateVector(sum,g[2],ProductList);
                
                z := [,];
                
                z[1] := [MAJORANA_ConjugateRow(row,g[1],unknowns,ProductList)];
                z[2] := [MAJORANA_RemoveNullSpace(y,ProductList[6])];

                MAJORANA_Append(z,mat,vec);
                
            fi;
        od;
    fi;

    return fail;
    
    end );
    
InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(GramMatrix, AlgebraProducts, EigenVectors, ProductList)
    
    local   dim,
            i,
            j,
            evals,
            mat,
            vec,
            unknowns,
            alpha,
            beta,
            gamma,
            bad,
            null,
            bad_mat,
            x;
    
    dim := Size(ProductList[1]);
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    # Find unknown algebra products from eigenvectors
    
    x := MAJORANA_EigenvectorsAlgebraUnknowns(AlgebraProducts, EigenVectors, ProductList);
    
    mat := ShallowCopy(x[1]);
    vec := ShallowCopy(x[2]);
    unknowns := ShallowCopy(x[3]);
    
    # Find unknown algebra products from the resurrection principle
    
    for i in ProductList[10][1] do     
        for evals in [[1,2],[2,1],[1,3],[2,3]] do         
            for beta in EigenVectors[i][evals[2]] do
                for gamma in EigenVectors[i][evals[1]] do  
                
                    if MAJORANA_AlgebraProduct(beta,gamma,AlgebraProducts,ProductList) = false then
                    
                        bad_mat := [];
                                            
                        for alpha in EigenVectors[i][evals[1]] do 
                        
                             MAJORANA_Resurrection(alpha, beta, gamma, i, evals, mat,vec, bad_mat, GramMatrix, AlgebraProducts, EigenVectors, ProductList);

                        od;
                        
                        if bad_mat <> [] then
                            
                            bad := MAJORANA_FindBadIndices(gamma,unknowns);
                            
                            null := NullspaceMat(List(bad_mat, x -> x{bad}));
                            
                            for j in [1..Size(null)] do 
                            
                                alpha := null[j]*bad_mat + Sum(null[j])*beta;
                                
                                MAJORANA_Resurrection(alpha, Sum(null[j])*beta, gamma, i, evals, mat,vec, [], GramMatrix, AlgebraProducts, EigenVectors, ProductList);
                                
                            od;
                        fi;
                    fi;
                    
                    if mat <> [] and Size(mat) > Size(mat[1]) then 
                    
                        Display(["Resurrection", evals]);
            
                        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                                
                        mat := ShallowCopy(x[1]);
                        vec := ShallowCopy(x[2]);
                        unknowns := ShallowCopy(x[3]);
                    fi;
                    
                od;
            od;                        
        od;
    od;
    
    if mat <> [] and Size(mat) <= Size(mat[1]) then 
                    
        Display(["Resurrection final"]);

        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                
        mat := ShallowCopy(x[1]);
        vec := ShallowCopy(x[2]);
        unknowns := ShallowCopy(x[3]);
    fi;
    
    if ProductList[6] <> [] and ProductList[6] <> false then                  
                            
        x := MAJORANA_NullSpaceAlgebraProducts(unknowns, AlgebraProducts, ProductList);
        
        MAJORANA_Append(x,mat,vec);
        
        if mat <> [] then 
            
            Display("Nullspace");
        
            MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
        fi;            
    fi;
    
    end );
    
InstallGlobalFunction(MAJORANA_NullSpaceAlgebraProducts,

    function(UnknownAlgebraProducts, AlgebraProducts, ProductList)
    
    local i, j, m, k, row, sum, dim, y, mat, vec, a, x, record, pos;
    
    dim := Size(ProductList[1]);
    
    mat := [];
    vec := [];
    record := [];
    
    for j in Union(ProductList[10]{[1,2]}) do 
        a := [1..dim]*0; a[j] := 1;        

        for k in [1..Size(ProductList[6][2])] do
            
            x := MAJORANA_SeparateAlgebraProduct(a,ProductList[6][2][k],UnknownAlgebraProducts,AlgebraProducts,ProductList);

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
    
    
InstallGlobalFunction( MAJORANA_LongCoordinates,
    
    function(t,ProductList)
    
    local   i,      # loop over coordinates
            dim,    # size of coordinates
            x;      # coordinate;
            
    dim := Size(ProductList[1]);

    for i in [t+1..dim] do
        
        x := ProductList[1][i];
    
        if Order(x) = 3 then 
    
            Append(ProductList[5],[i,i]);
            Append(ProductList[2],[x,x^2]);
            
        elif Order(x) = 4 then 

            Append(ProductList[5],[i,i]);
            Append(ProductList[2],[x,x^3]);
            
        elif Order(x) = 5 then 
            Append(ProductList[5],[i,-i,-i,i]);
            Append(ProductList[2],[x,x^2,x^3,x^4]); 
        fi;
    od;
    
    ProductList[2] := Flat(ProductList[2]);
    
    end );
    
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
    
    local   sol,        # solution of system
            sign,       # correct sign of 5A axes
            i,          # loop over <UnknownAlgebraProducts>
            x,          # element of <UnknownAlgebraProducts>
            y,          # orbit of x
            g,          # conj element of x
            unsolved,
            prod,
            perm,
            nonzero,
            j;
    
    if mat <> [] then
    
        Display([Size(mat),Size(mat[1])]);
    
        sol := MAJORANA_SolutionMatVecs(mat,vec);
        
        Display("Solved it!");

        if sol <> false then
            for i in [1..Size(UnknownAlgebraProducts)] do
            
                if not i in sol[2] then
                
                    x := UnknownAlgebraProducts[i]; 
                    
                    y := ProductList[3][x[1]][x[2]];
                    g := ProductList[4][x[1]][x[2]][2];
                
                    
                    if y > 0 then 
                        sign := 1;
                    else
                        sign := -1;
                        y := -y;
                    fi;
                    
                    if AlgebraProducts[y] = false then 

                        AlgebraProducts[y] := sign*MAJORANA_ConjugateVector(sol[1][i],g,ProductList);

                        AlgebraProducts[y] := MAJORANA_RemoveNullSpace(AlgebraProducts[y],ProductList[6]);
                        
                    fi;                   
                fi;                
            od;
        fi;
        
        # take the remaining relations and remove any products which are now known
    
        unsolved := [];
        
        for i in [1..Size(UnknownAlgebraProducts)] do 
        
            x := UnknownAlgebraProducts[i]; 
                        
            y := ProductList[3][x[1]][x[2]];
            
            if y > 0 then 
                sign := 1;
            else
                sign := -1;
                y := -y;
            fi;
            
            prod := AlgebraProducts[y];
                            
            if prod <> false then 
                
                g := ProductList[4][x[1]][x[2]][1];
                    
                for j in [1..Size(sol[3][1])] do 
                    sol[3][2][j] := sol[3][2][j] - sign*sol[3][1][j][i]*MAJORANA_ConjugateVector(prod,g,ProductList);
                od; 
            else
                Add(unsolved,i);
            fi;
        od;
           
        nonzero := [];
                            
        for j in [1..Size(sol[3][1])] do  
            sol[3][1][j] := sol[3][1][j]{unsolved};
            
            if ForAny(sol[3][1][j], x -> x <> 0) then 
                sol[3][2][j] := MAJORANA_RemoveNullSpace(sol[3][2][j],ProductList[6]);
                Add(nonzero,j);
            fi;
        od;
        
        sol[3][1] := sol[3][1]{nonzero};
        sol[3][2] := sol[3][2]{nonzero};
        
        return [sol[3][1],sol[3][2],UnknownAlgebraProducts{unsolved}];
    else
        return [[],[],UnknownAlgebraProducts];
    fi;
    
    end );
    
InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, UnknownInnerProducts, GramMatrix)
    
    local   Solution,   # solution of system
            i,          # loop over <UnknownInnerProducts>
            x;          # element of <UnknownInnerProducts>    
    
    Solution := MAJORANA_SolutionMatVecs(mat,vec);

    if Solution <> false then                    
        
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
    fi;
    
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
                elif x > 0 then
                    ProductList[6] := [[],[]];
                fi; 
                
                for i in [1..Size(ProductList[6][2])] do
                    Add(ProductList[6][1], MAJORANA_PositionLastOne(ProductList[6][2][i]));
                od;
                              
            fi;

            if ProductList[6] <> false and ProductList[6][2] <> [] then

                # Change alg products to get rid of any axes not in the basis
                
                for i in [1..Size(ProductList[9])] do
                    if AlgebraProducts[i] <> false then
                        AlgebraProducts[i]:= MAJORANA_RemoveNullSpace(AlgebraProducts[i], ProductList[6]);
                    fi;
                od;

                # Change evecs to get rid of any axes not in the basis

                for j in ProductList[10][1] do
                    for k in [1..3] do                        
                        for x in [1..Size(EigenVectors[j][k])] do
                            EigenVectors[j][k][x] := MAJORANA_RemoveNullSpace(EigenVectors[j][k][x],ProductList[6]);
                        od;                                                   
                    od;                    
                od;
                
               # Append(EigenVectors[j][1],ProductList[6][2]);
            fi;
        fi;
        
        return true;
    
    end );
        
InstallGlobalFunction(MAJORANA_MoreEigenvectors,

    function(AlgebraProducts,EigenVectors,ProductList)
    
    local   i,
            j,
            dim,
            a,
            b,
            mat,
            x,
            table,
            ev;
            
    dim := Size(ProductList[1]);
    
    table := [0,1/4,1/32];

    for i in ProductList[10][1] do
                    
        a := [1..dim]*0; a[i] := 1;
    
        if Size(EigenVectors[i][1])+Size(EigenVectors[i][2])+Size(EigenVectors[i][3]) + 1 <> dim then
        
            mat:=[];

            for j in [1..dim] do
            
                b := [1..dim]*0; b[j] := 1;
                
                x := MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                
                if x <> false then
                    Add(mat,x);
                else
                    mat := [];
                    break;
                fi;
            od;

            if mat <> [] then 

                for ev in [1..3] do 
                    Display(["More eigenvectors",i,ev]); 
                    Append(EigenVectors[i][ev], ShallowCopy(NullspaceMat(mat - IdentityMat(dim)*table[ev])));
                    EigenVectors[i][ev] := ShallowCopy(BaseMat(EigenVectors[i][ev]));                    
                od;
            fi;
        fi; 
    od;
    
    return [true];
    
    end);

InstallGlobalFunction(ShapesOfMajoranaRepresentation,
    
    function(G,T)
    
    local   t,              # size of T
            i,              # indexes
            j,              #
            k,              #
            x,              # result of orbitals
            OrbitalsT,      # orbitals on T
            shape,          # one shape
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            Unknowns3X,     # orbitals which may be 3A or 3C
            Binaries,       # used to loop through options for shapes
            shapeslist,     # final list of shapes
            result;         #
    
    t := Size(T);

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;
    
    # Construct orbitals of G on T x T
    
    if IsAbelian(G) then 
        x := Cartesian(T,T);
    
        x := List(x,y -> [y]);

    else
        x := OrbitsDomain(G,Cartesian(T,T),OnPairs);
    fi;
    
    
    
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

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(OrbitalsT))[1];

    RepsSquares6A := [];

    for i in [1..Size(OrbitalsT)] do
    
        x := Representative(OrbitalsT[i]);
        
        if Order(x[1]*x[2]) = 1 then
            shape[i]:="1A";
        fi;
        if Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            shape[i]:="2A";
        fi;
        if Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            shape[i]:="2B";
        fi;
        if Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
        fi;
        if Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            shape[i]:="4A";
        fi;
        if Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            shape[i]:="4B";
        fi;
        if Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        fi;
        if Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(RepsSquares6A,(x[1]*x[2])^2);
        fi;
    od;

    # Check for inclusions of 3A in 6A

    for i in [1..Size(OrbitalsT)] do
        if shape[i][1] = '3' then
            j := 0;
            while j < Size(OrbitalsT[i]) do
                j := j + 1;
                if OrbitalsT[i][j][1]*OrbitalsT[i][j][2] in RepsSquares6A then
                    shape[i]:="3A";;
                    j:=Size(OrbitalsT[i])+1;;
                fi;
            od;
        fi;
    od;

    Unknowns3X:=[];

    for i in [1..Size(OrbitalsT)] do
        if shape[i] = ['3','X'] then
            Add(Unknowns3X,i);
        fi;
    od;
    
    Binaries:=AsList(FullRowSpace(GF(2),Size(Unknowns3X)));
    
    shapeslist := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do
        
        for j in [1..Size(Unknowns3X)] do
            k:=Unknowns3X[j];
            if Binaries[i][j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;            
        od;
        
        Add(shapeslist,ShallowCopy(shape));
        
    od;
    
    result := rec(  group := G,
                    involutions := T,
                    orbitals := OrbitalsT,
                    shapes := shapeslist);
    
    return result;

    end );

InstallGlobalFunction(MAJORANA_MainLoop,

    function(GramMatrix,AlgebraProducts,EigenVectors,ProductList)
    
    local   x, 
            dim,
            maindimensions,
            j,
            k;
            
    dim := Size(ProductList[1]);
    
                                ## STEP 5: INNER PRODUCTS M1 ##
                                
    MAJORANA_UnknownsAxiomM1(GramMatrix,AlgebraProducts,ProductList);
    
    x := MAJORANA_CheckNullSpace(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
    
    if x = false then
        return  MAJORANA_OutputError("The inner product is not positive definite"
                            , []
                            , ["Shape",GramMatrix,AlgebraProducts,EigenVectors,ProductList]);                                     
    fi;
                                                
                                ## STEP 6: FUSION ##                                        
                            

    # Use these eigenvectors and the fusion rules to find more
    
    maindimensions:=[];

    for j in ProductList[10][1] do
        for k in [1..3] do
            if Size(EigenVectors[j][k]) > 0 then
                EigenVectors[j][k] := ShallowCopy(BaseMat(EigenVectors[j][k]));
            fi;
        od;
        Add(maindimensions,   Size(EigenVectors[j][1])
                            + Size(EigenVectors[j][2])
                            + Size(EigenVectors[j][3])+1);
    od;

    if ForAny(maindimensions, x -> x < dim - 1) then                
    
        x := MAJORANA_Fusion(GramMatrix, AlgebraProducts,EigenVectors,ProductList);
        
        if not x[1] and ProductList[6] <> false then 
            return MAJORANA_OutputError(x[2],
                            x[3],
                            [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
        fi;
        
    fi;
    
                        ## STEP 8: RESURRECTION PRINCIPLE I ##
            
    MAJORANA_UnknownAlgebraProducts(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
    
                                ## STEP 9: MORE EVECS II ##

    # Check if we have full espace decomp, if not find it

    x := MAJORANA_MoreEigenvectors(AlgebraProducts,EigenVectors,ProductList);
    
    
                        ## STEP 10: INNER PRODUCTS FROM ORTHOGONALITY ##
        
                                
    # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products

    x := MAJORANA_FullOrthogonality(EigenVectors,GramMatrix, AlgebraProducts,ProductList);
    
    if not x[1] then 
        return MAJORANA_OutputError( x[2]
                        , x[3]
                        , [,GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
    fi;
    
    end);
    
        
InstallGlobalFunction(MajoranaRepresentation,

function(input,index)

    local   # Seress
            ProductList,  

            # indexing and temporary variables
            i, j, k, x, y, 

            # Step 0 - Set Up
            t, SizeOrbitals, OrbitalsT, G, T,
            
            # Step 1 - Shape
            Shape, 

            # Step 3 - Products and evecs I
            GramMatrix, GramMatrixFull, AlgebraProducts, EigenVectors, sign,

            # Step 4 - More products and evecs
            h, s, dim,
            
            vals, pos, OutputList, record, 

            falsecount, newfalsecount, maindimensions, newdimensions, switchmain;     

    OutputList :=  MAJORANA_SetUp(input,index);
    
    if Size(OutputList) <> 5 then 
        Error("Some error");
    fi;
    
    dim := Size(OutputList[5][1]);
    
    maindimensions:=[];

    for j in OutputList[5][10][1] do
        for k in [1..3] do
            if Size(OutputList[4][j][k]) > 0 then
                OutputList[4][j][k]:=ShallowCopy(BaseMat(OutputList[4][j][k]));
            fi;
        od;
        Add(maindimensions,   Size(OutputList[4][j][1])
                            + Size(OutputList[4][j][2])
                            + Size(OutputList[4][j][3])+1);
    od;
    
    falsecount := [0,0];
    
    if false in OutputList[2] then
        falsecount[1] := Size(Positions(OutputList[2],false));
    fi;
    
    if false in OutputList[3] then
        falsecount[2] := Size(Positions(OutputList[3],false));
    fi;
    
    if ForAll(maindimensions, x -> x = dim) and falsecount = [0,0] then 
        switchmain := 1;
    else
        switchmain := 0;
    fi;
    
    while switchmain = 0 do
                                
        MAJORANA_MainLoop(OutputList[2],OutputList[3],OutputList[4],OutputList[5]);
        
        newdimensions := [];
        
        for j in OutputList[5][10][1] do 
            Add(newdimensions, Size(OutputList[4][j][1]) 
                                + Size(OutputList[4][j][2]) 
                                + Size(OutputList[4][j][3]) + 1);
        od;
        
        newfalsecount := [0,0];
        
        if false in OutputList[2] then
            newfalsecount[1] := Size(Positions(OutputList[2],false));
        fi;
        
        if false in OutputList[3] then
            newfalsecount[2] := Size(Positions(OutputList[3],false));
        fi;
        
        Display([newfalsecount,falsecount]);
        
        if newfalsecount = [0,0] then
            break;
        elif newdimensions = maindimensions and newfalsecount = falsecount then

            return StructuralCopy(["Fail"
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
    od;

    x := StructuralCopy(["Success"
                ,
                ,
                , OutputList[1]
                , OutputList[2]
                , OutputList[3]
                , OutputList[4]
                , OutputList[5]]);
    
    return x;

end );

InstallGlobalFunction(MajoranaAlgebraTest,
    
    function(res)
    
    local   error,
            GramMatrixFull;
                            
    # Check bilinear form is positive definite
    
    GramMatrixFull := MAJORANA_FillGramMatrix(res[5], res[8]);

    if not false in res[5] and MAJORANA_PositiveDefinite(GramMatrixFull) <0 then
        return "Gram Matrix is not positive definite";
    fi;

    # Check that all triples obey axiom M1

    error := MAJORANA_AxiomM1(res[5],res[6],res[8]);

    if Size(error)>0 then
        return ["Algebra does not obey axiom M1", error];
    fi;

    # Check that eigenvectors obey the fusion rules

    error := MAJORANA_TestFusion(res[5],res[6],res[7],res[8]);

    if ForAny(error,x->Size(x)>0) then
        return ["Algebra does not obey fusion rules", error];
    fi;

    # Check that the eigenspaces are orthogonal

    error := MAJORANA_TestOrthogonality(res[5],res[6],res[7],res[8]);

    if Size(error) > 0 then
        return ["Eigenspaces are not orthogonal with respect to the inner product", error];
    fi;

    # Check M2

    # error := MAJORANA_AxiomM2(res[5],res[6],res[8]);

    # if error = -1 then
    #    return "Algebra does not obey axiom M2";
    # fi;
    
    return true;
    
    end );

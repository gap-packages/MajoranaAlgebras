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
    
    function(v, AlgebraProducts, ProductList)
    
    local   i,
            j,
            k,
            dim,
            bad;            
    
    bad := [];
    dim := Size(ProductList[1]);
                    
    for i in [1..dim] do
        if v[i] <> 0 then 
            for j in [1..dim] do 
               k :=  ProductList[3][i][j];
               
               if k > 0 then 
                    if AlgebraProducts[k] = false then 
                        Add(bad,j);
                    fi;
                else
                    if AlgebraProducts[-k] = false then 
                        Add(bad,j);
                    fi;
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

    function(a,b,i,evals,other_mat, new, GramMatrix, AlgebraProducts, ProductList)
    
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
            Add(new[1], x + (3/32)*u*y - 4*z);
        elif y <> false then 
            Add(other_mat, x - (1/32)*u*y);
        fi;
            
    else
        Add(new[pos],x);
    fi;
    
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
            
    # looks like we should get rid of this - check on a big example
            
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
                
                bad := MAJORANA_FindBadIndices(a,AlgebraProducts,ProductList);
                null := NullspaceMat(List(EigenVectors[i][evals[2]], x -> x{bad}));
                
                for j in [1..Size(null)] do 
                    
                    b := null[j]*EigenVectors[i][evals[2]];
                    
                    MAJORANA_FuseEigenvectors(a, b, i, evals, other_mat, new, GramMatrix, AlgebraProducts, ProductList);
                    
                od;
            od;
        od;
        
        if other_mat <> [] then 
            u := [1..dim]*0; u[i] := 1;
        
            bad := MAJORANA_FindBadIndices(u,AlgebraProducts, ProductList );
            null := NullspaceMat(List(other_mat, x -> x{bad}));
            
            for j in [1..Size(null)] do 
                z := MAJORANA_AlgebraProduct(u,null[j]*other_mat,AlgebraProducts, ProductList);
                Add(new[2], z);
                Add(new[1], null[j]*other_mat - 4*z);
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
            pos_1,          # position of conjugated element in longcoords
            pos_2;          # position of conjugated element in coords
    
    if g[1] <> () then 
        
        dim := Size(v);
        
        vec := [1..dim]*0;
        
        for i in [1..dim] do 
        
            if v[i] <> 0 then 
            
                vec[i^g[1]] := g[2][i^g[1]]*v[i];
                
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
                vecs,
                elts,
                pos,
                dim;    # size of vectors 

        dim:=Size(u);
        vec:=[1..dim]*0;

        if ForAll(u, x -> x = 0) or ForAll(v, x -> x = 0) then
            return u*0;
        fi;

        elts := [];
        vecs := [];

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
                            
                            pos := Position(elts,g);
                            
                            if pos <> fail then 
                                vecs[pos] := vecs[pos] + sign*u[dim - i + 1]*v[dim - j + 1]*x;
                            else
                                Add(elts,g);
                                Add(vecs,sign*u[dim - i + 1]*v[dim - j + 1]*x);
                            fi;
                        else
                            # cannot calculate product
                            return false;
                        fi;
                    fi;
                od;
            fi;
        od;
        
        for i in [1..Size(elts)] do 
            vec := vec + MAJORANA_ConjugateVector(vecs[i],elts[i], list);
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
            z,          # to be added to mat vec system
            g,          # conjugating element
            dim;        # size of ProductList[1]
    
    dim := Size(ProductList[1]);
    
    table := [0, 1/4, 1/32];
    
    mat := [];
    vec := [];
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
    
    for i in ProductList[10][1] do 
        for ev in [1..3] do 
            
            u := [1..dim]*0; u[i] := 1;
            
            for v in EigenVectors[i][ev] do
                
                x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,AlgebraProducts,ProductList);
                
                if ForAll(x[1], y -> y = 0) then 
                    if ForAny((x[2] + table[ev]*v), y -> y <> 0 ) then 
                        return [false,1,v];
                    fi;
                else                          
                    x[1] := [x[1]];
                    x[2] := [x[2] + table[ev]*v];
                    MAJORANA_Append(x,mat,vec);
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
            elts,
            vecs,
            x,          # vector with 1 in the ith position
            y,
            pos,        # position of unknown product 
            dim;        # dimension
    
    dim := Size(ProductList[1]);
    
    row := [1..Size(UnknownAlgebraProducts)]*0;
    sum := [1..dim]*0;
    
    elts := [];
    vecs := [];
    
    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then
                
                    k := ProductList[3][i][j];
                    
                    if k > 0 then 
                        sign := 1;
                    else
                        sign := -1;
                        k := -k;
                    fi;
                    
                    x := AlgebraProducts[k];
                    
                    if x <> false then 
                                                
                        g := ProductList[4][i][j][1];
                        
                        pos := Position(elts,g);
                        
                        if pos <> fail then 
                            vecs[pos] := vecs[pos] - sign*u[i]*v[j]*x;
                        else
                            Add(elts,g);
                            Add(vecs,- sign*u[i]*v[j]*x);
                        fi;
                    else
                    
                        if i < j then
                            l := [i,j];
                        else
                            l := [j,i];
                        fi;
                        
                        pos := Position(UnknownAlgebraProducts,l);
                        row[pos] := row[pos] + u[i]*v[j]; 
                    fi;
                fi;
            od;
        fi;
    od;
    
    for i in [1..Size(elts)] do 
        sum := sum + MAJORANA_ConjugateVector(vecs[i],elts[i],ProductList);
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
    
InstallGlobalFunction(MAJORANA_AddConjugates,

    function(row, sum, mat, vec, unknowns, ProductList)   
    
    local   g,
            x;
    
    if ForAny(row, x -> x <> 0) then 
        for g in ProductList[12] do
                
            x := [,];
            
            x[1] := [MAJORANA_ConjugateRow(row,g[1],unknowns,ProductList)];
            x[2] := MAJORANA_ConjugateVector(sum,g[2],ProductList);
            x[2] := [MAJORANA_RemoveNullSpace(x[2],ProductList[6])];

            MAJORANA_Append(x,mat,vec);
        od;
    fi;
    
    end );
    
InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(GramMatrix, AlgebraProducts, EigenVectors, ProductList)
    
    local   dim,
            i,
            j,
            u,
            evals,
            mat,
            vec,
            unknowns,
            alpha,
            beta,
            gamma,
            null,
            alpha_mat,
            bad,
            n,
            ev,
            row,
            sum,
            old_mat,
            x,
            y,
            z,
            w;
    
    dim := Size(ProductList[1]);;
    
    # Find unknown algebra products from eigenvectors
    
    x := MAJORANA_EigenvectorsAlgebraUnknowns(AlgebraProducts, EigenVectors, ProductList);
    
    mat := ShallowCopy(x[1]);
    vec := ShallowCopy(x[2]);
    unknowns := ShallowCopy(x[3]);
    
    # Find unknown algebra products from the resurrection principle
    
    for i in ProductList[10][1] do     
        
        u := [1..dim]*0;; u[i] := 1;;
    
        for evals in [[1,2],[2,1],[1,3],[2,3]] do  
        
            ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
            
            for beta in EigenVectors[i][evals[2]] do
                for gamma in EigenVectors[i][evals[1]] do  
                
                    x := MAJORANA_SeparateAlgebraProduct(beta, gamma, unknowns, AlgebraProducts, ProductList); 
                
                    if ForAny(x[1], y -> y <> 0) then
                        
                        alpha_mat := [];
                                            
                        for alpha in EigenVectors[i][evals[1]] do                             
                            Add(alpha_mat, alpha - beta);                        
                        od;
                        
                        bad := MAJORANA_FindBadIndices(gamma, AlgebraProducts, ProductList);
                        
                        null := NullspaceMat(List(alpha_mat, x-> x{bad}));
                        
                        for j in [1..Size(null)] do 
                        
                            n := Sum(null[j]);
                        
                            if n <> 0 then
                            
                                row := [];
                                sum := [];
                                
                                y := MAJORANA_AlgebraProduct(null[j]*alpha_mat, gamma, AlgebraProducts, ProductList);
                                
                                z := MAJORANA_SeparateAlgebraProduct(u, y, unknowns, AlgebraProducts, ProductList);
                            
                                row := row + z[1] + n*ev*x[1];
                                sum := sum + z[2] + n*ev*x[2];
                                
                                if evals[1] = 2 then 
                                    w := MAJORANA_InnerProduct(null[j]*alpha_mat, gamma, GramMatrix, ProductList);
                                    
                                    if w <> false then 
                                        sum := sum + (1/4)*w*u;
                                    else
                                        row := [];
                                    fi;
                                fi;
                                
                                if row <> [] then 
                                    MAJORANA_Append([[row],[sum]],mat,vec);
                                fi;
                            fi;
                        od;
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
    
    
    if false then                  
                            
        x := MAJORANA_NullSpaceAlgebraProducts(unknowns, AlgebraProducts, ProductList);
        
        MAJORANA_Append(x,mat,vec);
        
        if mat <> [] then 
            
            Display("Nullspace");
        
            x := MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
            
            mat := ShallowCopy(x[1]);
            vec := ShallowCopy(x[2]);
            unknowns := ShallowCopy(x[3]);
        fi;            
    fi;
    
    if mat <> [] then 
    
        old_mat := ShallowCopy(mat);
    
        for i in [1..Size(old_mat)] do 
            
            if unknowns <> [] then 
                MAJORANA_AddConjugates(mat[i],vec[i],mat,vec,unknowns,ProductList);
                
                if Size(mat) > Size(mat[1]) then 
                
                    Display(["All conjugates",i]);
                
                    x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                        
                    mat := ShallowCopy(x[1]);
                    vec := ShallowCopy(x[2]);
                    unknowns := ShallowCopy(x[3]);
                fi;
            fi;
        od;
        
        Display("All conjugates");
        
        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                    
        mat := ShallowCopy(x[1]);
        vec := ShallowCopy(x[2]);
        unknowns := ShallowCopy(x[3]);
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
                
                prod := MAJORANA_ConjugateVector(prod,g,ProductList);
                    
                for j in [1..Size(sol[3][1])] do 
                    sol[3][2][j] := sol[3][2][j] - sign*sol[3][1][j][i]*prod;
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

                ProductList[6] := MAJORANA_NullSpace(GramMatrixFull); 
                
                for x in ProductList[6][2] do
                    Add(ProductList[6][1], MAJORANA_PositionLastOne(x));
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
                
               Append(EigenVectors[j][1],ProductList[6][2]);
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

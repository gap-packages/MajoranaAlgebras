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
        for j in [i..dim] do 
            
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

# This creates new eigenvectors through fusion rules

InstallGlobalFunction( MAJORANA_Fusion,

function(a, b, j, AlgebraProducts, EigenVectors, GramMatrix, ProductList)

    local   u,                  # vector with 1 in the j th position
            x,                  # new eigenvector
            y,                  # inner product of u and x
            z,                  # algebra product of u and x
            w,                  # further product in 1/32, 1/32 case
            NewEigenVectors,    # list of new eigenvectors
            k,                  # run through ev_a-eigenvectors
            l,                  # run through ev_b-eigenvectors
            ev_a,               # first eigenvector
            ev_b,               # second eigenvector
            ev,                 # new eigenvector
            pos,                # index of new eigenvector
            dim,                # size of coordinates              
            FusionError,        # list of indexes which do not obey fusion
            unknowns,
            alpha,
            list;
            
    dim := Size(ProductList[1]);

    ev := MAJORANA_FusionTable[a+1][b+1];
    
    NewEigenVectors := [];
    FusionError := [];

    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);

    u := [1..dim]*0;; u[j] := 1;

    ev_a := EigenVectors[j][a];
    ev_b := EigenVectors[j][b];

    # the 1/4,1/4 case is special
    if (a=2) and (b=2) then 
        for k in [1..Size(ev_a)] do
            for l in [k..Size(ev_b)] do

                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );

                if x = false then 
                    alpha := MAJORANA_FindAlpha(j,a,b,ev_a[k],ev_b[l],unknowns,EigenVectors,ProductList);
                    
                    if alpha <> fail then 
                        x := MAJORANA_AlgebraProduct(alpha - ev_a[k],ev_b[l],AlgebraProducts,ProductList);
                    fi;
                fi;
                
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
            for l in [k..Size(ev_b)] do
                
                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );

                if x = false then 
                    alpha := MAJORANA_FindAlpha(j,a,b,ev_a[k],ev_b[l],unknowns,EigenVectors,ProductList);
                    
                    if alpha <> fail then 
                        x := MAJORANA_AlgebraProduct(alpha - ev_a[k],ev_b[l],AlgebraProducts,ProductList);
                    fi;
                fi;
                
                if x <> false then 
                    
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                    
                    if y <> false then
                        
                        z := MAJORANA_AlgebraProduct( u, x, AlgebraProducts, ProductList );
                        
                        if z <> false then 
                            z := z - y*u;
                            
                            w := MAJORANA_AlgebraProduct( u, z, AlgebraProducts, ProductList);
                        
                            if (w <> false) and ( w <> z/4) then
                                Add(FusionError,[j,k,l]);
                            else
                                Add(NewEigenVectors,z);
                            fi;
                        fi;
                    fi;
                    
                fi;
                
            od;
        od;
    elif not [a,b] in [[2,2],[3,3]] then 
        for k in [1..Size(ev_a)] do
            
            if a = b then 
                list := [k..Size(ev_b)];
            else
                list := [1..Size(ev_b)];
            fi;
            
            for l in list do

                x := MAJORANA_AlgebraProduct( ev_a[k], ev_b[l], AlgebraProducts, ProductList );
                
                if x = false then 
                    alpha := MAJORANA_FindAlpha(j,b,a,ev_a[k],ev_b[l],unknowns,EigenVectors,ProductList);
                    
                    if alpha <> fail then 
                        x := MAJORANA_AlgebraProduct(alpha - ev_a[k],ev_b[l],AlgebraProducts,ProductList);
                    fi;
                fi;                
                
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
            new;                # new eigenvectors
            
    for j in ProductList[10][1] do 
        for k in [1..3] do 
            if EigenVectors[j][k] <> [] then
                MAJORANA_ReversedEchelonForm(EigenVectors[j][k]);
            fi;
        od;
    od;

    for j in ProductList[10][1] do
        
        new := [ [], [], [] ];
        
        for k in [[1,1],[1,2],[1,3],[2,3],[2,2],[3,3]] do
        
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
        
        for k in [1..3] do
            Append(EigenVectors[j][k],new[k]);
            if Size(EigenVectors[j][k]) > 0 then
                EigenVectors[j][k]:=ShallowCopy(BaseMat(EigenVectors[j][k]));
            fi;
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
 
    
InstallGlobalFunction( MAJORANA_ConjugateVector, # can we store this instead of g?

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
    
    unknowns:=[];

    for j in [1..Size(ProductList[9])] do
        if GramMatrix[j] = false then
            Add(unknowns,j);
        fi;
    od;
    
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
    
    for i in ProductList[10][1] do
    
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

            MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
        fi;
    od;
            
    return [true];
    
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
    
    function(j,k, GramMatrix, AlgebraProducts, ProductList)
    
    local   i,                      # loop over representatives of orbitals
            dim,                    # size of coordinates
            pos,                    # pair representatives of orbital i
            res,                    # inner product of j th and k th basis vectors
            x,                      # intermediate products
            y,                      #
            z,                      #
            u,                      # basis vectors
            v,                      #
            w;                      #
    
    i := 1;
    
    dim := Size(ProductList[1]);
    
    res := false;
    
    while i < Size(AlgebraProducts) + 1 do 
    
        pos := ProductList[7][i];
    
        if AlgebraProducts[i] <> false and AlgebraProducts[i][j] <> 0 then 
        
            u := [1..dim]*0; u[k] := 1;
            v := [1..dim]*0; v[pos[1]] := 1;
        
            x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
            if x <> false then
                
                u := [1..dim]*0; u[pos[2]] := 1;
            
                y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                
                if y <> false then 
                    
                    v := StructuralCopy(AlgebraProducts[i]); v[j] := 0;
                    
                    w := [1..dim]*0; w[k] := 1;
                    
                    z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList);
                    
                    if z <> false then 
                                           
                        res := (y - z)/AlgebraProducts[i][j]; 
                                          
                    fi;
                fi;
                
            fi;
            
            if res <> false then 
                
                i := Size(AlgebraProducts) +  1;
            
            else
            
                u := [1..dim]*0; u[k] := 1;
                v := [1..dim]*0; v[pos[2]] := 1;
            
                x := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
                if x <> false then
                    
                    u := [1..dim]*0; u[pos[1]] := 1;
                
                    y := MAJORANA_InnerProduct(u, x, GramMatrix, ProductList);
                    
                    if y <> false then 
                        
                        v := StructuralCopy(AlgebraProducts[i]); v[j] := 0;
                        
                        w := [1..dim]*0; w[k] := 1;
                        
                        z := MAJORANA_InnerProduct(w,v,GramMatrix, ProductList);
                        
                        if z <> false then 
                                               
                            res := (y - z)/AlgebraProducts[i][j];  
                                              
                        fi;
                    fi;
                    
                fi;
            fi;
            
            if res <> false then 
                i := Size(AlgebraProducts) + 1;
            else
                i := i + 1;
            fi;
            
        else
        
            i := i + 1;
            
        fi;
    od;
    
    return res;
    
    end );
    
InstallGlobalFunction(MAJORANA_FullUnknownsAxiomM1,

    function(GramMatrix,AlgebraProducts,ProductList)

    local   dim,            # 
            i,              # loop over orbitals
            j,              # loop over an orbital 
            k,              # position of first elementz
            l,              # position of second element
            res,            # result of unknowns axiom M1
            sign;           # correct sign of 5A axes
        
    dim := Size(ProductList[1]);
        
    for i in [1..dim] do 
        for j in [1..dim] do 
        
            k := ProductList[3][i][j];
            
            if k < 0 then 
                sign := -1; 
            else
                sign := 1;
            fi;
            
            if GramMatrix[sign*k] = false then 
            
                res := MAJORANA_UnknownsAxiomM1(i,j, GramMatrix, AlgebraProducts, ProductList);
                
                if res <> false then 
                    GramMatrix[sign*k] := sign*res;
                fi;
            fi;
        od;
    od;
    
    end );
    
    
InstallGlobalFunction(MAJORANA_FullUnknownsAxiomM1Old,

    function(GramMatrix,AlgebraProducts,ProductList)

    local   switch,         # switch for while loop
            count,          # counts any newly found products
            pos_1,          # positions of orbitals rep in longcoordinates
            pos_2,          # 
            i,              # loop over orbitals
            j,              # loop over an orbital 
            k,              # position of first elementz
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
                
                    pos_1 := Position(ProductList[2],ProductList[9][i][j][1]);
                    pos_2 := Position(ProductList[2],ProductList[9][i][j][2]);
                    
                    k := ProductList[5][pos_1];
                    l := ProductList[5][pos_2];
                    
                    if k*l > 0 then 
                        sign := 1;
                        
                        if k < 0 then 
                            k := -k;
                            l := -l;
                        fi;
                        
                    else
                        sign := -1;
                        
                        if k < 0 then 
                            k := -k;
                        else
                            l := -l;
                        fi;
                    fi;
                    
                    res := MAJORANA_UnknownsAxiomM1(k,l,GramMatrix,AlgebraProducts,ProductList);
                    
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
    
InstallGlobalFunction(MAJORANA_UnknownAlgebraProductsAxiomM1, 

    function(j,k,T,GramMatrix,AlgebraProducts,ProductList)
    
    local   unknowns,           # unknown alg products
            dim,                # size of coordinates
            t,                  # size of T 
            mat,                # matrix of inner product values
            vec,                # vector of inner product values
            u,                  # basis vectors
            v,                  #
            w,                  #
            x,                  # unknown algebra product
            GramMatrixFull,     # part of the full gram matrix
            i,                  # loop over coordinates
            y,                  # intermediate products
            z;                  #
            
    t := Size(T);
    dim := Size(ProductList[1]);
    
    GramMatrixFull := MAJORANA_FillGramMatrix(GramMatrix,ProductList);
    
    v := [1..dim]*0; v[j] := 1;
    w := [1..dim]*0; w[k] := 1;

    mat := [];
    vec := [];
    
    for i in [1..dim] do 
    
        u := [1..dim]*0; u[i] := 1;
        
        y := MAJORANA_AlgebraProduct(u,v,AlgebraProducts,ProductList);
        
        if y <> false then 
            z := MAJORANA_InnerProduct(y,w,GramMatrix,ProductList);
            
            Add(mat,GramMatrixFull[i]);
            Add(vec,z);
        fi;
    od;            

    return([mat,vec]);
    
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
            k,
            g,
            sign,
            x,          # vector with 1 in the ith position
            y,
            pos,        # position of unknown product 
            dim;        # dimension
    
    dim := Size(AlgebraProducts[1]);
    
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
                    
                        if true then 
                        
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
                        
                            x := [1..dim]*0; x[i] := 1;
                            y := [1..dim]*0; y[j] := 1;
                            
                            sum := sum - u[i]*v[j]*MAJORANA_AlgebraProduct(x,y,AlgebraProducts,ProductList);
                        fi;
                        
                        
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
    
InstallGlobalFunction(MAJORANA_FindAlpha,

    function(i,ev_a,ev_b,beta,gamma,UnknownAlgebraProducts,EigenVectors,ProductList)
    
    local   j, 
            k,
            dim,
            bad,
            list,
            x,
            alpha,
            pos;
        
    dim := Size(ProductList[1]);
            
    list := [];
    bad := [];
                   
    for j in [1..dim] do
        if gamma[j] <> 0 then 
            for k in [1 .. j - 1] do 
                if [k,j] in UnknownAlgebraProducts then 
                    Add(bad,k);
                fi;
            od;
            for k in [j + 1 .. dim] do 
                if [j,k] in UnknownAlgebraProducts then 
                    Add(bad,k);
                fi;
            od;
        fi;
    od;
    
    bad := DuplicateFreeList(bad);
    
    Sort(bad);

    for j in [1..Size(EigenVectors[i][ev_b])] do 
        Add(list, StructuralCopy(EigenVectors[i][ev_b][j]{bad}));
    od;

    if ev_a = ev_b and beta in list then 
        pos := Position(list,beta);
        list[pos] := [1..Size(bad)]*0;
    else
        pos := 0;
    fi;
    
    x := SolutionMat(list,beta{bad});
    
    if x <> fail then 
    
        alpha := [];
    
        for j in [1..Size(x)] do
            if ev_a <> ev_b or j <> pos then
                alpha := alpha + x[j]*EigenVectors[i][ev_b][j];
            fi;
        od;
        
        return alpha;
        
    else
        return fail;
    fi;
    
    end );   
        
    
InstallGlobalFunction(MAJORANA_Resurrection,

    function(i,ev_a, ev_b, EigenVectors,UnknownAlgebraProducts,AlgebraProducts,ProductList,GramMatrix)
    
    local   a, 
            b, 
            sum, 
            row, 
            alpha, 
            beta, 
            gamma, 
            mat, 
            vec, 
            bad, 
            list, 
            j, 
            m, 
            k, 
            l, 
            dim, 
            u, 
            v, 
            x, 
            y, 
            z, 
            g, 
            sign, 
            record;
    
    dim := Size(AlgebraProducts[1]);
    
    mat := [];
    vec := [];
    record := [];
    
    for a in [1..Size(EigenVectors[i][ev_a])] do 
        for b in [1..Size(EigenVectors[i][ev_b])] do 
        
            if MAJORANA_AlgebraProduct(EigenVectors[i][ev_a][a],EigenVectors[i][ev_b][b],AlgebraProducts,ProductList) = false then 
        
                row := [1..Size(UnknownAlgebraProducts)]*0;
                sum := [1..dim]*0;
                
                beta := EigenVectors[i][ev_a][a];
                gamma := EigenVectors[i][ev_b][b];
                
                list := [];
                bad := [];
                
                x := MAJORANA_SeparateAlgebraProduct(MAJORANA_FusionTable[ev_a + 1][ev_b + 1]*beta,gamma,UnknownAlgebraProducts,AlgebraProducts,ProductList);
                
                row := row + x[1];
                sum := sum + x[2];
                   
                alpha := MAJORANA_FindAlpha(i,ev_a,ev_b,beta,gamma,UnknownAlgebraProducts,EigenVectors,ProductList);
                
                if alpha <> fail then 
                
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
                                if ProductList[6] <> false then 
                                    Error("Resurrection error");
                                fi;
                            fi;
                        else
                            
                            for g in DuplicateFreeList(ProductList[12]) do
                                if g <> false then 
                                
                                    if sum = [] then 
                                        Error("empty");
                                    fi;
                                    
                                    y := MAJORANA_ConjugateVector(sum,g[2],ProductList);
                                    Add(mat,MAJORANA_ConjugateRow(row,g[1],UnknownAlgebraProducts,ProductList));
                                    Add(vec,MAJORANA_RemoveNullSpace(y,ProductList[6]));
                                    
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

    function(EigenVectors,AlgebraProducts,ProductList,GramMatrix)
    
    local mat, vec, record, j, k, l, x, g, t,unknowns,y,u,v,dim, table;
    
    mat := [];
    vec := [];
    record := [];
    
    table := [0, 1/4, 1/32];
    
    dim := Size(ProductList[1]);
    
    if false in AlgebraProducts then 
    
        t := Size(EigenVectors);
        
        # put eigenvectors into reversed echelon form 
        
        unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
        
        for j in ProductList[10][1] do 
            for k in [1..3] do 
            
                if EigenVectors[j][k] <> [] then
                    MAJORANA_ReversedEchelonForm(EigenVectors[j][k]);
                fi;
                
                u := [1..dim]*0; u[j] := 1;
                
                for v in EigenVectors[j][k] do
                    
                    x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,AlgebraProducts,ProductList);
                    
                    if ForAll(x[1], y -> y = 0) then 
                        if ForAny((x[2] + table[k]*v), y -> y <> 0 ) then 
                            return [false,1,v];
                        fi;
                    elif not x[1] in mat then         
                        Add(mat, x[1]);
                        Add(vec, x[2] + table[k]*v); 
                        Add(record, [j,k,v]);
                    fi;
                od;
            od;
        od;
        
        Display("EigenVectors unknowns");
        
        y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                
        mat := ShallowCopy(y[1]);
        vec := ShallowCopy(y[2]);
        unknowns := ShallowCopy(y[3]);
        record := [];
        
        for j in ProductList[10][1] do 
                
            for k in [1..3] do
            
                for l in [1..2] do 
                    if not [k,l] in [[1,1],[2,2]] then 
                
                        x := MAJORANA_Resurrection(j,k,l,EigenVectors,unknowns,AlgebraProducts,ProductList,GramMatrix);
                        
                        if x[1] <> [] then 
                            MAJORANA_Append(x,mat,vec);
                            Append(record, x[3]);
                            
                            Display([k,l,"Fusion"]);
                        
                            y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                            
                            mat := ShallowCopy(y[1]);
                            vec := ShallowCopy(y[2]);
                            unknowns := ShallowCopy(y[3]);
                        fi;
                        
                    fi;
                od;
            od;
        od;
        
        if ProductList[6] <> [] and ProductList[6] <> false then                  
                                
            x := MAJORANA_NullSpaceAlgebraProducts(unknowns, AlgebraProducts, ProductList);
            
            MAJORANA_Append(x,mat,vec);
            
            if mat <> [] then 
                MAJORANA_SolutionAlgProducts(mat,vec, unknowns, AlgebraProducts, ProductList);
            fi;            
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
            perm,
            gp,
            perms,
            dim,    # size of coordinates
            g;      # conjugating element
            
    table := [[],[1],[1,2],[1,3],[1,2,3,4]];
    
    dim := Size(ProductList[1]);
    
    gp := AsList(ProductList[8]);
    perms := List(gp, x -> MAJORANA_FindVectorPermutation(x,ProductList));
    
    for i in [1..dim] do
        for j in [1..dim] do
            if ProductList[4][i][j] = 0 then 
            
                x := [ProductList[1][i],ProductList[1][j]];
                
                # create list of other indices which are going to have this elt
                
                list := [];
                
                for k in table[Order(x[1])] do 
                    for l in table[Order(x[2])] do
                        Add(list,[x[1]^k,x[2]^l]);
                    od;
                od;
                
                # find orbit
                
                k := ProductList[3][i][j];
            
                if k < 0 then 
                    k := -k;
                fi; 
            
                # find rep of orbit
            
                y := ProductList[9][k][1];
                
                l := 1;
                
                while l < Size(list) + 1 do 
            
                    g := RepresentativeAction(ProductList[8],y,list[l],OnPairs);
                
                    if g <> fail then
                    
                        perm := [perms[Position(gp,g)],perms[Position(gp,Inverse(g))]];
                    
                        for z in list do
                            
                            pos_1 := ProductList[5][Position(ProductList[2],z[1])];
                            pos_2 := ProductList[5][Position(ProductList[2],z[2])];
                            
                            if pos_1 < 0 then 
                                pos_1 := -pos_1;
                            fi;
                            
                            if pos_2 < 0 then 
                                pos_2 := -pos_2;
                            fi;
                            
                            ProductList[4][pos_1][pos_2] := perm;
                            ProductList[4][pos_2][pos_1] := perm;
                        od;
                        
                        l := Size(list) + 1;
                        
                    else 
                             
                        g := RepresentativeAction(ProductList[8],y,Reversed(list[l]),OnPairs);
                    
                        if g <> fail then 
                        
                            perm := [perms[Position(gp,g)],perms[Position(gp,Inverse(g))]];
                        
                            for z in list do
                            
                                pos_1 := ProductList[5][Position(ProductList[2],z[1])];
                                pos_2 := ProductList[5][Position(ProductList[2],z[2])];
                                
                                if pos_1 < 0 then 
                                    pos_1 := -pos_1;
                                fi;
                                
                                if pos_2 < 0 then 
                                    pos_2 := -pos_2;
                                fi;
                                
                                ProductList[4][pos_1][pos_2] := perm;
                                ProductList[4][pos_2][pos_1] := perm;
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
                            
            pos_1 := ProductList[5][Position(ProductList[2],z[1])];
            pos_2 := ProductList[5][Position(ProductList[2],z[2])];
            
            if pos_1 < 0 then 
                pos_1 := -pos_1;
            fi;
            
            if pos_2 < 0 then 
                pos_2 := -pos_2;
            fi;
                            
            
            ProductList[4][pos_1][pos_2] := [[(),[]],[(),[]]];
            ProductList[4][pos_2][pos_1] := [[(),[]],[(),[]]];
            
            ProductList[3][pos_1][pos_2] := i;
            ProductList[3][pos_2][pos_1] := i;
        od;
    od;
    
    end);
    
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
    
InstallGlobalFunction( MAJORANA_SetupOrbitals,

    function(ProductList, OrbitalsT)
    
    local   G,      # the group
            i,      # loop over orbitals
            x,      # representative of orbital
            y,      # orders of representatives
            table,  # table of orders of axes
            j,      # loop over orders of first axis
            k,      # loop over orders of second axis
            l;      # loop over orbitals
            
    G := ProductList[8];
            
    x := Cartesian(ProductList[1],ProductList[1]);

    ProductList[9] := List(Orbits(G,x,OnPairs),ShallowCopy);

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
    
InstallGlobalFunction(MAJORANA_SetupOrbits,

    function(T,ProductList)
    
    local   G,      # the group
            x,      # orbit    
            i,      # loop over T
            j,      # loop over orbits 
            k,      # representative
            t,      # size of T
            orbits, # orbits of G on C
            pos;    # position in longcoordinates
    
    # Construct orbits of G on T
    
    G := ProductList[8];
    t := Size(T);
    
    ProductList[11] := OrbitsDomain(G,T);
    
    for x in ProductList[11] do
        Add(ProductList[10][1], Position( T, Representative(x)));
    od;
    
    for i in [1..t] do
        
        j := 1;
        
        while j < Size(ProductList[11]) + 1 do 
            if T[i] in ProductList[11][j] then 
            
                ProductList[13][i] := j;
                
                k := ProductList[10][1][j];
                
                ProductList[12][i][1] := RepresentativeAction(G,T[k],T[i]);
                ProductList[12][i][2] := MAJORANA_FindVectorPermutation(ProductList[12][i][1],ProductList);
                
                j := Size(ProductList[11]) + 1;;
                
            else
                j := j + 1;
            fi;
        od;
    od;
    
    orbits := Orbits(G, ProductList[1]{[t+1 .. Size(ProductList[1])]});
    
    for i in [1..Size(orbits)] do 
        pos := Position(ProductList[2], Representative(orbits[i]));
        
        x := ProductList[5][pos];
        
        if x < 0 then 
            x := -x;
        fi;
        
        Add(ProductList[10][2],x);
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
                    
                    # Append(EigenVectors[j][1],ProductList[6][2]); 
                    
                od;
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
    
InstallGlobalFunction(MAJORANA_ConnectedComponents,

    function(T)
    
    local   graph,      # graph with vertices elements of T
            t,          # size of T
            i,          # loop over T
            j;          # loop over T
            
    LoadPackage("automata");
    
    t := Size(T);
    graph := NullMat(t,0);
    
    for i in [1..t] do 
        for j in [1..t] do 
            if Order(T[i]*T[j]) <> 2 or not T[i]*T[j] in T then 
                Add(graph[i],j);
            fi;
        od;
    od;
    
    return AutoConnectedComponents(graph);
    
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
    
InstallGlobalFunction(MAJORANA_DihedralProducts,

function(T,Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList)

    local   j,      # loop over reps of T
            k,      # loop over T
            x,      # orbital of i,j
            y,      #
            pos,    # position 
            vals,   # values    
            dim,    # size of coordinates
            s,      # elt of T
            h,      # axis
            t,      # size of T
            sign;   # corrects sign of 5A axes
            
    dim := Size(ProductList[1]);
    t   := Size(T);
    
    # (2,2) products and eigenvectors from IPSS10

    # Add eigenvectors from IPSS10

    for j in ProductList[10][1] do
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
            
                x := Position(ProductList[2],T[j]*T[k]);
            
                pos := [j, k, 0, 0, 0, 0];
                pos[3] := Position(T, T[j]*T[k]*T[j]);
                pos[4] := Position(T, T[k]*T[j]*T[k]);
                pos[5] := Position(T, T[j]*T[k]*T[j]*T[k]*T[j]);
                pos[6] := ProductList[5][x]; 

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

    for j in [1..Size(ProductList[9])] do

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
    
    end );
    
InstallGlobalFunction( MAJORANA_SetUp,
    
    function(input, index)
    
    local   G,              # group
            T,              # involutions
            Shape,          # chosen shape
            OrbitalsT,      # orbits of G on T
            t,              # size of T
            ProductList,    # lists needed to calculate products
            x,              # temporary variables
            i,              # index variables
            j,
            k,
            dim,            # size of coordinates
            SizeOrbitals,   # number of orbits of G on T x T
            GramMatrix,
            AlgebraProducts,
            EigenVectors;

    G           := input.group;
    T           := input.involutions;
    Shape       := input.shapes[index];
    OrbitalsT   := input.orbitals;  
    
    t := Size(T);  
    
    # Set up ProductList
    
    ProductList := [1..13]*0;
    
    ProductList[1] := StructuralCopy(T);

    # Create coordinates list

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
    
    ProductList[2]  := StructuralCopy(T);
    ProductList[3]  := NullMat(dim,dim);
    ProductList[4]  := NullMat(dim,dim);
    ProductList[5]  := [1..t];
    ProductList[6]  := false;
    ProductList[7]  := [];    
    ProductList[8]  := G;
    ProductList[10] := [[],[]];
    ProductList[12] := NullMat(t,2);
    ProductList[13] := [1..t]*0;
    ProductList[14] := [];

    MAJORANA_LongCoordinates(t, ProductList);
    
    MAJORANA_SetupOrbits(T, ProductList);
    
    MAJORANA_SetupOrbitals(ProductList, OrbitalsT);

    MAJORANA_PairRepresentatives(ProductList);

    MAJORANA_PairOrbits(ProductList);

    MAJORANA_PairConjElements(ProductList);
    
                                ## STEP 3: PRODUCTS AND EVECS I ##
                                
    SizeOrbitals := Size(ProductList[9]);

    # Set up algebra product and gram matrices

    AlgebraProducts := NullMat(1,SizeOrbitals)[1];
    GramMatrix := NullMat(1,SizeOrbitals)[1];

    for j in [1..SizeOrbitals] do
        AlgebraProducts[j]:=false;
        GramMatrix[j]:=false;
    od;

    # Set up eigenvector matrix

    EigenVectors := NullMat(t,3);

    for j in [1..t] do
        if j in ProductList[10][1] then
            for k in [1..3] do
                EigenVectors[j][k] := [];
            od;
        else
            for k in [1..3] do
                EigenVectors[j][k] := false;
            od;
        fi;
    od;
    
    # Start filling in values and products!

    MAJORANA_DihedralProducts(T,Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList);
    
    # Check algebra obeys axiom M1 at this stage

    x := MAJORANA_AxiomM1(GramMatrix,AlgebraProducts,ProductList);;
    
    if Size(x) > 0 then 
        return MAJORANA_OutputError( "Algebra does not obey axiom M1"
                            , x
                            , [Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList]);
    fi;
    
    return [Shape, GramMatrix, AlgebraProducts, EigenVectors, ProductList];
    
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
                                
    MAJORANA_FullUnknownsAxiomM1(GramMatrix,AlgebraProducts,ProductList);
    
    x := MAJORANA_CheckNullSpace(GramMatrix,AlgebraProducts,EigenVectors,ProductList);
    
    if x = false then
        return  MAJORANA_OutputError("The inner product is not positive definite"
                            , []
                            , [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);                                     
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
    
        x := MAJORANA_FullFusion(AlgebraProducts,EigenVectors, GramMatrix, ProductList);
        
        if not x[1] and ProductList[6] <> false then 
            return MAJORANA_OutputError(x[2],
                            x[3],
                            [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
        fi;
    fi;
                                ## STEP 7: PRODUCTS FROM EIGENVECTORS ##

    # Use eigenvectors to find more products
    
    if false then 
        x := MAJORANA_EigenvectorsAlgebraUnknowns(EigenVectors,AlgebraProducts,ProductList);
        
        if not x[1] then 
            if x[2] = 1 then 
                return MAJORANA_OutputError( "Error eigenvectors algebra unknowns"
                        , x[3]
                        , [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
            elif x[2] = 2 then 
                return MAJORANA_OutputError("Inconsistent system of unknown algebra products step 7"
                        , x[3]
                        , [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
            fi;
        fi;
    fi;
    
                        ## STEP 8: RESURRECTION PRINCIPLE I ##
            
    MAJORANA_FullResurrection(EigenVectors,AlgebraProducts,ProductList,GramMatrix);
    
                                ## STEP 9: MORE EVECS II ##

    # Check if we have full espace decomp, if not find it

    #x := MAJORANA_MoreEigenvectors(AlgebraProducts,EigenVectors,ProductList);
    
    #if not x[1] then 
    #     return MAJORANA_OutputError(x[2],[], OutputList);
    #fi;
    
                        ## STEP 10: INNER PRODUCTS FROM ORTHOGONALITY ##
        
                                
    # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products

    x := MAJORANA_FullOrthogonality(EigenVectors,GramMatrix, AlgebraProducts,ProductList);
    
    if not x[1] then 
        return MAJORANA_OutputError( x[2]
                        , x[3]
                        , [GramMatrix,AlgebraProducts,EigenVectors,ProductList]);
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
        
        if newfalsecount = [0,0] then
            break;
        elif newdimensions = maindimensions and newfalsecount = falsecount then
        
        Error("fail"); 

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

    if MAJORANA_PositiveDefinite(GramMatrixFull) <0 then
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

     error := MAJORANA_AxiomM2(res[5],res[6],res[8]);

    if error = -1 then
        return "Algebra does not obey axiom M2";
    fi;
    
    return true;
    
    end );

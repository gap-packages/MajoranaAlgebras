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
    dim := Size(ProductList.coords);
    
    for i in [1..dim] do
        for j in [i + 1..dim] do 
            
            k := ProductList.pairorbit[i][j];
        
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
            list,
            bad;            
    
    bad := [];
    dim := Size(ProductList.coords);
    list := [1..dim];
                    
    for i in [1..dim] do
        if v[i] <> 0 then 
            for j in list do 
                k :=  ProductList.pairorbit[i][j];
               
                if k > 0 then 
                    if AlgebraProducts[k] = false then 
                        Add(bad,j);
                        list := Difference(list,[j]);
                    fi;
                else
                    if AlgebraProducts[-k] = false then 
                        Add(bad,j);
                        list := Difference(list,[j]);
                    fi;
                fi;                
            od;
        fi;
    od;

    Sort(bad);
    
    return bad;
    
    end );        
    
# given two eigenvectors, if possible, finds product and adds it to appropriate set of evecs

InstallGlobalFunction( MAJORANA_FuseEigenvectors,

    function(a,b,i,evals,other_mat, new, GramMatrix, AlgebraProducts, ProductList)
    
    local   dim,
            u,
            test,
            new_ev,
            pos,
            x,
            y,
            z;
         
    dim := Size(ProductList.coords);
    u := [1..dim]*0; u[i] := 1;
    
    new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
    
    x := MAJORANA_AlgebraProduct(a,b,AlgebraProducts,ProductList);
                    
    if evals = [2,2] then 
        y := MAJORANA_InnerProduct(a,b,GramMatrix,ProductList);
        
        if y <> false then 
            Add(new[1], x - (1/4)*u*y);
            
            test := MAJORANA_AlgebraProduct( u, x - (1/4)*u*y, AlgebraProducts, ProductList); 
            
            if test <> false and not MAJORANA_InnerProduct(test,test,GramMatrix, ProductList) in [0, false] then 
                Error("Fusion error");
            fi;
            if not MAJORANA_InnerProduct( u, x - (1/4)*u*y, GramMatrix, ProductList) in [0, false] then 
                Error("Orthog error");
            fi;
        fi;
    elif evals = [3,3] then 
        y := MAJORANA_InnerProduct(a,b,GramMatrix,ProductList);
        z := MAJORANA_AlgebraProduct(u,x,AlgebraProducts, ProductList);
        
        if y <> false and z <> false then 
            Add(new[2], z - (1/32)*u*y);
            Add(new[1], x + (3/32)*u*y - 4*z);
            
            test := MAJORANA_AlgebraProduct( u, x + (3/32)*u*y - 4*z, AlgebraProducts, ProductList);
            
            if test <> false and not MAJORANA_InnerProduct(test,test,GramMatrix, ProductList) in [0, false] then 
                Error("Fusion error");
            fi;
            #if not MAJORANA_AlgebraProduct( u, z - (1/32)*u*y, AlgebraProducts, ProductList) in [(z - (1/32)*u*y)/4, false] then 
            #    Error("Fusion error");
            #fi;
            
            if not MAJORANA_InnerProduct( u, x + (3/32)*u*y - 4*z, GramMatrix, ProductList) in [0, false] then 
                Error("Orthog error");
            fi;
            if not MAJORANA_InnerProduct( u, z - (1/32)*u*y, GramMatrix, ProductList) in [0, false] then 
                Error("Orthog error");
            fi;
            
        elif y <> false then 
            Add(other_mat, x - (1/32)*u*y);
        fi;
            
    else
        Add(new[pos],x);
        
        test := MAJORANA_AlgebraProduct(u, x, AlgebraProducts, ProductList);
        
        if test <> false and not MAJORANA_InnerProduct(test - x*new_ev,test - x*new_ev,GramMatrix,ProductList) in [0, false] then 
            Error("Fusion error");
        fi;
        if not MAJORANA_InnerProduct(u, x, GramMatrix, ProductList) in [0, false] then 
            Error("Orthog error");
        fi;
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
            evecs,
            new_ev,
            pos,
            mat,
            x,
            y,
            z,
            null,
            bad,
            other_mat;
    
    dim := Size(ProductList.coords);
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts, ProductList);
    
    for i in ProductList.orbitreps do 
    
        if IsMutable(EigenVectors[i][1]) and Sum(List(EigenVectors[i], x -> Size(x))) < dim - 1 then 

            new := [ [], [], [] ];
            other_mat := [];
        
            for evals in [[1,1],[1,2],[2,2],[1,3],[2,3],[3,3]] do

                evecs := EigenVectors[i][evals[2]];
                
                for a in EigenVectors[i][evals[1]] do
                    if Size(evecs) <> 0 then 
                        bad := MAJORANA_FindBadIndices(a,AlgebraProducts,ProductList);
                        
                        if bad <> [] then  
                            null := NullspaceMat(List(evecs, x -> x{bad}));
                        else
                            null := IdentityMat(Size(evecs));
                        fi;
                        
                        for j in [1..Size(null)] do 
                            
                            b := null[j]*evecs;
                            
                            MAJORANA_FuseEigenvectors(a, b, i, evals, other_mat, new, GramMatrix, AlgebraProducts, ProductList);
                            
                        od;
                    fi;
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
                if EigenVectors[i][j] <> [] then 
                    EigenVectors[i][j] := ShallowCopy(BaseMat(EigenVectors[i][j]));
                fi;
            od;
        fi;
    od;
    
    return [true];    
    
    end );     
         
InstallGlobalFunction(MAJORANA_Append,

    function(x,mat,vec)

    local   pos;        # position of first non zero elt of row
        
    pos := PositionNonZero(x[1]);
    x[2] := x[2]/x[1][pos];
    x[1] := x[1]/x[1][pos];

    if not x[1] in mat then             
        Add(mat,x[1]);
        Add(vec,x[2]);
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
                    
                        #if i = j and Order(list.coords[dim - i + 1]) <> 5 then 
                        #    vec[dim - i + 1] := u[dim - i + 1]*v[dim - j + 1];
                        # fi;
                    
                        k := list.pairorbit[dim - i + 1][dim - j + 1];
                        
                        if k > 0 then 
                            sign := 1;
                        else
                            sign := -1;
                            k := -k;
                        fi;

                        x := AlgebraProducts[k];
                        
                        if x <> false then
                            
                            g := list.pairconj[dim - i + 1][dim - j + 1][1];
                            
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
                    
                        k := ProductList.pairorbit[i][j];
                        
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
                                return x[2][1]; 
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

    local   i,                  # loop over ProductList.coords
            j,                  # loop over ProductList.coords
            k,                  # pair orbit index
            dim,                # size of ProductList.coords
            GramMatrixFull;     # output matrix
    
    dim := Size(ProductList.coords);
    
    GramMatrixFull := NullMat(dim,dim);
    
    for i in [1..dim] do 
        for j in [1..dim] do
            
            k := ProductList.pairorbit[i][j];
            
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
            
    dim := Size(ProductList.coords);
            
    sum := 0;
    row := [1..Size(UnknownInnerProducts)]*0;

    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then 
                
                    m := ProductList.pairorbit[i][j];
                    
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

    return [row,[sum]];

    end );
    
InstallGlobalFunction(MAJORANA_Orthogonality,

    function(EigenVectors,GramMatrix, ProductList)
    
    local   i,          # loop over T
            ev,         # loop over eigenvalues
            k,          # loop over eigenvalues
            evecs_a,    #
            evecs_b,    #
            u,
            v,
            dim,
            x,          # res of orthogonality
            mat,        # matrix of unknown values
            vec,        # vector of known values   
            unknowns;     
            
    mat := [];
    vec := [];
    
    dim := Size(ProductList.coords);
    
    unknowns := Positions(GramMatrix,false);
    
    if Size(unknowns) = 0 then
        return;
    fi; 
    
    for i in ProductList.orbitreps do        
        for ev in Combinations([0..3],2) do  
            if ev[1] = 0 then 
                u := [1..dim]*0; u[i] := 1;
                evecs_a := [u];  
            else
                evecs_a := EigenVectors[i][ev[1]];
            fi;
            
            evecs_b := EigenVectors[i][ev[2]];
                
            for u in evecs_a do
                for v in evecs_b do

                    x := MAJORANA_SeparateInnerProduct( u, v,
                                                        unknowns,
                                                        GramMatrix,
                                                        ProductList);

                    if ForAny(x[1], y -> y <> 0) then
                        MAJORANA_Append(x, mat, vec);
                    fi;
                    
                od;
            od;
        od;
    od;
    
    if mat <> [] then 
        MAJORANA_SolutionInnerProducts(mat,vec, unknowns, GramMatrix);
    fi;        

    end );    
    
InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(GramMatrix, AlgebraProducts, EigenVectors, ProductList)

    local   i,          # loop over representatives
            t,
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
            dim;        # size of ProductList.coords
    
    dim := Size(ProductList.coords);
    t := Size(EigenVectors);
    
    table := [0, 1/4, 1/32];
    
    mat := [];
    vec := [];
    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(AlgebraProducts,ProductList);
    
    if ForAny(unknowns, x -> x[1] <= t) then 
    
    Info( InfoMajorana, 50, "Building eigenvector unknowns");
    
        for i in ProductList.orbitreps do 
            for ev in [1..3] do 
                
                u := [1..dim]*0; u[i] := 1;
                
                for v in EigenVectors[i][ev] do
                    
                    x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,AlgebraProducts,ProductList);
                    
                    x[2] := x[2] + table[ev]*v;
                    
                    if ForAll(x[1], y -> y = 0) then 
                        
                        if ForAny(x[2], y -> y <> 0 ) then
                            
                            y := MAJORANA_InnerProduct(x[2],x[2],GramMatrix, ProductList);
                            
                            if not y in [0, false] then 
                                Error("EigenVectors alg unknowns");
                            fi;
                        fi;
                    else 
                        MAJORANA_Append(x,mat,vec);
                    fi;                
                od;
            od;
        od;
        
        y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                
        return y;
    else
        return rec( mat := [],
                    vec := [],
                    unknowns := unknowns) ;
    fi;
    
    end);
    
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
    
    dim := Size(ProductList.coords);
    
    mat := [];
    vec := [];
    
    unknowns := Positions(GramMatrix, false);
    
    if unknowns = [] then 
        return;
    fi;
    
    for i in [1..dim] do 
        
        u := [1..dim]*0; u[i] := 1;
    
        for j in [1..Size(AlgebraProducts)] do 
            
            if AlgebraProducts[j] <> false then 
            
                pos := ProductList.pairreps[j];
                
                for k in [pos,Reversed(pos)] do
                
                    v := [1..dim]*0; v[k[1]] := 1;
                    w := [1..dim]*0; w[k[2]] := 1;
                
                    row := [];
                    sum := [0];
                
                    x := MAJORANA_SeparateInnerProduct(u, AlgebraProducts[j], unknowns, GramMatrix, ProductList);
                    
                    row := row + x[1];
                    sum := sum + x[2];
                
                    y := MAJORANA_AlgebraProduct(u, v, AlgebraProducts, ProductList);
                    
                    if y <> false then 
                        z := MAJORANA_SeparateInnerProduct(y, w, unknowns, GramMatrix, ProductList);
                        
                        row := row - z[1];
                        sum := sum - z[2];
                        
                        if ForAny(row, x -> x <> 0) then 
                            MAJORANA_Append([row,sum], mat, vec);
                        fi;
                    fi;     
                od;
            fi;
        od;
    od;
    
    if mat <> [] then 
        MAJORANA_SolutionInnerProducts(mat,vec,unknowns,GramMatrix);
    fi;
    
    end );

InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,unknowns,AlgebraProducts,ProductList)
    
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
    
    dim := Size(ProductList.coords);
    
    row := [1..Size(unknowns)]*0;
    sum := [1..dim]*0;
    
    elts := [];
    vecs := [];
    
    for i in [1..dim] do
        if u[i] <> 0 then
            for j in [1..dim] do
                if v[j] <> 0 then
                
                    k := ProductList.pairorbit[i][j];
                    
                    if k > 0 then 
                        sign := 1;
                    else
                        sign := -1;
                        k := -k;
                    fi;
                    
                    x := AlgebraProducts[k];
                    
                    if x <> false then 
                                                
                        g := ProductList.pairconj[i][j][1];
                        
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
                        
                        pos := Position(unknowns,l);
                        row[pos] := row[pos] + u[i]*v[j]; 
                    fi;
                fi;
            od;
        fi;
    od;
    
    for i in [1..Size(elts)] do 
        sum := sum + MAJORANA_ConjugateVector(vecs[i],elts[i],ProductList);
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
    
    if g <> () then 
    
        len     := Size(row);
        output  := [1..len]*0;
        
        for i in [1..Size(row)] do
            if row[i] <> 0 then 
        
                j := unknowns[i][1];
                k := unknowns[i][2];
                
                x := [0,0];
                
                pos_1 := Position(ProductList.longcoords,(ProductList.coords[j])^g);
                pos_2 := Position(ProductList.longcoords,(ProductList.coords[k])^g);
                
                x[1] := ProductList.poslist[pos_1];
                x[2] := ProductList.poslist[pos_2];
                
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

    function(mat,vec,unknowns,AlgebraProducts,ProductList) 
    
    local   new_mat,
            new_vec,
            nonzero,
            i,
            j,
            g,
            y,
            x;
        
    new_mat := [];
    new_vec := [];
            
    for i in [1..Size(mat)] do
        
        if ForAny(mat[i], x -> x <> 0) then 
            for g in ProductList.conjelts do
                            
                x := [,];
                
                x[1] := MAJORANA_ConjugateRow(mat[i],g[1],unknowns,ProductList);
                x[2] := MAJORANA_ConjugateVector(vec[i],g[2],ProductList);

                MAJORANA_Append(x,new_mat,new_vec);
            od;
        fi;
    od;
        
    if new_mat <> [] then 
        
        x := MAJORANA_SolutionAlgProducts(new_mat, new_vec, unknowns, AlgebraProducts, ProductList);
        
        if x.unknowns = [] then 
            return rec( mat := [], vec := [], unknowns := []);
        fi; 
        
        y := MAJORANA_RemoveKnownAlgProducts(   mat, 
                                                vec, 
                                                unknowns, 
                                                AlgebraProducts, 
                                                ProductList);
    
        
        new_mat := x.mat;
        new_vec := x.vec;
        
        nonzero := [];
        
        for j in [1..Size(new_mat)] do 
            if ForAll(new_mat, x ->  x = 0) then 
                Add(nonzero, j);
            fi;
        od;
        
        new_mat := new_mat{nonzero};
        new_vec := new_mat{nonzero};
        
        mat := y.mat;
        vec := y.vec;
        
        unknowns := ShallowCopy(x.unknowns);
    fi;
    
    return rec( mat := mat,
                vec := vec,
                unknowns := unknowns    );
    
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
            evecs,
            row,
            sum,
            old_mat,
            x,
            y,
            z,
            w;
    
    dim := Size(ProductList.coords);;
    
    # Find unknown algebra products from eigenvectors
    
    x := MAJORANA_EigenvectorsAlgebraUnknowns(GramMatrix, AlgebraProducts, EigenVectors, ProductList);

    mat := x.mat;
    vec := x.vec;
    unknowns := x.unknowns;
    
    if Size(unknowns) = 0 then return; fi;
    
    # Find unknown algebra products from the resurrection principle
    
    Info(   InfoMajorana, 50, "Building resurrection");

    for i in ProductList.orbitreps do     
        
        u := [1..dim]*0;; u[i] := 1;;
    
        for evals in [[1,2],[2,1],[1,3],[2,3]] do  
        
            ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
            
            evecs := EigenVectors[i][evals[2]];
            
            for beta in evecs do
                for gamma in EigenVectors[i][evals[1]] do  
                
                    bad := MAJORANA_FindBadIndices(gamma, AlgebraProducts, ProductList);
                    
                    if ForAny(beta{bad}, x -> x <> 0) then 
                        
                        x := MAJORANA_SeparateAlgebraProduct(beta, gamma, unknowns, AlgebraProducts, ProductList);
                        
                        alpha_mat := [];
                                            
                        for alpha in EigenVectors[i][evals[1]] do                             
                            Add(alpha_mat, alpha - beta);                        
                        od;
                        
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
                                
                                if row <> [] and ForAny(row, x -> x <> 0) then 
                                    MAJORANA_Append([row,sum],mat,vec);
                                fi;
                            fi;
                        od;
                    fi; 
                    
                    if mat <> [] and Size(mat) > Size(mat[1]) then 
            
                        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
                                
                        mat := x.mat;
                        vec := x.vec;
                        unknowns := x.unknowns;
                        
                        if Size(unknowns) = 0 then return; fi;
                        
                    fi;
                                       
                od;
            od;                                               
        od;
    od;
    
    if mat <> [] and Size(mat) <= Size(mat[1]) then 

        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, AlgebraProducts, ProductList);
        
        if Size(x.unknowns) = 0 then return; fi;
        
        MAJORANA_AddConjugates(x.mat, x.vec, x.unknowns, AlgebraProducts, ProductList);
        
    fi;
           
        
    
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

    function( mat, vec, unknowns, AlgebraProducts, ProductList)
    
    local   sol,        # solution of system
            sign,       # correct sign of 5A axes
            i,          # loop over <unknowns>
            x,          # element of <unknowns>
            y,          # orbit of x
            g,          # conj element of x
            unsolved,
            prod,
            perm,
            nonzero,
            j;
    
    if mat = [] then
        return rec( mat := [], vec := [], unknowns := unknowns    );
    fi;
    
    Info(   InfoMajorana, 40, 
            STRINGIFY("Solving a ", Size(mat), " x ", Size(mat[1]), " matrix") );

    sol := MAJORANA_SolutionMatVecs(mat,vec);
    
    Info(   InfoMajorana, 40, "Solved it!" );
    
    for i in [1..Size(unknowns)] do
    
        if sol.solutions[i] <> fail then  
        
            x := unknowns[i]; 
            
            y := ProductList.pairorbit[x[1]][x[2]];
            g := ProductList.pairconj[x[1]][x[2]][2];
            
            if y > 0 then 
                sign := 1;
            else
                sign := -1;
                y := -y;
            fi;
            
            if AlgebraProducts[y] = false then 
                AlgebraProducts[y] := sign*MAJORANA_ConjugateVector(sol.solutions[i],g,ProductList);              
            fi;                   
        fi;                
    od;
    
    Unbind(sol.solutions);
    
    x := MAJORANA_RemoveKnownAlgProducts(   sol.mat,
                                            sol.vec,
                                            unknowns,
                                            AlgebraProducts,
                                            ProductList         );
                                                    
    nonzero := [];
                    
    for j in [1..Size(x.mat)] do              
        if ForAny(x.mat[j], x -> x <> 0) then 
            Add(nonzero,j);
        fi;
    od;
    
    x.mat := x.mat{nonzero};
    x.vec := x.vec{nonzero};
                                        
    return rec( mat := x.mat,
                vec := x.vec,
                unknowns := x.unknowns    );
    
    end );
    
InstallGlobalFunction( MAJORANA_RemoveKnownAlgProducts,
    
    # Takes a system [mat, vec] of unknown algebra products and removes 
    # from the system any variables which have already been found 
    
    function( mat, vec, unknowns, AlgebraProducts, ProductList)
    
    local   unsolved,
            i,
            j,
            x,
            y,
            sign,
            g,
            prod,
            nonzero;

    unsolved := [];
    
    for i in [1..Size(unknowns)] do 
    
        x := unknowns[i]; 
                    
        y := ProductList.pairorbit[x[1]][x[2]];
        
        if y > 0 then 
            sign := 1;
        else
            sign := -1;
            y := -y;
        fi;
        
        prod := AlgebraProducts[y];
                        
        if prod <> false then 
            
            g := ProductList.pairconj[x[1]][x[2]][1];
            
            prod := MAJORANA_ConjugateVector(prod,g,ProductList);
                
            for j in [1..Size(mat)] do 
                vec[j] := vec[j] - sign*mat[j][i]*prod;
            od; 
        else
            Add(unsolved,i);
        fi;
    od;
       
    mat := List(mat, x -> x{unsolved});
    
    return rec( mat := mat, 
                vec := vec, 
                unknowns := unknowns{unsolved});
        
    end );
    
InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, UnknownInnerProducts, GramMatrix)
    
    local   sol,    # solution of system
            i,      # loop over <UnknownInnerProducts>
            x;      # element of <UnknownInnerProducts>    
    
    sol := MAJORANA_SolutionMatVecs(mat,vec);                   
        
    for i in [1..Size(sol.solutions)] do
        if sol.solutions[i] <> fail then
            x := UnknownInnerProducts[i]; 
            GramMatrix[x] := sol.solutions[i][1];
        fi;
    od;
    
    return [true];
    
    end );
    
InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(GramMatrix,ProductList)
    
    local   gram,     # full gram matrix
            null;     # nullspace of gram matrix
    
        if ForAll(GramMatrix, x -> x <> false) then 
            gram := MAJORANA_FillGramMatrix(GramMatrix, ProductList);
            null := NullspaceMat(gram);; 
        fi;
        
        return null;
    
    end );
        
InstallGlobalFunction(MAJORANA_MoreEigenvectors,

    function(AlgebraProducts,EigenVectors,ProductList,nullspace)
    
    local   i,
            j,
            dim,
            list,
            a,
            b,
            d,
            mat,
            x,
            table,
            ev;
            
    dim := Size(ProductList.coords);
    
    table := [0,1/4,1/32];

    for i in ProductList.orbitreps do
                    
        a := [1..dim]*0; a[i] := 1;
    
        if IsMutable(EigenVectors[i][1]) then 
            
            if Size(BaseMat(Union(EigenVectors[i][1],EigenVectors[i][2],EigenVectors[i][3],nullspace))) < dim - 1 then
        
                mat := [];

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
                
                    list := List(mat[dim], x -> DenominatorRat(x));
                    d := Maximum(list);
                
                    for ev in [1..3] do
                     
                        Info(   InfoMajorana, 50, 
                                STRINGIFY( "Finding ", table[ev], " eigenvectors for axis ", i) ); 
                                
                        EigenVectors[i][ev] := NullspaceMat(d*(mat - IdentityMat(dim)*table[ev]));
                    od;
                fi;
            fi;
        fi; 
    od;
    
    return [true];
    
    end);

InstallGlobalFunction(MAJORANA_MainLoop,

    function(rep)
    
    local   x, 
            dim,
            maindimensions,
            j,
            k;
            
    dim := Size(rep.setup.coords);
    
                                ## STEP 5: INNER PRODUCTS M1 ##
                     
    MAJORANA_UnknownsAxiomM1(rep.innerproducts,rep.algebraproducts,rep.setup);
    
    if IsMutable(rep.innerproducts) and not false in rep.innerproducts then 
        rep.nullspace := MAJORANA_CheckNullSpace(rep.innerproducts, rep.setup);
        MakeImmutable(rep.innerproducts);
    fi;
    
                        ## STEP 8: RESURRECTION PRINCIPLE I ##
            
    MAJORANA_UnknownAlgebraProducts(rep.innerproducts,rep.algebraproducts,rep.evecs,rep.setup);
    
                                ## STEP 9: MORE EVECS II ##

    # Check if we have full espace decomp, if not find it

    x := MAJORANA_MoreEigenvectors(rep.algebraproducts,rep.evecs,rep.setup, rep.nullspace);
    
    ## STEP 6: FUSION ## 
    
    MAJORANA_Fusion(rep.innerproducts, rep.algebraproducts,rep.evecs,rep.setup);   
    
                        ## STEP 10: INNER PRODUCTS FROM ORTHOGONALITY ##
       
    # Use orthogonality of eigenspaces to write system of unknown variables for missing inner products
    
    MAJORANA_Orthogonality(rep.evecs,rep.innerproducts, rep.setup);
        
    if IsMutable(rep.innerproducts) and not false in rep.innerproducts then 
        rep.nullspace := MAJORANA_CheckNullSpace(rep.innerproducts, rep.setup);
        MakeImmutable(rep.innerproducts);
    fi;
    
    end);
    
InstallGlobalFunction(MajoranaRepresentation,

function(input,index)

    local   i,
            j,
            rep,
            dim,
            maindimensions,
            newdimensions,
            falsecount,
            newfalsecount,
            switchmain;    

    rep :=  MAJORANA_SetUp(input,index);
    
    if Size(rep.group) > 120 then 
        MAJORANA_AllEmbeddings(rep);
    fi;
    
    dim := Size(rep.setup.coords);
    
    maindimensions:=[];

    for i in rep.setup.orbitreps do
        for j in [1..3] do
            if Size(rep.evecs[i][j]) > 0 then
                rep.evecs[i][j]:=ShallowCopy(BaseMat(rep.evecs[i][j]));
            fi;
        od;
        Add(maindimensions,   Size(rep.evecs[i][1])
                            + Size(rep.evecs[i][2])
                            + Size(rep.evecs[i][3]) + 1);
    od;
    
    falsecount := [0,0];
    
    if false in rep.algebraproducts then
        falsecount[1] := Size(Positions(rep.algebraproducts,false));
    fi;
    
    if false in rep.innerproducts then
        falsecount[2] := Size(Positions(rep.innerproducts,false));
    fi;
    
    if ForAll(maindimensions, x -> x = dim) and falsecount = [0,0] then 
        switchmain := 1;
    else
        switchmain := 0;
    fi;
    
    while switchmain = 0 do
                                
        MAJORANA_MainLoop(rep);
        
        newdimensions := [];
        
        for i in rep.setup.orbitreps do 
            Add(newdimensions,   Size(rep.evecs[i][1])
                               + Size(rep.evecs[i][2])
                               + Size(rep.evecs[i][3]) + 1);
        od;
        
        newfalsecount := [0,0];
        
        if false in rep.algebraproducts then
            newfalsecount[1] := Size(Positions(rep.algebraproducts,false));
        fi;
        
        if false in rep.innerproducts then
            newfalsecount[2] := Size(Positions(rep.innerproducts,false));
        fi;

        Info(InfoMajorana, 20,
                STRINGIFY( "There are ", newfalsecount[1], " unknown algebra products ") );
        Info(InfoMajorana, 20,
                STRINGIFY( "There are ", newfalsecount[2], " unknown inner products ") );

        if newfalsecount = [0,0] then
            break;
        elif newdimensions = maindimensions and newfalsecount = falsecount then
            
            Info( InfoMajorana, 10, "Fail" );
            
            return rep;
            break;
        else
            maindimensions := StructuralCopy(newdimensions);
            falsecount := StructuralCopy(newfalsecount);
        fi;
    od;

    Info( InfoMajorana, 10, "Success" );
    
    return rep;

end );

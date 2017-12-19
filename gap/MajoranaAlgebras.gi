#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Implementations
#

# Creates list of indexes [i,j] where product of i th and j th coordinate vectors is not known

BindGlobal( "MAJORANA_ExtractUnknownAlgebraProducts",

function(algebraproducts, setup)

    local   unknowns,       # list of unknown algebra products
            i,              # loop over coordinates
            j,              # loop over coordinates
            k,              # pair orbit index
            dim;            # size of coordinates
    
    unknowns := [];
    dim := Size(setup.coords);
    
    for i in [1..dim] do
        for j in [i + 1..dim] do 
            
            k := setup.pairorbit[i][j];
        
            if k < 0 then 
                k := -k;
            fi;
        
            if algebraproducts[k] = false then 
                Add(unknowns,[i,j]);
            fi;
        od;
    od;

    return AsSet(unknowns);
end);

# Finds the indices i such that v_i*v is not known

InstallGlobalFunction(MAJORANA_FindBadIndices,
    
    function(v, algebraproducts, setup)
    
    local   i,
            j,
            k,
            dim,
            list,
            bad;            
    
    bad := [];
    dim := Size(setup.coords);
    list := [1..dim];
                    
    for i in v!.indices[1] do
        for j in list do 
            k :=  setup.pairorbit[i][j];
           
            if k > 0 then 
                if algebraproducts[k] = false then 
                    Add(bad,j);
                    list := Difference(list,[j]);
                fi;
            else
                if algebraproducts[-k] = false then 
                    Add(bad,j);
                    list := Difference(list,[j]);
                fi;
            fi;                
        od;
    od;

    Sort(bad);
    
    return bad;
    
    end );        
    
# given two eigenvectors, if possible, finds product and adds it to appropriate set of evecs

InstallGlobalFunction( MAJORANA_FuseEigenvectors,

    function(a, b, i, evals, other_mat, new, innerproducts, algebraproducts, setup)
    
    local   dim,
            u,
            test,
            new_ev,
            pos,
            x,
            y,
            z;
         
    dim := Size(setup.coords);
    u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
    
    new_ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    pos := Position(MAJORANA_FusionTable[1], new_ev) - 1 ;
    
    x := MAJORANA_AlgebraProduct(a,b,algebraproducts,setup);
                    
    if evals = [2,2] then 
        y := MAJORANA_InnerProduct(a,b,innerproducts,setup);
        
        if y <> false then 
            new[1] := UnionOfRows(new[1], x - (1/4)*u*y);
        fi;
    elif evals = [3,3] then 
        y := MAJORANA_InnerProduct(a,b,innerproducts,setup);
        z := MAJORANA_AlgebraProduct(u,x,algebraproducts, setup);
        
        if y <> false and z <> false then 
            new[2] := UnionOfRows(new[2], z - (1/32)*u*y);
            new[1] := UnionOfRows(new[1], x + (3/32)*u*y - 4*z);            
        elif y <> false then 
            UnionOfRows(other_mat, x - (1/32)*u*y);
        fi;  
    else
        new[pos] := UnionOfRows(new[pos],x);
    fi;
    
    end );

# finds new eigenvectors using the fusion rules 

InstallGlobalFunction( MAJORANA_Fusion,

function(innerproducts, algebraproducts, evecs, setup)

    local   i,
            j,
            a,
            b,
            dim,
            unknowns,
            new,
            u,
            ev_a,
            ev_b,
            new_ev,
            new_dim,
            pos,
            z,
            null,
            bad,
            other_mat;
    
    dim := Size(setup.coords);
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts, setup);
    
    for i in setup.orbitreps do 
    
        if true then 

            new := [ [], [], [] ];
            
            for j in [1..3] do 
                new[j] := SparseZeroMatrix(1, dim, Rationals);
            od;
            
            other_mat := SparseZeroMatrix(1, dim, Rationals);
        
            for ev_a in [1..3] do
                for ev_b in [ev_a..3] do 
                    for j in [1..Nrows(evecs[i][ev_a])] do
                        
                        a := CertainRows(evecs[i][ev_a], [j]);

                        bad := MAJORANA_FindBadIndices(a,algebraproducts,setup);
                        
                        if bad <> [] then  
                            null := KernelMat(CertainColumns(evecs[i][ev_b], bad)).relations;
                        else
                            null := SparseIdentityMatrix(Nrows(evecs[i][ev_b]));
                        fi;
                        
                        for j in [1..Nrows(null)] do 
                            
                            b := CertainRows(null, [j])*evecs[i][ev_b];
                            
                            if MAJORANA_AlgebraProduct(a,b,algebraproducts,setup) = false then
                                Error("Pause");
                            fi;
                            
                            MAJORANA_FuseEigenvectors(a, b, i, [ev_a, ev_b], other_mat, new, innerproducts, algebraproducts, setup);
                            
                        od;
                    od;
                od;
            od;
            
            if other_mat <> [] then 
                u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
            
                bad := MAJORANA_FindBadIndices(u,algebraproducts, setup );
                null := KernelMat(CertainColumns(other_mat, bad)).relations;
                
                for j in [1..Nrows(null)] do 
                    
                    b := CertainRows(null, [j])*other_mat;
                
                    z := MAJORANA_AlgebraProduct(u,b,algebraproducts, setup);
                    new[2] := UnionOfRows(new[2], z);
                    new[1] := UnionOfRows(new[1], b - 4*z);
                od;
            fi;
        
            for j in [1..3] do 
                evecs[i][j] := UnionOfRows(evecs[i][j], new[j]);
                if evecs[i][j] <> [] then 
                    evecs[i][j] := EchelonMatDestructive(evecs[i][j]).vectors;
                fi;
            od;
        fi;        
    od;
    
    return [true];    
    
    end );     
         
InstallGlobalFunction(MAJORANA_Append,

    function(x,mat,vec)

    local   lcm;
    
    if not x[1]!.entries[1] in mat!.entries then  
    
        lcm := _FoldList2(x[1]!.entries[1], DenominatorRat, LcmInt);
        
        x := x*lcm;
        
        UnionOfRows(mat,x[1]);
        UnionOfRows(vec,x[2]);
    fi;
    
    end); 
    
InstallGlobalFunction( MAJORANA_ConjugateVector, 

    function(mat,g,setup)
    
    local   i,
            j,
            nrows,
            ncols,
            indices,
            entries,
            res,            
            pos;
    
    if g <> () then 
        
        nrows := Nrows(mat);
        ncols := Ncols(mat);
        
        indices := IndicesOfSparseMatrix(mat);
        entries := EntriesOfSparseMatrix(mat);
        
        res := SparseZeroMatrix(nrows, ncols, Rationals);
        
        for i in [1..nrows] do 
            for j in [1..Size(indices[i])] do 
                
                pos := g[indices[i][j]];
        
                if pos < 0 then 
                    SetEntry(res, i, -pos, -entries[i][j]); 
                else
                    SetEntry(res, i, pos, entries[i][j]);
                fi;
            od;
        od;
        
        return res;
    else
        return mat;
    fi;
    
    end );

InstallGlobalFunction(  MAJORANA_AlgebraProduct,

        function(u,v,algebraproducts,list) # If all the relevant products are known, returns the algebra product of u and v. If not, returns 0

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

        dim := Nrows(u);
        
        if ForAll(u!.indices, x -> x = []) or ForAll(u!.indices, x -> x = []) then 
            return SparseZeroMatrix(1, dim, Rationals);
        fi;
        
        vec := SparseZeroMatrix(1, dim, Rationals);

        elts := [];
        vecs := [];

        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do
                
                k := list.pairorbit[u!.indices[1][i]][v!.indices[1][j]];
                
                if k > 0 then 
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;

                x := algebraproducts[k];
                
                if x <> false then
                    
                    g := list.pairconj[u!.indices[1][i]][v!.indices[1][j]][1];
                    
                    pos := Position(elts,g);
                    
                    if pos <> fail then 
                        vecs[pos] := vecs[pos] + sign*u!.entries[1][i]*v!.entries[1][j]*x;
                    else
                        Add(elts,g);
                        Add(vecs,sign*u!.entries[1][i]*v!.entries[1][j]*x);
                    fi;
                else
                    # cannot calculate product
                    return false;
                fi;
            od;
        od;
        
        for i in [1..Size(elts)] do 
            vec := vec + MAJORANA_ConjugateVector(vecs[i],elts[i], list);
        od;
                
        return vec;
        
        end

        );

InstallGlobalFunction(  MAJORANA_InnerProduct,

    function(u, v, innerproducts, setup) # If all the relevant products are known, returns the algebra product of u and v. If not, returns [0]

        local   i,              # loop over u 
                j,              # loop over v
                k,              # pair orbit index
                sign,           # correct for 5A axes
                sum,            # output value
                unknowns,
                x;
        
        sum := 0;

        for i in Reversed([1..Size(u!.indices[1])]) do
            for j in Reversed([1..Size(v!.indices[1])]) do
                k := setup.pairorbit[u!.indices[1][i]][v!.indices[1][j]];
                
                if k > 0 then 
                    sign := 1;
                else
                    sign := -1;
                    k := -k;
                fi;
                
                if innerproducts[k] <> false then
                    sum := sum + sign*u!.entries[1][i]*v!.entries[1][j]*innerproducts[k];
                else
                    return false;
                fi;
            od;
        od;
        
        return sum;
        
        end 
        
        );

InstallGlobalFunction(MAJORANA_FillGramMatrix,

function(innerproducts, setup)

    local   i,                  # loop over setup.coords
            j,                  # loop over setup.coords
            k,                  # pair orbit index
            dim,                # size of setup.coords
            GramMatrixFull;     # output matrix
    
    dim := Size(setup.coords);
    
    GramMatrixFull := SparseZeroMatrix(dim, dim, Rationals);
    
    for i in [1..dim] do 
        for j in [1..dim] do
            
            k := setup.pairorbit[i][j];
            
            if k > 0 then 
                SetEntry(GramMatrixFull, i, j, innerproducts[k]);
            else
                SetEntry(GramMatrixFull, i, j, -innerproducts[-k]);
            fi;
        od;
    od;

    return GramMatrixFull;

    end

    );
    
InstallGlobalFunction(MAJORANA_SeparateInnerProduct,

    function(u,v,UnknownInnerProducts,innerproducts,setup)

    local   row,            # record values of unknowns 
            sum,            # record values of knowns
            dim,            # size of coordinates
            i,              # index for dim of u
            j,              # index for dim of v
            m,              # orbit of i,j
            pos,            # position of m in unknowns
            sign;           # correct sign of 5A axes
            
    dim := Size(setup.coords);
            
    sum := SparseZeroMatrix(1, 1, Rationals);
    row := SparseZeroMatrix(1, Size(UnknownInnerProducts), Rationals);

    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do
            
            m := setup.pairorbit[u!.indices[1][i]][v!.indices[1][j]];
            
            if m > 0 then 
                sign := 1;
            else
                sign := -1;
                m := -m;
            fi;

            if innerproducts[m] <> false then
                AddToEntry(sum, 1, 1, - sign*u!.entries[1][i]*v!.entries[1][j]*innerproducts[m]);
            else
                pos := Position(UnknownInnerProducts,m);
                AddToEntry(row, 1, pos, sign*u!.entries[1][i]*v!.entries[1][j]);
            fi;
        od;
    od;

    return [row,sum];

    end );
    
InstallGlobalFunction(MAJORANA_Orthogonality,

    function(evecs,innerproducts, setup)
    
    local   i,          # loop over T
            j, 
            k,
            ev,         # loop over eigenvalues
            evecs_a,    #
            evecs_b,    #
            u,
            v,
            dim,
            x,          # res of orthogonality
            mat,        # matrix of unknown values
            vec,        # vector of known values   
            unknowns;     
    
    dim := Size(setup.coords);
    
    unknowns := Positions(innerproducts,false);
    
    if Size(unknowns) = 0 then
        return;
    fi; 
    
    mat := SparseZeroMatrix(1, Size(unknowns), Rationals);
    vec := SparseZeroMatrix(1, 1, Rationals);
    
    for i in setup.orbitreps do        
        for ev in Combinations([0..3],2) do  
            if ev[1] = 0 then 
                u := SparseMatrix(1, dim, [[i]], [[1]], Rationals); 
                evecs_a := u;  
            else
                evecs_a := evecs[i][ev[1]];
            fi;
            
            evecs_b := evecs[i][ev[2]];
                
            for j in [1..Nrows(evecs_a)] do
                for k in [1..Nrows(evecs_b)] do

                    u := CertainRows(evecs_a, [j]);
                    v := CertainRows(evecs_b, [k]);

                    x := MAJORANA_SeparateInnerProduct( u, v,
                                                        unknowns,
                                                        innerproducts,
                                                        setup);

                    if x[1]!.indices[1] <> [] then
                        MAJORANA_Append(x, mat, vec);
                    fi;
                    
                od;
            od;
        od;
    od;
    
    if Nrows(mat) > 1 then 
        MAJORANA_SolutionInnerProducts(mat,vec, unknowns, innerproducts);
    fi;        

    end );    
    
InstallGlobalFunction(MAJORANA_EigenvectorsAlgebraUnknowns,

function(innerproducts, algebraproducts, evecs, setup)

    local   i,          # loop over representatives
            j,
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
            list,
            pos,
            dim;        # size of setup.coords
    
    dim := Size(setup.coords);
    t := Size(evecs);
    
    table := [0, 1/4, 1/32];    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts,setup);
    
    mat := SparseZeroMatrix(1, Size(unknowns), Rationals);
    vec := SparseZeroMatrix(1, dim, Rationals);
    
    Info( InfoMajorana, 50, "Building eigenvector unknowns");
    
    for i in setup.orbitreps do 
        
        list := Filtered(unknowns, x -> i in x);
    
        if list <> [] then 
         
            for ev in [1..3] do 
                
                u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
                
                for j in [1..Nrows(evecs[i][ev])] do
                    
                    v := CertainRows(evecs[i][ev], [j]);
                    
                    x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,algebraproducts,setup);
                    
                    x[2] := x[2] + table[ev]*v;
                    
                    if Size(x[1]!.indices[1]) = 1 then 
                        y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns, 
                                                        algebraproducts,
                                                        setup);
                                                        
                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;
                                                        
                        if unknowns = [] then return; fi;
                        
                    elif x[1]!.indices <> [] then 
                        MAJORANA_Append(x, mat, vec);                    
                    fi;                
                od;
            od;
        fi;
    od;
    
    if mat = [] then return rec( mat := [], vec := [], unknowns := unknowns ); fi;   
        
    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
            
    return y;
    
    end);
    
InstallGlobalFunction(MAJORANA_UnknownsAxiomM1,

    function(innerproducts, algebraproducts, setup)
    
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
    
    dim := Size(setup.coords);
    unknowns := Positions(innerproducts, false);
    
    if unknowns = [] then 
        return;
    fi;
    
    mat := SparseZeroMatrix(1, Size(unknowns), Rationals);
    vec := SparseZeroMatrix(1, dim, Rationals);
    
    for i in [1..dim] do 
        
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
    
        for j in [1..Size(algebraproducts)] do 
            
            if algebraproducts[j] <> false then 
            
                pos := setup.pairreps[j];
                
                for k in [pos,Reversed(pos)] do
                
                    v := SparseMatrix(1, dim, [[k[1]]], [[1]], Rationals); 
                    w := SparseMatrix(1, dim, [[k[2]]], [[1]], Rationals); 
                
                    row := SparseZeroMatrix(1, Size(unknowns), Rationals);
                    sum := SparseZeroMatrix(1, 1, Rationals);
                
                    x := MAJORANA_SeparateInnerProduct(u, algebraproducts[j], unknowns, innerproducts, setup);
                    
                    row := row + x[1];
                    sum := sum + x[2];
                
                    y := MAJORANA_AlgebraProduct(u, v, algebraproducts, setup);
                    
                    if y <> false then 
                        z := MAJORANA_SeparateInnerProduct(y, w, unknowns, innerproducts, setup);
                        
                        row := row - z[1];
                        sum := sum - z[2];
                        
                        if row!.indices <> [] then 
                            MAJORANA_Append([row,sum], mat, vec);
                        fi;
                    fi;     
                od;
            fi;
        od;
    od;
    
    if mat <> [] then 
        MAJORANA_SolutionInnerProducts(mat,vec,unknowns,innerproducts);
    fi;
    
    end );

InstallGlobalFunction(MAJORANA_SeparateAlgebraProduct,

    function(u,v,unknowns,algebraproducts,setup)
    
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
    
    dim := Size(setup.coords);
    
    row := SparseZeroMatrix(1, Size(unknowns), Rationals);
    sum := SparseZeroMatrix(1, dim, Rationals);
    
    elts := [];
    vecs := [];
    
    for i in [1..Size(u!.indices[1])] do
        for j in [1..Size(v!.indices[1])] do
            
            k := setup.pairorbit[u!.indices[1][i]][v!.indices[1][j]];
            
            if k > 0 then 
                sign := 1;
            else
                sign := -1;
                k := -k;
            fi;
            
            x := algebraproducts[k];
            
            if x <> false then 
                                        
                g := setup.pairconj[u!.indices[1][i]][v!.indices[1][j]];
                
                pos := Position(elts,g);
                
                if pos <> fail then 
                    vecs[pos] := vecs[pos] - sign*u!.entries[1][i]*v!.entries[1][j]*x;
                else
                    Add(elts,g);
                    Add(vecs,- sign*u!.entries[1][i]*v!.entries[1][j]*x);
                fi;
            else
            
                if i < j then
                    l := [i,j];
                else
                    l := [j,i];
                fi;
                
                pos := Position(unknowns,l);
                AddToEntry(row, 1, pos, u!.entries[1][i]*v!.entries[1][j]); 
            fi;
        od;
    od;
    
    for i in [1..Size(elts)] do 
        sum := sum + MAJORANA_ConjugateVector(vecs[i],elts[i],setup);
    od;
       
    return [row,sum];
    
    end);
    
InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns, setup)
    
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
    
        len     := Ncols(row);
        output  := SparseZeroMatrix(1, len, Rationals);
        
        for i in [1..Size(row!.indices[1])] do
        
            j := unknowns[row!.indices[1][i]][1];
            k := unknowns[row!.indices[1][i]][2];
            
            x := [0,0];
            
            pos_1 := Position(setup.longcoords,(setup.coords[j])^g);
            pos_2 := Position(setup.longcoords,(setup.coords[k])^g);
            
            x[1] := setup.poslist[pos_1];
            x[2] := setup.poslist[pos_2];
            
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
            SetEntry(output, 1, pos, sign*row!.entries[1][i]);
        od;
    
        return output;
    else
        return row;
    fi;
    
    end);     
    
InstallGlobalFunction(MAJORANA_SingleSolutions,
    
    function(evals, innerproducts, algebraproducts, evecs, setup)
    
    local   dim,
            x,
            i, j, k, l, m,
            alpha, beta, gamma,
            unknowns,
            indices,
            orb,
            alpha_mat;

    dim := Size(setup.coords);
    
    for i in setup.orbitreps do 
        for j in [1..Nrows(evecs[i][evals[2]])] do 
            for k in [1..Nrows(evecs[i][evals[1]])] do 
                
                beta := CertainRows(evecs[i][evals[2]], [j]);
                gamma := CertainRows(evecs[i][evals[1]], [k]);
                
                indices := [];
                
                for l in beta!.indices[1] do 
                    for m in gamma!.indices[1] do 
                        orb := setup.pairorbit[l][m];
                    
                        if orb < 0 then orb := -orb; fi;
                    
                        if algebraproducts[orb] = false then 
                            Add(indices, [l,m]);
                        fi;
                        
                        if Size(indices) > 1 then 
                            break;
                        fi;                
                    od;
                    
                    if Size(indices) > 1 then 
                        break;
                    fi;    
                od;
                
                if Size(indices) = 1 then 
                
                    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts, setup);
                
                    alpha_mat := SparseZeroMatrix(1, dim, Rationals);
                    
                    for alpha in evecs[i][evals[1]] do                             
                        alpha_mat := UnionOfRows(alpha_mat, alpha - beta);                        
                    od;
                
                    MAJORANA_Resurrection(i, alpha_mat, beta, gamma, evals, [], [], unknowns, innerproducts, algebraproducts, setup, true);
                fi;
                
            od;
        od;
    od;
    
    end );

InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(innerproducts, algebraproducts, evecs, setup, nullspace, group)
    
    local   dim,
            i,
            j,
            k,
            u,
            evals,
            mat,
            vec,
            unknowns,
            alpha,
            beta,
            gamma,
            pos,
            alpha_mat,
            nonzero,
            list,
            x,
            y,
            v;
    
    dim := Size(setup.coords);;
    
    # Find unknown algebra products from eigenvectors
    
    x := MAJORANA_EigenvectorsAlgebraUnknowns(innerproducts, algebraproducts, evecs, setup);

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;
    
    if unknowns = [] then return; fi;
    
    # Find unknown algebra products from the resurrection principle
    
    Info(   InfoMajorana, 50, "Building resurrection");

    for evals in [[1,2],[2,1],[1,3],[2,3]] do
        MAJORANA_SingleSolutions(evals, innerproducts, algebraproducts, evecs, setup);
    od;
    
    x := MAJORANA_RemoveKnownAlgProducts(    mat, vec, unknowns, 
                                        algebraproducts,
                                        setup );
    nonzero := [];
                    
    for i in [1..Size(x.mat)] do              
        if x.mat!.indices[i] = [] then 
            Add(nonzero,i);
        fi;
    od;
    
    mat := CertainRows(x.mat, nonzero);
    vec := CertainRows(x.vec, nonzero);
    unknowns := x.unknowns;

    for i in setup.orbitreps do        
        for evals in [[1,2],[2,1],[1,3],[2,3]] do  
            for j in [1..Nrows(evecs[i][evals[2]])] do
            
                beta := CertainRows(evecs[i][evals[2]], [j]);
            
                alpha_mat := SparseZeroMatrix(1, dim, Rationals);
                                            
                for alpha in evecs[i][evals[1]] do                             
                    alpha_mat := UnionOfRows(alpha_mat, alpha - beta);                        
                od;
            
                for k in [1..Nrows(evecs[i][evals[1]])] do  
                        
                    gamma := CertainRows(evecs[i][evals[1]]);
                    
                    x := MAJORANA_Resurrection(  i, alpha_mat, beta, gamma,  
                                            evals, mat, vec, unknowns,
                                            innerproducts, 
                                            algebraproducts, 
                                            setup, false);
                    
                    mat := x.mat; vec := x.vec; unknowns := x.unknowns;
                    
                    if unknowns = [] then return; fi;
                    
                    if mat <> [] and (Nrows(mat) > 8000 or Nrows(mat) > Ncols(mat)) then 
            
                        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
                                
                        mat := x.mat; vec := x.vec; unknowns := x.unknowns;
                        
                        if Size(unknowns) = 0 then return; fi;
                        
                        if false then 
                        for evals in [[1,2],[2,1],[1,3],[2,3]] do
                            MAJORANA_SingleSolutions(evals, innerproducts, algebraproducts, evecs, setup);
                        od;
                        
                        y := MAJORANA_RemoveKnownAlgProducts(    mat, vec, unknowns, 
                                                            algebraproducts,
                                                            setup );
                        nonzero := [];
                                        
                        for i in [1..Size(y.mat)] do              
                            if ForAny(y.mat[i], z -> z <> 0) then 
                                Add(nonzero,i);
                            fi;
                        od;
                        
                        mat := y.mat{nonzero};
                        vec := y.vec{nonzero};
                        unknowns := y.unknowns;
                        fi;
                        
                    fi;
                                       
                od;
            od;                                               
        od;
    od;
    
    if mat <> [] then 

        x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
        
        mat := x.mat; vec := x.vec; unknowns := x.unknowns;
    fi;
    
    if unknowns = [] then return; fi;
    
    if Nrows(nullspace) > 1 then return; fi;
    
    MAJORANA_NullspaceUnknowns(mat, vec, unknowns, algebraproducts, setup, nullspace, group);

    end );
    
InstallGlobalFunction( MAJORANA_NullspaceUnknowns,

    function(mat, vec, unknowns, algebraproducts, setup, nullspace, group)
    
    local   i,j,
            list,
            u,
            v,
            x,
            y,
            dim;
    
    Info( InfoMajorana, 50, "Building nullspace unknowns" );
    
    dim := Size(setup.coords);
    
    x := MAJORANA_Orbits(group, dim, setup);

    for i in x.orbitreps do 
    
        list := Filtered(unknowns, x -> i in x);
    
        if list <> [] then 
        
            u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
            
            for j in [1..Nrows(nullspace)] do
                
                v := CertainRows(nullspace, [j]);
             
                x := MAJORANA_SeparateAlgebraProduct(u,v,unknowns,algebraproducts,setup);
                
                if Size(x[1]!.indices[1]) = 1 then 
                    
                    y := MAJORANA_SolveSingleSolution(  x, mat, vec, unknowns, 
                                                        algebraproducts,
                                                        setup);
                                                        
                    mat := y.mat; vec := y.vec; unknowns := y.unknowns;
                                                        
                    if unknowns = [] then return; fi;
                    
                elif x[1].indices[1] = [] then 
                    MAJORANA_Append(x, mat, vec);                    
                fi;               
            od;
        fi;
    od;

    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
    
    end );
    
InstallGlobalFunction( MAJORANA_Resurrection,

    function(i, alpha_mat, beta, gamma, evals, mat, vec, unknowns, innerproducts, algebraproducts, setup, single)
    
    local   bad,
            dim,
            null,
            ev,
            u,
            j,
            n,
            nsum,
            row, 
            sum, 
            v,
            x,
            y,
            z,
            w,
            g,
            conj;
            
    bad := MAJORANA_FindBadIndices(gamma, algebraproducts, setup);
    dim := Size(setup.coords);
                    
    if Intersection(bad, beta!.indices[1]) <> [] then 
                        
        x := MAJORANA_SeparateAlgebraProduct(beta, gamma, unknowns, algebraproducts, setup);
                            
        null := KernelMat(CertainColumns(alpha_mat, bad)).relations;
                
        ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        for j in [1..Nrows(null)] do 
            
            n := CertainRows(null, [j]);
            
            nsum := Sum(n!.entries[1]);
        
            if nsum <> 0 then
            
                row := SparseZeroMatrix(1, Size(unknowns), Rationals);
                sum := SparseZeroMatrix(1, dim, Rationals);
                
                v := n*alpha_mat;
                
                y := MAJORANA_AlgebraProduct(v, gamma, algebraproducts, setup);
                
                z := MAJORANA_SeparateAlgebraProduct(u, y, unknowns, algebraproducts, setup);
            
                row := row + z[1] + nsum*ev*x[1];
                sum := sum + z[2] + nsum*ev*x[2];
                
                if evals[1] = 2 then 
                    w := MAJORANA_InnerProduct(n*alpha_mat, gamma, innerproducts, setup);
                    
                    if w <> false then 
                        sum := sum + (1/4)*w*u;
                    else
                        row := [];
                    fi;
                fi;
                
                if Size(row!.indices[1]) = 1 then 
                    
                    y := MAJORANA_SolveSingleSolution(  [row,sum], mat, vec, unknowns, 
                                                        algebraproducts,
                                                        setup);
                                                        
                    mat := y.mat; vec := y.vec; unknowns := y.unknowns;
                    
                    if unknowns = [] then 
                        return rec( mat := SparseZeroMatrix(1, 1, Rationals),
                                    vec := SparseZeroMatrix(1, 1, Rationals),
                                    unknowns := []  );
                    fi;
                    
                    x := MAJORANA_SeparateAlgebraProduct(beta, gamma, unknowns, algebraproducts, setup);
                    
                    if x[1]!.indices[1] = [] then 
                        return rec( mat := mat, vec := vec, unknowns := unknowns);
                    fi;
                    
                elif row <> [] and row!.indices[1] <> [] and not single then 
                    for g in setup.conjelts do

                        conj := [,];
                        
                        conj[1] := MAJORANA_ConjugateRow(   row, g[1],
                                                            unknowns,
                                                            setup );
                                                            
                        conj[2] := MAJORANA_ConjugateVector(    sum,g[2],
                                                                setup );
                        MAJORANA_Append(conj, mat, vec);
                    od;
                fi;
            fi;
        od;
    fi;
    
    return rec( mat := mat, vec := vec, unknowns := unknowns);
    
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

    function( mat, vec, unknowns, algebraproducts, setup)
    
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
    
    if ForAll(mat!.indices, x -> x = []) then
        return rec( mat := [], vec := [], unknowns := unknowns    );
    fi;
    
    Info(   InfoMajorana, 40, 
            STRINGIFY("Solving a ", Nrows(mat), " x ", Ncols(mat), " matrix") );

    sol := MAJORANA_SolutionMatVecs(mat,vec);
    
    Info(   InfoMajorana, 40, "Solved it!" );
    
    for i in [1..Size(unknowns)] do
    
        if sol.solutions[i] <> fail then  
        
            x := unknowns[i]; 
            
            MAJORANA_RecordSolution(    sol.solutions[i], x,
                                        algebraproducts,
                                        setup );
        fi;
    od;
    
    Unbind(sol.solutions);
    
    x := MAJORANA_RemoveKnownAlgProducts(   sol.mat,
                                            sol.vec,
                                            unknowns,
                                            algebraproducts,
                                            setup         );
                                                    
    nonzero := [];
                    
    for j in [1..Nrows(x.mat)] do              
        if x.mat!.indices[j] <> [] then 
            Add(nonzero,j);
        fi;
    od;
    
    x.mat := CertainRows(x.mat, nonzero);
    x.vec := CertainRows(x.vec, nonzero);
                                        
    return rec( mat := x.mat,
                vec := x.vec,
                unknowns := x.unknowns    );
    
    end );
    
InstallGlobalFunction( MAJORANA_SolveSingleSolution,

    function(x, mat, vec, unknowns, algebraproducts, setup) 
    
    local   elm, 
            y,
            nonzero,
            i;
            
    Info( InfoMajorana, 60, "Solved a single solution");
    
    elm := x[1]!.entries[1][1]; 
    x := x/elm;
    
    MAJORANA_RecordSolution(    x[2], unknowns[x[1]!.indices[1][1]],
                                algebraproducts,
                                setup );
    
    y := MAJORANA_RemoveKnownAlgProducts(    mat, vec, unknowns, 
                                        algebraproducts,
                                        setup );
    nonzero := [];
                    
    for i in [1..Nrows(y.mat)] do              
        if y.mat!.indices[i] = [] then 
            Add(nonzero,i);
        fi;
    od;
    
    y.mat := CertainRows(y.mat, nonzero);
    y.vec := CertainRows(y.vec, nonzero);
                                        
    return rec( mat := y.mat,
                vec := y.vec,
                unknowns := y.unknowns    );
    
    end );
    
    
InstallGlobalFunction( MAJORANA_RecordSolution,

    function( v, x, algebraproducts, setup)
    
    local   y,
            g,
            sign;
    
    y := setup.pairorbit[x[1]][x[2]];
    g := setup.pairconj[x[1]][x[2]][2];
    
    if y > 0 then 
        sign := 1;
    else
        sign := -1;
        y := -y;
    fi;
    
    if algebraproducts[y] = false then 
        algebraproducts[y] := sign*MAJORANA_ConjugateVector(v,g,setup);              
    fi;      
    
    end );
    
InstallGlobalFunction( MAJORANA_RemoveKnownAlgProducts,
    
    # Takes a system [mat, vec] of unknown algebra products and removes 
    # from the system any variables which have already been found 
    
    function( mat, vec, unknowns, algebraproducts, setup)
    
    local   unsolved,
            i,
            j,
            elm,
            x,
            y,
            sign,
            g,
            new,
            prod;

    unsolved := [];
    
    for i in [1..Size(unknowns)] do 
    
        x := unknowns[i]; 
                    
        y := setup.pairorbit[x[1]][x[2]];
        
        if y > 0 then 
            sign := 1;
        else
            sign := -1;
            y := -y;
        fi;
        
        prod := algebraproducts[y];
                        
        if prod <> false then 
            
            g := setup.pairconj[x[1]][x[2]][1];
            
            prod := MAJORANA_ConjugateVector(prod,g,setup);
            
            new := SparseZeroMatrix(Nrows(vec), Ncols(vec), Rationals);
            
            for j in [1..Nrows(vec)] do 
                elm := GetEntry(mat, j, i);
                if elm <> 0 then
                    vec!.entries := -sign*elm*prod!.entries;
                fi;
            od;
            
            vec := vec + new;
            
        else
            Add(unsolved,i);
        fi;
    od;
       
    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};
    
    for i in [1..Nrows(mat)] do
        if Size(mat!.indices[i]) = 1 then
            y := MAJORANA_SolveSingleSolution(   [mat[i],vec[i]], 
                                            mat, vec, 
                                            unknowns, 
                                            algebraproducts, 
                                            setup);
                                            
            mat := y.mat; vec := y.vec; unknowns := y.unknowns;
        fi;
    od;
    
    return rec( mat := mat, 
                vec := vec, 
                unknowns := unknowns);
        
    end );
    
InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, UnknownInnerProducts, innerproducts)
    
    local   sol,    # solution of system
            i,      # loop over <UnknownInnerProducts>
            x;      # element of <UnknownInnerProducts>    
    
    sol := MAJORANA_SolutionMatVecs(mat,vec);                   
        
    for i in [1..Size(sol.solutions)] do
        if sol.solutions[i] <> fail then
            x := UnknownInnerProducts[i]; 
            innerproducts[x] := sol.solutions[i]!.entries[1][1];
        fi;
    od;
    
    return [true];
    
    end );
    
InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(innerproducts,setup)
    
    local   gram,     # full gram matrix
            null;     # nullspace of gram matrix
    
        if ForAll(innerproducts, x -> x <> false) then 
            gram := MAJORANA_FillGramMatrix(innerproducts, setup);
            null := KernelMatDestructive(gram).relations;; 
        fi;
        
        return null;
    
    end );
        
InstallGlobalFunction(MAJORANA_MoreEigenvectors,

    function(algebraproducts,evecs,setup,nullspace)
    
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
            
    dim := Size(setup.coords);
    
    table := [0,1/4,1/32];

    for i in setup.orbitreps do
        
        if true then 
        
            a := [1..dim]*0; a[i] := 1;
        
            mat := [];

            for j in [1..dim] do
            
                b := [1..dim]*0; b[j] := 1;
                
                x := MAJORANA_AlgebraProduct(a,b,algebraproducts,setup);
                
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
                            
                    evecs[i][ev] := NullspaceMat(d*(mat - IdentityMat(dim)*table[ev]));
                od;
            fi;
        fi;
    od;
    
    end);

InstallGlobalFunction(MAJORANA_MainLoop,

    function(rep)
    
    local   x, 
            dim,
            maindimensions,
            i,
            evecs;
            
    dim := Size(rep.setup.coords);
    
                                ## STEP 5: INNER PRODUCTS M1 ##
                                
    MAJORANA_UnknownsAxiomM1(rep.innerproducts,rep.algebraproducts,rep.setup);
    
    if IsMutable(rep.innerproducts) and not false in rep.innerproducts then 
        rep.nullspace := MAJORANA_CheckNullSpace(rep.innerproducts, rep.setup);
        MakeImmutable(rep.innerproducts);
    fi;
                                
                                    ## STEP 6: FUSION ## 
    
    
    
                        ## STEP 8: RESURRECTION PRINCIPLE I ##
            
    MAJORANA_UnknownAlgebraProducts(rep.innerproducts,rep.algebraproducts,rep.evecs,rep.setup, rep.nullspace, rep.group);
    
                                ## STEP 9: MORE EVECS II ##

    # Check if we have full espace decomp, if not find it

    MAJORANA_MoreEigenvectors(rep.algebraproducts,rep.evecs,rep.setup, rep.nullspace);
    
    MAJORANA_Fusion(rep.innerproducts, rep.algebraproducts,rep.evecs,rep.setup); 
    
    for i in rep.setup.orbitreps do 
        
        evecs := Union( rep.evecs[i][1], rep.evecs[i][2], rep.evecs[i][3], rep.nullspace);
    
        if IsMutable(rep.evecs[i][1]) and Size( BaseMat( evecs )) = dim - 1 then 
            MakeImmutable(rep.evecs[i][1]);
        fi;
    od;
    
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
            falsecount,
            newfalsecount;  

    rep :=  MAJORANA_SetUp(input,index);
    
    if Size(rep.group) > 120 then 
        MAJORANA_AllEmbeddings(rep);
    fi;
    
    falsecount := [0,0];
    
    falsecount[1] := Size(Positions(rep.algebraproducts,false));
    falsecount[2] := Size(Positions(rep.innerproducts,false));
    
    while true do
                                
        MAJORANA_MainLoop(rep);
        
        newfalsecount := [0,0];

        newfalsecount[1] := Size(Positions(rep.algebraproducts,false));
        newfalsecount[2] := Size(Positions(rep.innerproducts,false));

        Info(InfoMajorana, 20,
            STRINGIFY( "There are ", newfalsecount[1], " unknown algebra products ") );
        Info(InfoMajorana, 20,
            STRINGIFY( "There are ", newfalsecount[2], " unknown inner products ") );

        if newfalsecount = [0,0] then
            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif newfalsecount = falsecount then
            Info( InfoMajorana, 10, "Fail" );
            return rep;
        else
            falsecount := StructuralCopy(newfalsecount);
        fi;
    od;
    
    end );

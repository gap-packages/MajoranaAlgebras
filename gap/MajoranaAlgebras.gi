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
        for j in [i..dim] do 
            
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
    
    if x <> false then
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
            elif other_mat <> false and y <> false then 
                other_mat := UnionOfRows(other_mat, x - (1/32)*u*y);
            fi;  
        else
            new[pos] := UnionOfRows(new[pos],x);
        fi;
    fi;
    
    end );

# finds new eigenvectors using the fusion rules 

InstallGlobalFunction( MAJORANA_Fusion,

function(rep)

    local   i,
            j,
            k, l,
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
    
    dim := Size(rep.setup.coords);
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(rep.algebraproducts, rep.setup);
    
    for i in rep.setup.orbitreps do 
    
        if IsMutable(rep.evecs[i]) then 

            new := [0,0,0];
            
            for j in [1..3] do 
                new[j] := CopyMat(rep.evecs[i][j]);
            od;
            
            other_mat := SparseMatrix(0, dim, [], [], Rationals);
        
            for ev_a in [1..3] do
                for ev_b in [ev_a..3] do 
                    for j in [1..Nrows(rep.evecs[i][ev_a])] do
                        
                        a := CertainRows(rep.evecs[i][ev_a], [j]);

                        bad := MAJORANA_FindBadIndices(a,rep.algebraproducts,rep.setup);
                        
                        if bad <> [] then  
                            null := KernelMat(CertainColumns(rep.evecs[i][ev_b], bad)).relations;
                        else
                            null := SparseIdentityMatrix(Nrows(rep.evecs[i][ev_b]));
                        fi;
                        
                        for k in [1..Nrows(null)] do 
                            
                            b := CertainRows(null, [k])*rep.evecs[i][ev_b];
                            
                            MAJORANA_FuseEigenvectors(a, b, i, [ev_a, ev_b], other_mat, new, rep.innerproducts, rep.algebraproducts, rep.setup);
                            
                            
                        od;                        
                    od;
                
                    for k in [1..3] do 
                        if Nrows(new[k]) > dim then                             
                            new[k] := EchelonMatDestructive(new[k]).vectors;
                        fi;
                    od;
                    
                    if Sum(List(new, x -> Nrows(x))) + Nrows(rep.nullspace) >= dim - 1 then 
                        break;
                    fi;
                od;
                
                if Sum(List(new, x -> Nrows(x))) + Nrows(rep.nullspace) >= dim - 1 then 
                    break;
                fi;
                                
            od;
            
            if Nrows(other_mat) > 0 then 
                u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
            
                bad := MAJORANA_FindBadIndices(u,rep.algebraproducts, rep.setup );
                null := KernelMat(CertainColumns(other_mat, bad)).relations;
                
                for j in [1..Nrows(null)] do 
                    
                    b := CertainRows(null, [j])*other_mat;
                
                    z := MAJORANA_AlgebraProduct(u,b,rep.algebraproducts, rep.setup);
                    new[2] := UnionOfRows(new[2], z);
                    new[1] := UnionOfRows(new[1], b - 4*z);
                od;
            fi;
        
            if ForAll(new, x -> x!.indices = []) then return; fi;
        
            for j in [1..3] do 
                rep.evecs[i][j] := new[j];
                rep.evecs[i][j] := MAJORANA_BasisOfEvecs(rep.evecs[i][j]);
            od;
        fi;        
    od;
    
    end );     
    
InstallGlobalFunction( MAJORANA_ConjugateVec, 

    function(mat,g,setup)
    
    local   i,
            j,
            nrows,
            ncols,
            indices,
            entries,
            res,            
            pos;
    
    if g <> [] then 
        
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

        dim := Ncols(u);
        
        if u!.indices[1] = [] or v!.indices[1] = [] then 
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
                    
                    g := list.pairconj[u!.indices[1][i]][v!.indices[1][j]];
                    
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
            x := MAJORANA_ConjugateVec(vecs[i],list.pairconjelts[elts[i]],list);
            AddRow(x!.indices[1],x!.entries[1],vec!.indices,vec!.entries,1);
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

    function(u,v,unknowns,innerproducts,setup)

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
    row := SparseZeroMatrix(1, Size(unknowns), Rationals);

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
                pos := Position(unknowns,m);
                AddToEntry(row, 1, pos, sign*u!.entries[1][i]*v!.entries[1][j]);
            fi;
        od;
    od;

    return [row,sum];

    end );
    
InstallGlobalFunction(MAJORANA_Orthogonality,

    function(rep)
    
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
    
    if not false in rep.innerproducts then 
        return;
    fi;
    
    dim := Size(rep.setup.coords);
    
    unknowns := Positions(rep.innerproducts,false);
    
    if Size(unknowns) = 0 then
        return;
    fi; 
    
    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);
    
    for i in rep.setup.orbitreps do        
        for ev in Combinations([0..3],2) do  
            if ev[1] = 0 then 
                u := SparseMatrix(1, dim, [[i]], [[1]], Rationals); 
                evecs_a := u;  
            else
                evecs_a := rep.evecs[i][ev[1]];
            fi;
            
            evecs_b := rep.evecs[i][ev[2]];
                
            for j in [1..Nrows(evecs_a)] do
                for k in [1..Nrows(evecs_b)] do

                    u := CertainRows(evecs_a, [j]);
                    v := CertainRows(evecs_b, [k]);

                    x := MAJORANA_SeparateInnerProduct( u, v,
                                                        unknowns,
                                                        rep.innerproducts,
                                                        rep.setup);

                    if x[1]!.indices[1] <> [] then
                        if not _IsRowOfSparseMatrix(mat, x[1]) then 
                            mat := UnionOfRows(mat, x[1]);
                            vec := UnionOfRows(vec, x[2]);
                        fi;
                    fi;
                    
                od;
            od;
        od;
    od;
    
    if Nrows(mat) > 1 then 
        MAJORANA_SolutionInnerProducts(mat,vec, unknowns, rep.innerproducts);
    fi;      
    
    if not false in rep.innerproducts then 
        rep.nullspace := MAJORANA_CheckNullSpace(rep.innerproducts, rep.setup);
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
            conj,
            list,
            pos,
            dim;        # size of setup.coords
    
    dim := Size(setup.coords);
    t := Size(evecs);
    
    table := [0, 1/4, 1/32];    
    unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts,setup);
    
    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, dim, [], [], Rationals);
    
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
                        
                    elif x[1]!.indices[1] <> [] and not _IsRowOfSparseMatrix(mat, x[1]) then
                        mat := UnionOfRows(mat, x[1]);
                        vec := UnionOfRows(vec, x[2]);               
                    fi;                
                od;
            od;
        fi;
    od;

    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
            
    return y;
    
    end);
    
InstallGlobalFunction(MAJORANA_AxiomM1,

    function(rep)
    
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
            
    if not false in rep.innerproducts then 
        return;
    fi;
    
    dim := Size(rep.setup.coords);
    unknowns := Positions(rep.innerproducts, false);
    
    mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
    vec := SparseMatrix(0, 1, [], [], Rationals);
    
    for i in [1..dim] do 
        
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
    
        for j in [1..Size(rep.algebraproducts)] do 
            
            if rep.algebraproducts[j] <> false then 
            
                pos := rep.setup.pairreps[j];
                
                for k in [pos,Reversed(pos)] do
                
                    v := SparseMatrix(1, dim, [[k[1]]], [[1]], Rationals); 
                    w := SparseMatrix(1, dim, [[k[2]]], [[1]], Rationals); 
                
                    row := SparseZeroMatrix(1, Size(unknowns), Rationals);;
                    sum := SparseZeroMatrix(1, 1, Rationals);
                
                    x := MAJORANA_SeparateInnerProduct(u, rep.algebraproducts[j], unknowns, rep.innerproducts, rep.setup);
                    
                    row := row + x[1];
                    sum := sum + x[2];
                
                    y := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
                    
                    if y <> false then 
                        z := MAJORANA_SeparateInnerProduct(y, w, unknowns, rep.innerproducts, rep.setup);
                        
                        row := row - z[1];
                        sum := sum - z[2];
                        
                        if row!.indices[1] <> [] then 
                            if not _IsRowOfSparseMatrix(mat, row) then
                                mat := UnionOfRows(mat, row);
                                vec := UnionOfRows(vec, sum);
                            fi;
                        fi;
                    fi;     
                od;
            fi;
        od;
        
        if Nrows(mat) > Ncols(mat) then 
            x := MAJORANA_SolutionInnerProducts(mat, vec, unknowns, rep.innerproducts);
            
            mat := x.mat; vec := x.vec; unknowns := x.unknowns;
            
            if unknowns = [] then return; fi;
        fi;
        
    od;
    
    MAJORANA_SolutionInnerProducts(mat,vec,unknowns,rep.innerproducts);
    
    if not false in rep.innerproducts then 
        rep.nullspace := MAJORANA_CheckNullSpace(rep.innerproducts, rep.setup);
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

                l := [u!.indices[1][i], v!.indices[1][j]]; 
                Sort(l);
                
                pos := Position(unknowns,l);
                AddToEntry(row, 1, pos, u!.entries[1][i]*v!.entries[1][j]); 
            fi;
        od;
    od;
    
    for i in [1..Size(elts)] do 
        sum := sum + MAJORANA_ConjugateVec(vecs[i],setup.pairconjelts[elts[i]],setup);
    od;
       
    return [row,sum];
    
    end);
    
InstallGlobalFunction(MAJORANA_ConjugateRow,

    function(row, g, unknowns, setup)
    
    local   output,     # output row
            len,        # length of row
            i,          # loop over length of row
            x,y,
            sign,       # corrects sign of 5A axis
            pos;        # position of new product
    
    if g <> [] then 
    
        len     := Ncols(row);
        output  := SparseZeroMatrix(1, len, Rationals);
        
        for i in [1..Size(row!.indices[1])] do
        
            x := unknowns[row!.indices[1][i]];
            y := g{x};
            
            sign := 1;
            
            if y[1] < 0 then sign := -sign; y[1] := -y[1]; fi;
            if y[2] < 0 then sign := -sign; y[2] := -y[2]; fi;
            
            Sort(y);
            
            pos := Position(unknowns,y);
            SetEntry(output, 1, pos, sign*row!.entries[1][i]);
        od;
    
        return output;
    else
        return row;
    fi;
    
    end);     
    
InstallGlobalFunction(MAJORANA_BasisOfEvecs,

    function(mat)
    
    local ech, dim;
    
    dim := Ncols(mat);
    
    ech := EchelonMatTransformation(CertainColumns(mat, [dim, dim - 1..1]));
    
    return ech.coeffs*mat;
    
    end);
    
    
InstallGlobalFunction(MAJORANA_UnknownAlgebraProducts,

    function(rep)
    
    local   dim, t, x, y, i, j, k, l, evals, mat, vec, unknowns, u, a, b, c, bad, null, g, conj, list, evecs_a, evecs_b, index, new_mat, new_vec, nonzero; 
            
    t := Size(rep.involutions);
    dim := Size(rep.setup.coords);
    
    x := MAJORANA_EigenvectorsAlgebraUnknowns(rep.innerproducts, rep.algebraproducts, rep.evecs, rep.setup);

    mat := x.mat; vec := x.vec; unknowns := x.unknowns;
    
    if unknowns = [] then return; fi;
    
    if Nrows(rep.nullspace) > 0 then 
    
        x := MAJORANA_NullspaceUnknowns(    mat, vec, unknowns, 
                                            rep.algebraproducts, 
                                            rep.setup, rep.nullspace, 
                                            rep.group);
    
        mat := x.mat; vec := x.vec; unknowns := x.unknowns;
    
        if unknowns = [] then return; fi;
    fi;
    
    Info(   InfoMajorana, 50, "Building resurrection");
    
    for evals in [[1,2],[2,1],[1,3],[2,3]] do     
        for i in rep.setup.orbitreps do 
        
            evecs_a := rep.evecs[i][evals[1]];
            evecs_b := rep.evecs[i][evals[2]];
        
            u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
            
            list := [,];
            
            list[1] := List([1..Nrows(evecs_a)], i -> []);
            list[2] := List([1..Nrows(evecs_a)], i -> []);
             
            for j in [1..Nrows(evecs_a)] do 
                bad := MAJORANA_FindBadIndices(CertainRows(evecs_a,[j]), rep.algebraproducts, rep.setup);
                for k in [1..Nrows(evecs_b)] do 
                    x := Size(Intersection(bad, evecs_b!.indices[k]));
                    if x = 1 then 
                        Add(list[1][j], k);
                    elif x > 1 then 
                        Add(list[2][j], k);
                    fi;
                od;
            od;
            
            for index in [1,2] do 
                for j in [1..Nrows(evecs_a)] do 
                    c := CertainRows(evecs_a, [j]);
                    for k in list[index][j] do 
                        b := CertainRows(evecs_b, [k]);
                        
                        bad := MAJORANA_FindBadIndices(c, rep.algebraproducts, rep.setup);
                        
                        for l in [1..Nrows(evecs_a)] do 
                        
                            a := CertainRows(evecs_a, [l]);
                            
                            if CertainColumns(a, bad) = CertainColumns(b, bad) then 
                            
                                x := MAJORANA_Resurrection(  u, a, b, c, evals, 
                                                        unknowns,
                                                        rep.innerproducts,
                                                        rep.algebraproducts,
                                                        rep.setup);
                                                        
                                # Error("pause");
                                
                                if x <> false and x[1]!.indices[1] <> [] then 
                                    if Size(x[1]!.indices[1]) = 1 then 
                                    
                                        y := MAJORANA_SolveSingleSolution( x, 
                                                            mat, vec, unknowns, 
                                                            rep.algebraproducts,
                                                            rep.setup);
                                        
                                        mat := y.mat; vec := y.vec; unknowns := y.unknowns;
                                        
                                        if unknowns = [] then return; fi;
                                    elif not _IsRowOfSparseMatrix(mat, x[1]) then
                                        mat := UnionOfRows(mat, x[1]);
                                        vec := UnionOfRows(vec, x[2]);
                                    fi;
                                fi;                                
                            fi;
                        od;
                        
                        if Nrows(mat) > Ncols(mat) then 
                            x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);

                            mat := x.mat; vec := x.vec; unknowns := x.unknowns;
                            
                            if unknowns = [] then return; fi;                            
                        fi;
                        
                    od;
                od;
            od;
        od;
    od;
    
    x := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, rep.algebraproducts, rep.setup);
    
    mat := x.mat; vec := x.vec; unknowns := x.unknowns;
                            
    if unknowns = [] then return; fi;
    
    Info(   InfoMajorana, 50, "All conjugates") ;
    
    new_mat := CopyMat(mat);
    new_vec := CopyMat(vec);
    
    for i in [1..Nrows(mat)] do 
        if mat!.indices[i] <> [] then 
            for g in rep.setup.conjelts do
                conj := [,];
                
                conj[1] := MAJORANA_ConjugateRow(   CertainRows(mat, [i]), 
                                                    g, unknowns,
                                                    rep.setup );
                                                    
                conj[2] := MAJORANA_ConjugateVec(   CertainRows(vec, [i]), 
                                                    g, rep.setup );
                                                    
                new_mat := UnionOfRows(new_mat, conj[1]);
                new_vec := UnionOfRows(new_vec, conj[2]);
            od;
            
            if Nrows(new_mat) > Ncols(new_mat)/2 or Nrows(new_mat) > 8000 then 
                x := MAJORANA_SolutionAlgProducts(new_mat, new_vec, unknowns, rep.algebraproducts, rep.setup);
                
                if x.unknowns = [] then 
                    return;
                fi;
                
                nonzero := [];
    
                for j in [1..Nrows(x.mat)] do              
                    if x.mat!.indices[j] <> [] then 
                        Add(nonzero,j);
                    fi;
                od;
                
                new_mat := CertainRows(x.mat, nonzero);
                new_vec := CertainRows(x.vec, nonzero);
                
                y := MAJORANA_RemoveKnownAlgProducts(mat, vec, unknowns, rep.algebraproducts, rep.setup);
                
                mat := y.mat; vec := y.vec; unknowns := y.unknowns;
                
            fi;
        fi;
    od;
    
    MAJORANA_SolutionAlgProducts(new_mat,new_vec,unknowns, rep.algebraproducts, rep.setup);
    
    end );

InstallGlobalFunction(MAJORANA_Resurrection, 

    function(u, a, b, c, evals, unknowns, innerproducts, algebraproducts, setup)
    
    local   x, y, ev, res; 
    
    res := MAJORANA_SeparateAlgebraProduct(b, c, unknowns, algebraproducts, setup);
    
    if res[1]!.indices[1] = [] then 
        return false;
    fi;
    
    ev := MAJORANA_FusionTable[evals[1] + 1][evals[2] + 1];
    
    res := ev*res;
    
    x := MAJORANA_AlgebraProduct(c, a - b, algebraproducts, setup);
    
    res := res + MAJORANA_SeparateAlgebraProduct(u, x, unknowns, algebraproducts, setup);
    
    if evals[1] = 2 then 
        y := MAJORANA_InnerProduct(a, c, innerproducts, setup);
        
        if y <> false then 
            res[2] := res[2] + (1/4)*y*u;
        else
            return false;
        fi;
    fi;
    
    return res;
    
    end );

InstallGlobalFunction( MAJORANA_NullspaceUnknowns,

    function(mat, vec, unknowns, algebraproducts, setup, nullspace, group)
    
    local   i,j, gens,
            list,
            u,
            v,
            x,
            y,
            dim;
    
    Info( InfoMajorana, 50, "Building nullspace unknowns" );
    
    dim := Size(setup.coords);
    
    gens := GeneratorsOfGroup(group);
    gens := List(gens, x -> MAJORANA_FindVectorPermutation(x, setup));
    
    x := MAJORANA_Orbits(gens, dim, setup);

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
                                                        
                    if unknowns = [] then 
                        return rec( mat := SparseMatrix(0, 0, [], [], Rationals),
                                    vec := SparseMatrix(0, 0, [], [], Rationals),
                                    unknowns := []); 
                    fi;
                    
                elif x[1]!.indices[1] <> [] then 
                    if not _IsRowOfSparseMatrix(mat, x[1]) then
                        mat := UnionOfRows(mat, x[1]);
                        vec := UnionOfRows(vec, x[2]);
                    fi;
                fi;               
            od;
        fi;
    od;

    y := MAJORANA_SolutionAlgProducts(mat,vec,unknowns, algebraproducts, setup);
    
    return rec( mat := y.mat, vec := y.vec, unknowns := y.unknowns);
    
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
            elm,
            switch,
            j;
    
    if ForAll(mat!.indices, x -> x = []) then
        return rec( mat := SparseMatrix(0, Ncols(mat), [], [], Rationals), 
                    vec := SparseMatrix(0, Ncols(vec), [], [], Rationals),
                    unknowns := unknowns    );
    fi;
    
    Info(   InfoMajorana, 40, 
            STRINGIFY("Solving a ", Nrows(mat), " x ", Ncols(mat), " matrix") );
            
    for i in [1..Nrows(mat)] do 
        x := _FoldList2(mat!.entries[i], DenominatorRat, LcmInt);
        mat!.entries[i] := mat!.entries[i]*x;
        vec!.entries[i] := vec!.entries[i]*x;
    od;
    
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
    
    if Size(unknowns) > Size(x.unknowns) then 
        switch := true;
    else
        switch := false;
    fi;
    
    mat := CertainRows(x.mat, nonzero);
    vec := CertainRows(x.vec, nonzero);
    unknowns := x.unknowns;

    if switch = true then
        x := MAJORANA_SolutionAlgProducts(mat, vec, unknowns, algebraproducts, setup);
        
        mat := x.mat;
        vec := x.vec; 
        unknowns := x.unknowns;
    fi;
                                        
    return rec( mat := mat,
                vec := vec,
                unknowns := unknowns    );
    
    end );
    
InstallGlobalFunction( MAJORANA_SolveSingleSolution,

    function(x, mat, vec, unknowns, algebraproducts, setup) 
    
    local   elm, 
            y,
            switch,
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
    if Nrows(y.mat) > 0 then 
        nonzero := [];
                        
        for i in [1..Nrows(y.mat)] do              
            if y.mat!.indices[i] <> [] then 
                Add(nonzero,i);
            fi;
        od;
        
        mat := CertainRows(y.mat, nonzero);
        vec := CertainRows(y.vec, nonzero);
        unknowns := y.unknowns;
        
        switch := true;
    
        while switch = true do 
        
            switch := false;

            for i in [1..Nrows(mat)] do 
                if Size(mat!.indices[i]) = 1 then 
                    switch := true;
                    elm := mat!.entries[i][1];
                    MAJORANA_RecordSolution(    CertainRows(vec, [i])*(1/elm), 
                                                unknowns[mat!.indices[i][1]], 
                                                algebraproducts, setup);
                fi;;
            od;
            
            if switch = true then 
                Info( InfoMajorana, 60, "Solved a new single solution"); 
            fi;
            
            x := MAJORANA_RemoveKnownAlgProducts(   mat,
                                                    vec,
                                                    unknowns,
                                                    algebraproducts,
                                                    setup         );
                                                            
            nonzero := [];
            
            for i in [1..Nrows(x.mat)] do              
                if x.mat!.indices[i] <> [] then 
                    Add(nonzero,i);
                fi;
            od;
            
            mat := CertainRows(x.mat, nonzero);
            vec := CertainRows(x.vec, nonzero);
            unknowns := x.unknowns;
        od;
        
    fi;
                                        
    return rec( mat := mat,
                vec := vec,
                unknowns := unknowns    );
    
    end );
    
    
InstallGlobalFunction( MAJORANA_RecordSolution,

    function( v, x, algebraproducts, setup)
    
    local   y,
            g,
            sign;
    
    y := setup.pairorbit[x[1]][x[2]];
    g := SP_Inverse(setup.pairconjelts[setup.pairconj[x[1]][x[2]]]);
    
    if y > 0 then 
        sign := 1;
    else
        sign := -1;
        y := -y;
    fi;
    
    if algebraproducts[y] = false then 
        algebraproducts[y] := sign*MAJORANA_ConjugateVec(v,g,setup);              
    fi;      
    
    end );
    
InstallGlobalFunction( MAJORANA_RemoveKnownInnProducts,

    function(mat, vec, unknowns, innerproducts)

    local   unsolved, 
            i, j, 
            elm,
            prod;

    unsolved := [];
    
    for i in [1..Size(unknowns)] do 
        prod := innerproducts[unknowns[i]];
        
        if prod <> false then 
            for j in [1..Nrows(vec)] do 
                elm := GetEntry(mat, j, i);
                
                if elm <> 0 then 
                    AddToEntry(vec, j, 1, -elm*prod);
                fi;
            od;
        else
            Add(unsolved, i);
        fi;
    od;
        
    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};
    
    return rec( mat := mat, 
                vec := vec, 
                unknowns := unknowns);
                
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
            switch,
            pos,
            prod;
            
    if Nrows(mat) = 0 then 
    
        unknowns := MAJORANA_ExtractUnknownAlgebraProducts(algebraproducts, setup);
        mat := SparseMatrix(0, Size(unknowns), [], [], Rationals);
        vec := SparseMatrix(0, Size(setup.coords), [], [], Rationals); 
        
        return rec( mat := mat, vec := vec, unknowns := unknowns);
    fi;

    unsolved := [];
    
    switch := false;
    
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
            
            switch := true;
            
            g := setup.pairconjelts[setup.pairconj[x[1]][x[2]]];
            
            prod := MAJORANA_ConjugateVec(prod,g,setup);
            
            for j in [1..Nrows(vec)] do  
                pos := Position(mat!.indices[j], i);
                if pos <> fail then
                    elm := mat!.entries[j][pos];
                    AddRow( prod!.indices[1],-sign*elm*prod!.entries[1], 
                            vec!.indices, vec!.entries, j);
                fi;
            od;            
        else
            Add(unsolved,i);
        fi;
    od;
    
    mat := CertainColumns(mat, unsolved);
    unknowns := unknowns{unsolved};
    
    return rec( mat := mat, 
                vec := vec, 
                unknowns := unknowns);
        
    end );
    
InstallGlobalFunction( MAJORANA_SolutionInnerProducts,

    function( mat, vec, unknowns, innerproducts)
    
    local   sol,    # solution of system
            i,      # loop over <unknowns>
            x,      # element of <unknowns> 
            nonzero;   

    sol := MAJORANA_SolutionMatVecs(mat,vec);                   
        
    for i in [1..Size(sol.solutions)] do
        if sol.solutions[i] <> fail then
            x := unknowns[i]; 
            if sol.solutions[i]!.entries[1] = [] then 
                innerproducts[x] := 0;
            else
                innerproducts[x] := sol.solutions[i]!.entries[1][1];
            fi;
        fi;
    od;
    
    x := MAJORANA_RemoveKnownInnProducts(   sol.mat, sol.vec,
                                            unknowns, innerproducts );
                                                    
    nonzero := [];
    
    for i in [1..Nrows(x.mat)] do              
        if x.mat!.indices[i] <> [] then 
            Add(nonzero,i);
        fi;
    od;
    
    x.mat := CertainRows(x.mat, nonzero);
    x.vec := CertainRows(x.vec, nonzero);
                                        
    return rec( mat := x.mat,
                vec := x.vec,
                unknowns := x.unknowns    );
    
    end );
    
InstallGlobalFunction(MAJORANA_CheckNullSpace,

    function(innerproducts,setup)
    
    local   gram,     # full gram matrix
            null;     # nullspace of gram matrix
    
        if ForAll(innerproducts, x -> x <> false) then 
            gram := MAJORANA_FillGramMatrix(innerproducts, setup);
            null := KernelMat(gram).relations;; 
        fi;
        
        return null;
    
    end );
        
InstallGlobalFunction(MAJORANA_MoreEigenvectors,

    function(rep)
    
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
            
    dim := Size(rep.setup.coords);
    
    table := [0,1/4,1/32];

    for i in rep.setup.orbitreps do
        
        if IsMutable(rep.evecs[i]) then 
        
            a := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
            mat := SparseZeroMatrix(dim, dim, Rationals);

            for j in [1..dim] do
            
                b := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
                
                x := MAJORANA_AlgebraProduct(a,b,rep.algebraproducts,rep.setup);
                
                if x <> false then
                    mat!.indices[j] := x!.indices[1];
                    mat!.entries[j] := x!.entries[1];
                else
                    mat := fail;
                    break;
                fi;
            od;

            if mat = fail then 
                return;
            fi;
            
            for ev in [1..3] do
             
                Info(   InfoMajorana, 50, 
                        STRINGIFY( "Finding ", table[ev], " eigenvectors for axis ", i) ); 
                        
                rep.evecs[i][ev] := KernelMat( mat - SparseIdentityMatrix(dim, Rationals)*table[ev]).relations;
                
                rep.evecs[i][ev] := MAJORANA_BasisOfEvecs(rep.evecs[i][ev]);
                
            od;
            
            MakeImmutable(rep.evecs[i]);
        fi;
    od;
    
    end);

InstallGlobalFunction(MAJORANA_MainLoop,

    function(rep)
                                
    MAJORANA_AxiomM1(rep);
                                    
    MAJORANA_Fusion(rep);
            
    MAJORANA_UnknownAlgebraProducts(rep);

    MAJORANA_MoreEigenvectors(rep);
    
    MAJORANA_Orthogonality(rep);

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


    

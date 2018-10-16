# Table of fusion rules

BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4,     0, 1/32]
            ,[ 1/32, 1/32,  1/32, 1/4 ] ]);

InstallGlobalFunction(MAJORANA_TestEvecs,

    function(rep)
    
    local   u, v, i, j, k, ev, x, y;
    
    for i in rep.setup.orbitreps do 
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);
        
        for j in [1..3] do
            ev := MAJORANA_FusionTable[1, j + 1];
            for k in [1..Nrows(rep.evecs[i, j])] do 
                v := CertainRows(rep.evecs[i, j], [k]);
                
                x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
                
                if not x in [ev*v, false, fail] then
                    Error("evecs");
                fi;
            od;
        od;
    od;
    
    return true;

    end );
                    
InstallGlobalFunction( MAJORANA_IsAxis,

    function(rep, v)
    
    local dim, i, mat, basis, ev, evecs, u;
    
    basis := Positions(rep.setup.nullspace.heads, 0);
    
    dim := Size(basis);
    
    mat := SparseMatrix(0, Size(basis), [], [], Rationals);

    for i in basis do
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        mat := UnionOfRows(mat, MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup));
    od;
    
    evecs := [];;
    
    # Test for primitivity
    
    if Nrows(KernelMat(mat - SparseIdentityMatrix(dim, Rationals)).relations) > 1 then 
        return false;
    fi;
    
    for ev in [0, 1/4, 1/32] do 
        Add(evecs, KernelMat(mat - ev*SparseIdentityMatrix(dim, Rationals)).relations);
    od;
    
    # Test fusion
    
    if MAJORANA_TestFusionAxis(v, evecs, rep.innerproducts, rep.algebraproducts, rep.setup) = false then 
        return false;
    fi;
    
    return evecs;
    
    end );
    
InstallGlobalFunction( MAJORANA_TestAxiomM8, # Not actually working, giving false negatives

    function(rep)

    local i, j, k, t, u, v, w, g, dim, pos, l, m, sign, evecs, T, S, im;

    dim := Size(rep.setup.coords);

    t := Size(rep.involutions);
    T := ShallowCopy(rep.involutions);
    S := [];

    # First check that 2A algebras give a third axis

    for i in [1..t] do
    
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        pos := [];
        pos[1] := Position(AsList(rep.group), rep.involutions[i]);
        
        for j in [1..t] do
        
            k := rep.setup.pairorbit[i, j];
        
            if rep.shape[k] = ['2','A'] then
            
                Add(T, Product(rep.involutions{[i,j]}));
            
                v := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
                pos[2] := Position(AsList(rep.group), rep.involutions[j]);
             
                w := (-1/8)*MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup) + u + v;
                
                # Check if w is an axis
                
                evecs := MAJORANA_IsAxis(rep, w);
                
                if  evecs = false then return false; fi;
            
                # Check that \tau(u)\tau(v) = \tau(w)
            
                g := SP_Product(rep.setup.pairconjelts[pos[1]], rep.setup.pairconjelts[pos[2]]);
                
                for l in [1..3] do 
                    if l = 3 then sign := -1; else sign := 1; fi;
                    
                    for m in [1..Nrows(evecs[l])] do 
                        im := MAJORANA_ConjugateVec( CertainRows(evecs, [m]), g);
                        im := RemoveMatWithHeads(im, rep.setup.nullspace);
                       
                        if MAJORANA_AlgebraProduct(w, im, rep.algebraproducts, rep.setup) <> sign*im then 
                            return false;
                        fi;
                    od;
                od;
            elif rep.shape[k] = ['2','B'] then 
                if Product(rep.involutions{[i,j]}) in T then 
                    return false;
                else
                    Add(S, Product(rep.involutions{[i,j]}));
                fi;
            fi;
        od;
    od;
    
    # Then check that the sets T and S are disjoint
    
    T := DuplicateFreeList(T);;
    S := DuplicateFreeList(S);;
    
    if Intersection(T, S) <> [] then
        return false;
    fi;
    
    return T;
    
    end );
    
InstallGlobalFunction( MAJORANA_TestFusionAxis, 

    function(u, evecs, innerproducts, algebraproducts, setup)
    
    local   dim,            # size of setup.coords
            evals,
            new,
            a,              # first eigenvalue
            b,              # second eigenvalue
            ev_a,           # a - eigenvectors
            ev_b,           # b - eigenvectors
            ev,             # new eigenvalue
            i,              # loop over T 
            j, 
            k,
            x,              # product of eigenvectors
            y;              # product of x with u
            
    dim := Size(setup.coords);
    
    for evals in [[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]] do 
    
        new := [0,0,0];
        
        for j in [1..3] do 
            new[j] := SparseMatrix(0, dim, [], [], Rationals);
        od;
                
        ev_a := evecs[evals[1]];
        ev_b := evecs[evals[2]];
        
        for j in [1..Nrows(ev_a)] do  
            a := CertainRows(ev_a, [j]);
            for k in [1..Nrows(ev_b)] do
                b := CertainRows(ev_b, [k]);
                MAJORANA_FuseEigenvectorsNoForm(  a, b, u, evals, new, 
                                            innerproducts,
                                            algebraproducts,
                                            setup );
            od;
        od;
    
        for j in [1..3] do
            ev := MAJORANA_FusionTable[1, j + 1];
            
            new[j] := EchelonMatDestructive(new[j]).vectors;
            
            for k in [1..Nrows(new[j])] do 
                a := CertainRows(new[j], [k]);
                x := MAJORANA_AlgebraProduct(u, a, algebraproducts, setup);
                if not x in [ev*a, false, fail] then 
                    Error("The algebra does not obey the fusion rules");
                fi;
            od;
        od;
    od;
    
    end );

# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(rep) 
        
    local   dim, u, i;
            
    dim := Size(rep.setup.coords);

    for i in rep.setup.orbitreps do        
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        MAJORANA_TestFusionAxis(u, rep.evecs[i], rep.innerproducts, rep.algebraproducts, rep.setup);
    od;
    
    return true;
    
    end );
        
InstallGlobalFunction(MajoranaAlgebraTest,
    
    function(rep)
    
    MAJORANA_TestAxiomM1(rep);

    MAJORANA_TestFusion(rep);

    MAJORANA_TestPrimitivity(rep);

    return true;
    
    end );
        
InstallGlobalFunction(MAJORANA_TestOrthogonality,

    function(innerproducts,evecs, setup) 
    
    # Tests that eigenspaces are orthogonal with respect to the inner product

        local   errorortho, # list of indices which do not obey orthogonality of eigenvectors
                u,          # vector with 1 in j th position
                a,          # first eigenvalue
                b,          # second eigenvalue
                ev_a,       # list of a - eigenvectors
                ev_b,       # list of b - eigenvectors
                i,
                j,          # loop over T
                k,
                v,          # a - eigenvector
                w,          # b - eigenvector
                x;          # inner product
        
        errorortho := [];

        for i in setup.orbitreps do

            u := SparseMatrix(1, Size(setup.coords), [[i]], [[1]], Rationals);
            
            for a in [1..3] do 
            
                # orthogonality with 1-eigenvectors
                
                ev_a := evecs[i, a];
                
                for j in [1..Nrows(ev_a)] do  
                    v := CertainRows(ev_a, [j]);
                    x := MAJORANA_InnerProduct(u, v, innerproducts, setup);
                    
                    if (x <> false) and (x <> 0) then 

                        Add(errorortho, [i,0,a,u,v]);
                    fi;
                od;
                
                # orthogonality with all other eigenvectors
                
                for b in [a+1..3] do 
                
                    ev_b := evecs[i, b];
                    
                    for j in [1..Nrows(ev_a)] do  
                        v := CertainRows(ev_a, [j]);
                        for k in [1..Nrows(ev_b)] do
                            w := CertainRows(ev_b, [k]);
                            
                            x := MAJORANA_InnerProduct(v, w, innerproducts, setup);
                            
                            if (x <> false) and (x <> 0) then 
                                Add(errorortho, [i,a,b,v,w]);
                            fi;
                        od;
                    od;
                od;
            od;
        od;
        
        if Size(errorortho) > 0 then Error("Orthog"); fi;
        
        return true;
        
        end );
        
# Checks if bilinear and algebra products obey Fusion, outputs a list which is empty if they do obey the axiom

InstallGlobalFunction(MAJORANA_TestAxiomM1,

    function(rep) 

    local   ErrorM1,    # setup of indices which do not obey Fusion
            j,          # loop over algebra products
            k,          # loop over setup.coords
            l,
            p,          # second product
            dim,        # size of setup.coords
            x,          # first inner product
            y,          # second inner product
            u,          # vectors
            w,          #
            v;          #

    dim := Size(rep.setup.coords);

    ErrorM1:=[];
    
    for j in Filtered([1..Size(rep.algebraproducts)], i -> rep.algebraproducts[i] <> fail) do
        if not rep.algebraproducts[j] in [false, fail] then
            for k in Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0) do 
                for l in [rep.setup.pairreps[j], Reversed(rep.setup.pairreps[j])] do  
                    
                    u := SparseMatrix(1, dim, [[l[1]]], [[1]], Rationals);
                    v := SparseMatrix(1, dim, [[l[2]]], [[1]], Rationals);
                    w := SparseMatrix(1, dim, [[k]], [[1]], Rationals);
                    
                    p := MAJORANA_AlgebraProduct(v,w,rep.algebraproducts,rep.setup);
                    
                    if not p in [fail, false] then
                        x := MAJORANA_InnerProduct(u,p,rep.innerproducts, rep.setup);
                        y := MAJORANA_InnerProduct(rep.algebraproducts[j],w,rep.innerproducts, rep.setup);
                        
                        if x <> false and y <> false and x <> y then 
                            # return false;
                            Error("Axiom M1");
                            Add(ErrorM1,[l[1], l[2] ,k]);
                        fi;
                        
                    fi;
                od;
            od;
        fi;
    od;
    
    if Size(ErrorM1) > 0 then Error("Axiom M1"); fi;
    
    return true;

    end );

InstallGlobalFunction( MAJORANA_TestPrimitivity,

    function(rep)
    
    local i, dim, u, v, j, x, mat, espace, basis;
    
    if false in rep.algebraproducts then return fail; fi;
    
    dim := Size(rep.setup.coords);
    
    basis := Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0);
    
    for i in rep.setup.orbitreps do 
        
        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        mat := SparseMatrix(0, Size(basis), [], [], Rationals);
        
        for j in basis do 
            v := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
        
            x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
            
            mat := UnionOfRows(mat, x);
        od;
        
        espace := KernelMat(mat - SparseIdentityMatrix(dim, Rationals));
        
        if Nrows(espace.relations) <> 1 then 
            Error("Primitivity");
        fi;
    od;
    
    return true;
    
    end );

InstallGlobalFunction( MAJORANA_IsComplete, 

    function(rep)
    
    if false in rep.algebraproducts then 
        return false; 
    else
        return true;
    fi;
    
    end );
    
InstallGlobalFunction( MAJORANA_TestPositiveDefiniteForm,

    function(rep)
    
    local dim, gram;
    
    dim := Size(rep.setup.coords);
    
    gram := MAJORANA_FillGramMatrix(Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0), rep.innerproducts, rep.setup);

    if MAJORANA_PositiveDefinite(ConvertSparseMatrixToMatrix(gram)) < 0 then 
        return false; 
    else 
        return true; 
    fi;
    
    end );
    
InstallGlobalFunction(MAJORANA_TestAxiomM2,

    function(rep) # Tests that the algebra obeys axiom M2

    local   B,      # matrix of inner products
            dim,    # size of setup.coords
            j,      # loop through setup.coords
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
            x3,     #
            basis;

    dim:=Size(rep.setup.coords);
    
    basis := Filtered([1..dim], i -> rep.setup.nullspace.heads[i] = 0);

    B:=NullMat(dim^2,dim^2);

    for j in basis do
        for k in basis do
            for l in basis do
                for m in basis do
                    
                    a := SparseMatrix(1, dim, [[j]], [[1]], Rationals);
                    b := SparseMatrix(1, dim, [[k]], [[1]], Rationals);
                    c := SparseMatrix(1, dim, [[l]], [[1]], Rationals);
                    d := SparseMatrix(1, dim, [[m]], [[1]], Rationals);

                    x0 := MAJORANA_AlgebraProduct(a,c,rep.algebraproducts,rep.setup);
                    x1 := MAJORANA_AlgebraProduct(b,d,rep.algebraproducts,rep.setup);
                    x2 := MAJORANA_AlgebraProduct(b,c,rep.algebraproducts,rep.setup);
                    x3 := MAJORANA_AlgebraProduct(a,d,rep.algebraproducts,rep.setup);
                
                    B[dim*(j-1) + k, dim*(l-1) +m]:=
                          MAJORANA_InnerProduct(x0,x1,rep.innerproducts, rep.setup)
                        - MAJORANA_InnerProduct(x2,x3,rep.innerproducts, rep.setup);
                od;
            od;
        od;
    od;
    
    if MAJORANA_PositiveDefinite(B) < 0 then 
        return false; 
    else 
        return true; 
    fi;

    end );
    
InstallGlobalFunction( MAJORANA_SubalgebraAxes,

    function(rep, t)
    
    local n, H, tau, axes, U, V, i, j, u, v, prod, field;
    
    if IsBound(rep.field) then 
        field := rep.field; 
    else
        field := Rationals;
    fi;
    
    n := Size(rep.setup.coords);
    
    H := Group(rep.involutions{t});;
        
    axes := DuplicateFreeList(Concatenation( List( t, i -> Orbit(H, i) ) ));
    
    tau := List(rep.involutions{axes}, x -> RestrictedPerm(x, axes));
    
    U := NullMat(Size(axes), n);
    
    for i in [1..Size(axes)] do 
        U[i, axes[i]] := One(field);
    od;
    
    while true do 
    
        V := StructuralCopy(U);
        
        for i in [1..Size(V)] do
            u := SparseMatrix( [V[i]], field);
            for j in [i + 1 .. Size(V)] do 
                v := SparseMatrix( [V[j]], field);
                prod := MAJORANA_AlgebraProduct( u, v, rep.algebraproducts, rep.setup);
                
                if prod in [false, fail] then return false; fi;
                
                Add(U, ConvertSparseMatrixToMatrix(prod)[1]);
            od;
        od;
    
        U := SemiEchelonMatDestructive(U).vectors;
        
        if Size(V) = Size(U) then break; fi;
        
    od;

    return [tau,U];
    
    end );
    
InstallGlobalFunction( IsMinimal3GeneratedAlgebra,

    function(rep)

    local n, dim, k, triples, U, t, type, dih, x, list, bad;

    k := Size(rep.involutions);;
    dim := MAJORANA_Dimension(rep);
    n := Size(rep.setup.coords);
    
    triples := Combinations(Filtered([1..k], i -> rep.setup.nullspace.heads[i] = 0), 3);;
    
    bad := [];
    
    for t in triples do 
        
        # Construct the full algebra generated by the axes
        
        U := MAJORANA_SubalgebraAxes(rep, t);
        
        if U <> false then  # in case of incomplete algebra
        
            # If the dimension of U is the same as the whole algebra then good.
            # Otherwise, if the dimension of U is the same as the dihedral algebra generated by two of the axes, good.
            # Otherwise, return false.
            
            if Size(U[2]) <> dim then 
            
                list := [];
            
                for x in Combinations(t, 2) do 
                    type := rep.shape[rep.setup.pairorbit[x[1], x[2]]];  
                    dih := MAJORANA_DihedralAlgebras.(type);
                
                    Add(list, Size(dih.setup.coords));
                od;
            
                if not Size(U[2]) in list then Add(bad, t); fi;
            
            fi;
        fi;
    od;
    
    if Size(bad) > 0 then return bad; fi;
    
    return true;
        
    end );
    
InstallGlobalFunction( MAJORANA_TestSetup,

    function(rep)
    
    local dim, i, j, k, g, sign, sign_k, im;
    
    dim := Size(rep.setup.coords);
    
    for i in [1 .. dim] do 
        for j in [i .. dim] do 
            k := rep.setup.pairorbit[i, j];
            g := rep.setup.pairconj[i, j];
            g := rep.setup.pairconjelts[g];
            
            sign_k := 1;
            
            if k < 0 then k := -k; sign_k := -1; fi;
            
            im := g{rep.setup.pairreps[k]};
            
            sign := 1;
            
            if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
            if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;
            
            if SortedList(im) <> [i,j] then Error("Does not conjugate to correct pair"); fi;
            
            if sign <> sign_k then Error("Sign error"); fi;
        od;
    od;
    
    end );
        
        
        
        
         

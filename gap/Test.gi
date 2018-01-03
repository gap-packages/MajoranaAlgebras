# Table of fusion rules

BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4,     0, 1/32]
            ,[ 1/32, 1/32,  1/32, 1/4 ] ]);

InstallGlobalFunction(MAJORANA_EvecsTest,

    function(rep)
    
    local   u, v, i, j, k, ev, x, y;
    
    for i in rep.setup.orbitreps do 
        u := SparseMatrix(1, Size(rep.setup.coords), [[i]], [[1]], Rationals);
        
        for j in [1..3] do
            ev := MAJORANA_FusionTable[1][j + 1];
            for k in [1..Nrows(rep.evecs[i][j])] do 
                v := CertainRows(rep.evecs[i][j], [k]);
                
                x := MAJORANA_AlgebraProduct(u, v, rep.algebraproducts, rep.setup);
                
                if not x in [ev*v, false] then
                    y := MAJORANA_InnerProduct(x - ev*v, x - ev*v, rep.innerproducts, rep.setup);
                    if not y in [0, false] then
                        Error("evecs");
                    fi;
                fi;
            od;
        od;
    od;

    end );
                    

# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(innerproducts,algebraproducts,evecs,setup) 
        
    local   errorfusion,    # list of indices which do not obey fusion rules
            dim,            # size of setup.coords
            a,              # first eigenvalue
            b,              # second eigenvalue
            ev_a,           # a - eigenvectors
            ev_b,           # b - eigenvectors
            ev,             # new eigenvalue
            u,              # vector with 1 in i th position
            i,              # loop over T 
            j, 
            k,
            v,              # a - eigenvector
            w,              # b - eigenvector
            x,              # product of eigenvectors
            y,              # product of x with u
            z,              # inner product where needed
            x0,             # further product in 1/32 case
            test;

    errorfusion:=[];

    dim := Size(setup.coords);

    for i in setup.orbitreps do

        u := SparseMatrix(1, dim, [[i]], [[1]], Rationals);
        
        for a in [1..3] do 
            for b in [a..3] do 
                if not [a,b] in [[2,2],[3,3]] then 
                    
                    ev_a := evecs[i][a];
                    ev_b := evecs[i][b];
                    
                    ev := MAJORANA_FusionTable[a + 1][b + 1];
                    
                    for j in [1..Nrows(ev_a)] do  
                        v := CertainRows(ev_a, [j]);
                        if MAJORANA_InnerProduct(v, v, innerproducts, setup) <> 0 then 
                            for k in [1..Nrows(ev_b)] do
                                w := CertainRows(ev_b, [k]);
                                if MAJORANA_InnerProduct(w, w, innerproducts, setup) <> 0 then 
                            
                                    x := MAJORANA_AlgebraProduct(v,w,algebraproducts,setup);
                                    
                                    if x <> false then
                                        y:=MAJORANA_AlgebraProduct(u,x,algebraproducts,setup);
                                        
                                        if y <> false and y <> ev * x then 
                                            
                                            test := MAJORANA_InnerProduct(y - ev*x, y - ev*x, innerproducts, setup);
                                            
                                            if not test in [0, false] then 
                                                Error("Fusion");
                                            
                                                Add(errorfusion,[i,a,b,v,w]);
                                            fi;
                                        fi;
                                    fi;
                                fi;
                            od;
                        fi;
                    od;
                    
                # 1/4,1/4 fusion is special
                    
                elif [a,b] = [2,2] then 

                    ev_a := evecs[i][a];
                    ev_b := evecs[i][b];
                    
                    ev := 0; 
                    
                    for j in [1..Nrows(ev_a)] do  
                        v := CertainRows(ev_a, [j]);
                        if MAJORANA_InnerProduct(v, v, innerproducts, setup) <> 0 then 
                            for k in [1..Nrows(ev_b)] do
                                if MAJORANA_InnerProduct(w, w, innerproducts, setup) <> 0 then 
                                    w := CertainRows(ev_b, [k]);
                                
                                    x := MAJORANA_AlgebraProduct(v,w,algebraproducts,setup);
                                    
                                    if x <> false then
                                        
                                        z := MAJORANA_InnerProduct(u,x,innerproducts, setup);
                                            
                                        if z <> false then 
                                            
                                            x := x - z*u;
                                            y := MAJORANA_AlgebraProduct(u,x,algebraproducts,setup);
                                            
                                            if y <> false and y <> ev * x then 
                                                test := MAJORANA_InnerProduct(y - ev*x, y - ev*x, innerproducts, setup);
                                            
                                                if not test in [0, false] then 
                                                    Error("Fusion");
                                                    Add(errorfusion,[i,a,b,v,w]);
                                                fi;
                                            fi;
                                        fi;
                                    fi;
                                fi;
                            od;
                        fi;
                    od;
                
                # 1/32,1/32 fusion is even more special
                    
                else
                
                    ev_a := evecs[i][a];
                    ev_b := evecs[i][b];
                    
                    ev := 0; 
                    
                    for j in [1..Nrows(ev_a)] do  
                        v := CertainRows(ev_a, [j]);
                        if MAJORANA_InnerProduct(v, v, innerproducts, setup) <> 0 then 
                            for k in [1..Nrows(ev_b)] do
                                w := CertainRows(ev_b, [k]);
                                if MAJORANA_InnerProduct(w,w, innerproducts, setup) <> 0 then 
                            
                                    x := MAJORANA_AlgebraProduct(v,w,algebraproducts,setup);
                                    
                                    if x <> false then 
                                    
                                        x0 := MAJORANA_AlgebraProduct(u,x,algebraproducts,setup);
                                        
                                        if x0 <> false then 
                                        
                                            y := MAJORANA_InnerProduct(u,x,innerproducts,setup);
                        
                                            if y <> false then 
                                                
                                                x0 := x0 - y*u;
                                                
                                                z := MAJORANA_AlgebraProduct(u,x0,algebraproducts,setup);
                                                
                                                if (z <> false) and (z <> x0*(1/4)) then  
                                                
                                                    test := MAJORANA_InnerProduct(z - x0*(1/4),z - x0*(1/4), innerproducts, setup);
                                            
                                                    if not test in [0, false] then 
                                                        Error("Fusion");
                                                        Add(errorfusion,[i,a,b,v,w]);
                                                    fi;
                                                
                                                fi;
                                            fi;
                                        fi;
                                    fi;
                                fi;
                            od;
                        fi;
                    od;
                fi;
            od;
        od;
    od;
    
    return errorfusion;
    
    end
    
    );
        
InstallGlobalFunction(MajoranaAlgebraTest,
    
    function(rep)
    
    local   error,
            GramMatrixFull;
    
    MAJORANA_EvecsTest(rep);
    
    # Check bilinear form is positive definite

    if false then
        GramMatrixFull := MAJORANA_FillGramMatrix(rep.innerproducts, rep.setup);
        
        if MAJORANA_PositiveDefinite(GramMatrixFull) <0 then
            return "Gram Matrix is not positive definite";
        fi;
 
    fi;

    # Check that all triples obey Fusion

    error := MAJORANA_AxiomM1(rep.innerproducts,rep.algebraproducts,rep.setup);

    if Size(error)>0 then
        return ["Algebra does not obey Fusion", error];
    fi;

    # Check that eigenvectors obey the fusion rules

    error := MAJORANA_TestFusion(rep.innerproducts,rep.algebraproducts,rep.evecs,rep.setup);

    if ForAny(error,x->Size(x)>0) then
        return ["Algebra does not obey fusion rules", error];
    fi;

    # Check that the eigenspaces are orthogonal

    error := MAJORANA_TestOrthogonality(rep.innerproducts,rep.evecs,rep.setup);

    if Size(error) > 0 then
        return ["Eigenspaces are not orthogonal with reppect to the inner product", error];
    fi;

    # Check M2

    # error := MAJORANA_AxiomM2(rep.innerproducts,rep.algebraproducts,rep.setup);

    # if error = -1 then
    #    return "Algebra does not obey axiom M2";
    # fi;
    
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
                
                ev_a := evecs[i][a];
                
                for j in [1..Nrows(ev_a)] do  
                    v := CertainRows(ev_a, [j]);
                    x := MAJORANA_InnerProduct(u, v, innerproducts, setup);
                    
                    if (x <> false) and (x <> 0) then 
                    
                        Error("Orthog");
                        Add(errorortho, [i,0,a,u,v]);
                    fi;
                od;
                
                # orthogonality with all other eigenvectors
                
                for b in [a+1..3] do 
                
                    ev_b := evecs[i][b];
                    
                    for j in [1..Nrows(ev_a)] do  
                        v := CertainRows(ev_a, [j]);
                        for k in [1..Nrows(ev_b)] do
                            w := CertainRows(ev_b, [k]);
                            
                            x := MAJORANA_InnerProduct(v, w, innerproducts, setup);
                            
                            if (x <> false) and (x <> 0) then 
                                Error("Orthog");
                                Add(errorortho, [i,a,b,v,w]);
                            fi;
                        od;
                    od;
                od;
            od;
        od;
        
        return errorortho;
        
        end
        
        );
        
# Checks if bilinear and algebra products obey Fusion, outputs a list which is empty if they do obey the axiom

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(innerproducts,algebraproducts,setup) 

    local   ErrorM1,    # setup of indices which do not obey Fusion
            j,          # loop over algebra products
            k,          # loop over setup.coords
            p,          # second product
            dim,        # size of setup.coords
            x,          # first inner product
            y,          # second inner product
            u,          # vectors
            w,          #
            v;          #

    dim := Size(setup.coords);

    ErrorM1:=[];
    
    for j in [1..Size(algebraproducts)] do
        if algebraproducts[j] <> false then
            for k in [1..dim] do 
            
                u := SparseMatrix(1, dim, [[setup.pairreps[j][1]]], [[1]], Rationals);
                v := SparseMatrix(1, dim, [[setup.pairreps[j][2]]], [[1]], Rationals);
                w := SparseMatrix(1, dim, [[k]], [[1]], Rationals);
                
                p := MAJORANA_AlgebraProduct(v,w,algebraproducts,setup);
                
                if p <> false then
                    x := MAJORANA_InnerProduct(u,p,innerproducts, setup);
                    y := MAJORANA_InnerProduct(algebraproducts[j],w,innerproducts, setup);
                    
                    if x <> false and y <> false and x <> y then 
                        Error("Axiom M1");
                        Add(ErrorM1,[j,k]);
                    fi;
                    
                fi;
            od;
        fi;
    od;
                
    return ErrorM1;

    end

    );

InstallGlobalFunction(MAJORANA_AxiomM2,

        function(innerproducts,algebraproducts,setup) # Tests that the algebra obeys axiom M2

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
                x3;     #

        dim:=Size(algebraproducts[1]);

        B:=NullMat(dim^2,dim^2);

        for j in [1..dim] do
            for k in [1..dim] do
                for l in [1..dim] do
                    for m in [1..dim] do
                        
                        a := [1..dim]*0; a[j] := 1; 
                        b := [1..dim]*0; b[k] := 1;
                        c := [1..dim]*0; c[l] := 1;
                        d := [1..dim]*0; d[m] := 1;
                        
                        x0 := MAJORANA_AlgebraProduct(a,c,algebraproducts,setup);
                        x1 := MAJORANA_AlgebraProduct(b,d,algebraproducts,setup);
                        x2 := MAJORANA_AlgebraProduct(b,c,algebraproducts,setup);
                        x3 := MAJORANA_AlgebraProduct(a,d,algebraproducts,setup);
                    
                        B[dim*(j-1) + k][dim*(l-1) +m]:=
                              MAJORANA_InnerProduct(x0,x1,innerproducts, setup)
                            - MAJORANA_InnerProduct(x2,x3,innerproducts, setup);
                    od;
                od;
            od;
        od;
        
        return MAJORANA_PositiveDefinite(B);

        end

        );

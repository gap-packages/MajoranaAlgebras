# Table of fusion rules

BindGlobal("MAJORANA_FusionTable",
           [ [    1,    0,   1/4, 1/32]
            ,[    0,    0,   1/4, 1/32]
            ,[  1/4,  1/4,     0, 1/32]
            ,[ 1/32, 1/32,  1/32, 1/4 ] ]);

# Checks if algebra obeys the fusion rules, outputs list which is empty if it does obey fusion rules

InstallGlobalFunction(MAJORANA_TestFusion,

    function(GramMatrix,AlgebraProducts,EigenVectors,ProductList) 
        
    # list should be of the form [ProductList[1],ProductList[2],ProductList[3],ProductList[4],ProductList[5],ProductList[6]]
        
        local   errorfusion,    # list of indices which do not obey fusion rules
                dim,            # size of ProductList[1]
                a,              # first eigenvalue
                b,              # second eigenvalue
                ev_a,           # a - eigenvectors
                ev_b,           # b - eigenvectors
                ev,             # new eigenvalue
                u,              # vector with 1 in i th position
                j,              # loop over T 
                v,              # a - eigenvector
                w,              # b - eigenvector
                x,              # product of eigenvectors
                y,              # product of x with u
                z,              # inner product where needed
                x0;             # further product in 1/32 case

        errorfusion:=[];

        dim := Size(AlgebraProducts[1]);

        for j in ProductList[10] do

            u := [1..dim]*0; u[j]:=1;
            
            for a in [1..3] do 
                for b in [a..3] do 
                    if not [a,b] in [[2,2],[3,3]] then 
                        
                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := MAJORANA_FusionTable[a + 1][b + 1];

                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then
                                
                                    y:=MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                    
                                    if y <> false and y <> ev * x then 
                                        Add(errorfusion,[j,a,b,v,w]);
                                    fi;
                                fi;
                            od;
                        od;
                        
                    # 1/4,1/4 fusion is special
                        
                    elif [a,b] = [2,2] then 

                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := 0; 
                        
                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then
                                    
                                    z := MAJORANA_InnerProduct(u,x,GramMatrix, ProductList);
                                        
                                    if z <> false then 
                                        
                                        x := x - z*u;
                                        y := MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                        
                                        if y <> false and y <> ev * x then 
                                            Add(errorfusion,[j,a,b,v,w]);
                                        fi;
                                        
                                    fi;
                                fi;
                            od;
                        od;
                    
                    # 1/32,1/32 fusion is even more special
                        
                    else
                    
                        ev_a := EigenVectors[j][a];
                        ev_b := EigenVectors[j][b];
                        
                        ev := 0; 
                        
                        for v in ev_a do
                            for w in ev_b do
                            
                                x := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,ProductList);
                                
                                if x <> false then 
                                
                                    x0 := MAJORANA_AlgebraProduct(u,x,AlgebraProducts,ProductList);
                                    
                                    if x0 <> false then 
                                    
                                        y := MAJORANA_InnerProduct(u,x,GramMatrix,ProductList);
                    
                                        if y <> 0 then 
                                            
                                            x0 := x0 - y*u;
                                            
                                            z := MAJORANA_AlgebraProduct(u,x0,AlgebraProducts,ProductList);
                                            
                                            if (z <> false) and (z <> x0/4) then  
                                            
                                                Add(errorfusion,[j,a,b,v,w]);
                                            
                                            fi;
                                        fi;
                                    fi;
                                fi;
                            od;
                        od;
                    fi;
                od;
            od;
        od;
        
        return errorfusion;
        
        end
        
        );
        
InstallGlobalFunction(MAJORANA_TestOrthogonality,

    function(GramMatrix,AlgebraProducts,EigenVectors, ProductList) # Tests that eigenspaces are orthogonal with respect to the inner product

        local   errorortho, # list of indices which do not obey orthogonality of eigenvectors
                u,          # vector with 1 in j th position
                a,          # first eigenvalue
                b,          # second eigenvalue
                ev_a,       # list of a - eigenvectors
                ev_b,       # list of b - eigenvectors
                j,          # loop over T
                v,          # a - eigenvector
                w,          # b - eigenvector
                x;          # inner product
        
        errorortho := [];

        for j in ProductList[10] do

            u := [1..Size(AlgebraProducts[1])]*0; u[j]:=1;
            
            for a in [1..3] do 
            
                # orthogonality with 1-eigenvectors
                
                ev_a := EigenVectors[j][a];
                
                for v in ev_a do
                    x := MAJORANA_InnerProduct(u, v, GramMatrix, ProductList);
                    
                    if (x <> false) and (x <> 0) then 
                        Add(errorortho, [j,1,a,u,v]);
                    fi;
                od;
                
                # orthogonality with all other eigenvectors
                
                for b in [a+1..3] do 
                
                    ev_b := EigenVectors[j][b];
                    
                    for v in ev_a do
                        for w in ev_b do
                            x := MAJORANA_InnerProduct(v, w, GramMatrix, ProductList);
                            
                            if (x <> false) and (x <> 0) then 
                                Add(errorortho, [j,a,b,u,v]);
                            fi;
                        od;
                    od;
                od;
            od;
        od;
        
        return errorortho;
        
        end
        
        );
        
# Checks if bilinear and algebra products obey axiom M1, outputs a list which is empty if they do obey the axiom

InstallGlobalFunction(MAJORANA_AxiomM1,

    function(GramMatrix,AlgebraProducts,list) 

    # list should be of the form [ProductList[1],ProductList[2],ProductList[3],ProductList[4],ProductList[5]]

        local   ErrorM1,    # list of indices which do not obey axiom M1
                j,          # loop over algebra products
                k,          # loop over ProductList[1]
                p,          # second product
                dim,        # size of ProductList[1]
                x,          # first inner product
                y,          # second inner product
                u,          # vectors
                w,          #
                v;          #

        dim:=Size(AlgebraProducts[1]);

        ErrorM1:=[];
        
        for j in [1..Size(AlgebraProducts)] do
            if AlgebraProducts[j] <> false then
                for k in [1..dim] do 
                    
                    u := NullMat(1,dim)[1];
                    u[list[7][j][1]] := 1;
                    
                    v := NullMat(1,dim)[1];
                    v[list[7][j][2]] := 1;
                    
                    w := NullMat(1,dim)[1];
                    w[k] := 1;
                    
                    p := MAJORANA_AlgebraProduct(v,w,AlgebraProducts,list);
                    
                    if p <> false then
                        x := MAJORANA_InnerProduct(u,p,GramMatrix, list);
                        y := MAJORANA_InnerProduct(AlgebraProducts[j],w,GramMatrix, list);
                        
                        if x <> false and y <> false and x <> y then 
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

        function(GramMatrix,AlgebraProducts,ProductList) # Tests that the algebra obeys axiom M2

        local   B,      # matrix of inner products
                dim,    # size of ProductList[1]
                j,      # loop through ProductList[1]
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
                              MAJORANA_InnerProduct(x0,x1,GramMatrix, ProductList)
                            - MAJORANA_InnerProduct(x2,x3,GramMatrix, ProductList);
                    od;
                od;
            od;
        od;
        
        return MAJORANA_PositiveDefinite(B);

        end

        );

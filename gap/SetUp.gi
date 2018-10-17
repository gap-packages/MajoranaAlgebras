##
## Finds all possible shapes of a Majorana representation of the form (G,T,V)
## that obey axiom M8.
##

InstallGlobalFunction(ShapesOfMajoranaRepresentationAxiomM8,

    function(G,T)

    local   gens, g, perm,
            t,              # size of T
            i,              # indices
            j,
            k,
            x,              # result of orbitals
            shape,          # one shape
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            unknowns,       # indices of 3X axes
            pos,            # positions
            Binaries,       # used to loop through options for shapes
            input;          #

    t := Size(T);

    # Check that T obeys axiom M8

    for i in [1..t] do
        for j in [1..t] do
            if Order(T[i]*T[j]) = 6 and not (T[i]*T[j])^3 in T then
                Error("The set T does not obey axiom M8");
            fi;
        od;
    od;

    # Construct orbitals of  on T x T

    gens := [];

    for g in GeneratorsOfGroup(G) do
        perm := [];
        for i in [1..t] do
            Add(perm, Position(T, T[i]^g));
        od;
        Add(gens, perm);
    od;

    input := rec();

    input.pairorbit := NullMat(t,t);
    input.pairconj  := NullMat(t,t);
    input.pairreps  := [];
    input.orbitals  := [];
    input.pairconjelts := [ [1..t] ];
    input.coords := [1..t];
    input.involutions := T;
    input.group       := G;

    MAJORANA_Orbitals(gens, 0, input);

    input.orbitals := List( input.orbitals, x -> List(x, y -> T{y}) );

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(input.pairreps))[1];

    RepsSquares6A := [];
    unknowns := [];;

    for i in [1..Size(input.pairreps)] do

        x := T{input.pairreps[i]};

        if Order(x[1]*x[2]) = 1 then
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 and x[1]*x[2] in T then
            shape[i]:="2A";
        elif Order(x[1]*x[2]) = 2 and not x[1]*x[2] in T then
            shape[i]:="2B";
        elif Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
            Add(unknowns,i);
        elif Order(x[1]*x[2]) = 4 and not (x[1]*x[2])^2 in T then
            shape[i]:="4A";
        elif Order(x[1]*x[2]) = 4 and (x[1]*x[2])^2 in T then
            shape[i]:="4B";
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(RepsSquares6A,(x[1]*x[2])^2);
        else
            Error("This is not a 6-transposition group");
        fi;
    od;

    # Check for inclusions of 2A and 3A in 6A

    for i in unknowns do
        if ForAny(input.orbitals[i], x -> x[1]*x[2] in RepsSquares6A) then
            shape[i]:="3A";;
            unknowns := Difference(unknowns, [i]);
        fi;
    od;

    Binaries := AsList(FullRowSpace(GF(2),Size(unknowns)));

    input.shapes := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do

        for j in [1..Size(unknowns)] do
            k := unknowns[j];
            if Binaries[i, j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;
        od;

        Add(input.shapes,ShallowCopy(shape));
    od;

    return input;

    end );

##
## Finds all possible shapes of a Majorana representation of the form (G,T,V)
##

InstallGlobalFunction(ShapesOfMajoranaRepresentation,

    function(G,T)

    local   gens,
            t,              # size of T
            i,              # indices
            j,
            k,
            perm,
            x,              # result of orbitals
            ind,            # list of indices
            orbs,           # orbitals on T
            shape,          # one shape
            RepsSquares4X,  # (ts)^2 where o(ts) = 4
            RepsSquares6A,  # (ts)^2 where o(ts) = 6
            RepsCubes6A,    # (ts)^3 where o(ts) = 6
            gph,            # digraph of 2X, 4X inclusions
            cc,             # connected components of gph
            pos,            # positions
            Binaries,       # used to loop through options for shapes
            input;          #

    t := Size(T);

    # Construct orbitals of  on T x T

    input := rec();

    input.pairorbit := NullMat(t,t);
    input.pairconj  := NullMat(t,t);
    input.pairreps  := [];
    input.orbitals  := [];
    input.pairconjelts := [ [1..t] ];
    input.coords := [1..t];
    input.involutions := T;
    input.group       := G;

    gens := [];

    for g in GeneratorsOfGroup(G) do
        perm := [];
        for i in [1..t] do
            Add(perm, Position(T, T[i]^g));
        od;
        Add(gens, perm);
    od;

    MAJORANA_Orbitals(gens, 0, input);

    input.orbitals := List( input.orbitals, x -> List(x, y -> T{y}) );

    # Determine occurances of 1A, 2A, 2B, 4A, 4B 5A, 6A in shape

    shape := NullMat(1,Size(input.pairreps))[1];

    RepsSquares4X := [];
    RepsSquares6A := [];
    RepsCubes6A := [];

    ind := NullMat(6,0);;

    for i in [1..Size(input.pairreps)] do

        x := T{input.pairreps[i]};

        if Order(x[1]*x[2]) = 1 then
            shape[i] := "1A";
        elif Order(x[1]*x[2]) = 2 then
            shape[i]:="2X";
            Add(ind[2],i);
        elif Order(x[1]*x[2]) = 3 then
            shape[i]:="3X";
            Add(ind[3],i);
        elif Order(x[1]*x[2]) = 4 then
            shape[i]:="4X";
            Add(ind[4],i);
            Add(RepsSquares4X, (x[1]*x[2])^2);
        elif Order(x[1]*x[2]) = 5 then
            shape[i]:="5A";
        elif Order(x[1]*x[2])=6 then
            shape[i]:="6A";
            Add(ind[6],i);
            Add(RepsSquares6A,(x[1]*x[2])^2);
            Add(RepsCubes6A,(x[1]*x[2])^3);
        fi;
    od;

    # Check for inclusions of 2X in 4X

    gph := NullMat(Size(input.orbitals), 0);

    for i in ind[2] do
        for x in input.orbitals[i] do
            pos := Positions(RepsSquares4X, x[1]*x[2]);

            if pos <> [] then
                Append(gph[i], ind[4]{pos} + Size(ind[2]));
            fi;
        od;
    od;

    gph := List(gph, DuplicateFreeList);

    cc := AutoConnectedComponents(gph);

    # Check for inclusions of 2A and 3A in 6A

    for i in ind[3] do
        if ForAny(input.orbitals[i], x -> x[1]*x[2] in RepsSquares6A) then
            shape[i]:="3A";;
            ind[3] := Difference(ind[3], [i]);
        fi;
    od;

    for i in ind[2] do
        if ForAny(input.orbitals[i], x -> x[1]*x[2] in RepsCubes6A) then

            shape[i]:="2A";;

            for x in cc do
                if i in x then
                    for j in Intersection(ind[2],x) do
                        shape[j] := "2A";
                    od;
                    for j in Intersection(ind[4],x - Size(ind[2])) do
                        shape[j] := "4B";
                    od;

                    cc := Difference(cc, [x]);

                fi;
            od;
        fi;
    od;

    cc := Filtered(cc, x -> Size(Intersection(x,ind[2])) > 0);

    Binaries := AsList(FullRowSpace(GF(2),Size(ind[3]) + Size(cc)));

    input.shapes := [];

    # Add new values in the shape

    for i in [1..Size(Binaries)] do

        for j in [1..Size(ind[3])] do
            k:=ind[3, j];
            if Binaries[i, j] = 1*Z(2) then
                shape[k]:="3A";
            else
                shape[k]:="3C";
            fi;
        od;

        for j in [1 .. Size(cc)] do

            if Binaries[i, j + Size(ind[3])] = 1*Z(2) then
                for k in Intersection(ind[2],cc[j]) do
                    shape[k] := "2A";
                od;
                for k in Intersection(ind[4],cc[j] - Size(ind[2])) do
                    shape[k] := "4B";
                od;
            else
                for k in Intersection(ind[2],cc[j]) do
                    shape[k] := "2B";
                od;
                for k in Intersection(ind[4],cc[j] - Size(ind[2])) do
                    shape[k] := "4A";
                od;
            fi;
        od;

        Add(input.shapes,ShallowCopy(shape));
    od;

    return input;

    end );

##
## The main setup function for the algorithm <MajoranaRepresentation>
##

InstallGlobalFunction( MAJORANA_SetUp,

    function(input, index, axioms)

    local rep, s, t, i, j, k, gens, orbs, dim, algebras;

    rep         := rec( group       := input.group,
                        involutions := input.involutions,
                        shape       := input.shapes[index],
                        axioms      := axioms   );

    t := Size(rep.involutions);

    rep.setup   := rec( coords          := [1..t],
                        coordmap        := HashMap( t*t ),
                        pairorbit       := StructuralCopy(input.pairorbit),
                        pairconj        := StructuralCopy(input.pairconj),
                        pairconjelts    := StructuralCopy(input.pairconjelts),
                        pairreps        := ShallowCopy(input.pairreps)       );

    for i in [1..t] do
        rep.setup.coordmap[i] := i;
        rep.setup.coordmap[rep.involutions[i]] := i;
    od;

    algebras := MAJORANA_DihedralAlgebras;

    ## Orbits on axes for eigenvectors

    gens := GeneratorsOfGroup(input.group);
    gens := List(gens, x -> MAJORANA_FindPerm(x, rep, rep));

    orbs := MAJORANA_Orbits(gens, t, rep.setup);

    rep.setup.conjelts := orbs.conjelts;
    rep.setup.orbitreps := orbs.orbitreps;

    ## Set up products and eigenvectors

    s := Size(rep.setup.pairreps);

    rep.algebraproducts := List([1..s], x -> false);
    rep.innerproducts   := List([1..s], x -> false);
    rep.evecs           := NullMat(t,3);

    for j in [1..t] do
        if j in rep.setup.orbitreps then
            for k in [1..3] do
                rep.evecs[j, k] := SparseMatrix(0, t, [], [], Rationals);
            od;
        else
            for k in [1..3] do
                rep.evecs[j, k] := false;
            od;
        fi;
    od;

    ## Embed dihedral algebras

    for i in Positions(rep.shape, "4B") do
        MAJORANA_EmbedDihedralAlgebra( i, rep, algebras.4B );
    od;

    for i in Positions(rep.shape, "6A") do
        MAJORANA_EmbedDihedralAlgebra( i, rep, algebras.6A );
    od;

    for i in PositionsProperty(rep.shape, x -> not x in [ "1A", "4B", "6A" ]) do
        MAJORANA_EmbedDihedralAlgebra( i, rep, algebras.(rep.shape[i]) );
    od;

    ## Finish off setup

    dim := Size(rep.setup.coords);

    rep.setup.nullspace := rec(     heads := [1 .. dim]*0,
                                    vectors := SparseMatrix( 0, dim, [], [], Rationals) );

    for i in rep.setup.orbitreps do
        for j in [1..3] do
            rep.evecs[i, j]!.ncols := dim;
            rep.evecs[i, j] := MAJORANA_BasisOfEvecs(rep.evecs[i, j]);
        od;
    od;

    for i in gens do MAJORANA_ExtendPerm( i, rep); od;

    MAJORANA_Orbitals( gens, t, rep.setup);

    for i in [1..Size(rep.setup.pairreps)] do

        if not IsBound(rep.algebraproducts[i]) then
            rep.algebraproducts[i] := false;
            rep.innerproducts[i] := false;
        elif rep.algebraproducts[i] <> false then
            rep.algebraproducts[i]!.ncols := dim;
        fi;
    od;

    return rep;

    end );

##
## Given the dihedral algebra generated by the axes in <rep.setup.pairreps[i]>,
## adds any new 2A, 3A, 4A or 5A axes to <rep.setup.coords> and <rep.setup.coordmap>
## and adds any new products and eigenvectors coming from the dihedral algebra.
##

InstallGlobalFunction( MAJORANA_EmbedDihedralAlgebra,

    function( i, rep, subrep )

    local   dim, t, gens, x, inv, elts, emb, j, im, orbit, y, k, sign;

    dim := Size(rep.setup.coords);
    t := Size(rep.involutions);

    x := rep.setup.pairreps[i];
    inv := rep.involutions{x};

    ## Add new basis vector(s) and their orbit(s) and extend pairconj and pairorbit matrices

    MAJORANA_AddNewVectors( rep, subrep, inv);

    ## Find the embedding of the subrep into the main algebra

    emb := MAJORANA_FindEmbedding( rep, subrep, inv );

    ## Add any new orbits

    gens := GeneratorsOfGroup(rep.group);
    gens := List( gens, g -> MAJORANA_FindPerm(g, rep, rep) );

    for j in [1 .. Size(subrep.setup.pairreps)] do

        im := emb{subrep.setup.pairreps[j]};

        if im[1] < 0 then im[1] := -im[1]; fi;
        if im[2] < 0 then im[2] := -im[2]; fi;

        orbit := rep.setup.pairorbit[im[1], im[2]];

        ## If need be, add a new orbit

        if orbit = 0 then
            MAJORANA_NewOrbital(im, gens, rep.setup);
        fi;
    od;

    ## Embed products and evecs

    MAJORANA_Embed( rep, subrep, emb );

    end );

##
## Finds the embedding of a dihedral algebra <subrep> into <rep>
##

InstallGlobalFunction( MAJORANA_FindEmbedding,

    function( rep, subrep, inv)

    local imgs, emb, pos, x, gens;

    ## Find the images of the embedding of the subrep into the main rep

    gens := GeneratorsOfGroup(subrep.group);

    imgs := List(subrep.setup.coords, w -> MAJORANA_MappedWord(rep, subrep, w, gens, inv) );

    emb := [];

    for x in imgs do
        pos := rep.setup.coordmap[x];
        if pos = fail then
            pos := rep.setup.coordmap[ Product( rep.involutions{x} )];
        fi;

        Add( emb, pos);
    od;

    return emb;

    end );

##
## If new vectors have been added to <setup.coords> then extend the action
## of an existing perm to these new vectors.
##

InstallGlobalFunction( MAJORANA_ExtendPerm,

    function(perm, rep)

    local dim, new_dim, i, im, sign, pos;

    new_dim := Size(rep.setup.coords);
    dim := Size(perm);

    for i in [dim + 1 .. new_dim] do
        im := perm{rep.setup.coords[i]};

        ## Keep track of sign changes

        sign := 1;

        if im[1] < 0 then im[1] := -im[1]; sign := -sign; fi;
        if im[2] < 0 then im[2] := -im[2]; sign := -sign; fi;

        if im[1] > im[2] then
            im := im{[2,1]};
        fi;

        ## Find the new vector in <setup.coords>

        pos := rep.setup.coordmap[im];

        if pos = fail then
            pos := rep.setup.coordmap[ Product( rep.involutions{im} ) ];
        fi;

        Add(perm, sign*pos);
    od;

    end);

##
## Adds any additional 2A, 3A, 4A or 5A vectors coming from the dihedral
## algebra <subrep>.
##

InstallGlobalFunction( MAJORANA_AddNewVectors,

    function(rep, subrep, inv)

    local i, list, list_5A, new, new_5A, x, vec, g, im, k, dim, gens;

    dim := Size(rep.setup.coords);
    gens := GeneratorsOfGroup(subrep.group);

    for i in [Size(subrep.involutions) + 1 .. Size(subrep.setup.coords)] do

        ## Find the new vectors to be added to <setup.coordmap>
        ## TODO - do we want to change this to hashmaps?

        list := Positions(subrep.setup.poslist, i);
        list_5A := Positions(subrep.setup.poslist, -i);

        new := []; new_5A := [];

        for x in subrep.setup.longcoords{list} do
            Add( new, MAJORANA_MappedWord(rep, subrep, x, gens, inv));
        od;

        for x in subrep.setup.longcoords{list_5A} do
            Add( new_5A, MAJORANA_MappedWord(rep, subrep, x, gens, inv));
        od;

        MAJORANA_AddConjugateVectors( rep, new, new_5A );
    od;

    for x in rep.setup.pairorbit do
        Append(x, [dim + 1 .. Size(rep.setup.coords)]*0 );
    od;

    for x in rep.setup.pairconj do
        Append(x, [dim + 1 .. Size(rep.setup.coords)]*0 );
    od;

    Append(rep.setup.pairorbit, NullMat( Size(rep.setup.coords) - dim , Size(rep.setup.coords) ));
    Append(rep.setup.pairconj, NullMat( Size(rep.setup.coords) - dim , Size(rep.setup.coords) ));

    for g in rep.setup.pairconjelts do  MAJORANA_ExtendPerm( g, rep); od;

    for g in rep.setup.conjelts do MAJORANA_ExtendPerm( g, rep); od;

    end );

##
## When adding a new set of vectors to <rep.setup.coordmap>, also adds
## all of their images under the group action.
##

InstallGlobalFunction( MAJORANA_AddConjugateVectors,

    function( rep, new, new_5A )

    local   vec, g, im, im_5A, k, elts, x;

    ## Find which (if any) new vectors are not yet in <setup.coordmap>

    vec := First(new, x -> x in rep.setup.coordmap);

    if vec <> fail then
        new := Filtered(new, x -> not x in rep.setup.coordmap);
        new_5A := Filtered(new_5A, x -> not x in rep.setup.coordmap);
    fi;

    if vec = fail and rep.axioms <> "NoAxioms" then
        vec := First(new, x -> Product( rep.involutions{x} ) in rep.setup.coordmap );

        if vec <> fail then
            new := Filtered(new, x -> not Product( rep.involutions{x} ) in rep.setup.coordmap);
            new_5A := Filtered(new_5A, x -> not Product( rep.involutions{x} ) in rep.setup.coordmap);
        fi;
    fi;

    if new = [] then return; fi;

    ## If vec = fail then all vectors are new and a new rep will be added
    ## to <setup.coords>. Otherwise, there are some new vectors but these
    ## are equal to an existing element of <setup.coords>

    for g in rep.setup.pairconjelts do

        im := List(new, x -> SortedList( g{ x } ));
        im := Filtered( im, x -> not x in rep.setup.coordmap );

        if rep.axioms = "AllAxioms" then
            im := Filtered( im, x -> not Product( rep.involutions{x} ) in rep.setup.coordmap);
        fi;

        if im <> [] then

            im_5A := List(new_5A, x -> SortedList( g{ x } ));

            ## If need be, add a new vector to coords, otherwise find index of existing vector

            if vec = fail then
                Add( rep.setup.coords, im[1] );
                k := Size(rep.setup.coords);
            else
                k := rep.setup.coordmap[ SortedList( g{ vec } )];
                if k = fail then
                    k := rep.setup.coordmap[ Product(rep.involutions{ g{ vec } }) ];
                fi;
            fi;

            ## Add new vectors to <setup.coordmap>

            for x in im do rep.setup.coordmap[ x ] := k; od;
            for x in im_5A do rep.setup.coordmap[ x ] := -k; od;

            if rep.axioms = "AllAxioms" then
                for x in im do
                    rep.setup.coordmap[ Product( rep.involutions{ x } )] := k;
                od;
                for x in im_5A do
                    rep.setup.coordmap[ Product( rep.involutions{ x } )] := -k;
                od;
            fi;
        fi;
    od;

    end );

##
## Here <g> is either a group element or a homomorphism from <subrep.group>
## to <rep.group>. In the first case, <subrep> = <rep> and the func
## return <g> as a permutation on <rep.setup.coords>. In the second case,
## returns <g> as a list sending <subrep.setup.coords> to <rep.setup.coords>.
##

InstallGlobalFunction( MAJORANA_FindPerm,

    function(g, rep, subrep)

    local   dim, i, list, im, sign, vec;

    dim := Size(subrep.setup.coords);
    list := [1..dim]*0;

    for i in [1..dim] do

        vec := subrep.setup.coords[i];

        if IsRowVector(vec) then

            im := list{vec};

            sign := 1;

            if im[1] < 0 then sign := -sign; im[1] := -im[1]; fi;
            if im[2] < 0 then sign := -sign; im[2] := -im[2]; fi;

            if im[1] > im[2] then im := im{[2,1]}; fi;

            list[i] := rep.setup.coordmap[ im ];

            if list[i] = fail then
                list[i] := rep.setup.coordmap[ Product( rep.involutions{im} ) ];
            fi;
        else
            list[i] := rep.setup.coordmap[ rep.involutions[i]^g ];
        fi;
    od;

    return list;

    end);

##
## Optional function for use by the user after calling <ShapesOfMajoranaRepresentation>
##

InstallGlobalFunction( MAJORANA_RemoveDuplicateShapes,

    function(input)

    local autgp, inner_autgp, outer_auts, perm, g, i, pos, im;

    autgp := AutomorphismGroup(input.group);
    inner_autgp := InnerAutomorphismsAutomorphismGroup(autgp);
    outer_auts := [];

    for g in RightTransversal(autgp, inner_autgp) do
        if AsSet(OnTuples(input.involutions, g)) = AsSet(input.involutions) then
            perm := [];

            for i in [1..Size(input.orbitals)] do

                im := OnPairs(Representative(input.orbitals[i]), g);

                pos := PositionProperty(input.orbitals, x -> im in x);

                if pos = fail then pos := PositionProperty(input.orbitals, x -> Reversed(im) in x); fi;

                if pos = fail then Error(); fi;

                Add(perm, pos);
            od;

            Add(outer_auts, perm);
        fi;
    od;

    for i in [1..Size(input.shapes)] do
        if IsBound(input.shapes[i]) then
            for g in outer_auts do
                pos := Position(input.shapes, input.shapes[i]{g});
                if pos <> fail and pos <> i then
                    Unbind(input.shapes[pos]);
                fi;
            od;
        fi;
    od;

    input.shapes := Compacted(input.shapes);

    end );

##
## <MappedWord> for indices or pairs of indices referring to elements of coords
##

InstallGlobalFunction(MAJORANA_MappedWord,

    function(rep, subrep, w, gens, imgs)

    local im;

    if IsRowVector(w) then
        im := List(w, i -> MappedWord(subrep.setup.coords[i], gens, imgs));

        return SortedList(List(im, x -> Position(rep.involutions, x )));
    else
        return Position(rep.involutions, MappedWord(w, gens, imgs) );
    fi;

    end );

##
## This is used only the the main algorithm to record algebra products once found
##

InstallGlobalFunction(SP_Inverse,

    function(perm)

    local l, inv, i;

    if perm = [] then return []; fi;

    l := Length(perm);

    inv := [1..l];

    for i in [1..l] do
        if perm[i] > 0 then
            inv[perm[i]] := i;
        else
            inv[-perm[i]] := -i;
        fi;
    od;

    return inv;

    end);

##
## Not currently in use, here for reference
##

InstallGlobalFunction(SP_Product,

    function( perm1, perm2) # Perms must be of same length!

    local prod, i;

    prod := [];

    for i in perm1 do
        if i > 0 then
            Add(prod, perm2[i]);
        else
            Add(prod, -perm2[-i]);
        fi;
    od;

    return prod;

    end );

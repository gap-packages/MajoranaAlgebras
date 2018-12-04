##
## Not currently used - maximal subgroups is generally faster
##

InstallGlobalFunction( "MAJORANA_AllEmbeddings",

    function(rep)

    local   x, gens, subgp, T, subrep, ex, i, ev;

    # Loop over the algebraproducts that are not yet known
    for x in Positions(rep.algebraproducts, false) do
        if rep.algebraproducts[x] = false then

            # Find the involutions that generate the subgroup that contains the vectors in this product
            gens := ShallowCopy(rep.setup.pairreps[x]);

            while ForAny(gens, i -> i > Size(rep.involutions)) do
                gens := Flat( rep.setup.coords{gens});
            od;

            subgp := Subgroup(rep.group, rep.involutions{gens});

            # If this subgroup is proper
            if Size(subgp) < Size(rep.group) then

                # If this subgroup is generated by involutions of the main rep
                T := Intersection(subgp, rep.involutions);
                if T <> [] and Size(Group(T)) = Size(subgp) then

                    ex := ShapesOfMajoranaRepresentationAxiomM8(subgp,T);

                    # If the shapes of the reps match then build the rep
                    for i in [1..Size(ex.shapes)] do
                        if IsSubsetSet(AsSet(rep.shape), AsSet(ex.shapes[i])) then

                            Info(   InfoMajorana, 10,
                                    STRINGIFY("Constructing subrep of ", StructureDescription(subgp) ) );

                            subrep := MajoranaRepresentation(ex,i,rep.axioms);

                            # Embed the constructed rep
                            MAJORANA_EmbedKnownRep(rep, subrep);
                        fi;
                    od;
                fi;
            fi;
        fi;
    od;

    # Find bases of the eigenvectors
    for i in rep.setup.orbitreps do
        for ev in RecNames(rep.evec[i]) do
            if Nrows(rep.evecs[i].(ev)) > 0 then
                rep.evecs[i].(ev) := ReversedEchelonMatDestructive(rep.evecs[i].(ev)).vectors;
            fi;
        od;
    od;

    end );

##
## Attempts to construct Majorana representations of the maximal subgroups of
## G and embeds the results into the incomplete algebra
##

InstallGlobalFunction( "MAJORANA_MaximalSubgps",

    function(rep,axioms)

    local   max, inv, i, j, ev, ex, subrep;

    # Find all maximal subgroups and their intersections with T
    max := MaximalSubgroupClassReps(rep.group);
    max := Filtered(max, x -> Size(x) > 12);

    inv := List(max, x -> Intersection(AsList(x), rep.involutions));
    inv := Filtered(DuplicateFreeList(inv), x -> x <> []);

    max := List(inv, Group);

    # Loop over these maximal subgroups
    for i in [1..Size(max)] do
        ex := ShapesOfMajoranaRepresentationAxiomM8(max[i], inv[i]);
        # Loop over all possible shapes
        for j in [1..Size(ex.shapes)] do
            # If the shapes match
            if IsSubsetSet(AsSet(rep.shape), AsSet(ex.shapes[j])) then

                Info(   InfoMajorana, 10,
                        STRINGIFY("Constructing subrep of ", StructureDescription(ex.group) ) );

                # Construct the subrep and embed it
                subrep := MajoranaRepresentation(ex,j,axioms);
                MAJORANA_EmbedKnownRep(rep, subrep);
            fi;
        od;
    od;

    # Find bases for the eigevectors and the nullspace
    for i in rep.setup.orbitreps do
        for ev in RecNames(rep.evecs[i]) do
            rep.evecs[i].(ev) := ReversedEchelonMatDestructive(rep.evecs[i].(ev)).vectors;
        od;
    od;

    if Nrows(rep.setup.nullspace.vectors) > 0 then
        rep.setup.nullspace := ReversedEchelonMatDestructive(rep.setup.nullspace.vectors);
    fi;

    end );

##
## For a given injective homomorphism <emb> from subrep.group to rep.group,
## if possible, find an automorphism of subrep.group such that the shapes and
## involutions of rep and subrep match. Otherwise, return false.
##

InstallGlobalFunction( "MAJORANA_CheckEmbedding",

    function(rep, subrep, emb)

    local   aut, g, aut_emb, im, i, x, pos, k;

    aut := AutomorphismGroup(subrep.group);

    # For each automorphism of subrep.group
    for g in aut do

        # Compose the automorphism with the embedding
        aut_emb := CompositionMapping2(emb, g);

        # Check if the image of the involutions of the subrep is a subset of the
        # involutions of the rep
        im := AsSet(Image(aut_emb, subrep.involutions));
        if not IsSubsetSet(AsSet(rep.involutions), im) then
            return false;
        fi;

        # Check if the shape of the image of the subrep matches with the shape
        # of the rep
        for i in [1..Size(subrep.shape)] do
            if subrep.shape[i, 1] in ['2','3','4'] then

                x := subrep.setup.pairreps[i];
                # TODO we might have to fix this if the subrep is of dim 0

                im := OnPairs(subrep.involutions{x}, aut_emb);

                pos := List(im, x -> Position(rep.involutions, x));

                k := UnorderedOrbitalRepresentative(rep.setup.orbitalstruct, pos);

                if subrep.shape[i] <> rep.shape[k] then
                    return false;
                fi;
            fi;
        od;

        # If both shape and involutions match, return the automorphism
        return g;
    od;

    return false;

    end);

InstallGlobalFunction( "MAJORANA_EmbedKnownRep",

    function( rep, subrep)

    local   embs,
            i,
            g,
            emb;

    embs := IsomorphicSubgroups(rep.group, subrep.group);

    # For each embedding of subrep.group into rep.group
    for i in [1..Size(embs)] do
        g := MAJORANA_CheckEmbedding(rep, subrep, embs[i]);

        # Use the group embedding to embed the algebras
        if g <> false then
            emb := CompositionMapping2(embs[i], g);

            MAJORANA_Embed(rep, subrep,  emb);
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

    function(emb, rep, subrep)

    local   dim, i, list, im, sign, vec;

    dim := Size(subrep.setup.coords);
    list := [1..dim]*0;

    # Loop over the spanning set of the subrepresentation
    for i in [1..dim] do
        vec := subrep.setup.coords[i];

        # If spanning set element is a row vector then use the images we have already calculated
        if IsRowVector(vec) then

            im := list{vec};

            # Adjust the sign
            sign := 1;

            if im[1] < 0 then sign := -sign; im[1] := -im[1]; fi;
            if im[2] < 0 then sign := -sign; im[2] := -im[2]; fi;

            if im[1] > im[2] then im := im{[2,1]}; fi;

            list[i] := rep.setup.coordmap[ im ];

            if list[i] = fail then
                list[i] := rep.setup.coordmap[ Product( rep.involutions{im} ) ];
            fi;
        else
            # Otherwise, use the images of the group elements
            list[i] := rep.setup.coordmap[ subrep.involutions[i]^emb ];
        fi;
    od;

    return list;

    end);

InstallGlobalFunction( "MAJORANA_Embed",

    function(rep, subrep, emb)

    local   i, im, ev, k, g, v, sign, x, l, temp;

    if not IsRowVector(emb) then
        emb := MAJORANA_FindPerm(emb, rep, subrep);
    fi;

    for i in [1..Size(subrep.algebraproducts)] do

        if subrep.setup.pairreps[i] <> fail then
            sign := 1;

            # Find the image of the pair representative and adjust the sign
            im := emb{subrep.setup.pairreps[i]};
            if im[1] < 0 then sign := -sign; im[1] := -im[1]; fi;
            if im[2] < 0 then sign := -sign; im[2] := -im[2]; fi;

            # Find the corresponding pair orbit in rep
            g := UnorderedOrbitalCanonizingElement(rep.setup.orbitalstruct, im);
            k := rep.setup.pairrepsmap[ OnPairs( im, g ) ];

            if k < 0 then sign := -sign; k := -k; fi;

            g := ListSignedPerm( Inverse(g), Size(rep.setup.coords));

            # Record the new algebraproduct
            if not IsBound(rep.algebraproducts[k]) or rep.algebraproducts[k] = false then
                if subrep.algebraproducts[i] <> false then
                    v := MAJORANA_ImageMat(subrep.algebraproducts[i], emb, rep, subrep);
                    rep.algebraproducts[k] := sign*MAJORANA_ConjugateVec(v,g);
                fi;
            fi;

            # Record the new inner product
            if IsBound(rep.innerproducts) then
                if not IsBound(rep.innerproducts[k]) or rep.innerproducts[k] = false then
                    if subrep.innerproducts[i] <> false then
                        rep.innerproducts[k] := sign*subrep.innerproducts[i];
                    fi;
                fi;
            fi;
        fi;
    od;

    # Record the new eigevectors
    for i in subrep.setup.orbitreps do

        k := emb[i];

        # Use the inverse of the conjugating element
        g := SP_Inverse(rep.setup.conjelts[k]);

        # Record new eigenvectors
        for ev in RecNames(subrep.evecs[i]) do
            if Nrows(subrep.evecs[i].(ev)) > 0 then
                im := MAJORANA_ImageMat(subrep.evecs[i].(ev), emb, rep, subrep);
                for v in Iterator(im) do
                    v := MAJORANA_ConjugateVec(v, g);
                    rep.evecs[g[k]].(ev) := UnionOfRows(rep.evecs[g[k]].(ev), v);
                od;
            fi;
        od;
    od;

    # Record new nullspace vectors
    if Nrows( subrep.setup.nullspace.vectors ) > 0 then
        im := MAJORANA_ImageMat(subrep.setup.nullspace.vectors, emb, rep, subrep);
        rep.setup.nullspace.vectors := UnionOfRows(rep.setup.nullspace.vectors, im);
    fi;

    end );

##
## Finds the image of the matrix <mat> under the embedding <emb>
##

InstallGlobalFunction( "MAJORANA_ImageMat",

    function(mat, emb, rep, subrep)

    local   res, i, j, sign, im;

    res := SparseZeroMatrix(mat!.nrows, mat!.ncols, Rationals);

    # Loop over all nonzero matrix elements
    for i in [1..mat!.nrows] do
        for j in [1..Size(mat!.indices[i])] do

            # Find the image
            im := emb[mat!.indices[i, j]];

            # Adjust the sign
            sign := 1;;
            if im < 0 then sign := -sign; im := -im; fi;

            # Record the image
            SetEntry(res, i, im, sign*mat!.entries[i, j]);
        od;
    od;

    return res;

    end );

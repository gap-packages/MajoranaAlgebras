# Compute an orbital structure which has the
# following properties
#
# Input: G acting on Omega via Act

# Output: Orbital Structure
#         - A set of representatives of orbitals
#         - For any pair (i,j) efficiently
#           determine which orbital it belongs to
#         - For any pair (i,j) efficiently
#           compute an element g in G that takes (i,j)
#           to a representative in the set of
#           representatives
#         - iterate over orbit given a rep
#         - iterate over a transversal for
#           a stabilizer of a pair

InstallGlobalFunction(MAJORANA_OrbitalStructure,
# gens  - generators of a group acting on Omega
# Omega - the domain
# Act   - action of Group(gens) on Omega
function(gens, Omega, Act)
    local o, so, i, res;

    # Result will be an orbital structure that allows
    # some stuff to be done with orbitals
    res := rec( gens := gens, group := Group(gens), act := Act );

    # Orbits. Currently we'll just choose the
    # first element in each orbit as orbit rep
    res.orbits := Orbits( res.group, Omega, Act );

    # This is so we're able to determine which orbit the
    # first point of a pair is in, and then get an element
    # of G that maps said point to the orbit rep
    # Really what we want is a schreier trees/vectors here
    res.orbreps := [];
    res.orbnums := HashMap(Size(Omega));
    res.orbstabs := [];
    for o in [1..Length(res.orbits)] do
        Add(res.orbreps, res.orbits[o][1]);
        res.orbstabs[o] := rec( );
        res.orbstabs[o].stab := Stabilizer(res.group, res.orbits[o][1], Act);
        res.orbstabs[o].orbs := Orbits(res.orbstabs[o].stab, Omega, Act);
        res.orbstabs[o].orbnums := HashMap(Size(Omega));
        res.orbstabs[o].orbreps := [];
        for so in [1..Length(res.orbstabs[o].orbs)] do
            Add(res.orbstabs[o].orbreps, res.orbstabs[o].orbs[so][1]);
            for i in res.orbstabs[o].orbs[so] do
                res.orbstabs[o].orbnums[i] := so;
            od;
        od;

        for i in res.orbits[o] do
            res.orbnums[i] := o;
        od;
    od;

    return res;
end);

InstallGlobalFunction(MAJORANA_OrbitalRep,
function(os, pair)

    # Returns a representative (as a pair of elements of the G set) for the
    # orbital that contains <pair>.

    local fo, so, p;
    fo := os.orbnums[pair[1]];
    p := RepresentativeAction(os.group, pair[1], os.orbreps[fo], os.act );
    so := os.orbstabs[fo].orbnums[os.act(pair[2], p)];
    return [ os.orbreps[fo], os.orbstabs[fo].orbreps[so] ];
end);

InstallGlobalFunction(MAJORANA_OrbitalReps,
function(os)
    return Union(List(os.orbreps, i -> [i,i])
                , ListX( [1..Length(os.orbreps)]
                       , k -> os.orbstabs[k].orbreps
                       , {x,y} -> [os.orbreps[x],y] ) );
end);


# TODO: fix this.
InstallGlobalFunction(MAJORANA_OrbitalCanonizingElement,
function(os, pair)

    # Returns a group elements that maps <pair> to its orbital representative
    # (as given by MAJORANA_OrbitalRep).

    local fo, so, p1, p2;

    fo := os.orbnums[pair[1]];
    p1 := RepresentativeAction(os.group, pair[1], os.orbreps[fo], os.act);
    so := os.orbstabs[fo].orbnums[os.act(pair[2], p1)];
    p2 := RepresentativeAction( os.orbstabs[fo].stab, os.act(pair[2], p1)
                               , os.orbstabs[fo].orbreps[so], os.act);

    return p1 * p2;
end);

InstallGlobalFunction(MAJORANA_OrbitalCanonizingElementInverse,
function(os, pair)
    return MAJORANA_OrbitalCanonizingElement(os, pair)^-1;
    # Returns a group elements that maps the orbital representative of <pair>
    # to <pair> itself. This will be the inverse of the output of
    # MAJORANA_OrbitalCanonizingElement( os, pair )
end);

# Acting on sets of size 2
InstallGlobalFunction(MAJORANA_UnorderedOrbitalRep,
function(os, p)
    local a, b, oa, ob, r1, r2, p1, p2, tmp, tmp2;

    a := p[1];
    b := p[2];

    oa := os.orbnums[a];
    ob := os.orbnums[b];

    r1 := os.orbreps[oa];
    r2 := os.orbreps[ob];

    if r1 = r2 then
    # a and b are in the same orbit
        p1 := RepresentativeAction(os.group, a, r1, os.act);
        p2 := RepresentativeAction(os.group, b, r1, os.act);
        tmp := [r1, os.act(b,p1)];
        tmp2 := [r1, os.act(a, p2)];
        return Minimum(tmp,tmp2);
    elif r2 < r1 then
    # b has the smaller orbrep, i.e. in the orbitalrep
    # has the smaller representative. We swap roles of a and b
            tmp := b; b := a; a := tmp;
        tmp := r2; r2 := r1; r1 := tmp;
        tmp := ob; ob := oa; oa := tmp;
    fi;

    # Move a to the smaller rep
    p1 := RepresentativeAction(os.group, a, r1, os.act);
    # Now look in the point stabiliser of r1 what
    # element we can map b^p1 to
    ob := os.orbstabs[oa].orbnums[os.act(b, p1)];
    return [ r1, os.orbstabs[oa].orbreps[ob]];
end);

InstallGlobalFunction(MAJORANA_UnorderedOrbitalCanonizingElement,
function(os, pair)

    # Returns a group elements that maps <pair> to its orbital representative
    # (as given by MAJORANA_OrbitalRep).

    local a, b, oa, ob, r1, r2, p1, p2, tmp, tmp2;

    a := pair[1];
    b := pair[2];

    oa := os.orbnums[a];
    ob := os.orbnums[b];

    r1 := os.orbreps[oa];
    r2 := os.orbreps[ob];

    if r1 = r2 then
        p1 := RepresentativeAction(os.group, a, r1, os.act);
        p2 := RepresentativeAction(os.group, b, r1, os.act);
        tmp := [r1, os.act(b,p1)];
        tmp2 := [r1, os.act(a, p2)];
        if tmp < tmp2 then
            return p1;
        else
            return p2;
        fi;
    elif r2 < r1 then
        tmp := b; b := a; a := tmp;
        tmp := r2; r2 := r1; r1 := tmp;
        tmp := ob; ob := oa; oa := tmp;
    fi;

    p1 := RepresentativeAction(os.group, a, r1, os.act);
    b := os.act(b, p1);
    ob := os.orbstabs[oa].orbnums[b];
    p2 := RepresentativeAction(os.orbstabs[oa].stab, b, os.orbstabs[oa].orbreps[ob], os.act);

    return p1 * p2;
end);

InstallGlobalFunction(MAJORANA_UnorderedOrbitalCanonizingElementInverse,
     {os, pair} -> MAJORANA_UnorderedOrbitalCanonizingElement(os, pair) ^ -1);

# TODO: fix this
InstallGlobalFunction(MAJORANA_UnorderedOrbitalReps,
function(os)
    local ordered_reps, unordered_reps;

    ordered_reps := MAJORANA_OrbitalReps(os);
    unordered_reps := List(ordered_reps, x -> MAJORANA_UnorderedOrbitalRep(os, x));

    return DuplicateFreeList(unordered_reps);
end);

InstallGlobalFunction(MAJORANA_OrbitalTransversalIterator,
function( os, rep )
    local r, fo, so;

    # Make sure we have *the* rep, not *a* rep
    rep := MAJORANA_OrbitalRep(os, rep);

    r := rec( orb := HashMap()
            , new := [ [ rep, [] ] ]
            , NextIterator := function(iter)
                local i, pntp, pnt, npntp, npnt;

                pntp := Remove(iter!.new, 1);
                pnt := pntp[1];

                for i in [1..Length(os.gens)] do
                    npnt := OnTuples(pnt, os.gens[i]);
                    if not npnt in iter!.orb then
                        npntp := [ npnt, Concatenation(pntp[2], [i]) ];
                        iter!.orb[npnt] := npntp[2];
                        Add(iter!.new, npntp);
                    fi;
                od;
                if Length(pntp[2]) = 0 then
                    return One(os.group);
                else
                    return Product(List(pntp[2], i -> os.gens[i]));
                fi;
            end
            , IsDoneIterator := iter -> iter!.new = []
            , ShallowCopy := iter -> rec( orb := StructuralCopy(iter!.orb)
                                        , new := ShallowCopy(iter!.new) )
            );
    r.orb[rep] := [];
    return IteratorByFunctions(r);
end);

# For now, a disguised orbit algorithm which is probably better
# than computing RepresentativeAction all the time!
InstallGlobalFunction(MAJORANA_UnorderedOrbitalTransversalIterator,
function( os, rep )
    local r, fo, so;

    # Make sure we have *the* rep, not *a* rep
    rep := MAJORANA_UnorderedOrbitalRep(os, rep);

    r := rec( orb := HashMap()
            , new := [ [ rep, [] ] ]
            , NextIterator := function(iter)
                local i, pntp, pnt, npntp, npnt;

                pntp := Remove(iter!.new, 1);
                pnt := pntp[1];

                for i in [1..Length(os.gens)] do
                    npnt := Set(pnt, x -> os.act(x, os.gens[i]));
                    if not npnt in iter!.orb then
                        npntp := [ npnt, Concatenation(pntp[2], [i]) ];
                        iter!.orb[npnt] := npntp[2];
                        Add(iter!.new, npntp);
                    fi;
                od;
                if Length(pntp[2]) = 0 then
                    return One(os.group);
                else
                    return Product(List(pntp[2], i -> os.gens[i]));
                fi;
            end
            , IsDoneIterator := iter -> iter!.new = []
            , ShallowCopy := iter -> rec( orb := StructuralCopy(iter!.orb)
                                        , new := ShallowCopy(iter!.new) )
            );
    r.orb[rep] := [];
    return IteratorByFunctions(r);
end);

BindGlobal("UnorderedOrbitalTest",
function(os, domain)
    local o, orbs, reps, muoreps;
    orbs := Orbits(os.group, Combinations(domain, 2), OnSets);
    reps := List(orbs, Minimum);

    muoreps := MAJORANA_UnorderedOrbitalReps(os);
    for o in domain do
        if o in os.orbreps then
            if not ([o,o] in muoreps) then
                Error("The element ", o, " is an orbit representative, but its diagonal is not a representative");
            fi;
            Remove(muoreps, Position(muoreps, [o,o]));
        else
            if [o,o] in muoreps then
                Error("The element ", o, " is not an orbit representative, but its diagonal is a representative");
            fi;
        fi;
    od;
    if not IsSet(muoreps) then
        Error("Representatives do not form a set");
    fi;
    if not Set(muoreps) = Set(reps) then
        Print("Elements in muoreps that are not in reps: ", Difference(muoreps, reps), "\n");
        Print("Elements in reps that are not in muoreps: ", Difference(reps, muoreps), "\n");
        Error("differences between orbital reps and reps");
    fi;
    return true;
end);

BindGlobal("UnorderedOrbitalTransversalTest",
function(os, domain)
    local o, orbs, r, reps, i, iter, e, p;

    orbs := Orbits(os.group, Combinations(domain, 2), OnSets);
    reps := List(orbs, Minimum);

    for r in reps do
        o := ShallowCopy(Filtered(orbs, x -> r in x)[1]);
        iter := MAJORANA_UnorderedOrbitalTransversalIterator(os, r);
        for i in iter do
            e := OnSets(r, i);
            p := Position(o, e);
            if p = fail then
                Print("element ", e, " = ", r, "^", i, " not found\n");
            else
                Remove(o, p);
            fi;
        od;
        if not IsEmpty(o) then
            Error("leftover elements in orbit: ", o, " with rep ", r, "\n");
        fi;
    od;
    return true;
end);

BindGlobal("UnorderedOrbitalCanonizingTest",
function(os, domain)
    local o, orbs, r, rep, i, iter, e, p, can;

    orbs := Orbits(os.group, Combinations(domain, 2), OnSets);
    for o in orbs do
        for i in o do
            rep := MAJORANA_UnorderedOrbitalRep(os, i);
            can := MAJORANA_UnorderedOrbitalCanonizingElement(os, i);
            if rep <> OnSets(i, can) then;
                Error("element ", i, " is not canonized by ", can, "\n");
            fi;
        od;
    od;
    return true;
end);

BindGlobal("OrbitalTest",
function(os, domain)
    local o, orbs, reps, muoreps;
    orbs := Orbits(os.group, Tuples(domain, 2), OnTuples);
    reps := List(orbs, Minimum);

    muoreps := MAJORANA_OrbitalReps(os);
    if not IsSet(muoreps) then
        Error("Representatives do not form a set");
    fi;
    if not Set(muoreps) = Set(reps) then
        Print("Elements in muoreps that are not in reps: ", Difference(muoreps, reps), "\n");
        Print("Elements in reps that are not in muoreps: ", Difference(reps, muoreps), "\n");
        Error("differences between orbital reps and reps");
    fi;
    return true;
end);

BindGlobal("OrbitalTransversalTest",
function(os, domain)
    local o, orbs, r, reps, i, iter, e, p;

    orbs := Orbits(os.group, Arrangements(domain, 2), OnTuples);
    reps := List(orbs, Minimum);

    for r in reps do
        o := ShallowCopy(Filtered(orbs, x -> r in x)[1]);
        iter := MAJORANA_OrbitalTransversalIterator(os, r);
        for i in iter do
            e := OnTuples(r, i);
            p := Position(o, e);
            if p = fail then
                Print("element ", e, " = ", r, "^", i, " not found\n");
                Error("");
            else
                Remove(o, p);
            fi;
        od;
        if not IsEmpty(o) then
            Error("leftover elements in orbit: ", o, " with rep ", r, "\n");
        fi;
    od;
    return true;
end);

BindGlobal("OrbitalCanonizingTest",
function(os, domain)
    local o, orbs, r, rep, i, iter, e, p, can;

    orbs := Orbits(os.group, Combinations(domain, 2), OnTuples);
    for o in orbs do
        for i in o do
            rep := MAJORANA_OrbitalRep(os, i);
            can := MAJORANA_OrbitalCanonizingElement(os, i);
            if rep <> OnTuples(i, can) then;
                Print("element ", i, " is not canonized by ", can, "\n");
            fi;
        od;
    od;
    return true;
end);

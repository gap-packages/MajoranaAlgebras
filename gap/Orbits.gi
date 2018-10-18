InstallGlobalFunction(MAJORANA_Orbits,

    function(gens, t, setup)
    
    local   conjelts,
            orbitreps,
            i, j,
            orb,
            gen,
            elts,
            count,
            dim,
            p,
            q,
            h,
            g;

    conjelts := [1..t]*0;
    orbitreps := [];
    
    dim := Size(setup.coords);
    
    for i in [1..t] do 
        if conjelts[i] = 0 then 
            
            Add(orbitreps, i);
            conjelts[i] := [1..dim];
            
            orb := [i];
            elts := [[1..dim]];
            
            count := 0;
            
            for p in orb do 
            
                count := count + 1;
                h := elts[count];
                
                for gen in gens do 
                    q := gen[p];
                    
                    g := [];
    
                    for j in h do 
                        if j > 0 then 
                            Add(g, gen[j]);
                        else
                            Add(g, -gen[-j]);
                        fi;
                    od;
                    
                    if q < 0 then q := -q; fi;
                    
                    if conjelts[q] = 0 then 
                        Add(orb, q);
                        Add(elts, g);
                        conjelts[q] := g;
                    fi;
                od;
            od;
        fi;
    od;
    
    conjelts := DuplicateFreeList(conjelts);
    
    return rec( conjelts := conjelts,
                orbitreps := orbitreps  );
                        
    end ); 
   
InstallGlobalFunction(MAJORANA_Orbitals,

    function(gens,t,setup)
    
    local   dim, i, j, orb;
    
    dim := Size(setup.coords);

    for i in [1..dim] do 
        for j in [Maximum(i,t + 1)..dim] do 

            if setup.pairorbit[i][j] = 0 then 
                
                orb := MAJORANA_NewOrbital([i,j], gens, setup);
                
                if IsBound(setup.orbitals) then 
                    Add(setup.orbitals, Immutable(orb));
                fi;
            fi;
        od;
    od;
    
    end );


InstallGlobalFunction( MAJORANA_NewOrbital,

    function( pnt, gens, setup)

    local orb, elts, count, y, p, h, q, z, g, i, j, k, pos, im, new, old, dim, sign, gen;

    # New Orbital, record as a new representative
    Add(setup.pairreps, pnt);

    # Number of points we're acting on
    dim := Size(setup.coords);

    # Orbit contains the point we're starting with
    orb := [ pnt ];
    # What does elts do?
    elts := [ [1..dim] ];

    count := 0;

    # The number of the orbit we're enumerating
    y := Size(setup.pairreps);

    setup.pairorbit[pnt[1]][pnt[2]] := y;
    setup.pairorbit[pnt[2]][pnt[1]] := y;

    # The element that's conguating is the identity
    setup.pairconj[pnt[1]][pnt[2]] := 1;
    setup.pairconj[pnt[2]][pnt[1]] := 1;

    # Orbit enumeration?
    for p in orb do

        count := count + 1;
        h := elts[count];

        # Orbit enumeration, every generator
        for gen in gens do

            # Generator applied to p?
            q := gen{p};

            # Make a permutation out of a signed permutation
            if q[1] < 0 then q[1] := -q[1]; fi;
            if q[2] < 0 then q[2] := -q[2]; fi;


            # the orbit in which this point lies
            z := setup.pairorbit[q[1]][q[2]];

            # Not detected yet
            if z = 0 then

                g := [];

                for i in h do
                    if i > 0 then
                        Add(g, gen[i]);
                    else
                        Add(g, -gen[-i]);
                    fi;
                od;

                Add( orb, q );
                Add( elts, g);

                if Product(g{orb[1]}) < 0 then
                    sign := -1;
                else
                    sign := 1;
                fi;

                setup.pairorbit[q[1]][q[2]] := sign*y;
                setup.pairorbit[q[2]][q[1]] := sign*y;

                pos := Position(setup.pairconjelts, g);

                if pos = fail then
                    Add(setup.pairconjelts, g);
                    pos := Size(setup.pairconjelts);
                fi;

                # The element that takes this pair to
                # the representative
                setup.pairconj[q[1]][q[2]] := pos;
                setup.pairconj[q[2]][q[1]] := pos;
            fi;
        od;
    od;

    return orb;

    end );

# Compute an orbital structure which has the
# following properties
#
# Input: G acting on Omega via Act
#        (though, for the time being,
#         G <= S_n acting on [1..n] OnPoints)
#
# Output: Orbital Structure
#         - A set of representatives of orbitals
#         - For any pair (i,j) efficiently
#           determine which orbital it belongs to
#         - For any pair (i,j) efficiently
#           compute an element g in G that takes (i,j)
#           to a representative in the set of
#           representatives

InstallGlobalFunction(MAJORANA_OrbitalStructure,
# gens  - generators of the permrep on C
#         as lists
# t     - starting point for second component
# setup - the setup structure
function(gens, Omega, Act)
    local o, so, i, res;

    # Result will be an orbital structure that allows
    # some stuff to be done with orbitals
    res := rec( group := Group(gens) );

    # Orbits. Currently we'll just choose the
    # first element in each orbit as orbit rep
    res.orbits := Orbits( res.group, Omega, Act );

    # This is so we're able to determine which orbit the
    # first point of a pair is in, and then get an element
    # of G that maps said point to the orbit rep
    # Really what we want is a schreier trees/vectors here
    res.orbreps := [];
    res.orbnums := [];
    res.orbstabs := [];
    for o in [1..Length(res.orbits)] do
        Add(res.orbreps, res.orbits[o][1]);
        res.orbstabs[o] := rec( );
        res.orbstabs[o].stab := Stabilizer(res.group, res.orbits[o][1], Act);
        res.orbstabs[o].orbs := Orbits(res.orbstabs[o].stab, Omega, Act);
        res.orbstabs[o].orbnums := [];
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

InstallGlobalFunction(MAJORANA_OrbitalStructureSigned,
function(gens, Omega, Act)
    local res;

    res := MAJORANA_OrbitalStructure(List(gens, x->x![1]), Omega, Act);
    res.signedgens := gens;
    return res;
end);

InstallGlobalFunction(MAJORANA_OrbitalRep,
function(os, pair)
    local fo, so, p;
    fo := os.orbnums[pair[1]];
    p := RepresentativeAction(os.group, pair[1], os.orbreps[fo] );
    so := os.orbstabs[fo].orbnums[pair[2]^p];
    return [ os.orbreps[fo], os.orbstabs[fo].orbreps[so] ];
end);

InstallGlobalFunction(MAJORANA_OrbitalRepAct,
function(os, pair)
    local fo, so, p1, p2;

    fo := os.orbnums[pair[1]];
    p1 := RepresentativeAction(os.group, pair[1], os.orbreps[fo]);
    so := os.orbstabs[fo].orbnums[pair[2]^p1];
    p2 := RepresentativeAction(os.orbstabs[fo].stab, pair[2]^p1, os.orbstabs[fo].orbreps[so]);

    return p1 * p2;
end);

# getting [i,j], we union the orbits [i,j] and [j,i]
# -> compute orbreps of both, take the smaller one
InstallGlobalFunction(MAJORANA_OrbitalRepUnion,
function(os, p)
    local r1, r2;

    r1 := MAJORANA_OrbitalRep(os, p);
    r2 := MAJORANA_OrbitalRep(os, p{[2,1]});

    if r1 < r2 then
        return r1;
    else
        return r2;
    fi;
end);

InstallGlobalFunction(MAJORANA_OrbitalRepUnions,
function(os)
    local reps, reps2, p, q;

    reps := Set(Union( List( [1..Length(os.orbreps)]
                       , k -> ListX(os.orbreps, os.orbstabs[k].orbreps
                                    , {x,y} -> MAJORANA_OrbitalRepUnion(os, [x,y]) ) ) ) );
    return reps;
end);

# Still have to try both [i,j] and [j,i]
# And we might want to think about caching.
InstallGlobalFunction(MAJORANA_OrbitalRepActSigned,
function(os, pair)
    local ra, fact;

    ra := MAJORANA_OrbitalRepAct(os, pair);
    fact := Factorization(os.group, ra);

    return MappedWord(fact, GeneratorsOfGroup(FamilyObj(fact)!.freeGroup), os.signedgens);
end);

MAJORANA_SomeOrbTest :=
function()
    local ex, rep, gens, orbs, t, ra, fact, maresult;

    t := NanosecondsSinceEpoch();
    ex := A7();;
    t := NanosecondsSinceEpoch() - t;
    Print("ex  setup: ", t/1000000., "\n");

    t := NanosecondsSinceEpoch();
    rep := MAJORANA_SetUp(ex, 1, "AllAxioms");
    t := NanosecondsSinceEpoch() - t;
    Print("rep setup: ", t/1000000., "\n");

    t := NanosecondsSinceEpoch();
    gens := GeneratorsOfGroup(rep.group);
    gens := List(gens, x -> MAJORANA_FindPerm(x, rep, rep));
    gens := List(gens, SignedPermList);

    orbs := MAJORANA_OrbitalStructureSigned( gens
                                           , [1..Length(rep.setup.coords)]
                                           , OnPoints );
    t := NanosecondsSinceEpoch() - t;
    Print("orb setup: ", t/1000000., "\n");

    t := NanosecondsSinceEpoch();
    ra := MAJORANA_OrbitalRepActSigned(orbs, [216, 106]);
    t := NanosecondsSinceEpoch() - t;
    Print("repr calc: ", t/1000000., "\n");

    return [ra, rep, orbs];
end;

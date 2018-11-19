gap> os := MAJORANA_OrbitalStructure([(1,2), (1,2,3,4,5)], [1..5], OnPoints);;
gap> MAJORANA_OrbitalRep(os, [5,2]);
[ 1, 2 ]
gap> ForAll(Tuples([1..5], 2), t -> OnTuples(t, MAJORANA_OrbitalCanonizingElement(os, t)) = MAJORANA_OrbitalRep(os, t));
true
gap> ForAll(Tuples([1..5], 2), t -> OnTuples(t, MAJORANA_OrbitalCanonizingElementInverse(os, t)^-1) = MAJORANA_OrbitalRep(os, t));
true
gap> os := MAJORANA_OrbitalStructure([(1,2,3)(4,5,6), (1,10)(2,20)], [1..20], OnPoints);;
gap> ForAll(Tuples([1..20], 2), t -> OnTuples(t, MAJORANA_OrbitalCanonizingElement(os, t)) = MAJORANA_OrbitalRep(os, t));
true
gap> ForAll(Arrangements([1..20], 2), t -> OnSets(Set(t), MAJORANA_UnorderedOrbitalCanonizingElement(os, t)) = MAJORANA_UnorderedOrbitalRep(os, t));
true
gap> OrbitalTest(os, [1..20]);
true
gap> OrbitalCanonizingTest(os, [1..20]);
true
gap> OrbitalTransversalTest(os, [1..20]);
true
gap> UnorderedOrbitalTest(os, [1..20]);
true
gap> UnorderedOrbitalCanonizingTest(os, [1..20]);
true
gap> UnorderedOrbitalTransversalTest(os, [1..20]);
true

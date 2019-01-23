gap> os := OrbitalStructure([(1,2), (1,2,3,4,5)], [1..5], OnPoints);;
gap> OrbitalRepresentative(os, [5,2]);
[ 1, 2 ]
gap> ForAll(Tuples([1..5], 2), t -> OnTuples(t, OrbitalCanonizingElement(os, t)) = OrbitalRepresentative(os, t));
true
gap> ForAll(Tuples([1..5], 2), t -> OnTuples(t, OrbitalCanonizingElementInverse(os, t)^-1) = OrbitalRepresentative(os, t));
true
gap> os := OrbitalStructure([(1,2,3)(4,5,6), (1,10)(2,20)], [1..20], OnPoints);;
gap> ForAll(Tuples([1..20], 2), t -> OnTuples(t, OrbitalCanonizingElement(os, t)) = OrbitalRepresentative(os, t));
true
gap> ForAll(Arrangements([1..20], 2), t -> OnSets(Set(t), UnorderedOrbitalCanonizingElement(os, t)) = UnorderedOrbitalRepresentative(os, t));
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
gap> os := OrbitalStructure([ (1,13,4,14,5)(2,10,12,9,8)(3,7,15,6,11)(16,17,18,20,19), (1,2,3)(4,6,5)(7,10,13)(8,12,14)(9,11,15)(16,18,21)(17,19,20) ], [1..21], OnPoints);;
gap> OrbitalTest(os, [1..21]);
true
gap> OrbitalCanonizingTest(os, [1..21]);
true
gap> OrbitalTransversalTest(os, [1..21]);
true
gap> UnorderedOrbitalTest(os, [1..21]);
true
gap> UnorderedOrbitalCanonizingTest(os, [1..21]);
true
gap> UnorderedOrbitalTransversalTest(os, [1..21]);
true


#! @Chapter Orbital Structures
#!
#! The functions for orbital structures are based on recent work in permutation
#! group algorithms. An orbital structure contains information about orbits and
#! stabilisers of a group acting on a set for the purposes of quickly
#! determining representatives, canonising elements, and transversal elements
#! (directed) orbitals (orbits of ordered pairs of elements of the domain), and
#! undirected orbitals, i.e. orbits of sets of size two.
#!
#! @Section Examples
#! To create an orbital structure we need a group, a set, and an action
#! @BeginExample
#! gap> os := OrbitalStructure([
#! > (1,13,4,14,5)(2,10,12,9,8)(3,7,15,6,11)(16,17,18,20,19),
#! > (1,2,3)(4,6,5)(7,10,13)(8,12,14)(9,11,15)(16,18,21)(17,19,20) ],
#! > [1..21],
#! > OnPoints);;
#! gap> OrbitalRepresentative(os, [16,2]);
#! [ 16, 2 ]
#! gap> UnorderedOrbitalRepresentative(os, [16,2]);
#! [ 1, 20 ]
#! gap> AllOrbitalRepresentatives(os)
#! [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 1, 4 ], [ 1, 5 ], [ 1, 6 ], [ 1, 16 ],
#!   [ 1, 18 ], [ 1, 20 ], [ 16, 1 ], [ 16, 2 ], [ 16, 3 ], [ 16, 16 ], [ 16, 17 ] ]
#! gap> AllUnorderedOrbitalRepresentatives(os)
#! [ [ 1, 1 ], [ 1, 2 ], [ 1, 4 ], [ 1, 5 ], [ 1, 6 ], [ 1, 16 ], [ 1, 18 ],
#!   [ 1, 20 ], [ 16, 16 ], [ 16, 17 ] ]
#! @EndExample
#!
DeclareCategory( "IsOrbitalStructure", IsObject );
BindGlobal("OrbitalStructureFamily", NewFamily("OrbitalStructureFamily", IsObject ) );
DeclareRepresentation("IsOrbitalStructureRep", IsOrbitalStructure and IsComponentObjectRep, [] );
BindGlobal("OrbitalStructureType", NewType(OrbitalStructureFamily, IsOrbitalStructureRep ) );


#! @Description
#! Given generators, a set, and an action function create an orbital structure.
DeclareGlobalFunction( "OrbitalStructure" );
DeclareGlobalFunction( "OS_OrbitRepresentative" );
DeclareGlobalFunction( "OS_CanonisingElement" );
DeclareGlobalFunction( "OS_CanonisingElementAndRepresentative" );
DeclareGlobalFunction( "OS_StabilizerOf" );


#! @Description
#! Given an orbital structure and a pair of elements
#! of the domain of the orbital structure, return a
#! canonical representative of this pair in its orbit
#! of ordered pairs.
DeclareGlobalFunction( "OrbitalRepresentative" );

#! @Description
#! Return the set of all canonical representatives of orbits
#! of pairs under the action of the orbital structure.
DeclareGlobalFunction( "AllOrbitalRepresentatives" );

#! @Description
#! Given an orbital structure and a pair, return the element
#! <M>g</M> of the group that maps the given pair to the
#! canonical representative.
DeclareGlobalFunction( "OrbitalCanonizingElement" );
DeclareGlobalFunction( "OrbitalCanonizingElementInverse" );

#! @Description
#! Given an orbital structure and a pair, return an iterator
#! that produces a group element for every element in the orbit
#! of the pair that takes the canonical representative of this orbit
#! to that element.
DeclareGlobalFunction( "OrbitalTransversalIterator" );

#! @Description
#! Given an orbital structure and a pair of elements
#! of the domain of the orbital structure, return a
#! canonical representative of this pair in its orbit
#! of sets.
DeclareGlobalFunction( "UnorderedOrbitalRepresentative" );

#! @Description
#! Return the set of all canonical representatives of orbits
#! of sets under the action of the orbital structure.
DeclareGlobalFunction( "AllUnorderedOrbitalRepresentatives" );

#! @Description
#! Given an orbital structure and a set of size 2, return an iterator
#! that produces a group element for every element in the orbit
#! of the set that takes the canonical representative of this orbit
#! to that element.
DeclareGlobalFunction( "UnorderedOrbitalTransversalIterator" );

#! @Description
#! Given an orbital structure and a set of size 2, return the element
#! <M>g</M> of the group that maps the given set to the
#! canonical representative.
DeclareGlobalFunction( "UnorderedOrbitalCanonizingElement" );
DeclareGlobalFunction( "UnorderedOrbitalCanonizingElementInverse" );


DeclareCategory( "IsOrbitalStructure", IsObject );
BindGlobal("OrbitalStructureFamily", NewFamily("OrbitalStructureFamily", IsObject ) );
DeclareRepresentation("IsOrbitalStructureRep", IsOrbitalStructure and IsComponentObjectRep, [] );
BindGlobal("OrbitalStructureType", NewType(OrbitalStructureFamily, IsOrbitalStructureRep ) );

DeclareGlobalFunction( "OrbitalStructure" );
DeclareGlobalFunction( "OS_OrbitRepresentative" );
DeclareGlobalFunction( "OS_CanonisingElement" );
DeclareGlobalFunction( "OS_CanonisingElementAndRepresentative" );
DeclareGlobalFunction( "OS_StabilizerOf" );


DeclareGlobalFunction( "OrbitalRepresentative" );
DeclareGlobalFunction( "AllOrbitalRepresentatives" );
DeclareGlobalFunction( "OrbitalCanonizingElement" );
DeclareGlobalFunction( "OrbitalCanonizingElementInverse" );
DeclareGlobalFunction( "OrbitalTransversalIterator" );

DeclareGlobalFunction( "UnorderedOrbitalRepresentative" );
DeclareGlobalFunction( "AllUnorderedOrbitalRepresentatives" );
DeclareGlobalFunction( "UnorderedOrbitalTransversalIterator" );
DeclareGlobalFunction( "UnorderedOrbitalCanonizingElement" );
DeclareGlobalFunction( "UnorderedOrbitalCanonizingElementInverse" );

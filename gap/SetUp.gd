DeclareGlobalFunction( "MAJORANA_SetUp" );

#! @Chapter Majorana Representations
#! @Section Majorana Representations

#! @Arguments G,T
#! @Returns a record with a component <A>shapes</A>
#! @Description Takes a group <A>G</A> and a <A>G</A>-invariant set of generating involutions
#! <A>T</A>. Returns a list of possible shapes of a Majorana Representation of the form
#! <A>(G,T,V)</A> that is stored in the <A>shapes</A> component of the output.
DeclareGlobalFunction( "ShapesOfMajoranaRepresentation" );

#! @Arguments G,T
#! @Returns a record with a component <A>shapes</A>
#! @Description Performs exactly the same function as <Ref Func="ShapesOfMajoranaRepresentation"/>
#! but gives only those shapes that obey axiom M8. 
DeclareGlobalFunction( "ShapesOfMajoranaRepresentationAxiomM8" );

DeclareGlobalFunction( "MAJORANA_FindPerm" );

DeclareGlobalFunction( "SP_Inverse" );

DeclareGlobalFunction( "SP_Product" );

DeclareGlobalFunction( "MAJORANA_MappedWord" );

DeclareGlobalFunction( "MAJORANA_RemoveDuplicateShapes" );

DeclareGlobalFunction( "MAJORANA_EmbedDihedralAlgebra" );

DeclareGlobalFunction( "MAJORANA_FindTauMap" );

DeclareGlobalFunction( "MAJORANA_AddNewVectors" );

DeclareGlobalFunction( "MAJORANA_AddConjugateVectors" );

DeclareGlobalFunction( "MAJORANA_FindEmbedding" );


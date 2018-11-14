#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Declarations
#

DeclareInfoClass( "InfoMajorana" );

################################################################################
##
## The main Majorana representation function.
##
################################################################################

#! @Arguments input, index, [options]
#! @Returns a record giving a Majorana representation
#! @Description This takes two or three arguments, the first of which must be
#! the output of the function <Ref Func="ShapesOfMajoranaRepresentation"/> and
#! the second of which is the index of the desired shape in list <A>input.shapes</A>.
#!
#! If the optional argument <A>options</A> is given then it must be a record.
#! The following components of <A>options</A> are recognised:
#! <List>
#!  <Mark><C>axioms</C></Mark>
#!  <Item> This component must be bound to the string <A>"AllAxioms"</A> or
#!  <A>"NoAxioms"</A>. If bound to <A>"AllAxioms"</A> then the algorithm assumes the axioms
#!  2Aa, 2Ab, 3A, 4A and 5A as in Seress (2012). If bound to <A>"NoAxioms"</A> then
#!  the algorithm only assumes the Majorana axioms M1 - M7. The default value is
#!  <A>"AllAxioms"</A>. </Item>
#!  <Mark><C>form</C></Mark>
#!  <Item> If this is bound to <A>true</A> then the algorithm assume the existence
#!  of an inner product (as in the definition of a Majorana algebra). Otherwise, if
#!  bound to <A>false</A> then no inner product is assumed (and we are in fact
#!  constructing an axial algebra that satisfies the Majorana fusion law).
#!  The default value is <A>true</A>.</Item>
#!  <Mark><C>embedding</C></Mark>
#!  <Item> If this is bound to <A>true</A> then the algorithm first attempts to construct
#!  large subalgebras of the final representation before starting the main construction.
#!  The default value is <A>false</A>.</Item>
#! </List>
#! @ChapterInfo Majorana representations, The main function
DeclareGlobalFunction( "MajoranaRepresentation" );

################################################################################
##
## The main loop functions
##
################################################################################

DeclareGlobalFunction( "MAJORANA_MainLoop" );

DeclareGlobalFunction( "MAJORANA_FindInnerProducts" );

DeclareGlobalFunction( "MAJORANA_Fusion" );

DeclareGlobalFunction( "MAJORANA_FindAlgebraProducts" );

################################################################################
##
## Functions used in fusion
##
################################################################################

DeclareGlobalFunction( "MAJORANA_AddEvec" );

DeclareGlobalFunction( "MAJORANA_FuseEigenvectors" );

DeclareGlobalFunction( "MAJORANA_FuseEigenvectorsNoForm" );

DeclareGlobalFunction( "MAJORANA_CheckBasis" );

DeclareGlobalFunction( "MAJORANA_IntersectEigenspaces" );

################################################################################
##
## Functions used in MAJORANA_FindAlgebraProducts
##
################################################################################

DeclareGlobalFunction( "MAJORANA_EigenvectorsAlgebraUnknowns" );

DeclareGlobalFunction( "MAJORANA_NullspaceUnknowns" );

DeclareGlobalFunction( "MAJORANA_Resurrection" );

DeclareGlobalFunction( "MAJORANA_AllConjugates" );

################################################################################
##
## The product functions
##
################################################################################

#! @Arguments u, v, algebraproducts, setup
#! @Returns the algebra product of vectors <A>u</A> and <A>v</A>
#! @Description The arguments <A>u</A> and <A>v</A> must be row vectors in sparse
#! matrix format. The arguments <A>algebraproducts</A> and <A>setup</A> must be
#! the components with these names of a representation as outputted by
#! <Ref Func="MajoranaRepresentation"/>. The output is the algebra product of
#! <A>u</A> and <A>v</A>, also in sparse matrix format.
#! @ChapterInfo Functions for calculating with Majorana representations, Calculating products
DeclareGlobalFunction( "MAJORANA_AlgebraProduct" );

#! @Arguments u, v, innerproducts, setup
#! @Returns the inner product of vectors <A>u</A> and <A>v</A>
#! @Description The arguments <A>u</A> and <A>v</A> must be row vectors in sparse
#! matrix format. The arguments <A>innerproducts</A> and <A>setup</A> must be
#! the components with these names of a representation as outputted by
#! <Ref Func="MajoranaRepresentation"/>. The output is the inner product of
#! <A>u</A> and <A>v</A>.
#! @ChapterInfo Functions for calculating with Majorana representations, Calculating products
DeclareGlobalFunction( "MAJORANA_InnerProduct" );

################################################################################
##
## Functions for finding indices that give unknown algebra products
##
################################################################################

DeclareGlobalFunction( "MAJORANA_FindBadIndices" );

DeclareGlobalFunction( "MAJORANA_ListOfBadIndicesForResurrection" );

################################################################################
##
## The conjugating functions
##
################################################################################

DeclareGlobalFunction( "MAJORANA_ConjugateVec" );

DeclareGlobalFunction( "MAJORANA_ConjugateRow" );

################################################################################
##
## Ancilliary functions for finding unknown algebra products
##
################################################################################

DeclareGlobalFunction( "MAJORANA_SeparateAlgebraProduct" );

DeclareGlobalFunction( "MAJORANA_SolutionAlgProducts" );

DeclareGlobalFunction( "MAJORANA_SolveSingleSolution" );

DeclareGlobalFunction( "MAJORANA_RecordSolution" );

DeclareGlobalFunction( "MAJORANA_RemoveKnownAlgProducts" );

################################################################################
##
## Ancilliary functions for finding unknown inner products
##
################################################################################

DeclareGlobalFunction( "MAJORANA_SeparateInnerProduct" );

DeclareGlobalFunction( "MAJORANA_RemoveKnownInnProducts" );

DeclareGlobalFunction( "MAJORANA_SingleInnerSolution" );

DeclareGlobalFunction( "MAJORANA_SolutionInnerProducts" );

DeclareGlobalFunction( "MAJORANA_FillGramMatrix" );

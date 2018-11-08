#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Declarations
#

#! @Arguments input, index, [axioms]
#! @Returns a record giving a Majorana representation
#! @Description This takes two or three arguments, the first of which must be the output of the function
#! <Ref Func="ShapesOfMajoranaRepresentation"/> and the second of which is the index of the desired shape in list <A>input.shapes</A>.
#! The optional third variable is a string that may take one of the following two values
#! * <A>"NoAxioms"</A>: the algorithm assumes no axioms beyond the main axioms of Majorana theory;
#! * <A>"AllAxioms"</A>: the algorithm assumes the axioms 2Aa, 2Ab, 3A, 4A and 5A.
#! If no third argument is given then the default values is set to be <A>"AllAxioms"</A>.
#! @ChapterInfo Majorana representations, The main function
DeclareGlobalFunction( "MajoranaRepresentation" );

DeclareGlobalFunction( "MAJORANA_BasisOfEvecs" );

DeclareGlobalFunction( "MAJORANA_ConjugateVec" );

#! @Arguments u, v, algebraproducts, setup
#! @Returns the algebra product of vectors <A>u</A> and <A>v</A>
#! @Description The arguments <A>u</A> and <A>v</A> must be row vectors in sparse
#! matrix format. The arguments <A>algebraproducts</A> and <A>setup</A> must be
#! the components with these names of a representation as outputted by
#! <Ref Func="MajoranaRepresentation"/>. The output is the algebra product of
#! <A>u</A> and <A>v</A>, also in sparse matrix representation.
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

DeclareGlobalFunction( "MAJORANA_FindBadIndices" );

DeclareGlobalFunction( "MAJORANA_AddEvec" );

DeclareGlobalFunction( "MAJORANA_FuseEigenvectors" );

DeclareGlobalFunction( "MAJORANA_Fusion" );

DeclareGlobalFunction( "MAJORANA_CheckBasis" );

DeclareGlobalFunction( "MAJORANA_FillGramMatrix" );

DeclareGlobalFunction( "MAJORANA_SeparateInnerProduct" );

DeclareGlobalFunction( "MAJORANA_EigenvectorsAlgebraUnknowns" );

DeclareGlobalFunction( "MAJORANA_AxiomM1" );

DeclareGlobalFunction( "MAJORANA_Orthogonality" );

DeclareGlobalFunction( "MAJORANA_SeparateAlgebraProduct" );

DeclareGlobalFunction( "MAJORANA_RecordSolution" );

DeclareGlobalFunction( "MAJORANA_ConjugateRow" );

DeclareGlobalFunction( "MAJORANA_UnknownAlgebraProducts" );

DeclareGlobalFunction( "MAJORANA_ListOfBadIndicesForResurrection" );

DeclareGlobalFunction( "MAJORANA_AllConjugates" );

DeclareGlobalFunction( "MAJORANA_NullspaceUnknowns" );

DeclareGlobalFunction( "MAJORANA_Resurrection" );

DeclareGlobalFunction( "MAJORANA_SolutionAlgProducts" );

DeclareGlobalFunction( "MAJORANA_SolveSingleSolution" );

DeclareGlobalFunction( "MAJORANA_SingleInnerSolution" );

DeclareGlobalFunction( "MAJORANA_RemoveKnownAlgProducts" );

DeclareGlobalFunction( "MAJORANA_RemoveKnownInnProducts" );

DeclareGlobalFunction( "MAJORANA_SolutionInnerProducts" );

DeclareGlobalFunction( "MAJORANA_CheckNullSpace" );

DeclareGlobalFunction( "MAJORANA_MainLoop" );

DeclareInfoClass( "InfoMajorana" );

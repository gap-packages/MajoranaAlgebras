#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Declarations
#

#! @Description
#!   Constructs a Majorana representation given a finite group <A>G</A> 
#!   and a <A>G</A> invariant set of generating involutions <A>T</A>.
#! @Arguments G, T
DeclareGlobalFunction( "MajoranaRepresentation" );

DeclareGlobalFunction( "MAJORANA_ChangeFieldOfRep" );

DeclareGlobalFunction( "MAJORANA_SetValue" );

DeclareGlobalFunction( "MAJORANA_BasisOfEvecs" );

DeclareGlobalFunction( "MAJORANA_ConjugateVec" );

DeclareGlobalFunction( "MAJORANA_AlgebraProduct" );

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

DeclareGlobalFunction( "MAJORANA_SeparateAlgebraProduct" );

DeclareGlobalFunction( "MAJORANA_RecordSolution" );

DeclareGlobalFunction( "MAJORANA_ConjugateRow" );

DeclareGlobalFunction( "MAJORANA_UnknownAlgebraProducts" );

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

DeclareGlobalFunction( "MAJORANA_AdjointAction" );

DeclareGlobalFunction( "MAJORANA_NaiveProduct" );

DeclareGlobalFunction( "MAJORANA_FindFusionTable" );

DeclareGlobalFunction( "MAJORANA_MainLoop" );

DeclareInfoClass( "InfoMajorana" );

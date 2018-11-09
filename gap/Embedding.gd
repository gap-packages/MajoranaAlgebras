#! @Chapter Embedding representations
#! @Section Embedding a single representation

#! @Arguments rep, subrep
#! @Description Searches for possible embeddings of subrep into rep. If such an
#! embedding exists, then use the structure of <A>subrep</A> to record
#! new products and eigenvectors of <A>rep</A>.
DeclareGlobalFunction( "MAJORANA_EmbedKnownRep" );

DeclareGlobalFunction( "MAJORANA_Embed" );

DeclareGlobalFunction( "MAJORANA_ImageVector" );

DeclareGlobalFunction( "MAJORANA_CheckEmbedding" );

#! @Section Embedding multiple representations

DeclareGlobalFunction( "MAJORANA_MaximalSubgps" );

DeclareGlobalFunction( "MAJORANA_AllEmbeddings" );

DeclareGlobalFunction( "MAJORANA_Image" );

DeclareGlobalFunction( "MAJORANA_EmbedDihedral" );

DeclareGlobalFunction( "MAJORANA_FindPerm" );

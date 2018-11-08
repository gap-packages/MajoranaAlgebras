#! @Chapter Shapes of a Majorana representation
#! @Section The shape functions

#! @Arguments G,T
#! @Returns a record with a component <A>shapes</A>
#! @Description Takes a group <A>G</A> and a <A>G</A>-invariant set of generating involutions
#! <A>T</A>. Returns a list of possible shapes of a Majorana Representation of the form
#! <A>(G,T,V)</A> that is stored in the <A>shapes</A> component of the output.
DeclareGlobalFunction( "ShapesOfMajoranaRepresentation" );

#! @Arguments G,T
#! @Returns a record with a component <A>shapes</A>
#! @Description Performs exactly the same function as <Ref Func="ShapesOfMajoranaRepresentation"/>
#! but gives only those shapes at obey axiom M8. That is to say, we additionally assume
#! that if $t,s \in T$ such that $|ts| = 2$ then the dihedral subalgebra $\langle \langle a_t, a_s \rangle \rangle$
#! is of type $2A$ if and only if $ts \in T$ (and otherwise is of type $2B$).
DeclareGlobalFunction( "ShapesOfMajoranaRepresentationAxiomM8" );

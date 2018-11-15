#! @Chapter Shapes of a Majorana representation
#! @Section The shapes functions

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

#! @Arguments G,T
#! @Returns true if <A>(G,T)</A> is a 6-transposition group, otherwise returns false
#! @Description For a group <A>G</A> and a subset <A>T</A> of <A>G</A>, returns true
#! if all of the following conditions are satisfied:
#! *<A>T</A> is a set of involutions that generate <A>G</A>;
#! *<A>T</A> is closed under conjugation by <A>G</A>;
#! *the order of the product of two elements of <A>T</A> is at most 6.
DeclareGlobalFunction( "MAJORANA_IsSixTranspositionGroup" );

#! @Arguments input
#! @Description If an automorphism of the group <A>G</A> stabilises the set
#! <A>T</A> then it induces an action on the pairs of elements of <A>T</A> and
#! therefore on the shapes of a possible Majorana representation of the form
#! <A>(G,T,V)</A>. If one shape is mapped to another in this way then their
#! corresponding algebras must be isomorphic.
#!
#! This function takes the record <A>input</A> as produced by the function
#! <Ref Func="ShapesOfMajoranaRepresentation"/> or <Ref Func="ShapesOfMajoranaRepresentationAxiomM8"/>
#! and replaces <A>input.shapes</A> with a list of shapes such that no two
#! can be mapped to each other by an automorphism of <A>G</A>.
DeclareGlobalFunction( "MAJORANA_RemoveDuplicateShapes" );

DeclareGlobalFunction( "MAJORANA_RecordSubalgebras" );

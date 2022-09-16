# Signed permutations

#! @Chapter Signed Permutations
#! We provide <E>signed permutations</E>, that is permutations
#! that can additionally change the sign of their result.
#!
#! Assume <M>n \in \mathbb{N}</M>, then a signed permutation on <M>n</M> points
#! is a permutation <M>\pi</M> on <M>\{ 1 \ldots n \}</M> together with signs
#! <M>sgn : \{ 1 .. n \} \rightarrow \{-1,1\}</M>.
#
#! A signed permutation on <M>n</M> points acts on the set
#! <M>\{ -n \ldots 1, 1 \ldots n \}</M> by
#! <M> \omega ^ { (\pi, sgn) } = sgn(\omega)\cdot sgn(|\omega|^\pi) \cdot (|\omega|^\pi) </M>.
#!
#! We provide two representations of signed permutations, one as a list of images
#! <Ref Filt="IsSignedPermListRep" Label="for IsSignedPerm and IsPositionalObjectRep"/> and one formed as pair of a permutation and
#! a sign map <Ref Filt="IsSignedPermRep" Label="for IsSignedPerm and IsPositionalObjectRep"/>. Our benchmarks indicate that a list of
#! images is the better representation, and hence this is the default.
#!
#! To get started with signed permutations consider the following example
#! @BeginExample
#! gap> s := SignedPerm([2,-1]);
#! <signed permutation in list rep>
#! gap> 1 ^ s;
#! 2
#! gap> 2 ^ s;
#! -1
#! gap> OnPoints(2, s);
#! -1
#! @EndExample
#!
#! One can form groups out of signed permutations
#!
#! @BeginExample
#! gap> r := SignedPerm([-1,3,-2,4]);; t := SignedPerm([3,1,4,2]);;
#! gap> G := Group(r,t);
#! <group with 2 generators>
#! gap> Size(G);
#! 32
#! gap> Orbit(G, 1, OnPoints);
#! [ 1, -1, 3, -3, -2, 4, 2, -4 ]
#! gap> Stabilizer(G, 1, OnPoints);
#! <group of size 4 with 9 generators>
#! @EndExample
#!
#! Note that currently the package does not make an effort to exploit the special
#! structure of signed permutation groups as permutation groups.
#!
#! @Section Different Representations
#!
#! To create signed permutations in the different representations, we provide a constructor.
#! @BeginExample
#! gap> r := NewSignedPerm(IsSignedPermRep, [-1,3,-2,4]);;
#! gap> t := SignedPerm(IsSignedPermRep, [3,1,4,2]);;
#! gap> G := Group(r,t);
#! <group with 2 generators>
#! gap> Size(G);
#! 32
#! gap> r := NewSignedPerm(IsSignedPermListRep, [-1,3,-2,4]);;
#! gap> t := SignedPerm(IsSignedPermListRep, [3,1,4,2]);;
#! gap> G := Group(r,t);
#! <group with 2 generators>
#! gap> Size(G);
#! 32
#! @EndExample
#!
#! @Section Low-Level Descriptions
#! @Description
#! Category of signed permutations
DeclareCategory("IsSignedPerm",
                IsAssociativeElement and
                IsExtLElement and
                IsExtRElement and
                IsMultiplicativeElement and
                IsMultiplicativeElementWithOne and
                IsMultiplicativeElementWithInverse and
                IsFiniteOrderElement );

BindGlobal("SignedPermFamily", NewFamily("SignedPermFamily", IsSignedPerm));

DeclareCategoryCollections( "IsSignedPerm" );
DeclareCategoryCollections( "IsSignedPermCollection" );
InstallTrueMethod( IsGeneratorsOfMagmaWithInverses
                 , IsSignedPermCollection );


#! @Description
#! Convert a signed permutation into a list of images, equivalent
#! to List([1..LargestMovedPoint(s)], x -> x^s);
#! @Arguments perm
DeclareOperation( "ListSignedPerm", [ IsSignedPerm ] );
#! @Description
#! Convert a signed permutation to a list of images of length <A>len</A>.
#! Arguments perm, len
DeclareOperation( "ListSignedPerm", [ IsSignedPerm, IsPosInt] );

#! @Description
#! Given a list of signed images create a signed permutation object
#! in <Ref Filt="IsSignedPermListRep" Label="for IsSignedPerm and IsPositionalObjectRep"/>.
DeclareGlobalFunction("SignedPerm");
DeclareConstructor( "NewSignedPerm", [ IsSignedPerm, IsList ] );
DeclareConstructor( "NewSignedPerm", [ IsSignedPerm, IsPerm, IsList ] );

#! @Description
#! Representation of signed permutations as a permutation and a vector of signs.
DeclareRepresentation("IsSignedPermRep", IsSignedPerm and IsPositionalObjectRep, []);
BindGlobal("SignedPermType", NewType(SignedPermFamily, IsSignedPermRep));

#! @Description
#! Representation of signed permutations as a list of signed images
DeclareRepresentation("IsSignedPermListRep", IsSignedPerm and IsPositionalObjectRep, []);
BindGlobal("SignedPermListType", NewType(SignedPermFamily, IsSignedPermListRep));

#! @Description
#! Only act as a permutation on <M>\{ 1\ldots n\}</M>
DeclareGlobalFunction("OnPosPoints");

#! @Description
#! The largest point that is moved by the signed permutation, where moving includes
#! changing the sign.
DeclareAttribute("LargestMovedPoint", IsSignedPerm );

#! @Description
#! Create a random list of images that can be
#! used to create a signed permutation.
DeclareGlobalFunction("RandomSignedPermList");

#! @Description
#! Create a random signed permutation
DeclareGlobalFunction("RandomSignedPerm");


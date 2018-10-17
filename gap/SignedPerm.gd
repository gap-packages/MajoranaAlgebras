# Signed permutations

DeclareCategory("IsSignedPerm",
                IsMultiplicativeElement and IsMultiplicativeElementWithOne and IsMultiplicativeElementWithInverse );
BindGlobal("SignedPermFamily", NewFamily("SignedPermFamily"));
DeclareRepresentation("IsSignedPermRep", IsSignedPerm and IsPositionalObjectRep, []);
BindGlobal("SignedPermType", NewType(SignedPermFamily, IsSignedPermRep));

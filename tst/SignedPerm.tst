# Signed permutations
gap> sp := SignedPerm([1,2,3,4,5,6,7,8]);
<signed permutation in list rep>
gap> sp := SignedPerm([1,1,1,1,1]);
Error, list does not define a signed permutation
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := SignedPerm(ls);;
gap> t := SignedPerm(lt);;
gap> ListSignedPerm(SignedPerm(ls)) = ls;
true
gap> ListSignedPerm(SignedPerm(lt)) = lt;
true
gap> ListSignedPerm(s, 10);
[ 3, -2, 1, 4, 6, -5, 7, 8, 9, 10 ]
gap> s * t;;
gap> t * s;;
gap> One(s);; One(t);;
gap> Inverse(s);;
gap> Inverse(t);;
gap> Inverse(s) * s = One(s);
true
gap> IsOne(s * Inverse(s));
true
gap> IsOne(Inverse(s) * s);
true
gap> IsOne(t * Inverse(t));
true
gap> IsOne(Inverse(t) * t);
true
gap> s := SignedPerm([-3,2,1]);;
gap> 1 ^ s;
-3
gap> (-1)^s;
3
gap> s = t;
false
gap> s < t;
true
gap> t < s;
false
gap> t > s;
true
gap> G := Group(s,t);; Size(G);
8

# Signed permutations in (perm, sign) rep
gap> sp := NewSignedPerm(IsSignedPermRep, [1,2,3,4,5,6,7,8]);
<signed permutation (), [  ]>
gap> sp := NewSignedPerm(IsSignedPermRep, [1,1,1,1,1]);
Error, l does not define a permutation
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := NewSignedPerm(IsSignedPermRep, ls);;
gap> t := NewSignedPerm(IsSignedPermRep, lt);;
gap> ListSignedPerm(SignedPerm(ls)) = ls;
true
gap> ListSignedPerm(SignedPerm(lt)) = lt;
true
gap> s * t;;
gap> t * s;;
gap> One(s);; One(t);;
gap> Inverse(s);;
gap> Inverse(t);;
gap> Inverse(s) * s = One(s);
true
gap> IsOne(s * Inverse(s));
true
gap> IsOne(Inverse(s) * s);
true
gap> IsOne(t * Inverse(t));
true
gap> IsOne(Inverse(t) * t);
true
gap> s := SignedPermList([-3,2,1]);;
gap> 1 ^ s;
-3
gap> (-1)^s;
3
gap> s = t;
false
gap> s < t;
false
gap> t < s;
true
gap> t > s;
false

# Constructors
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := SignedPerm(IsSignedPermListRep, ls);;
gap> t := SignedPerm(IsSignedPermListRep, lt);;
gap> ListSignedPerm(s) = ls;
true
gap> ListSignedPerm(t) = lt;
true
gap> IsOne(Inverse(s*t)*s*t);
true
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := SignedPerm(IsSignedPermRep, ls);;
gap> t := SignedPerm(IsSignedPermRep, lt);;
gap> ListSignedPerm(s) = ls;
true
gap> ListSignedPerm(t) = lt;
true
gap> IsOne(Inverse(s*t)*s*t);
true
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := SignedPerm(IsSignedPermRep, ls);;
gap> t := SignedPerm(IsSignedPermListRep, lt);;
gap> ListSignedPerm(s) = ls;
true
gap> ListSignedPerm(t) = lt;
true
gap> IsOne(Inverse(s*t)*s*t);
Error, no method found! For debugging hints type ?Recovery from NoMethodFound
Error, no 1st choice method found for `*' on 2 arguments

#
gap> TestSomeRandomPerms();;

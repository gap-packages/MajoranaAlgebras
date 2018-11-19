gap> sp := SignedPermList([1,2,3,4,5,6,7,8]);
<signed permutation (), [  ]>
gap> sp := SignedPermList([1,1,1,1,1]);
Error, l does not define a permutation
gap> ls := [3,-2,1,4,6,-5];; lt := [-1,-2,-3,5,4,10,7,8,-9,6];;
gap> s := SignedPermList(ls);;
gap> t := SignedPermList(lt);;
gap> ListSignedPerm(SignedPermList(ls)) = ls;
true
gap> ListSignedPerm(SignedPermList(lt)) = lt;
true
gap> s * t;
<signed permutation ( 1, 3)( 4, 5,10, 6), [ Z(2)^0, 0*Z(2), Z(2)^0, Z(2)^0,
  0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
gap> t * s;
<signed permutation ( 1, 3)( 4, 6,10, 5), [ Z(2)^0, 0*Z(2), Z(2)^0, 0*Z(2),
  Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
gap> One(s); One(t);
<signed permutation (), [  ]>
<signed permutation (), [  ]>
gap> Inverse(s);
<signed permutation (1,3)(5,6), [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2),
  Z(2)^0 ]>
gap> Inverse(t);
<signed permutation ( 4, 5)( 6,10), [ Z(2)^0, Z(2)^0, Z(2)^0, 0*Z(2), 0*Z(2),
  0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
gap> Inverse(s) * s = One(s);
true
gap> Inverse(t) * s;
<signed permutation ( 1, 3)( 4, 6,10, 5), [ Z(2)^0, 0*Z(2), Z(2)^0, 0*Z(2),
  Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
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
gap> G := Group(s,t);; Size(G);
8
gap> TestSomeRandomPerms();;

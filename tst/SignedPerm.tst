gap> sp := SignedPermList([1,2,3,4,5,6,7,8]);
<signed permutation (), [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 
0*Z(2), 0*Z(2) ]>
gap> sp := SignedPermList([1,1,1,1,1]);
Error, l does not define a permutation
gap> s := SignedPermList([3,-2,1,4,6,-5]);;
gap> t := SignedPermList([-1,-2,-3,5,4,10,7,8,-9,6]);;
gap> s * t;
<signed permutation ( 1, 3)( 4, 5,10, 6), [ Z(2)^0, 0*Z(2), Z(2)^0, 0*Z(2), 
0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
gap> t * s;
<signed permutation ( 1, 3)( 4, 6,10, 5), [ Z(2)^0, 0*Z(2), Z(2)^0, 0*Z(2), 
0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ]>
gap> One(s); One(t);
<signed permutation (), [ 0*Z(2) ]>
<signed permutation (), [ 0*Z(2) ]>
gap> Inverse(s);
<signed permutation (1,3)(5,6), [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), Z(2)^0, 
0*Z(2) ]>
gap> Inverse(t);
<signed permutation ( 4, 5)( 6,10), [ Z(2)^0, Z(2)^0, Z(2)^0, 0*Z(2), 0*Z(2), 
0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]>
gap> Inverse(s) * s = One(s);
true
gap> Inverse(t) * s;
<signed permutation ( 1, 3)( 4, 6,10, 5), [ Z(2)^0, 0*Z(2), Z(2)^0, 0*Z(2), 
0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ]>
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

gap> mat := [ [ 1, 0, -1, -1, -1, 1 ], [ -1, 1, -2, -1, -2, 0 ], [ -1, 2, -2, -1, 2, -3 ], [ -1, 3, 0, -2, 1, -4 ], [ 0, 1, 0, -2, 3, 1 ] ];;
gap> MAJORANA_ReversedEchelonForm(mat);
gap> Display(mat);
[ [  -1,   0,   0,   0,   0,   1 ],
  [   0,   0,   0,   0,   1,   0 ],
  [  -2,   0,   0,   1,   0,   0 ],
  [   0,   0,   1,   0,   0,   0 ],
  [  -3,   1,   0,   0,   0,   0 ] ]
gap> mat := IdentityMat(5,5);;
gap> MAJORANA_ReversedEchelonForm(mat);
gap> Display(mat);
[ [  0,  0,  0,  0,  1 ],
  [  0,  0,  0,  1,  0 ],
  [  0,  0,  1,  0,  0 ],
  [  0,  1,  0,  0,  0 ],
  [  1,  0,  0,  0,  0 ] ]
gap> mat := NullMat(5,6);;
gap> MAJORANA_ReversedEchelonForm(mat);
gap> Display(mat);
[ [  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0 ] ]
gap> mat := [[1,-1,0],[2,-2,0]];;
gap> MAJORANA_NullSpace(mat);
[ [  ], [ [ 0, 0, 1 ], [ 1, 1, 0 ] ] ]

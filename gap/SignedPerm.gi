BindGlobal( "TrimSignedPerm",
function(p)
    local ls;
    ls := Length(p![2]);
    while ls > 0 and IsZero(p![2][ls]) do ls := ls - 1; od;
    p![3] := Maximum(LargestMovedPoint(p![1]), ls);
    p![2] := (p![2]){[1..p![3]]} + ListWithIdenticalEntries(p![3], Z(2) * 0);
    return p;
end);


BindGlobal( "SignedPermList",
function(l)
    local sp, s, ts;

    sp := [];

    # Signs
    s := PositionsProperty(l, x -> x < 0);
    l{s} := -l{s};

    sp![1] := PermList(l);
    if sp![1] = fail then
        ErrorNoReturn("l does not define a permutation");
    fi;

    sp![2] := ListWithIdenticalEntries(Length(l), Z(2) * 0);
    sp![2]{s} := List(s, x -> Z(2)^0);
    sp![2] := Permuted(sp![2], sp![1]);

    TrimSignedPerm(sp);
    ConvertToVectorRep(sp![2]);

    return Objectify( SignedPermType, sp );
end);

BindGlobal( "ListSignedPerm",
function( arg )
    local sp, len, p, s, l, i;

    sp := arg[1]; p := sp![1]; s := sp![2];

    if Length(arg) = 2 then
        len := arg[2];
    else
        len := Length(s);
    fi;

    l := ListPerm(p, len);

    for i in [ 1 .. len ] do
        if IsBound(s[l[i]]) and IsOne(s[l[i]]) then
            l[i] := -l[i];
        fi;
    od;

    return l;

end );

BindGlobal( "SignedPerm",
function(p, s)
    return Objectify( SignedPermType, TrimSignedPerm([ p, s ]));
end);

InstallMethod(ViewObj, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    Print("<signed permutation ", sp![1], ", ", sp![2], ">");
end);

InstallMethod(PrintObj, "for signed permutations",
[ IsSignedPermRep ],
function(sp)
    Print("<signed permutation ", sp![1], ", ", sp![2], ">");
end);

InstallMethod(\*, "for signed permutations",
  [ IsSignedPermRep, IsSignedPermRep ],
function(l, r)
    local sp;
    return Objectify( SignedPermType
                    , TrimSignedPerm( [ l![1] * r![1]
                                      , Permuted(Zero(r![2]) + l![2], r![1]) + r![2] ] ) );
end);

InstallMethod(InverseOp, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, TrimSignedPerm([ sp![1]^-1, Permuted(sp![2], sp![1]^-1) ]) );
end);

InstallMethod(OneImmutable, "for signed permutations",
              [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, [ (), [] ]);
end);

InstallMethod(OneMutable, "for signed permutations",
              [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, [ (), [] ]);
end);

InstallMethod(IsOne, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    return IsOne(sp![1]) and IsZero(sp![2]);
end);

InstallMethod(\^, "for an integer and a signed permutation",
  [ IsInt, IsSignedPermRep ],
function(pt, sp)
    local sign;

    sign := SignInt(pt);
    pt := AbsInt(pt);

    pt := pt ^ sp![1];
    if IsOne(sp![2][pt]) then
        sign := -sign;
    fi;
    return sign * pt;
end);

InstallMethod( \=, "for a signed permutation and a signed permutation",
               [ IsSignedPermRep, IsSignedPermRep ],
function(l,r)
    # Trim?
    if l![1] = r![1] then
        return l![2] = r![2];
    fi;
    return false;
end);

InstallMethod( \<, "for a signed permutation and a signed permutation",
               [ IsSignedPermRep, IsSignedPermRep ],
function(l,r)
    # Trim?
    if l![1] < r![1] then
        return true;
    elif l![1] = r![1] then
        return l![2] < r![2];
    fi;
    return false;
end);

# for an int and a signed perm
InstallGlobalFunction(OnPosPoints,
    { pnt, elm } -> pnt ^ elm![1]);

InstallTrueMethod( IsGeneratorsOfMagmaWithInverses
                 , IsSignedPermCollection );



#### bla
BindGlobal( "TrimSignedPermList",
function(p)
    local ls;
    ls := Length(p);
    while ls > 0 and p[ls] = ls do
        Unbind(p[ls]);
        ls := ls - 1;
    od;
    return p;
end);

BindGlobal( "SignedPermL",
function(list)
    return Objectify( SignedPermListType, [ TrimSignedPermList(list) ]);
end);

InstallMethod(ViewObj, "for signed permutations (list rep)",
  [ IsSignedPermListRep ],
function(sp)
    Print("<signed permutation in list rep>");
end);

InstallMethod(PrintObj, "for signed permutations",
[ IsSignedPermListRep ],
function(sp)
    Print("<signed permutation in list rep>");
end);

InstallMethod(\*, "for signed permutations",
  [ IsSignedPermListRep, IsSignedPermListRep ],
function(l, r)
    local degree, res, i, ll, rr;
    degree := Maximum(Length(l![1]), Length(r![1]));
    res := [1..degree];
    ll := [1..degree];
    ll{[1..Length(l![1])]} := l![1];
    rr := [1..degree];
    rr{[1..Length(r![1])]} := r![1];
    for i in [1..degree] do
        res[i] := SignInt(ll[i]) * rr[AbsInt(ll[i])];
    od;
    return Objectify(SignedPermListType, [ TrimSignedPermList( res ) ] );
end);

InstallMethod(InverseOp, "for signed permutations",
  [ IsSignedPermListRep ],
function(sp)
    local id, rhs, i;

    id := [1..Length(sp![1])];
    rhs := ShallowCopy(sp![1]);

    for i in [1..Length(id)] do
        id[i] := SignInt(rhs[i]) * i;
        rhs[i] := AbsInt(rhs[i]);
    od;

    SortParallel(rhs, id);

    return Objectify(SignedPermListType, [ TrimSignedPermList( id ) ] );
end);

InstallMethod(OneImmutable, "for signed permutations",
              [ IsSignedPermListRep ],
function(sp)
    return Objectify(SignedPermListType, [ [] ]);
end);

InstallMethod(OneMutable, "for signed permutations",
              [ IsSignedPermListRep ],
function(sp)
    return Objectify(SignedPermListType, [ [] ]);
end);

InstallMethod(IsOne, "for signed permutations",
  [ IsSignedPermListRep ],
function(sp)
    return IsEmpty(sp![1]);
end);

InstallMethod(\^, "for an integer and a signed permutation",
  [ IsInt, IsSignedPermListRep ],
function(pt, sp)
    local apt;

    apt := AbsInt(pt);
    if IsBound(sp![1][apt]) then
        return SignInt(pt) * sp![1][apt];
    else
        return pt;
    fi;
end);

InstallMethod( \=, "for a signed permutation and a signed permutation",
               [ IsSignedPermListRep, IsSignedPermListRep ],
function(l,r)
    # Trim?
    return l![1] = r![1];
end);

InstallMethod( \<, "for a signed permutation and a signed permutation",
               [ IsSignedPermListRep, IsSignedPermListRep ],
function(l,r)
    return l![1] < r![1];
end);


InstallMethod( NewSignedPerm, "for perm, vec rep",
               [ IsSignedPermRep, IsList ],
function(filt, list)
    return SignedPermList(ShallowCopy(list));
end);

InstallMethod( NewSignedPerm, "for list rep",
               [ IsSignedPermListRep, IsList ],
function(filt, list)
    return SignedPermL(ShallowCopy(list));
end);

BindGlobal("RandomSignedPermList",
function(degree)
    local res, i;

    res := Permuted([1..degree], Random(SymmetricGroup(degree)));
    res := List(res, x -> Random([-1,1]) * x);

    return res;
end);


BindGlobal("RandomSignedPerm",
function(filt, degree)
    return NewSignedPerm(filt, RandomSignedPermList(degree));
end);

BindGlobal("TestSomeRandomPerms",
function()
    local perms, i, p, q, r, t1, t2, t3, t4, t;

    perms := List([1..100], x->RandomSignedPermList(40));

    t1 := 0;
    t2 := 0;
    t3 := 0;
    t4 := 0;
 
    for i in [1..500] do
        perms := List([1..Random([1..100])], x->RandomSignedPermList(1000));
        t := NanosecondsSinceEpoch();
        p := ListSignedPerm(Product(perms, x -> NewSignedPerm(IsSignedPermRep, x)));
        t1 := t1 + (NanosecondsSinceEpoch() - t);
        t := NanosecondsSinceEpoch();
        q := Product(perms, x->NewSignedPerm(IsSignedPermListRep, x))![1];
        t2 := t2 + (NanosecondsSinceEpoch() - t);
        if p <> q then
            Error("Products do not match\n");
        fi;
        for p in perms do
            q := NewSignedPerm(IsSignedPermRep, p);
            t := NanosecondsSinceEpoch();
            r := Inverse(q);
            t3 := t3 + (NanosecondsSinceEpoch() - t);
            if not IsOne(r * q) then
                Error("inverse is not inverse");
            fi;
            q := NewSignedPerm(IsSignedPermListRep, p);
            t := NanosecondsSinceEpoch();
            r := Inverse(q);
            t4 := t4 + (NanosecondsSinceEpoch() - t);
            if not IsOne(r * q) then
                Error("inverse is not inverse");
            fi;
        od;
    od;
    return [t1,t2,t3,t4];
end);


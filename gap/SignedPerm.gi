BindGlobal( "SignedPermList",
function(l)
    local p, s, ts;

    s := PositionsProperty(l, x -> x < 0);
    l{s} := -l{s};
    p := PermList(l);
    if p = fail then
        ErrorNoReturn("l does not define a permutation");
    fi;

    ts := ListWithIdenticalEntries(Length(l), Z(2) * 0);
    ts{s} := List(s, x -> Z(2)^0);
    ConvertToVectorRep(ts);

    return Objectify( SignedPermType, [ p, ts ] );
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
        if IsBound(s[i]) and IsOne(s[i]) then
            l[i] := -l[i];
        fi;
    od;

    return l;

end );

BindGlobal( "SignedPerm",
function(p, s)
    return Objectify( SignedPermType, [ p, s ]);
end);

InstallMethod(ViewObj, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    Print("<signed permutation ", sp![1], ", ", sp![2], ">");
end);

InstallMethod(\*, "for signed permutations",
  [ IsSignedPermRep, IsSignedPermRep ],
function(l, r)
    return Objectify(SignedPermType, [ l![1] * r![1]
                                     # This is really ugly!
                                     , l![2] + Permuted(Zero(l![2]) + r![2], l![1]) ] );
end);

InstallMethod(InverseOp, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, [ sp![1]^-1, Permuted(sp![2], sp![1]^-1) ] );
end);

InstallMethod(OneImmutable, "for signed permutations",
              [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, [ (), [0 * Z(2)] ]);
end);

InstallMethod(OneMutable, "for signed permutations",
              [ IsSignedPermRep ],
function(sp)
    return Objectify(SignedPermType, [ (), [0 * Z(2)] ]);
end);

InstallMethod(IsOne, "for signed permutations",
  [ IsSignedPermRep ],
function(sp)
    return IsOne(sp![1]) and IsZero(sp![2]);
end);

InstallMethod(\^, "for an integer and a signed permutation",
  [ IsInt, IsSignedPermRep ],
function(pt, sp)
    local spt, sign;

    if pt < 0 then
        spt := -pt;
        sign := -1;
    else
        spt := pt;
        sign := 1;
    fi;

    if IsBound(sp![2][spt]) and IsOne(sp![2][spt]) then
        sign := -sign;
    fi;
    return sign * (spt^sp![1]);
end);

InstallMethod( \=, "for a signed permutation and a signed permutation",
               [ IsSignedPermRep, IsSignedPermRep ],
function(l,r)
    local lmp;

    if l![1] = r![1] then
        lmp := LargestMovedPoint(l![1]);
        return l![2]{ [1..lmp ]} = r![2]{ [1..lmp] };
    fi;
    return false;
end);

InstallMethod( \<, "for a signed permutation and a signed permutation",
               [ IsSignedPermRep, IsSignedPermRep ],
function(l,r)
    local lmp;

    if l![1] < r![1] then
        return true;
    elif l![1] = r![1] then
        lmp := LargestMovedPoint(l![1]);
        return l![2]{ [1..lmp ]} < r![2]{ [1..lmp] };
    fi;
    return false;
end);

# for an int and a signed perm
InstallGlobalFunction(OnPosPoints,
    { pnt, elm } -> pnt ^ elm![1]);

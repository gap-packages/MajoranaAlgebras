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
    local spt, sign;

    if pt < 0 then
        spt := -pt;
        sign := -1;
    else
        spt := pt;
        sign := 1;
    fi;

    if IsOne(sp![2][spt]) then
        sign := -sign;
    fi;
    return sign * (spt^sp![1]);
end);



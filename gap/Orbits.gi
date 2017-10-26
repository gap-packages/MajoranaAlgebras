InstallGlobalFunction( MAJORANA_Orbits,

    function(G, D, act)
    
    local       gens,
                orbs,
                elements, 
                orb,
                sort,
                plist,
                pos,
                use,
                o,
                result;
  
    gens := GeneratorsOfGroup(G);
    
    sort := Length(D)>0 and CanEasilySortElements(D[1]);
    plist := IsPlistRep(D);
  
    if not plist then
        use := BlistList([1..Length(D)],[]);
    fi;
    
    orbs := [];
    elements := [];
    pos := 1;
    
    while Length(D) > 0 and pos <= Length(D) do
    
        orb := MAJORANA_OrbitOp( G, D[pos], gens, act );
        
        Add( orbs, orb[1] );
        Add(elements, orb[2]);
            
        if plist then
            if sort then
                D:=Difference(D,orb[1]);
                MakeImmutable(D); # to remember sortedness
            else
                D:=Filtered(D,i-> not i in orb[1]);
            fi;
        else
            for o in orb[1] do
                use[PositionCanonical(D,o)]:=true;
            od;
            
            # not plist -- do not take difference as there may be special
            # `PositionCanonical' method.
            
            while pos <= Length(D) and use[pos] do 
                pos := pos + 1;
            od;
        fi;
    od;
    
    result := rec(  orbits := orbs,
                    elts := elements );
    
    return result;
    
    end );
    
InstallGlobalFunction( MAJORANA_OrbitOp, 

    function( G, pnt, gens, act )
    
    local   D,
            orb,
            orb_elts,
            d,
            gen,
            i,
            g,
            h,
            p,
            count,
            result;
    
    D := DomainForAction(pnt,gens,act);
    
    pnt := Immutable(pnt);
    
    d := NewDictionary(pnt, false, D);
    
    orb := [ pnt ];
    orb_elts := [ Identity(G) ];
    
    AddDictionary(d,pnt);
    
    g := Identity(G);
    
    count := 0;
    
    for p in orb do
        
        count := count + 1;
        h := orb_elts[count];
        
        for gen in gens do
    
            i := act(p,gen);
            g := h*gen;
            
            MakeImmutable(i);
            
            if not KnowsDictionary(d,i) then
                Add( orb, i );
                AddDictionary(d,i);
                Add( orb_elts, g);
            fi;
        od;
    od;
    
    result := [Immutable(orb),Immutable(orb_elts)];
    
    return result;
    
    end );

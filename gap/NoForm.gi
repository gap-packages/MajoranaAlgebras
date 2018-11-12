
InstallGlobalFunction( MAJORANA_MainLoopNoForm,

    function(rep)

    MAJORANA_Fusion(rep, false);

    return MAJORANA_UnknownAlgebraProducts(rep, false);

    end);

InstallGlobalFunction( MajoranaRepresentationNoForm,

    function(arg)

    local   rep, unknowns, algebras, main;

    if Size(arg) = 2 then
        arg[3] := "AllAxioms";
    fi;

    rep :=  MAJORANA_SetUp(arg[1], arg[2], arg[3]);

    rep.innerproducts := false;

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        main := MAJORANA_MainLoopNoForm(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );

        if not false in rep.algebraproducts then

            MAJORANA_Fusion(rep, false);
            MAJORANA_IntersectEigenspaces(rep);

            Info( InfoMajorana, 10, "Success" );
            return rep;
        elif ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            Info( InfoMajorana, 10, "Fail" );
            rep.system := main;
            return rep;
        fi;
    od;

    end );

InstallGlobalFunction( NClosedMajoranaRepresentationNoForm,

    function(rep)

    local products, unknowns;

    products := Positions(rep.algebraproducts, false);

    MAJORANA_NClosedSetUp(rep, products[1]);

    while true do

        unknowns := Positions(rep.algebraproducts, false);

        MAJORANA_MainLoopNoForm(rep);

        Info(InfoMajorana, 20, STRINGIFY( "There are ", Size(Positions(rep.algebraproducts, false)), " unknown algebra products ") );

        if not false in rep.algebraproducts then
            Info( InfoMajorana, 10, "Success" );
            return;
        fi;

        if ForAll(rep.algebraproducts{unknowns}, x -> x = false) then
            products := Filtered(products, x -> rep.algebraproducts[x] = false);

            if products = [] then
                Info( InfoMajorana, 10, "Fail" );
                return;
            else
                MAJORANA_NClosedSetUp(rep, products[1]);
            fi;
        fi;
    od;

    end );

InstallGlobalFunction( MajoranaAlgebraTestNoForm,

    function(rep)

    MAJORANA_TestPrimitivity(rep);

    MAJORANA_TestFusion(rep);

    end );

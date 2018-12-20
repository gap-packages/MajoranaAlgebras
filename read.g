#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# Reading the implementation part of the package.
#
#

#! @ChapterInfo Introduction, Info levels
#! @Description The default info level of <A>InfoMajorana</A> is 0. No information
#! is printed at this level. If the info level is at least 10 then <A>Success</A> is printed if the
#! algorithm has produced a complete Majorana algebra, otherwise <A>Fail</A>
#! is printed. If the info level is at least 20 then more information is printed
#! about the progress of the algorithm, up to a maximum info level of 100.
DeclareInfoClass( "InfoMajorana" );

ReadPackage( "MajoranaAlgebras", "gap/SignedPerm.gi");

ReadPackage( "MajoranaAlgebras", "gap/OrbitalStructure.gi");

ReadPackage( "MajoranaAlgebras", "gap/DihedralAlgebras.gi");

ReadPackage( "MajoranaAlgebras", "gap/LinAlg.gi");

ReadPackage( "MajoranaAlgebras", "gap/Test.gi");

ReadPackage( "MajoranaAlgebras", "gap/SetUp.gi" );

ReadPackage( "MajoranaAlgebras", "gap/Shapes.gi" );

ReadPackage( "MajoranaAlgebras", "gap/MajoranaAlgebras.gi");

ReadPackage( "MajoranaAlgebras", "gap/Examples.gi");

ReadPackage( "MajoranaAlgebras", "gap/ExamplesPaper.gi");

ReadPackage( "MajoranaAlgebras", "gap/NClosed.gi");

ReadPackage( "MajoranaAlgebras", "gap/NoForm.gi");

ReadPackage( "MajoranaAlgebras", "gap/Orbits.gi");

ReadPackage( "MajoranaAlgebras", "gap/Embedding.gi");

ReadPackage( "MajoranaAlgebras", "gap/Irreducibles.gi");

ReadPackage( "MajoranaAlgebras", "gap/AxialAlgebras.gi");

ReadPackage( "MajoranaAlgebras", "gap/TauMaps.gi");

ReadPackage( "MajoranaAlgebras", "gap/Miscellaneous.gi");

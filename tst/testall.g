#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
SetInfoLevel(InfoPackageLoading, 4);
LoadPackage( "MajoranaAlgebras" );

TestDirectory(DirectoriesPackageLibrary( "MajoranaAlgebras", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error

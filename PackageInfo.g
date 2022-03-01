#
# MajoranaAlgebras: A package for constructing Majorana algebras and representations.
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "MajoranaAlgebras",
Subtitle := "A package for constructing Majorana algebras and representations",
Version := "1.4",
Date := "06/12/2018", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Markus",
    LastName := "Pfeiffer",
    WWWHome := "https://markusp.morphism.de/",
    Email := "markus.pfeiffer@st-andrews.ac.uk",
    PostalAddress := "School of Computer Science, University of St Andrews, North Haugh, KY16 9SX",
    Place := "St Andrews, UK",
    Institution := "University of St Andrews",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Madeleine",
    LastName := "Whybrow",
    WWWHome := "https://www.madeleinewhybrow.wordpress.com",
    Email := "mlw10@ic.ac.uk",
    PostalAddress := "Department of Mathematics, Imperial College, South Kensington, SW7 2AZ",
    Place := "London, UK",
    Institution := "Imperial College London",
  ),
],

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/MWhybrow92/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
SupportEmail   := "mlw10@ic.ac.uk",
PackageWWWHome  := "https://MWhybrow92.github.io/MajoranaAlgebras/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "MajoranaAlgebras is a package for constructing Majorana representations of finite groups. It also offers some functions to calculate with a constructed Majorana representation. The main constructive functions use the algorithm described in the preprint Constructing Majorana Representations (https://arxiv.org/abs/1803.10723) by Markus Pfeiffer and Madeleine Whybrow.",

PackageDoc := rec(
  BookName  := "MajoranaAlgebras",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "A package for constructing Majorana algebras and representations.",
),

Dependencies := rec(
  GAP := ">= 4.8",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.5" ]
                         , [ "automata", ">= 1.13"]
                         , [ "Gauss", ">=0" ]
                         , [ "datastructures", ">=0.2.2" ] ],
  SuggestedOtherPackages := [ [ "Char0Gauss", ">=0" ] ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

Keywords := [ "Majorana algebras", "axial algebras" ],

AutoDoc := rec(
    TitlePage := rec(
        Copyright :=
"""&copyright; 2018 by Markus Pfeiffer and Madeleine Whybrow<P/>
&MajoranaAlgebras; package is free software;
you can redistribute it and/or modify it under the terms of the
<URL Text="GNU General Public License">http://www.fsf.org/licenses/gpl.html</URL>
as published by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.""",
    ),
),

));

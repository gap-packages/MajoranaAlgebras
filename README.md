[![Build Status](https://travis-ci.org/MWhybrow92/MajoranaAlgebras.svg?branch=master)](https://travis-ci.org/MWhybrow92/MajoranaAlgebras)
[![Code Coverage](https://codecov.io/github/gap-system/gap/coverage.svg?branch=master&token=)](https://codecov.io/gh/MWhybrow92/MajoranaAlgebras)

# MajoranaAlgebras

This is a GAP package for constructing Majorana algebras of finite groups. The package implements the algorithm described in the paper "Constructing Majorana representations" (https://arxiv.org/abs/1803.10723) by M. Pfeiffer and M. Whybrow.

## Getting Started

To get the latest version of the package, either download the archive file `MajoranaAlgebras-master.zip` and inside the `pkg` subdirectory of your GAP installation or clone the MajoranaAlgebras repository by doing 

    git clone https://github.com/MWhybrow92/MajoranaAlgebras.git
    
inside the `pkg` subdirectory of your GAP installation.

Finally, start GAP in the usual way and type

    LoadPackage("MajoranaAlgebras");

## Constructing Majorana representations

Take a group `G` and a set `T` of involutions in `G` which generate `G` and are closed under conjugation. The set `T` must be a duplicate-free mutable list of involutions. 

    G := AlternatingGroup(5);;
    T := Filtered(AsList(G), x -> Order(x) = 2);;
    
To find the shapes of all possible Majorana representations `(G,T,V)`, use the function `ShapesOfMajoranaRepresentation`, or if you are interested in only shapes which obey axiom M8, use the function `ShapesOfMajoranaRepresentationAxiomM8`. The output of each of these functions is a record whose component `shapes` gives a list of all possible shapes. It is now up to the user to decide which shape they want to use in the construction.

    input := ShapesOfMajoranaRepresentation(G,T);;
    input.shapes;
    
The next step is to use the function `MajoranaRepresentation`. This takes two arguments, the first of which is the record `input` which was returned by the function `ShapesOfMajoranaRepresentation` and the second of which is the index of the desired shape in the list `input.shapes`.

    rep := MajoranaRepresentation(input, 1);;
    
This function returns the (potentially incomplete) representation. The function `MAJORANA_IsComplete`returns `true` if all algebra product values are known and `false` otherwise. 

    MAJORANA_IsComplete(rep);

If the representation is not complete then the user may use the function `NClosedMajoranaRepresentation` in order to attempt construction of the 3-closed algebra. This function takes as its argument the output of the function `MajoranaRepresentation`. 
    
    NClosedMajoranaRepresentation(rep);;

This function should be called `n-2` times in order to attempt construction of the `n`-closed algebra. 

Once a completed algebra has been found, the function `MAJORANA_Dimension` returns the dimension of the algebra. If all inner products have been found, but there are some missing algebra products then this value will be a lower bound on the true value of the dimension. 

    MAJORANA_Dimension(rep);

## Choice of axioms 

The function `MajoranaAlgebras` takes an optional third argument which is a string which make take one of the following values: `"AllAxioms"`, `"AxiomM8"`, `"NoAxioms"`. 

The string `"All axioms"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7, as well as the additional axioms M8, 3A, 4A and 5A as defined on page 8 of the paper "Constructing Majorana representations" (https://arxiv.org/abs/1803.10723). 

The string `"AxiomM8"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7, as well as the additional axioms M8 but no others. 

The string `"NoAxioms"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7 but no others. 

## Testing the output

Once a representation (complete or otherwise) has been constructed, the user may use the following functions to test properties of the algebra. Each function takes as its input the record returned by `MajoranaRepresentation`. Each function enters a break loop if the property does not hold in the algebra, returns true if the property can be determined to be true on the complete or incomplete algebra and fail otherwise.
    
    MAJORANA_TestAxiomM1(rep);
    MAJORANA_TestAxiomM2(rep);
    MAJORANA_TestFusion(rep);
    MAJORANA_TestPrimitivity(rep);
    
Alternatively, the function `MajoranaAlgebraTest` takes the same argument and runs all tests, with the exception of `MAJORANA_TestAxiomM2` which is a particularly expensive test.

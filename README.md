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

## Choice of axioms 

The function `MajoranaAlgebras` takes an optional third argument which is a string which make take one of the following values: `"AllAxioms"`, `"AxiomM8"`, `"NoAxioms"`. 

The string `"All axioms"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7, as well as the additional axioms M8, 3A, 4A and 5A as defined on page 8 of the paper "Constructing Majorana representations" (https://arxiv.org/abs/1803.10723). 

The string `"AxiomM8"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7, as well as the additional axioms M8 but no others. 

The string `"NoAxioms"` chooses a version of the algorithm which uses the Majorana axioms M1 - M7 but no others. 

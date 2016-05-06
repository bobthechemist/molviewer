(* ::Package:: *)

(* Mathematica Package *)

(* Created by the Wolfram Workbench Jul 27, 2013 *)

BeginPackage["Molviewer`"]
(* Exported symbols added here with SymbolName::usage *) 

makemolecule::usage = "
	makemolecule[atoms, bonds, bondtypes, coords]
	Constructs the graphics elements for a molecule"
molview::usage = "
	molview[graphics]
	Displays a grpahics object with default settings to make
	molecular viewing more appealing"
getsampledata::usage = "
	getsampledata[] 
	returns the structure of ethyl crotonate in an array of 
	{atoms, bonds, types, coords}"

	
Begin["`Private`"]
(* Implementation of the package *)

Clear[colorrules, sizerules]; 

colorrules[] = ColorData["Atoms","ColorRules"];
sizerules[] = DeleteCases[Table[ElementData[z,"Abbreviation"]->
	QuantityMagnitude@(ElementData[z,"AtomicRadius"]/2),{z,118}],x_/;Head[x[[2]]]==QuantityMagnitude];


Clear[makeatomcomplex];
makeatomcomplex[atoms_, coords_] := GraphicsComplex[coords, {
   MapThread[{#1 /. colorrules[], 
      Sphere[#2, #1 /. sizerules[]]} &, {atoms, 
     Range[Length[atoms]]}]}]

Clear[makebondcomplex]
makebondcomplex[atoms_, bonds_, bondtypes_, coords_, i_] := 
 Module[{b1, hb},
  hb   = Mean@coords[[#]] & /@ ({#[[1]], #[[2]]} & /@ bonds);
  b1 = coords[[bonds[[All, i]]]];
  {EdgeForm[None],
   GraphicsComplex[Join[b1, hb], {
     Riffle[atoms[[bonds[[All, i]]]] /. colorrules[],
      Table[GeometricTransformation[
          Cylinder[{#, Length@b1 + #}, 10],
          TranslationTransform[ k Normalize[hb[[#]]\[Cross]b1[[#]]]]
          ], {k, 
          bondtypes[[#]] /. {"Single" -> {0}, "Double" -> {-15, 15}, 
            "Aromatic" -> {-15, 15}, "Triple" -> {-15, 0, 15}}}] & /@ 
       Range[Length@b1]
      ]
     }]}]

Clear[addhighlight];
addhighlight[coords_, hilight_] :=
 GraphicsComplex[
  coords, {Opacity[0.6], Green, Sphere[#, 35] & /@ hilight}]

Clear[addlabels];
addlabels[coords_, labels_, offset_] :=
 GraphicsComplex[coords,
  {Black, MapThread[
    Text[#1, offset + #2,
      BaseStyle -> {FonstSize -> 12, FontFamily -> "Helvetica", 
        Orange}] &, {labels, coords}]}]

Clear[makemolecule];
makemolecule::badcoord = 
  "2D coordinates suspected, use PadRight[#,3]&/@coords";
Options[makemolecule] =
  {hilight -> {}, labels -> {}, labeloffset -> {1, 1, 1}};
makemolecule[{a_,b_,c_,d_}]:=makemolecule[{a,b,c,d}];
makemolecule[atoms_, bonds_, bondtypes_, coords_, OptionsPattern[]] := 
	{
  		makeatomcomplex[atoms, coords],
  		makebondcomplex[atoms, bonds, bondtypes, coords, #] & /@ {1, 2},
  		If[Length@OptionValue[hilight] > 0,
  			addhighlight[coords, OptionValue[hilight]]],
  		If[Length@OptionValue[labels] == Length@atoms,
  				addlabels[coords, OptionValue[labels], OptionValue[labeloffset]]]
  	}

Clear[molview];
molview[molecule_] := Graphics3D[{
   Specularity[GrayLevel[1], 100], molecule
   },
  Boxed -> False,
  Lighting -> "Neutral",
  SphericalRegion -> True, Background -> Black]
  
  
  
Clear[getsampledata];
getsampledata[]:=First@Import["623-70-1.sdf", #] & /@ 
	{"VertexTypes", "EdgeRules", "EdgeTypes", "VertexCoordinates"};
    
End[]

EndPackage[]


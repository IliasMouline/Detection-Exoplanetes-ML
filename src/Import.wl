(* ::Package:: *)

ImportKeplerLightCurve[chemin_String] := Module[
  {donnees, temps, flux, qualite, masqueValide},

  donnees = Quiet[Import[chemin, {"TableData", 2}]];

  If[!MatrixQ[donnees, NumericQ] && !MatrixQ[donnees],
    Return[Missing["ImportEchoue"]]
  ];
  If[Dimensions[donnees][[2]] < 10,
    Return[Missing["ImportEchoue"]]
  ];

  temps = donnees[[All, 1]];
  flux = donnees[[All, 8]];
  qualite = donnees[[All, 10]];

  masqueValide = MapThread[
    (NumericQ[#1] && NumericQ[#2] && #3 == 0) &,
    {temps, flux, qualite}
  ];

  temps = Pick[temps, masqueValide];
  flux = Pick[flux, masqueValide];

  If[Length[temps] < 10,
    Missing["DonneesInsuffisantes"],
    <|"Times" -> temps, "Flux" -> flux|>
  ]
];

ImportEtoileComplete[listeFichiers_List] := Module[
  {courbes, courbesValides, tempsTotal, fluxTotal, ordre},

  courbes = ImportKeplerLightCurve /@ listeFichiers;
  courbesValides = Select[courbes, AssociationQ];

  If[courbesValides === {},
    Return[Missing["AucuneDonneeValide"]]
  ];

  tempsTotal = Flatten[#["Times"] & /@ courbesValides];
  fluxTotal = Flatten[#["Flux"] & /@ courbesValides];

  ordre = Ordering[tempsTotal];

  <|"Times" -> tempsTotal[[ordre]], "Flux" -> fluxTotal[[ordre]]|>
];

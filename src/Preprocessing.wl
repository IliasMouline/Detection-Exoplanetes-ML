detecterSegments[temps_List, seuil_: 1.] :=
  Split[Range[Length[temps]], temps[[#2]] - temps[[#1]] < seuil &];

fusionnerSegments[segments_List, flux_List, tailleMin_: 50] :=
  Module[{segs = segments, i, tailles, medianes, meilleurVoisin,
    compteur = 0, maxIterations},
   maxIterations = 3 Length[segments];
   While[
    Length[Select[segs, Length[#1] < tailleMin &]] > 0 &&
     Length[segs] > 1 && compteur < maxIterations,
    compteur++;
    tailles = Length /@ segs;
    i = First[Position[tailles, x_ /; x < tailleMin, {1}, 1]][[1]];
    medianes = (Median[flux[[#1]]] &) /@ segs;
    Which[
     i == 1,
     segs = Join[{Join[segs[[1]], segs[[2]]]}, segs[[3 ;; All]]],
     i == Length[segs],
     segs = Join[segs[[1 ;; -3]], {Join[segs[[-2]], segs[[-1]]]}],
     True,
     meilleurVoisin =
      If[Abs[medianes[[i]] - medianes[[i - 1]]] <
        Abs[medianes[[i]] - medianes[[i + 1]]], i - 1, i + 1];
     If[meilleurVoisin == i - 1,
      segs = Join[segs[[1 ;; i - 2]], {Join[segs[[i - 1]], segs[[i]]]},
        segs[[i + 1 ;; All]]],
      segs = Join[segs[[1 ;; i - 1]], {Join[segs[[i]], segs[[i + 1]]]},
        segs[[i + 2 ;; All]]]
      ]
     ];
    ];
   If[compteur >= maxIterations,
    Print["fusionnerSegments : nombre max d'iterations atteint"]];
   segs
   ];

normaliserParSegments[temps_List, flux_List, segments_List] :=
  Module[{fluxNormalise},
   fluxNormalise = flux;
   Do[
    fluxNormalise[[segments[[s]]]] =
      flux[[segments[[s]]]]/Median[flux[[segments[[s]]]]];
    , {s, Length[segments]}];
   fluxNormalise
   ];

sigmaClip[flux_List, nSigma_: 4] :=
  Module[{mediane, ecartType, seuil, ecarts, masque},
   mediane = Median[flux];
   ecartType = StandardDeviation[flux];
   seuil = nSigma ecartType;
   ecarts = Abs[flux - mediane];
   masque = (#1 < seuil &) /@ ecarts;
   masque
   ];

detecterRupturesNiveau[flux_List, tailleFenetre_: 50, seuilSaut_: 0.01] :=
  Module[{medianesGlissantes, differences, positionsRupture},
   If[Length[flux] < 4*tailleFenetre, Return[{}]];
   medianesGlissantes =
    Table[Median[flux[[Max[1, i - tailleFenetre] ;; i]]],
     {i, tailleFenetre, Length[flux], tailleFenetre}];
   differences = Abs[Differences[medianesGlissantes]];
   positionsRupture = Position[differences, x_ /; x > seuilSaut, {1}];
   (Flatten[positionsRupture]*tailleFenetre) + tailleFenetre
   ];

redecouperSegmentParIndices[indicesSegment_List, fluxComplet_List,
   tailleFenetre_: 50, seuilSaut_: 0.01] :=
  Module[{fluxSeg, ruptures, debuts, fins, sousListesIndices},
   fluxSeg = fluxComplet[[indicesSegment]];
   ruptures = detecterRupturesNiveau[fluxSeg, tailleFenetre, seuilSaut];
   If[ruptures === {}, Return[{indicesSegment}]];
   debuts = Prepend[ruptures, 1];
   fins = Append[ruptures - 1, Length[indicesSegment]];
   sousListesIndices =
    Table[indicesSegment[[debuts[[k]] ;; fins[[k]]]], {k, 1,
      Length[debuts]}];
   sousListesIndices
   ];

diagnostiquerEtoile[kepid_, tempsMin_: 50, seuilGap_: 2.5,
   tailleFenetre_: 50, seuilSaut_: 0.01] :=
  Module[{temps, flux, segBruts, segFusionnes, fluxNorm,
    segmentsSuspects, plots},
   temps = datasetFinal[kepid]["Times"];
   flux = datasetFinal[kepid]["Flux"];
   segBruts = detecterSegments[temps, seuilGap];
   segFusionnes = fusionnerSegments[segBruts, flux, tempsMin];
   fluxNorm = normaliserParSegments[temps, flux, segFusionnes];
   segmentsSuspects =
    Select[Range[Length[segFusionnes]],
     detecterRupturesNiveau[fluxNorm[[segFusionnes[[#1]]]],
        tailleFenetre, seuilSaut] != {} &];
   Print["KIC ", kepid, " : ", Length[segFusionnes], " segments, ",
    Length[segmentsSuspects], " suspects : ", segmentsSuspects];
   If[segmentsSuspects != {},
    plots =
     Table[ListLinePlot[fluxNorm[[segFusionnes[[i]]]],
       PlotRange -> All, ImageSize -> 380,
       PlotLabel -> "Segment " <> ToString[i],
       Epilog -> {Red, PointSize[Medium],
         Point[{detecterRupturesNiveau[fluxNorm[[segFusionnes[[i]]]],
             tailleFenetre, seuilSaut],
            fluxNorm[[segFusionnes[[i]]]][[
             detecterRupturesNiveau[fluxNorm[[segFusionnes[[i]]]],
              tailleFenetre, seuilSaut]]]}]}], {i, segmentsSuspects}];
    Print[GraphicsGrid[Partition[plots, UpTo[3]], ImageSize -> 1200]]
    ];
   Association[Temps -> temps, Flux -> flux,
    SegFusionnes -> segFusionnes, FluxNorm -> fluxNorm,
    SegmentsSuspects -> segmentsSuspects]
   ];

finaliserEtoile[diagnostic_Association, segmentsAvalider_List, nSigma_: 4] :=
  Module[{temps, flux, segFusionnes, fluxNorm, segmentsFinaux,
    fluxRenorm, masqueClip, tempsClip, fluxClip},
   temps = diagnostic[Temps];
   flux = diagnostic[Flux];
   segFusionnes = diagnostic[SegFusionnes];
   fluxNorm = diagnostic[FluxNorm];
   segmentsFinaux =
    Flatten[Table[
      If[MemberQ[segmentsAvalider, i],
       redecouperSegmentParIndices[segFusionnes[[i]], fluxNorm, 50,
        0.01], {segFusionnes[[i]]}], {i, Length[segFusionnes]}], 1];
   fluxRenorm = normaliserParSegments[temps, flux, segmentsFinaux];
   masqueClip = sigmaClip[fluxRenorm, nSigma];
   tempsClip = Pick[temps, masqueClip];
   fluxClip = Pick[fluxRenorm, masqueClip];
   Print["Segments : ", Length[segFusionnes], " -> ",
    Length[segmentsFinaux], " | Points retirés : ",
    Count[masqueClip, False], "/", Length[fluxRenorm]];
   Association[Times -> tempsClip, Flux -> fluxClip,
    SegmentsFinaux -> segmentsFinaux]
   ];

replierPhase[temps_List, flux_List, period_, time0bk_] := Module[
  {phase, phaseCentree, ordre, phaseTriee, fluxTrie},
  phase = Mod[temps - time0bk, period];
  phaseCentree = Map[If[# > period/2, # - period, #] &, phase];
  ordre = Ordering[phaseCentree];
  phaseTriee = phaseCentree[[ordre]];
  fluxTrie = flux[[ordre]];
  <|"Phase" -> phaseTriee, "Flux" -> fluxTrie|>
];

binnerVueLocale[phase_List, flux_List, duree_, periode_: Infinity, nbPoints_: 201, facteurFenetre_: 3] :=
  Module[{demiFenetre, masqueFenetre, phaseFenetre, fluxFenetre,
    largeurBin, indicesBins, groupes, fluxBinne, phaseBinCentres, bornesBins},
   demiFenetre = Min[facteurFenetre*duree, periode/2];
   masqueFenetre = Thread[Abs[phase] <= demiFenetre];
   phaseFenetre = Pick[phase, masqueFenetre];
   fluxFenetre = Pick[flux, masqueFenetre];
   If[Length[phaseFenetre] < nbPoints, Return[Missing["PasAssezDePointsDansLaFenetre"]]];
   largeurBin = (2*demiFenetre)/nbPoints;
   indicesBins = Clip[Floor[(phaseFenetre + demiFenetre)/largeurBin] + 1, {1, nbPoints}];
   groupes = GroupBy[Transpose[{indicesBins, fluxFenetre}], First -> Last];
   fluxBinne = Table[If[KeyExistsQ[groupes, i], Median[groupes[i]], Missing["BinVide"]], {i, nbPoints}];
   bornesBins = Subdivide[-demiFenetre, demiFenetre, nbPoints];
   phaseBinCentres = MovingAverage[bornesBins, 2];
   <|"PhaseBin" -> phaseBinCentres, "FluxBin" -> fluxBinne|>
   ];

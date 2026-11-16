(* UNTESTED: written without access to Wolfram Language in the authoring environment. *)
(* Cross-check of one rule: edge/vertex counts, beta1 and planarity of the final state. *)
rule = {{x, y}, {x, z}} -> {{x, y}, {y, z}, {y, w}};
init = EdgeList[GridGraph[{8, 8}]] /. UndirectedEdge[a_, b_] :> {Min[a, b], Max[a, b]};
beta1[edges_] := Module[{g = Graph[UndirectedEdge @@@ edges]},
  Length[edges] - VertexCount[g] + Length[ConnectedComponents[g]]];
Table[
  With[{s = ResourceFunction["WolframModel"][rule, init, t, "FinalState"]},
    {t, Length[s], Length[DeleteDuplicates[Flatten[s]]], beta1[s],
     PlanarGraphQ[Graph[DeleteDuplicates[UndirectedEdge @@@ Select[s, #[[1]] =!= #[[2]] &]]]]}],
  {t, 0, 6}] // TableForm
(* Note: WolframModel's default update order differs from the Python script's; compare invariants only. *)

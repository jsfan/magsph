open Vectorop;;
(* module Particle = struct *)
  type properties = {
      (* physical properties *)
      r: Vectorop.space;
      v: Vectorop.space;
      b: Vectorop.space;
      uth: float;
      e: float;
      rho: float;
      p: float;
      cs: float;
      m: float;
      h: float;
      alpha: float;
      (* rates of change *)
      f: Vectorop.space;
      dbdt: Vectorop.space;
      sij: Vectorop.tensor;
      drhodt: float;
      dudt: float;
      dedt: float;
      dalphadt: float;      
    }    

(* checks if 2 particles are neighbours, returns particle in a one element list if neighbour *)
  let check_neighbour part1 part2 =
    let dist = Vectorop.norm (Vectorop.diff part1.r part2.r) in
    if ((dist < part1.h) && (dist == 0.)) then
      part2 :: []
    else
      []
(* end *)

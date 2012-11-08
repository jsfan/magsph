open Vectorop;;
open Particle;;

let rec init_particles gamma pmass particles cpos distvec bound1 bound2 =
  let pi = 3.14159265 and
  cpos =
    if cpos = (zerovec cpos) then
      sum bound1 (scale 0.5 distvec)
    else
      cpos in
  let uthglob = 1. and bglob = r1 0. in
  let csloc = sqrt (gamma *. (gamma - 1.) * uthglob + bglob**2) in
  let vloc = r1 (0.025 *. csloc *. sin (scale (2. *. pi) cpos)) in
  let densloc = 1. +. 0.025 *. sin (scale (2. *. pi) cpos);
  let particle = { r = cpos;
		   v = vloc;
		   b = bglob;
		   uth = uthglob;
		   e = 0.5 *. dotprod vloc vloc +. uthglob +. 0.5 *. dotprod bglob bglob /. densloc;
		   rho = densloc;
		   p = 0.;
		   cs = csloc;
		   m = pmass;
		   h = hfac *.;
		   alpha = 0.;
		   f = r1 0.;
		   dbdt = r1 0.;
		   sij = tensor_r1 0.;
		   drhodt = 0.;
		   dudt = 0.;
		   dedt = 0.;
		   dalphadt =  0.} in
  let newparticles = particle :: particles in
  match bound1,bound2,distvec,cpos with
    Space1D b1,Space1D b2,Space1D dvec, Space1D cloc ->
      if ((cloc.x1 +. dvec.x1) < b2.x1) then
	let newloc = Space1D { x1 = cloc.x1 +. dvec.x1 } in
	init_particles newparticles newloc (Space1D dvec) bound1 bound2
      else
	newparticles
  | Space2D b1, Space2D b2, Space2D dvec, Space2D cloc ->
      if ((cloc.x2 +. dvec.x2) < b2.x2) then
	let newloc = Space2D { cloc with x2 = cloc.x2 +. dvec.x2 } in
	init_particles newparticles newloc (Space2D dvec) bound1 bound2
      else
	if ((cloc.y2 +. dvec.y2) < b2.y2) then
	  let newloc = Space2D { x2 = b1.x2; y2 = cloc.y2 +. dvec.y2 } in
	  init_particles newparticles newloc (Space2D dvec) bound1 bound2
	else
	  newparticles
(*  | Space3D zerovec, Space3D cpos -> *)
  | _ -> raise Not_implemented



let rec init gamma boundaries n =
  let cond = match boundaries with
    xl :: xr :: [] -> (* 1D *)
      let dist = (abs_float (xr -. xl)) /. float_of_int n and pmass = 1. /. float_of_int n in
      let pdist = r1 dist in
      let boundl = r1 xl in
      let boundr = r1 xr in
      (pdist,boundl,boundr)
  | xl :: xr :: yt :: yb :: [] -> (* 2D *)
      let dist = sqrt ((abs_float ((xr -. xl) *. (yt -. yb))) /. float_of_int n) and pmass = 1. /. float_of_int n in
      let pdist = r2 dist dist in
      let bound1 = r2 xl yb in
      let bound2 = r2 xr yt in
      (pdist,bound1,bound2)
(*  | xl :: xr :: yt :: yb :: zf :: zb :: [] -> (* 3D *)
      let dist = (xr -. xl) /. n in
       *)
  | _ -> raise Not_implemented
  in
  let particle_separation = (fun (psep,_,_) -> psep) cond in
  let boundary1 = (fun (_,bound,_) -> bound) cond in
  let boundary2 = (fun (_,_,bound) -> bound) cond in
  init_particles gamma pmass particle_separation [] (zerovec particle_separation) particle_separation boundary1 boundary2

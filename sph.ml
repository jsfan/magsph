open Vectorop;;
open Particle;;
(* module SPH = struct *)
(* calculate stress tensor for dv/dt *)
  let stresstensor ~pressure:p ~magneticfield:magfield =
    match magfield with
      (Space1D b) -> tensor_r1 (-.p +. b.x1 *. b.x1 -. 0.5 *. dotprod magfield magfield)
    | (Space2D b) -> tensor_r2 (-.p +. b.x2 *. b.x2 -. 0.5 *. dotprod magfield magfield) (b.x2 *. b.y2) (-.p +. b.y2 *. b.y2 -. 0.5 *. dotprod magfield magfield)
    | (Space3D b) -> tensor_r3 (-.p +. b.x3 *. b.x3 -. 0.5 *. dotprod magfield magfield) (b.x3 *. b.y3) (b.x3 *. b.z3) (-.p +. b.y3 *. b.y3 -. 0.5 *. dotprod magfield magfield) (b.y3 *. b.z3) (-.p +. b.z3 *. b.z3 -. 0.5 *. dotprod magfield magfield)
    | _ -> raise Not_implemented

(* gradient of the kernel *)
  let gradw ~vector1:rvec1 ~vector2:rvec2 ~average_h:hbar =
    let rrel = norm (diff rvec1 rvec2) in
    let unnormalised =
      match (rvec1,rvec2) with
	(Space1D ra, Space1D rb) -> r1 (-.3. *. rrel /. (hbar ** 2.) +. 1.25 *. (rrel ** 2.) /. (hbar ** 3.))
      | (Space2D ra, Space2D rb) -> r2 (-.3. /. (hbar ** 2.) *. (ra.x2 -. rb.x2) +. 1.25 *. rrel /. (hbar ** 3.) *. (ra.x2 -. rb.x2)) (-.3. /. (hbar ** 2.) *. (ra.y2 -. rb.y2) +. 1.25 *. rrel /. (hbar ** 3.) *. (ra.y2 -. rb.y2))
      | (Space3D ra, Space3D rb) -> r3 (-.3. /. (hbar ** 2.) *. (ra.x3 -. rb.x3) +. 1.25 *. rrel /. (hbar ** 3.) *. (ra.x3 -. rb.x3)) (-.3. /. (hbar ** 2.) *. (ra.y3 -. rb.y3) +. 1.25 *. rrel /. (hbar ** 2.) *. (ra.y3 -. rb.y3)) (-.3. /. (hbar ** 2.) *. (ra.z3 -. rb.z3) +. 1.25 *. rrel /. (hbar ** 2.) *. (ra.z3 -. rb.z3))
      | _ -> raise Not_implemented
    in let normfac = match unnormalised with
      (Space1D wu) -> 2. /. 3.
    | (Space2D wu) -> 10. /. (2. *. 3.14159265)
    | (Space3D wu) -> 1. /. 3.14159265
    | _ -> raise Not_implemented
    in scale normfac unnormalised

(* calculate the average of 2 floats *)
  let average val1 val2 = 0.5 *. (val1 +. val2)

(* calculate the signal velocity *)
  let vsig ~particle1:part1 ~particle2:part2 = 
    let beta = 1. in
    let rvec1 = part1.r in
    let rvec2 = part2.r in
    let rrel = norm (diff rvec1 rvec2) in
    sqrt (part1.cs ** 2. +. (norm part1.b) ** 2.) +. sqrt (part2.cs ** 2. +. (norm part2.b) ** 2.) -. beta *. norm (diff part1.v part2.v) /. rrel


(* ghost particle routines *)
(* determine distances to boundaries *)
  let boundary_distances particle boundaries =
    let r = particle.r in
    match (r, boundaries) with
      (Space1D vec, leftb :: rightb :: otherblist) -> 
	[vec.x1 -. leftb; rightb -. vec.x1]
    | (Space2D vec, leftb :: rightb :: bottomb :: topb :: otherblist) -> 
	[vec.x2 -. leftb; rightb -. vec.x2; vec.y2 -. bottomb; topb -. vec.y2]
    | (Space3D vec, leftb :: rightb :: bottomb :: topb :: frontb :: backb :: []) ->
	[vec.x3 -. leftb; rightb -. vec.x3; vec.y3 -. bottomb; topb -. vec.y3; vec.z3 -. frontb; backb -. vec.z3]
    | _,_ -> raise Not_implemented

(* determine ghostparticles of one particle *)
  let find_ghostparticles ~particle:particle maxh boundaries =
    let loc = particle.r
    (* 2*maxh is the limit for ghostparticles *)
    and maxh2 = 2. *. maxh in
    let bounddist = boundary_distances particle boundaries in
    match bounddist,loc,boundaries with
      (xl :: xr :: [],Space1D rvec,leftb :: rightb :: []) (* 1D *) ->
	let ghostp = if xl < maxh2 then [{ particle with r = Space1D { x1 = rightb +. xl } }] else [] in
	if xr < maxh2 then { particle with r = Space1D { x1 = leftb -. xr } } :: ghostp else ghostp
    | (xl :: xr :: yb :: yt :: [],Space2D rvec,leftb :: rightb :: bottomb :: topb :: []) (* 2D *) ->
	let ghostp = if xl < maxh2 then [{ particle with r = Space2D { x2 = rightb +. xl; y2 = rvec.y2 } }] else [] in
	let ghostp = if xr < maxh2 then { particle with r = Space2D { x2 = leftb -. xr; y2 = rvec.y2 } } :: ghostp else ghostp in
	let ghostp = if yb < maxh2 then { particle with r = Space2D { x2 = rvec.x2; y2 = topb +. yb } } :: ghostp else ghostp in
	let ghostp = if yt < maxh2 then { particle with r = Space2D { x2 = rvec.x2; y2 = bottomb -. yt } } :: ghostp else ghostp in
	let ghostp = if (norm (Space2D {x2 = xl; y2 = yb})) < maxh2 then { particle with r = Space2D { x2 = rightb +. xl; y2 = topb +. yb } } :: ghostp else ghostp in
	let ghostp = if (norm (Space2D {x2 = xr; y2 = yb})) < maxh2 then { particle with r = Space2D { x2 = leftb -. xr; y2 = topb +. yb } } :: ghostp else ghostp in
	let ghostp = if (norm (Space2D {x2 = xl; y2 = yt})) < maxh2 then { particle with r = Space2D { x2 = rightb +. xr; y2 = bottomb -. yt } } :: ghostp else ghostp in		
	if (norm (Space2D {x2 = xr; y2 = yt})) < maxh2 then { particle with r = Space2D { x2 = leftb -. xr; y2 = bottomb -. yt } } :: ghostp else ghostp
    (* | xl :: xr :: yb :: yt :: zf :: zb :: [] (* 3D *) -> not implemented, yet *)
    | _ -> raise Not_implemented

(* compile list of all ghostparticles *)
  let rec get_ghostparticles particlelist ?(ghostparticles = []) maxh boundaries =
    match particlelist with
      particle :: partlst ->
	let newghostparticles = find_ghostparticles particle maxh boundaries in
	get_ghostparticles partlst ~ghostparticles:(newghostparticles @ ghostparticles) maxh boundaries
    | [] -> ghostparticles

(* shift particles back into the computational domain *)
  let rec shift_particles ?(shifted = []) particles boundaries =
    match particles,boundaries with
      particle :: remparticles,bound ->
	begin
	  let newr = begin
	    match particle.r,bound with
	      Space1D loc, xl :: xr :: [] ->
		if loc.x1 < xl then
		  Space1D { x1 = xr -. (xl -. loc.x1) }
		else
		  if loc.x1 > xr then
		    Space1D { x1 = xl +. (loc.x1 -. xr) }
		else
		    Space1D loc

	    | Space2D loc2, xl :: xr :: yb :: yt :: [] ->
		if loc2.x2 < xl then
		  Space2D { loc2 with x2 = xr -. (xl -. loc2.x2) }
		else
		  let Space2D newloc = if loc2.x2 > xr then
		    Space2D { loc2 with x2 = xl +. (loc2.x2 -. xr) }
		  else
		    Space2D loc2 in
		  if loc2.y2 < yb then
		    Space2D { newloc with x2 = yt -. (yb -. loc2.x2) }
		  else
		    if loc2.y2 > yt then
		      Space2D { newloc with y2 = xl +. (loc2.x2 -. xr) }
		    else		  
		      Space2D newloc
(*	| Space3D loc, xl :: xr :: yb :: yt :: zf :: zb :: [] -> *)
	    | _,_ -> raise Not_implemented
	  end in
	  shift_particles ~shifted:({ particle with r = newr } :: shifted) remparticles boundaries
	end
     | [],bound -> shifted
  

(* routines to find neighbours *)
(* wrapper for check_neighbour *)
  let rec check_particles particle partlist neighbours =
    match partlist with
      [] -> neighbours
    | testpart :: rempart ->
	let newneighbour = check_neighbour particle testpart in
	check_particles particle rempart newneighbour

(* get the list of all neighbours of a particle *)
  let get_neighbours ~particle:particle ~all_particles:all_particles =
    check_particles particle all_particles []


(* rates of change *)

(* get the force one particle exerts on another *)
  let pairforce sigma ~particle_a:ptcl ~neighbour:neighbour =
    let first_term part1 part2 = tensorsum (scaletensor (1. /. (part1.rho ** 2.)) part1.sij) (scaletensor (1. /. (part2.rho ** 2.)) part2.sij) in
    let sigma = (average ptcl.h neighbour.h)/.(vsig ptcl neighbour) in
    (sigma,scale neighbour.m (tensorprod (first_term ptcl neighbour) (gradw ptcl.r neighbour.r (average ptcl.h neighbour.h))))
	  
(* get new force for one particle *)
(* parameters are one particle and list of its neighbours *)
  let get_force ~particle_a:particle ~list_of_neighbours:neighbours ~sigma:sigma = 
    let sumforces sigma part prelimforce neighbour = 
      let (sigma,force) = pairforce sigma part neighbour in
      (sigma,sum prelimforce force) in
    Helpers.sphfold sigma (fun sigma cur -> if cur > sigma then cur else sigma) (sumforces sigma particle) (zerovec particle.f) neighbours
      
(* get new density for one particle *)
(* parameters are one particle and list of its neighbours *)
  let get_rate_density ~particle_a:particle ~list_of_neighbours:neighbours = 
    let pair_drhodt ~particle_a:ptcl ~neighbour:neighbour =
      neighbour.m *. (dotprod (diff ptcl.v neighbour.v) (gradw ptcl.r neighbour.r (average ptcl.h neighbour.h))) in
    let sum_drhodt part prelim_drhodt neighbour = prelim_drhodt +. (pair_drhodt part neighbour) in
    List.fold_left (sum_drhodt particle) 0.0 neighbours

(* wrappers for rates routines *)
  let update_forces sigma allparticles maxh boundaries particlelist particle =
    let ghostparticles = get_ghostparticles allparticles maxh boundaries in
    let neighbours = get_neighbours particle (allparticles @ ghostparticles) in
    let (new_sigma,new_force) = get_force particle neighbours sigma in
    (new_sigma,{ particle with f = new_force } :: particlelist)

  let update_rate_density allparticles maxh boundaries particlelist particle =
    let ghostparticles = get_ghostparticles allparticles maxh boundaries in
    let neighbours = get_neighbours particle (allparticles @ ghostparticles) in    
    let new_drhodt = get_rate_density particle neighbours in
    { particle with drhodt = new_drhodt } :: particlelist


(* integrator *)

(* integrate forward by half a time step *)
let halfstep dt vec =
  scale (0.5 *. dt) vec

(* integrate forward by a full time step *)
let fullstep dt vec =
  scale dt vec

(* predictor step *)
  let predict dt gamma particlelist particle =
    let rnew = sum particle.r (halfstep dt particle.v) in
    let vnew = sum particle.v (halfstep dt particle.f) in
    let uthnew = particle.uth *. (1. +. 0.5 *. dt *. (gamma -. 1.) *. particle.drhodt /. particle.rho) in
    let hnew = particle.h *. (1. -. 0.5 *. dt *. particle.drhodt /. particle.rho) in
    let rhonew = particle.rho +. 0.5 *. dt *. particle.drhodt in
    let pnew = (gamma -. 1.) *. uthnew *. rhonew in
    let csnew = sqrt (gamma *. pnew /. rhonew) in
    (hnew,{ particle with r = rnew; v = vnew; p = pnew; cs = csnew; h = hnew; f = zerovec particle.f; dbdt = zerovec particle.dbdt; drhodt = 0.; dudt = 0.; dedt = 0.; dalphadt = 0. } :: particlelist)

(* corrector step *)
  let correct_locvel dt particlelist particle =
    let vnew = sum particle.v (fullstep dt particle.f) in
    let rnew = sum particle.r (fullstep dt vnew) in
    { particle with r = rnew; v = vnew; dbdt = zerovec particle.dbdt; drhodt = 0.; dudt = 0.; dedt = 0.; dalphadt = 0. } :: particlelist

  let correct dt gamma particlelist particle =
    let rhonew = particle.rho +. 0.5 *. dt *. particle.drhodt in
    let uthnew = particle.uth /. (1. -. 0.5 *. dt *. (gamma -. 1.) *. particle.drhodt /. rhonew) in
    let hnew = particle.h /. (1. +. 0.5 *. dt *. particle.drhodt /. rhonew) in
    (hnew,{ particle with uth = uthnew; rho = rhonew; h = hnew; dbdt = zerovec particle.dbdt; drhodt = 0.; dudt = 0.; dedt = 0.; dalphadt = 0. } :: particlelist)

(* end *)

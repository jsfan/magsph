open Vectorop;;
open Particle;;
open Sph;;
open Init;;
open Io;;

let rec timeloop maxstep boundaries t dt step gamma sigma maxh particles =
  let (maxh,particles_predicted) = Helpers.sphfold maxh ((fun maxh cur -> if cur > maxh then cur else maxh)) (predict dt gamma) [] particles in
  let (new_sigma,particles_newforce) = Helpers.sphfold sigma (fun sigma cur -> if cur > sigma then cur else sigma) (update_forces sigma particles maxh boundaries) [] particles_predicted in
  let particles_newlocvel = List.fold_left (correct_locvel dt) [] particles_newforce in
  let particles_shifted = shift_particles particles_newlocvel boundaries in
  let particles_newdrhodt = List.fold_left (update_rate_density particles_shifted maxh boundaries) [] particles_shifted in
  let (maxh,newparticles) = Helpers.sphfold maxh ((fun maxh cur -> if cur > maxh then cur else maxh)) (correct dt gamma) [] particles_newdrhodt
  and time = t +. dt and newstep = step + 1 and dt = 0.5 *. new_sigma in
  if step > maxstep then
    (newparticles,time)
  else
    let _ = gnuplot_output newparticles step in
    timeloop maxstep boundaries time dt newstep gamma sigma maxh newparticles;;

let sigma = 1000. in
let maxh = 0. in
let step = 0 in
let gamma = 4./.3. and n = 150 in
let particles = init gamma [0.;1.] n and maxstep = 100 in
timeloop maxstep [0.;1.] 0. 0. step gamma sigma maxh particles

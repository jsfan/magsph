(* I/O routines *)
open Printf;;
open Vectorop;;
open Particle;;

let print_vector chan vector =
  match vector with
    Space1D vec -> fprintf chan "%.10f" vec.x1
  | Space2D vec -> fprintf chan "%.10f %.10f" vec.x2 vec.y2
  | Space3D vec -> fprintf chan "%.10f %.10f %.10f" vec.x3 vec.y3 vec.z3

let rec write_line out particles step =
  match particles with
    particle :: moreparticles -> (* print line of values *)
      print_vector out particle.r;
      fprintf out " ";
      print_vector out particle.v;      
      fprintf out " ";
      print_vector out particle.b;
      fprintf out " ";
      fprintf out "%.10f" particle.uth;
      fprintf out " ";
      fprintf out "%.10f" particle.p;
      fprintf out " ";
      fprintf out "%.10f" particle.rho;
      fprintf out " ";
      fprintf out "%.10f" particle.h;
      fprintf out " ";
      fprintf out "%.10f" particle.alpha;
      fprintf out " ";
      fprintf out "%.10f" particle.cs;
      fprintf out "\n";
      write_line out moreparticles step
  | [] -> close_out out


let rec gnuplot_output particles step =
  let filename = String.concat "-" ["sphout";sprintf "%0.6d" step] in
  let out = open_out filename in
  fprintf out "# r v b uth p rho h alpha cs step\n";
  write_line out particles step

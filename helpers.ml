exception Unknown_type

(* take a pair and return first element *)
let firstel tup =
  match tup with
    (fst,_) -> fst
(*  | _ -> raise Unknown_type *)

(* take a pair and return first element *)
let secondel tup =
  match tup with
    (_,snd) -> snd
(*  | _ -> raise Unknown_type *)


(* folds left while updating a certain parameter according to sidefun *)
let rec sphfold side sidefun predfun dest start =
  match start with
    [] -> (side,dest)
  | lsthead :: newstart ->
      let newvals = predfun dest lsthead in
      let comp =  firstel newvals and newdest = secondel newvals in
      let newside = sidefun side comp in
      sphfold newside sidefun predfun newdest newstart

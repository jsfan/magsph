(* module Vectorop = struct *)
  exception Not_implemented
  type space1d = { x1: float }
  type space2d = { x2: float; y2: float }
  type space3d = { x3: float; y3: float; z3: float }
  type space = Space1D of space1d | Space2D of space2d | Space3D of space3d
  let r1 a = Space1D { x1 = a }
  let r2 a b = Space2D { x2 = a; y2 = b }
  let r3 a b c = Space3D { x3 = a; y3 = b; z3 = c }
  type tensor1d = { xx1: float }
  type tensor2d = { xx2: float; xy2: float; yy2: float }
  type tensor3d = { xx3: float; xy3: float; xz3: float; yy3: float; yz3: float; zz3: float }
  type tensor = Tensor1D of tensor1d | Tensor2D of tensor2d | Tensor3D of tensor3d
  let tensor_r1 a = Tensor1D { xx1 = a }
  let tensor_r2 a b c = Tensor2D { xx2 = a; xy2 = b; yy2 = c}
  let tensor_r3 a b c d e f = Tensor3D { xx3 = a; xy3 = b; xz3 = c; yy3 = d; yz3 = e; zz3 = f }
  let sum v1 v2 = 
    match (v1,v2) with
      (* vector sums *)
      (Space1D vec1, Space1D vec2) -> r1 (vec1.x1 +. vec2.x1)
    | (Space2D vec1, Space2D vec2) -> r2 (vec1.x2 +. vec2.x2) (vec1.y2 +. vec2.y2)
    | (Space3D vec1, Space3D vec2) -> r3 (vec1.x3 +. vec2.x3) (vec1.y3 +. vec2.y3) (vec1.z3 +. vec2.z3)
    | _ -> raise Not_implemented

  let tensorsum t1 t2 =
        match (t1,t2) with
      (* tensor sums *)
      (Tensor1D tensor1, Tensor1D tensor2) -> tensor_r1 (tensor1.xx1 +. tensor2.xx1)
    | (Tensor2D tensor1, Tensor2D tensor2) -> tensor_r2 (tensor1.xx2 +. tensor2.xx2) (tensor1.xy2 +. tensor2.xy2) (tensor1.yy2 +. tensor2.yy2)
    | (Tensor3D tensor1, Tensor3D tensor2) -> tensor_r3 (tensor1.xx3 +. tensor2.xx3) (tensor1.xy3 +. tensor2.xy3) (tensor1.xz3 +. tensor2.xz3) (tensor1.yy3 +. tensor2.yy3) (tensor1.yz3 +. tensor2.yz3) (tensor1.zz3 +. tensor2.zz3)
    | _ -> raise Not_implemented

  let dotprod v1 v2 = 
    match (v1,v2) with
      (* scalar product *)
      (Space1D vec1, Space1D vec2) -> vec1.x1 *. vec2.x1
    | (Space2D vec1, Space2D vec2) -> vec1.x2 *. vec2.x2 +. vec1.y2 *. vec2.y2
    | (Space3D vec1, Space3D vec2) -> vec1.x3 *. vec2.x3 +. vec1.y3 *. vec2.y3 +. vec1.z3 *. vec2.z3
    | _ -> raise Not_implemented

  let tensorprod t v =
    match (t,v) with
      (* tensor product with vector *)
      (Tensor1D tens, Space1D vec) -> r1 (tens.xx1 *. vec.x1)
    | (Tensor2D tens, Space2D vec) -> r2 (tens.xx2 *. vec.x2 +. tens.xy2 *. vec.y2) (tens.xy2 *. vec.x2 +. tens.yy2 *. vec.y2)
    | (Tensor3D tens, Space3D vec) -> r3 (tens.xx3 *. vec.x3 +. tens.xy3 *. vec.y3 +. tens.xz3 *. vec.z3) (tens.xy3 *. vec.x3 +. tens.yy3 *. vec.y3 +. tens.yz3 *. vec.z3) (tens.xz3 *. vec.x3 +. tens.yz3 *. vec.y3 +. tens.zz3 *. vec.z3)
    | _ -> raise Not_implemented

  let scale (f: float) (v1: space) =
    match v1 with
      (* scale 1D vector with a scalar *)
      (Space1D vec1) -> r1 (vec1.x1 *. f)
      (* scale 2D vector with a scalar *)
    | (Space2D vec1) -> r2 (vec1.x2 *. f) (vec1.y2 *. f)
      (* scale 3D vec1tor with a scalar *)
    | (Space3D vec1) -> r3 (vec1.x3 *. f) (vec1.y3 *. f) (vec1.z3 *. f)

  let scaletensor (f: float) (tensor1: tensor) =
    match tensor1 with
      (* scale a 1D tensor *)
      (Tensor1D t1) -> tensor_r1 (t1.xx1 *. f)  
      (* scale a 2D tensor *)
    | (Tensor2D t1) -> tensor_r2 (t1.xx2 *. f) (t1.xy2 *. f) (t1.yy2 *. f)	  
      (* scale a 3D tensor *)
    | (Tensor3D t1) -> tensor_r3 (t1.xx3 *. f) (t1.xy3 *. f) (t1.xz3 *. f) (t1.yy3 *. f) (t1.yz3 *. f) (t1.zz3 *. f)

  let diff v1 v2 = 
    sum v1 (scale (-1.) v2)

  let crossprod v1 v2 =
    match (v1,v2) with
      (* vector product *)
      (Space3D vec1, Space3D vec2) -> r3 (vec1.y3 *. vec2.z3 -. vec2.y3 *. vec1.z3) (vec1.z3 *. vec2.x3 -. vec1.x3 *. vec2.z3) (vec1.x3 *. vec2.y3 -. vec1.y3 *. vec2.x3)
    | _ -> raise Not_implemented
  
  let norm v1 =
    sqrt (dotprod v1 v1)

  let zerovec sample =
    match sample with
      (* scale 1D vector with a scalar *)
      (Space1D sample) -> r1 (0.)
      (* scale 2D vector with a scalar *)
    | (Space2D sample) -> r2 (0.) (0.)
      (* scale 3D vec1tor with a scalar *)
    | (Space3D sample) -> r3 (0.) (0.) (0.)

(* end *)

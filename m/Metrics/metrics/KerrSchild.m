(* ::Package:: *)

Module[
  {ld, lu, h, eta, regAdd, reg0, reg2, reg4, reg},
  ld={1,
      (r[x,y,z] x+a y)/(r[x,y,z]^2+a^2),
      (r[x,y,z] y-a x)/(r[x,y,z]^2+a^2),
      z/r[x,y,z]};
  lu={-1,
      (r[x,y,z] x+a y)/(r[x,y,z]^2+a^2),
      (r[x,y,z] y-a x)/(r[x,y,z]^2+a^2),
      z/r[x,y,z]};
  h=(M r[x,y,z])/(r[x,y,z]^2+a^2 z^2/r[x,y,z]^2) // Simplify;
  
  eta={{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
  (* Various regularisation methods *)
  (* Note: Don't use an explicit case distinction (If, Piecewise)
     here; this leads to complex code that Kranc is not able to handle
     well. Instead, use Max or similar functions. *)
  regAdd[r_] := Sqrt[epsilon^2 + r^2];
  reg0[r_] := Max[r, epsilon];
  reg2[r_] := Max[r, (epsilon^2 + r^2)/(2 epsilon)];
  reg4[r_] := Max[r, (3 epsilon^4 + 6 epsilon^2 r^2 - r^4)/(8 epsilon^3)];
  reg = reg4;   (* Choose *)
  
  {
    "Name" -> "KerrSchild",
    "Description" -> "Kerr spacetime in Kerr-Schild coordinates",
    "Dimensions" -> 4,
    "Coordinates" -> {t, x, y, z},
    "Parameters" -> {M, a, epsilon},
    "Metric" -> eta + 2 h Outer[Times, ld, ld],
    "InverseMetric" -> eta - 2 h Outer[Times, lu, lu],
    "SignDet" -> -1,
    "Shorthands" -> {
      r[x,y,z] -> reg[Sqrt[(x^2+y^2+z^2-a^2+
                            Sqrt[(x^2+y^2+z^2-a^2)^2+4 a^2 z^2])/2]]},
    "SimplifyHints" -> {
      (4 a^2 z^2 + (-a^2 + x^2 + y^2 + z^2)^2)^n_ :>
      (2 r[x,y,z]^2 + a^2 - x^2 - y^2 - z^2)^(2 n),
      (r[x,y,z]^m_)^n_ :> r[x,y,z]^(m n),
      (r[x,y,z]^n_ /; n>=4 :>
       r[x,y,z]^(n-4) (a^2 z^2 - r[x,y,z]^2 (a^2 - x^2 - y^2 - z^2))),
      (2 a^2 Z^2 + (-a^2 + x^2 + y^2 + z^2) r[x, y, z]^2 ->
       a^2 z^2 + r[x, y, z]^4)}
  }]

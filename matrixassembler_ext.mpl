# Matrices assembling for noncommutative geometry-inspired
# dirty black-holes' quasi-normal modes computation
# The extreme case
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local f::function, g::function, xe::numeric, mu::numeric, a::numeric, fp::function, gp::function, V::function, F::function, Fp::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Determination of the event horizon:
  f := (xe, mu) -> 1 - 1/xe*erf(mu*xe) + 2*mu/sqrt(Pi)*exp(-mu^2*xe^2):
  sol := fsolve({f(xe,mu) = 0, diff(f(xe,mu),xe) = 0}, {xe = 0.79, mu = 1.90}):
  assign(sol); # we assign the values of critical parameters
  printf("Found BH's event horizon      xe = %f\n", xe);
  printf("Found BH's dimensionless mass mu = %f\n", mu);
  # Physical parameters:
  a := 1 + 1/6*xe*(4*mu^4*xe^4 - 3*mu^3*xe^3 - 4*mu^2*xe^2 + 3*mu*xe - 4)*exp(-1/2*mu*(xe - 1))/(mu^2*xe^2 - 1)^2:
  # Other parameters appearing in the model:
  f := y -> 1 - 1/(2*xe)*(1 - y)*erf(2*mu*xe/(1 - y)) + 2*mu*exp(-4*mu^2*xe^2/(1 - y)^2)/sqrt(Pi):
  fp := y -> erf(2*mu*xe/(1 - y))/(2*xe) - 2*exp(-4*mu^2*xe^2/(1 - y)^2)*mu/((1 - y)*sqrt(Pi)) - 16*mu^3*xe^2*exp(-4*mu^2*xe^2/(1 - y)^2)/((1 - y)^3*sqrt(Pi)):
  g := y -> (1/2*mu*(1 - erf(2*mu*xe/(1 - y)) + 4*mu*xe*exp(-4*mu^2*xe^2/(1 - y)^2)/((1 - y)*sqrt(Pi)))):
  gp := y -> (-16*mu^4*xe^3*exp(-4*mu^2*xe^2/(1 - y)^2)/((1 - y)^4*sqrt(Pi))):
  F := y -> f(y)*exp(-g(y)):
  Fp := y -> fp(y)*exp(-g(y)) - f(y)*gp(y)*exp(-g(y)):
  V := y -> (1 - y)^2/16*f(y)*exp(-2*g(y))*((-s^2 + 1)*(1 - y)*fp(y) - (1 - s)*(1 - y)*f(y)*gp(y) + L*(L + 1)):
  eta := y -> (1 + y)/(1 - y) + (1 - y)*exp(-1/2*mu*(xe - 1))/((mu^2*xe^2 - 1)*(1 + y)):
  etap := y -> 1/(1 - y) + (1 + y)/(1 - y)^2 - exp(-mu*(xe - 1)/2)/((mu^2*xe^2 - 1)*(1 + y)) - (1 - y)*exp(-mu*(xe - 1)/2)/((mu^2*xe^2 - 1)*(1 + y)^2):
  etapp := y -> 2/(1 - y)^2 + 2*(1 + y)/(1 - y)^3 + 2*exp(-mu*(xe - 1)/2)/((mu^2*xe^2 - 1)*(1 + y)^2) + 2*(1 - y)*exp(-mu*(xe - 1)/2)/((mu^2*xe^2 - 1)*(1 + y)^3):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*V(y)/(-y^2 + 1)^2:
  L01 := y -> F(y)*(-2*(1 - y)*F(y) + (1 - y)^2*Fp(y))/(4*(1 + y)^2):
  L02 := y -> (1 - y)^2/(4*(1 + y)^2)*F(y)^2:
  L10 := y -> xe/2*(1 - y)^2*F(y)*Fp(y)*etap(y)/(1 + y)^2 + xe/2*(1 - y)*F(y)^2*((1 - y)*etapp(y) - 2*etap(y))/(1 + y)^2 + 1/2*(1 - y)*(2 - a*(1 - y))*F(y)*Fp(y)/(1 + y)^3 - F(y)^2*(2 - 2*a*(1 - y) + 1/2*a*(1 - y)^2)/(1 + y)^4:
  L11 := y -> (1 - y)*F(y)^2*(xe*(-y^2 + 1)*etap(y) + 2 - a*(1 - y))/(1 + y)^3:
  L12 := y -> 0:
  L20 := y -> 4*xe^2/(-y^2 + 1)^2 - F(y)^2*(xe*(-y^2 + 1)*etap(y) + 2 - a*(1 - y))^2/(1 + y)^4:
  L21 := y -> 0:
  L22 := y -> 0:
  P := y -> simplify(add(a[j]*ChebyshevT(j, y), j=0..n-1)):
  M0 := Matrix(n):
  M1 := Matrix(n):
  M2 := Matrix(n):
  for i from 1 to n do
    xi := cos((2.0*i-1.0)*Pi/(2.0*n)); # Chebyshev roots collocation points
    expr0 := evalf(L00(xi)*P(xi) + L01(xi)*subs(x=xi, diff(P(x),x)) + L02(xi)*subs(x=xi, diff(P(x),x$2))):
    expr1 := evalf(L10(xi)*P(xi) + L11(xi)*subs(x=xi, diff(P(x),x)) + L12(xi)*subs(x=xi, diff(P(x),x$2))):
    expr2 := evalf(L20(xi)*P(xi) + L21(xi)*subs(x=xi, diff(P(x),x)) + L22(xi)*subs(x=xi, diff(P(x),x$2))):
    for j from 1 to n do
      M0[i,j] := coeff(expr0, a[j-1]):
      M1[i,j] := coeff(expr1, a[j-1]):
      M2[i,j] := coeff(expr2, a[j-1]):
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(n, string);
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc:
# Matrices assembling for noncommutative geometry-inspired
# dirty black-holes' quasi-normal modes computation
# The non-extreme case
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  mu::numeric, # mu : dimensionless BH's mass
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local f::function, g::function, xh::numeric, alpha::numeric, fp::function, gp::function, V::function, F::function, Fp::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Determination of the event horizon:
  f := x -> 1 - erf(mu*x)/x + 2*mu/sqrt(Pi)*exp(-mu^2*x^2):
  xh := fsolve(f(x) = 0, x = 0.98):
  printf("Found BH's event horizon xh = %f\n", xh);
  # Physical parameters:
  alpha := exp(-1/2*mu*(xh - 1))/(1 - 4*mu^3*xh^2*exp(-mu^2*xh^2)/sqrt(Pi)):
  # Other parameters appearing in the model:
  f := y -> 1 - 1/(2*xh)*(1 - y)*erf(2*mu*xh/(1 - y)) + 2*mu*exp(-4*mu^2*xh^2/(1 - y)^2)/sqrt(Pi):
  fp := y -> erf(2*mu*xh/(1 - y))/(2*xh) - 2*exp(-4*mu^2*xh^2/(1 - y)^2)*mu/((1 - y)*sqrt(Pi)) - 16*mu^3*xh^2*exp(-4*mu^2*xh^2/(1 - y)^2)/((1 - y)^3*sqrt(Pi)):
  g := y -> (1/2*mu*(1 - erf(2*mu*xh/(1 - y)) + 4*mu*xh*exp(-4*mu^2*xh^2/(1 - y)^2)/((1 - y)*sqrt(Pi)))):
  gp := y -> (-16*mu^4*xh^3*exp(-4*mu^2*xh^2/(1 - y)^2)/((1 - y)^4*sqrt(Pi))):
  F := y -> f(y)*exp(-g(y)):
  Fp := y -> fp(y)*exp(-g(y)) - f(y)*gp(y)*exp(-g(y)):
  V := y -> (1 - y)^2/16*f(y)*exp(-2*g(y))*((-s^2 + 1)*(1 - y)*fp(y) - (1 - s)*(1 - y)*f(y)*gp(y) + L*(L + 1)):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*V(y)/((1 + y)*(1 - y)^2):
  L01 := y -> F(y)*(-2*(1 - y)*F(y) + (1 - y)^2*Fp(y))/(4 + 4*y):
  L02 := y -> (1 - y)^2*F(y)^2/(4*(1 + y)):
  L10 := y -> F(y)*(Fp(y)*(2*xh + 1 - y) - F(y) + xh*alpha*(1 - y)*((3 + y)*F(y) - (-y^2 + 1)*Fp(y))/(1 + y)^2)/(2*(1 + y)):
  L11 := y -> F(y)^2*((2*xh + 1 - y)/(1 + y) - xh*alpha*(1 - y)^2/(1 + y)^2):
  L12 := y -> 0:
  L20 := y -> 4*xh^2/((1 + y)*(1 - y)^2) - ((1 + y)*(2*xh + 1 - y) - xh*alpha*(1 - y)^2)^2*F(y)^2/((1 + y)^3*(1 - y)^2):
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
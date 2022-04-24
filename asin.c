#define POLY2(x, c1, c0) MLA(x, C2V(c1), C2V(c0))
#define POLY4(x, x2, c3, c2, c1, c0) MLA(x2, MLA(x, C2V(c3), C2V(c2)), MLA(x, C2V(c1), C2V(c0)))
#define POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0) MLA(x4, POLY4(x, x2, c7, c6, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY12(x, x2, x4, x8, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY4(x, x2, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
//3.1415926218032836914
double xasin(double d) {
  int o = fabsk(d) < 0.5;
  double x2 = o ? (d*d) : ((1-fabsk(d))*0.5), x = o ? fabsk(d) : SQRT(x2), u;

  double x4 = x2 * x2, x8 = x4 * x4, x16 = x8 * x8;
  u = POLY12(x2, x4, x8, x16,
	     +0.3161587650653934620e-1,
	     -0.1581918243329996649e-1,
	     +0.1929045477267910672e-1,
	     +0.6606077476277170612e-2,
	     +0.1215360525577377332e-1,
	     +0.1388715184501609213e-1,
	     +0.1735956991223614600e-1,
	     +0.2237176181932048349e-1,
	     +0.3038195928038132230e-1,
	     +0.4464285681377102431e-1,
	     +0.7500000000378581610e-1,
	     +0.1666666666666497542e+0);

  u = mla(u, x * x2, x);
  
  double r = o ? u : (M_PI/2 - 2*u);
  r = mulsign(r, d);

  return r;
}

#include <armadillo>

using namespace std;
using namespace arma;

double val(vec x)
{
	int row = 3, n = x.n_elem, i, j;
	int col = n / 3;

	mat geo(n, 1); geo = x;
	geo.reshape(row, col);
	

	double value, a = 0, b = 0, temp;

	for (i = 0; i < col; i++)
	{
		for (j = 0; j < i; j++) {

			if (j != i)
			{
				temp = norm(geo.col(j) - geo.col(i));
				a = a + pow(temp, -6 );
				b = b + pow(temp, -12);
			}
		}
	}
	value = b - 2 * a;

	return value;
}
vec Diff(vec x)
{
	int row = 3, n = x.n_elem,i,j;
	int col = n / 3;

	double r;

	vec out,ddf=zeros<vec>(3),temp;

	mat geo(n, 1); geo = x;
	
	geo.reshape(row, col);

	
	for (j = 0; j < col; j++) {
		if (j != 0)
		{
			temp = geo.col(j) - geo.col(0);
			r = norm(temp);
			ddf = ddf + 12 * (temp*pow(r, -8) - temp*pow(r, -14));
		}
	}
	out = ddf;

	for (i = 1; i < col; i++)
	{   
		ddf.zeros();
		for (j = 0; j < col; j++) {
			if (j != i)
			{   
				temp = geo.col(j) - geo.col(i);
				r = norm(temp);
				ddf = ddf + 12 * (temp*pow(r, -8) - temp*pow(r, -14));
			}
		}
		out = join_cols(out, ddf);
	}
	return out;
}

double cubeminpoint(double x1,double x2, double f1,double f2,double g1,double g2)
{
	mat A(4, 4);
	vec B(4);
	vec X(4);

	double a, b, c, d, p1, p2, ff1, ff2, min, minp;
	
	A(0, 0) = x1*x1*x1; A(0, 1) = x1*x1; A(0, 2) = x1; A(0, 3) = 1;
	A(1, 0) = x2*x2*x2; A(1, 1) = x2*x2; A(1, 2) = x2; A(1, 3) = 1;
	A(2, 0) =  3*x1*x1; A(2, 1) =  2*x1; A(2, 2) =  1; A(2, 3) = 0;
	A(3, 0) =  3*x2*x2; A(3, 1) =  2*x2; A(3, 2) =  1; A(3, 3) = 0;

	B(0) = f1; B(1) = f2; B(2) = g1; B(3) = g2;

	X = solve(A, B);
	
	a = X(0); b = X(1); c = X(2); d = X(3);

	p1 = ( -b + sqrt(b*b - 3*a*c)) / (3*a);
	p2 = ( -b - sqrt(b*b - 3*a*c)) / (3*a);

	ff1 = a*p1*p1*p1 + b*p1*p1 + c*p1 +d;
	ff2 = a*p2*p2*p2 + b*p2*p2 + c*p2 +d;

	

	min = f1 < f2 ? f1 : f2;
	minp= f1 < f2 ? x1 : x2;

	if (ff1 < min) {
		if ((x1 < p1) && (p1 < x2) || (x2 < p1) && (p1 < x1)) {
			min = ff1;
			minp= p1;
		}
	}
	if (ff2 < min) {
		if ((x1 < p2) && (p2 < x2) || (x2 < p2) && (p2 < x1)) {
			min = ff2;
			minp= p2;
		}
	}

	return minp;
}
double cubemin(double x1, double x2, double f1, double f2, double g1, double g2)
{
	mat A(4, 4);
	vec B(4);
	vec X(4);

	double a, b, c, d, p1, p2, ff1, ff2, min, minp;

	A(0, 0) = x1*x1*x1; A(0, 1) = x1*x1; A(0, 2) = x1; A(0, 3) = 1;
	A(1, 0) = x2*x2*x2; A(1, 1) = x2*x2; A(1, 2) = x2; A(1, 3) = 1;
	A(2, 0) = 3 * x1*x1; A(2, 1) = 2 * x1; A(2, 2) = 1; A(2, 3) = 0;
	A(3, 0) = 3 * x2*x2; A(3, 1) = 2 * x2; A(3, 2) = 1; A(3, 3) = 0;

	B(0) = f1; B(1) = f2; B(2) = g1; B(3) = g2;

	X = solve(A, B);

	a = X(0); b = X(1); c = X(2); d = X(3);

	p1 = (-b + sqrt(b*b - 3 * a*c)) / (3 * a);
	p2 = (-b - sqrt(b*b - 3 * a*c)) / (3 * a);

	ff1 = a*p1*p1*p1 + b*p1*p1 + c*p1 + d;
	ff2 = a*p2*p2*p2 + b*p2*p2 + c*p2 + d;



	min = f1 < f2 ? f1 : f2;
	minp = f1 < f2 ? x1 : x2;

	if (ff1 < min) {
		if ((x1 < p1) && (p1 < x2) || (x2 < p1) && (p1 < x1)) {
			min = ff1;
			minp = p1;
		}
	}
	if (ff2 < min) {
		if ((x1 < p2) && (p2 < x2) || (x2 < p2) && (p2 < x1)) {
			min = ff2;
			minp = p2;
		}
	}

	return min;
}
double intpolfgf(double x1, double x2, double f1,  double g1, double f2 )
{
	mat A(3, 3);
	vec B(3);
	vec X(3);

	double a, b, c, p, ff, min, minp;

	A(0, 0) = x1*x1; A(0, 1) = x1; A(0, 2) = 1; 
	A(1, 0) =  2*x1; A(1, 1) =  1; A(1, 2) = 0; 
	A(2, 0) = x2*x2; A(2, 1) = x2; A(2, 2) = 1; 
	

	B(0) = f1; B(1) = g1; B(2) = f2;

	X = solve(A, B);

	a = X(0); b = X(1); c = X(2); 

	p = -b/(2*a);
	ff = a*p*p + b*p + c;

	 min = f1 < f2 ? f1 : f2;
	minp = f1 < f2 ? x1 : x2;

	if (ff < min) {
		if ((x1 < p) && (p < x2) || (x2 < p) && (p < x1)) {
			min = ff;
			minp = p;
		}
	}
	
	return minp;
}
double intpolfgg(double x1, double x2, double f1, double g1, double g2)
{
	mat A(3, 3);
	vec B(3);
	vec X(3);

	double a, b, c,  p, ff, min, minp ,f2;

	A(0, 0) = x1*x1; A(0, 1) = x1; A(0, 2) = 1;
	A(1, 0) = 2 *x1; A(1, 1) = 1; A(1, 2) = 0;
	A(2, 0) = 2 *x2; A(2, 1) = 1; A(2, 2) = 0;


	B(0) = f1; B(1) = g1; B(2) = g2;

	X = solve(A, B);

	a = X(0); b = X(1); c = X(2);

	p = -b / (2 * a);
	ff = a*p*p + b*p + c;
	f2 = a*x2*x2 + b*x2 + c;

	min = f1 < f2 ? f1 : f2;
	minp = f1 < f2 ? x1 : x2;

	if (ff < min) {
		if ((x1 < p) && (p < x2) || (x2 < p) && (p < x1)) {
			min = ff;
			minp = p;
		}
	}

	return minp;
}

double phi(vec x, vec dir, double a)
{
	return val(x + a*dir);
}
double dphi(vec x, vec dir, double a)
{
	return dot(dir, Diff(x + a*dir));
}
double bphi(vec x, vec dir, double a,double u)
{
	double out=val(x + a*dir) - val(x) - u*a*dot(dir,Diff(x));
	return out;
}
double dbphi(vec x, vec dir, double a,double u)
{
	return dphi(x,dir,a)-u*dphi(x,dir,0);
}
double getminp(vec x, vec dir, double al, double au)
{
	double minv, delta;
	int i, xu;
	vec vvv(9);
	delta = (au - al) / 10;
	while (delta>0.0001)
	{
		for (i = 0; i < 9; i++) { vvv(i) = phi(x, dir, al + (i + 0.5)*delta); }
		minv = min(vvv);
		xu = index_min(vvv);
		al = al + (xu + 0.5)*delta - 0.5*delta;
		au = al + (xu + 0.5)*delta + 0.5*delta;
		delta = (au - al) / 10;

	}
	return (al + au) / 2;
}

vec Linesearch(vec x, vec d, double u, double mu, double amax)
{
	vec x_next;
	
	double al, al_next, au, au_next;
	double at, at_next, atemp;
	double fl, ft, fu, gl, gt, gu;
	al = 0;
	au = amax;
	at = amax / 2;
	while((phi(x, d, at)>(phi(x, d, 0) +u*at*dphi(x, d, 0)))||(abs(dphi(x, d, at))>mu*abs(dphi(x, d,0))))
	{
		    if ((bphi(x, d, at,u)<=0)&&(dphi(x, d, at)>=0))
			{
				if (phi(x, d, at) > phi(x, d, al)) { al_next = al; au_next = at; }
				else if ((dphi(x, d, at)*(al - at)) >= 0) { al_next = at; au_next = au; }
				else { al_next = at; au_next = al; }
			    gl = dphi(x, d, al_next);gt = dphi(x, d, at); gu = dphi(x, d, au_next);
			}
			else
			{ 
				if (bphi(x, d, at,u) > bphi(x, d, al,u)) { al_next = al; au_next = at; }
				else if ((dbphi(x, d, at,u)*(al - at)) >= 0) { al_next = at; au_next = au; }
				else { al_next = at; au_next = al; }
			    gl = dbphi(x, d, al_next,u); gt = dbphi(x, d, at,u); gu = dbphi(x, d, au_next,u);
			}

			fl = phi(x, d, al_next); ft = phi(x, d, at); fu = phi(x, d, au_next);

	        double ac = cubeminpoint(al, at, fl, ft, gl, gt);
	        double aq = intpolfgf(al, at, fl, ft, gl);
	        double as = intpolfgf(al, at, fl, gl, gt);
	        
	        if (ft > fl)
	        {
	        	if (abs(ac - al) < abs(aq - al)) { at_next = ac; }
	        	else { at_next = (aq + ac) / 2; }
	        }
	        else if (gt*gl < 0)
	        {
	        	if (abs(ac - at) >= abs(as - at)) { at_next = ac; }
	        	else { at_next = as; }
	        }
	        else if (abs(gt)<abs(gl))
	        {
	        	if (abs(ac - at) < abs(as - at)) { at_next = ac; }
	        	else { at_next = as; }
	        
	        	atemp = at_next + 0.66*(au - at);
	        	if (at > al) { at_next = atemp < at_next ? atemp : at_next; }
	        	else { at_next = atemp > at_next ? atemp : at_next; }
	        }
	        else
	        {
	        	at_next = cubeminpoint(au, at, fu, ft, gu, gt);
	        }
			al = al_next;
			at = at_next;
			au = au_next;
	}
	at = getminp(x,d, al, au);
	x_next = x + at*d;
	return x_next;
}

field<vec> cauchy(vec x, vec l, vec u, vec g, double theta, mat W, mat M)
{
	field<vec> Fout(2, 1);

	int N = W.n_rows;
	int m2 = W.n_cols;
	int i, b;
	vec t(N), d(N), eb(N), xout(N);
	vec p(m2), c = zeros<vec>(m2), omiga(m2);
	uvec F, btemp;

	double f, ff, deltatmin, deltat, told, tt;
	double xb, zb, gb;
	for (i = 0; i < N; i++)
	{
		if (g[i]<0) { t[i] = (x[i] - u[i]) / g[i]; }
		else if (g[i]>0) { t[i] = (x[i] - l[i]) / g[i]; }
		else { t[i] = INFINITY; }

		if (t[i] == 0) { d[i] = 0; }
		else { d[i] = -g[i]; }
	}

	//===========Initialize========

	p = W.t()*d;
	f = -dot(d, d);
	ff = -theta*f - dot(p, M*p);
	deltatmin = -f / ff;
	told = 0;
	F = find(t>0);
	tt = min(t.elem(F));
	btemp = find(t == tt);
	b = btemp[0];
	deltat = tt;

	//==========Processing============
	while (deltatmin>deltat) {
		if (d[b] >= 0) { xb = u[b]; }
		else { xb = l[b]; }
		zb = xb - x[b];
		x[b] = xb;
		gb = g[b];
		c = c + deltat*p;
		eb[b] = 1;
		omiga = W.t()*eb;
		eb[b] = 0;
		f = f + deltat*ff + gb*gb + theta*gb*zb - gb*dot(omiga, M*c);//vector multiply
		ff = ff - theta*gb*gb - 2 * gb*dot(omiga, M*p) - gb*gb*dot(omiga, M*omiga);//vector multiply
		p = p + gb*omiga;
		d[b] = 0;
		deltatmin = -f / ff;
		told = tt;
		t[b] = 0;
		F = find(t>0);
		tt = min(t.elem(F));
		btemp = find(t == tt);
		b = btemp[0];
		deltat = tt - told;
	}
	deltatmin = deltatmin > 0 ? deltatmin : 0;
	told = told + deltatmin;
	xout = x + told*d;
	c = c + deltatmin*p;

	Fout(0, 0) = xout;
	Fout(1, 0) = c;

	return Fout;
}
vec cauchypoint(vec x, vec l, vec u, vec g, double theta, mat W, mat M)
{
	int N = W.n_rows;
	int m2 = W.n_cols;
	int i, b;
	vec t(N), d(N), eb(N), xout(N);
	vec p(m2), c = zeros<vec>(m2), omiga(m2);
	uvec F, btemp;

	double f, ff, deltatmin, deltat, told, tt;
	double xb, zb, gb;
	for (i = 0; i < N; i++)
	{
		if (g[i]<0) { t[i] = (x[i] - u[i]) / g[i]; }
		else if (g[i]>0) { t[i] = (x[i] - l[i]) / g[i]; }
		else { t[i] = INFINITY; }

		if (t[i] == 0) { d[i] = 0; }
		else { d[i] = -g[i]; }
	}

	//===========Initialize========

	p = W.t()*d;
	f = -dot(d, d);
	ff = -theta*f - dot(p, M*p);
	deltatmin = -f / ff;
	told = 0;
	F = find(t);
	tt = min(t.elem(F));
	btemp = find(t == tt);
	b = btemp[0];
	deltat = tt;

	//==========Processing============
	while (deltatmin>deltat) {
		if (d[b] >= 0) { xb = u[b]; }
		else { xb = l[b]; }
		zb = xb - x[b];
		x[b] = xb;
		gb = g[b];
		c = c + deltat*p;
		eb[b] = 1;
		omiga = W.t()*eb;
		eb[b] = 0;
		f = f + deltat*ff + gb*gb + theta*gb*zb - gb*dot(omiga, M*c);//vector multiply
		ff = ff - theta*gb*gb - 2 * gb*dot(omiga, M*p) - gb*gb*dot(omiga, M*omiga);//vector multiply
		p = p + gb*omiga;
		d[b] = 0;
		deltatmin = -f / ff;
		told = tt;
		t[b] = 0;
		F = find(t);
		tt = min(t.elem(F));
		btemp = find(t == tt);
		b = btemp[0];
		deltat = tt - told;
	}
	deltatmin = deltatmin > 0 ? deltatmin : 0;
	told = told + deltatmin;
	xout = x + told*d;

	return xout;
}
vec Cauchydual(mat Y, mat S, vec x, vec l, vec u, vec g, double theta)
{
	int n = Y.n_rows, m = Y.n_cols, m2; m2 = 2 * m;
	int i, j, b,ta;

	double f, ff, deltatmin, deltat, told, tt;
	double xb, zb, gb,a;

	vec t(n), d(n), eb = zeros<vec>(n);
	vec p(m2), c = zeros<vec>(m2), omiga(m2);

	mat W(n, m2), M(m2, m2);
	mat ZERO = zeros<mat>(m, m), D = zeros<mat>(m, m), R = zeros<mat>(m, m), L = zeros<mat>(m, m);
	mat H = eye<mat>(m2, m2), I = eye<mat>(n, n),A;

	uvec F, btemp,FA;

	vec dstar, xout,lstar,xc;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (i <= j) {
				R(i, j) = dot(S.col(i), Y.col(j));
				if (i == j) { D(i, j) = dot(S.col(i), Y.col(j)); }
			}
			else { L(i, j) = dot(S.col(i), Y.col(j)); }
		}
	}


	W = join_rows(Y, theta*S);
	M = inv_sympd(join_cols(join_rows(-D, L.t()), join_rows(L, theta*S.t()*S)));

	for (i = 0; i < n; i++)
	{
		if (g[i]<0) { t[i] = (x[i] - u[i]) / g[i]; }
		else if (g[i]>0) { t[i] = (x[i] - l[i]) / g[i]; }
		else { t[i] = INFINITY; }

		if (t[i] == 0) { d[i] = 0; }
		else { d[i] = -g[i]; }
	}

	//===========Initialize========

	p = W.t()*d;
	f = -dot(d, d);
	ff = -theta*f - dot(p, M*p);
	deltatmin = -f / ff;
	told = 0;
	F = find(t>0);
	tt = min(t.elem(F));
	btemp = find(t == tt);
	b = btemp[0];
	deltat = tt;

	//==========Processing============
	while (deltatmin>deltat) {
		if (d[b] >= 0) { xb = u[b]; }
		else { xb = l[b]; }
		zb = xb - x[b];
		x[b] = xb;
		gb = g[b];
		c = c + deltat*p;
		eb[b] = 1;
		omiga = W.t()*eb;
		eb[b] = 0;
		f = f + deltat*ff + gb*gb + theta*gb*zb - gb*dot(omiga, M*c);//vector multiply
		ff = ff - theta*gb*gb - 2 * gb*dot(omiga, M*p) - gb*gb*dot(omiga, M*omiga);//vector multiply
		p = p + gb*omiga;
		d[b] = 0;
		deltatmin = -f / ff;
		told = tt;
		t[b] = 0;
		F = find(t>0);
		tt = min(t.elem(F));
		btemp = find(t == tt);
		b = btemp[0];
		deltat = tt - told;
	}
	deltatmin = deltatmin > 0 ? deltatmin : 0;
	told = told + deltatmin;
	xc = x + told*d;
	c = c + deltatmin*p;

	ta = n - F.n_elem;

	R = -inv(R);
	W = join_rows(Y / theta, S);
	M = join_cols(join_rows(ZERO, R), join_rows(R.t(), R.t()*(D + Y.t()*Y / theta)*R));
	H = H / theta + W*M*W.t();

	if (ta == 0) {
		dstar = -H*g;
	}
	else {
		FA = find(t <= 0);
		A = I.cols(FA);
		lstar = inv(A.t()*H*A)*A*(H*g + xc - x);
		dstar = -H*(A*lstar + g);
	}
	vec gg, aa;
	gg = x + dstar - xc;
	for (i = 0; i < n; i++)
	{
		if (gg[i]<0) { aa[i] = (l[i] - xc[i]) / gg[i]; }
		else if (gg[i]>0) { aa[i] = (u[i] - xc[i]) / gg[i]; }
		else { aa[i] = INFINITY; }

	} 
	a = min(aa.elem(F));
	xout = xc + a*gg;
	return xout;
}

double Amax(vec xc, vec l, vec u, vec gg)
{
	int n = xc.n_elem,i;
	vec aa;
	for (i = 0; i < n; i++)
	{
		if (gg[i]<0) { aa[i] = (l[i] - xc[i]) / gg[i]; }
		else if (gg[i]>0) { aa[i] = (u[i] - xc[i]) / gg[i]; }
		else { aa[i] = INFINITY; }

	}
	return min(aa);
}
vec Proj(vec x, vec l, vec u)
{
	int n = x.n_elem, i;
	vec out;
	for (i = 0; i < n; i++)
	{
		if (x[i] < l[i]) { out[i] = l[i]; }
		else if (x[i] > u[i]) { out[i] = u[i]; }
		else { out[i] = x[i]; }
	}
	return out;
}

vec twoloop(mat Y, mat S, vec g)
{
	int n = Y.n_cols,i;
	vec q; vec a; vec z; vec b; vec R;
	q = g;
	for (i = n - 1; i >= 0; i--)
	{
		R(i) = 1 / dot(Y.col(i), S.col(i));
		a(i) = R(i)*dot(S.col(i),q);
		q = q - a(i)*Y.col(i);
	}
	z = dot(Y.col(n - 1), S.col(n - 1)) / dot(Y.col(n - 1), Y.col(n - 1))*q;
	for (i = 0; i < n; i++)
	{
		b(i) = R(i)*dot(Y.col(i), z);
		z = z - (a(i)-b(i))*S.col(i);
	}
	return -z;
}
vec oneloop(mat Y, mat S, vec g)
{
	int n = Y.n_rows, m = Y.n_cols, m2; m2 = 2 * m;
	double theta;
	int i, j;
	mat W(n, m2), M(m2, m2);
	mat ZERO = zeros<mat>(m, m), D = zeros<mat>(m, m), R = zeros<mat>(m, m), L = zeros<mat>(m, m);
	mat H = eye<mat>(m2, m2), I = eye<mat>(n, n), A;

	theta = dot(Y.col(m-1), Y.col(m - 1)) / dot(S.col(m - 1), Y.col(m - 1));

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (i <= j) {
				R(i, j) = dot(S.col(i), Y.col(j));
				if (i == j) { D(i, j) = dot(S.col(i), Y.col(j)); }
			}
			else { L(i, j) = dot(S.col(i), Y.col(j)); }
		}
	}

	R = -inv(R);
	W = join_rows(Y / theta, S);
	M = join_cols(join_rows(ZERO, R), join_rows(R.t(), R.t()*(D + Y.t()*Y / theta)*R));
	H = H / theta + W*M*W.t();

	
	return -H*g;
}


mat Activeset(vec x, vec u, vec l)
{
	int n = x.n_elem;
	mat out(n, n, fill::zeros);
	int i;

	for (i = 0; i < n; i++)
	{
		if (x(i) > u(i) || x(i) < l(i))
		{
			out(i, i) = 1;
		}
	}
	return out;
}
mat Activeset(vec x, vec u, vec l,double eps)
{
	int n = x.n_elem;
	mat out(n, n, fill::zeros);
	int i;

	for (i = 0; i < n; i++)
	{
		if (x(i) > u(i)-eps || x(i) < l(i) + eps)
		{
			out(i, i) = 1;
		}
	}
	return out;
}
mat Inactive(vec x, vec u, vec l)
{
	int n = x.n_elem;
	mat out(n, n, fill::eye);
	int i;

	for (i = 0; i < n; i++)
	{
		if (x(i) > u(i) || x(i) < l(i))
		{
			out(i, i) = 0;
		}
	}
	return out;
}
mat Inactive(vec x, vec u, vec l,double eps)
{
	int n = x.n_elem;
	mat out(n, n, fill::eye);
	int i;

	for (i = 0; i < n; i++)
	{
		if (x(i) > u(i)-eps|| x(i) < l(i)+eps)
		{
			out(i, i) = 0;
		}
	}
	return out;
}

vec BFGS(vec x)
{
	//m is the store vector number;
	//n is the para number 
	//x is the start point ;

	int n = x.n_elem;

	mat B(n, n, fill::eye);
	mat I(n, n, fill::eye);
	vec x_next = zeros<vec>(n), x_last = zeros<vec>(n), g_next = zeros<vec>(n), g_last = zeros<vec>(n), d;
	vec sk, yk;

	double epsilon = pow(10, -6);
	double un = 0.0001, mu = 0.9, amax;

	mat y(3, 1), s(3, 1);

	x_next = x;
	g_next = Diff(x_next);

	while (norm(g_next) > epsilon)
	{
		sk = x_next - x_last;
		yk = g_next - g_last;
		y = yk; s = sk;
		B = B + kron(y, y.t()) / dot(y, s) - B*kron(s, s.t())*B / dot(s, B*s);
		d = -inv(B)*g_next;
		//---no inverse version------------------------------
		//B = (I - kron(s, y.t()) / dot(y, s))*B*(I - kron(y, s.t()) / dot(y, s)) + kron(s, s.t()) / dot(y, s);
		//d = -B*g_next;
		//---------------------------------------------------
		x_last = x_next;
		g_last = g_next;
		amax = (0.01 / abs(dphi(x_last, d, 0))) / un;//0.01 is the precesion.
		x_next = Linesearch(x_last, d, un, mu, amax);
		g_next = Diff(x_next);
	}

	return x_next;

}

vec BFGSB(vec x, vec l, vec u)
{
	//m is the store vector number;
	//n is the para number 
	//x is the start point ;

	int n = x.n_elem;

	mat B,A,In;
	vec x_next = zeros<vec>(n), x_last = zeros<vec>(n), g_next = zeros<vec>(n), g_last = zeros<vec>(n), d;
	vec sk, yk;
	mat I(n, n, fill::eye);
	double epsilon = pow(10, -6);
	double un = 0.0001, mu = 0.9, amax,rho;

	mat y(3, 1), s(3, 1);

	x_next = x;
	g_next = Diff(x_next);
	A = Activeset(x_next, u, l);
	In = Inactive(x_next, u, l);
	B = In;

	while (norm(g_next) > epsilon)
	{
		A = Activeset(x_next, u, l);
		In = Inactive(x_next, u, l);
		sk = x_next - x_last;
		yk = g_next - g_last;
		yk = In*yk;
		rho = dot(sk, yk);
		sk = In*sk;
		y = yk; s = sk;
		if (rho > 0) { B = B + (I - kron(s, y.t()) / dot(y, s))*In*B*In*(I - kron(y, s.t()) / dot(y, s)) + kron(s, s.t()) / dot(y, s); }
		else { B = In; }
		d = -(A + B)*g_next;
		amax = (0.01 / abs(dphi(x_next, d, 0))) / un;//0.01 is the precesion.
		x_last = x_next;
		g_last = g_next;
		x_next = Linesearch(x_last, d, un, mu, amax);
		g_next = Diff(x_next);
	}

	return x_next;

}

vec BFGSB_eps(vec x, vec l, vec u)
{
	//m is the store vector number;
	//n is the para number 
	//x is the start point ;

	int n = x.n_elem;

	mat B, A, In;
	vec x_next = zeros<vec>(n), x_last = zeros<vec>(n), g_next = zeros<vec>(n), g_last = zeros<vec>(n), d;
	vec sk, yk ;
	mat I(n, n, fill::eye);
	double epsilon = pow(10, -6);
	double un = 0.0001, mu = 0.9, amax, rho, eps,pg;

	mat y(3, 1), s(3, 1);

	x_next = x;
	g_next = Diff(x_next);
	pg = norm(x - Proj(x_next - g_next,l,u));
	eps = min((u - l) / 2)< pg ? min((u - l) / 2): pg;
	A = Activeset(x_next, u, l, eps);
	In = Inactive(x_next, u, l, eps);
	B = In;

	while (norm(g_next) > epsilon)
	{	
		A = Activeset(x_next, u, l, eps);
		In = Inactive(x_next, u, l, eps);
		sk = x_next - x_last;
		yk = g_next - g_last;
		yk = In*yk;
		rho = dot(sk, yk);
		sk = In*sk;
		y = yk; s = sk;
		if (rho > 0) { B = B + (I - kron(s, y.t()) / dot(y, s))*In*B*In*(I - kron(y, s.t()) / dot(y, s)) + kron(s, s.t()) / dot(y, s); }
		else { B = In; }
		d = -(A + B)*g_next;
		amax = (0.01 / abs(dphi(x_next, d, 0))) / un;//0.01 is the precesion.
		x_last = x_next;
		g_last = g_next;
		x_next = Linesearch(x_last, d, un, mu, amax);
		g_next = Diff(x_next);
		pg = norm(x - Proj(x_next - g_next, l, u));
		eps = min((u - l) / 2)< pg ? min((u - l) / 2) : pg;
	}

	return x_next;

}

vec LBFGS(int m, vec x)
{
	//m is the store vector number;
	//n is the para number 
	//x is the start point ;
	
	int n = x.n_elem;
	mat Y(n, 1);
	mat S(n, 1);

	vec x_next, x_last, g_next, g_last, s_next, y_next, dir;
	vec R;

	double delta, epsilon, eps = pow(10, -16);
	double un = 0.0001, mu = 0.9;

	x_last = x;
	g_last = Diff(x_last);
	Y = x_last;
	S = g_last;
	delta = 1;
	double amax;
	while (delta>0.00001)
	{
		dir = twoloop(Y, S, g_last);
		amax = (0.01/abs(dphi(x_last,dir,0))) / un;//0.01 is the precesion.
		x_next = Linesearch(x_last, dir, un, mu, amax);
		g_next = Diff(x_next);
		s_next = x_next - x_last;
		y_next = g_next - g_last;
		epsilon = norm(y_next, "inf"); epsilon = eps*epsilon*epsilon;
		if (dot(s_next, y_next) > epsilon)
		{
			S = join_rows(S, s_next);
			Y = join_rows(Y, y_next);			
			if ((int)S.n_cols > m) { S.shed_col(0); Y.shed_col(0); }
		}
	    delta = norm(x_next - x_last, "inf");
		x_last = x_next;
		g_last = g_next;
	}

	return x_last;

}

vec LBFGSB(int m, vec x, vec l, vec u)
{
	//m is the store vector number;
	//n is the para number 
	//x is the start point ;
	//u is the simple uppon boundary;
	//l is the simple lower boundary;

	int n = x.n_elem;
	mat Y(n, 1);
	mat S(n, 1);

	vec x_next, x_last, g_next, g_last, s_next, y_next, dir;
	double theta, delta, epsilon, eps = pow(10, -16);
	double un=0.0001, mu=0.9, amax;

	x_last = x;
	g_last = Diff(x_next);
	Y = x_last;
	S = g_last;
	delta = norm(Proj(x_last, l, u) - x_last, "inf");
	theta = dot(g_last, g_last) / dot(x_last, g_last);
	while (delta>0.00001)
	{
		dir = Cauchydual(Y, S, x_last, l, u, g_last, theta)-x_last;
		amax = Amax(x_last, l, u, dir);
		x_next = Linesearch(x_last,dir,un,mu,amax);
		g_next = Diff(x_last);
		s_next = x_next - x_last;
		y_next = g_next - g_last;
		epsilon = norm(y_next, "inf"); epsilon = eps*epsilon*epsilon;
		if (dot(s_next, y_next) > epsilon)
		{
			S = join_rows(S, s_next);
			Y = join_rows(Y, y_next);
			if ((int)S.n_cols > m) { S.shed_col(0); Y.shed_col(0); }
		}
		theta = dot(y_next, y_next) / dot(s_next, y_next);
		x_last = x_next;
		g_last = g_next;
		delta = norm(Proj(x_last, l, u) - x_last, "inf");
	}

	return x_last;

}
















//=============================TEST PART=======================

double tphi(double x, double a)
{
	double beta = 2;
	return -a / (a*a + beta);
}

double tdphi(double x, double a)
{
	double beta = 2;
	return (a*a - beta) / ((a*a + beta)*(a*a + beta));
}

double tbphi(double x, double a, double u)
{
	return tphi(x, a) - tphi(x, 0) - u*a*tdphi(x, 0);
}

double tdbphi(double x, double a, double u)
{
	return tdphi(x, a) - u*tdphi(x, 0);
}

double tgetminp(double x, double al, double au)
{
	double minv, delta;
	int i, xu;
	vec vvv(9);
	delta = (au - al) / 10;
	while (delta>0.0001)
	{
		for (i = 0; i < 9; i++) { vvv(i) = tphi(x, al + (i + 0.5)*delta); }
		minv = min(vvv);
		xu = index_min(vvv);
		al = al + (xu + 0.5)*delta - 0.5*delta;
		au = al + (xu + 0.5)*delta + 0.5*delta;
		delta = (au - al) / 10;

	}
	return (al + au) / 2;
}

double testline(double x)
{
	double x_next, amax = 4;
	double u = 0.0001, mu = 0.1;
	double al, al_next, au, au_next;
	double at, at_next, atemp;
	double fl, ft, fu, gl, gt, gu;
	al = 0;
	au = amax;
	at = 0.01;
	while ((tphi(x, at)>(tphi(x, 0) + u*at*tdphi(x, 0))) || (abs(tdphi(x, at))>mu*abs(tdphi(x, 0))))
	{
		if ((tbphi(x, at, u) <= 0) && (tdphi(x, at) >= 0))
		{
			if (tphi(x, at) > tphi(x, al)) { al_next = al; au_next = at; }
			else if ((tdphi(x, at)*(al - at)) >= 0) { al_next = at; au_next = au; }
			else { al_next = at; au_next = al; }
			gl = tdphi(x, al_next); gt = tdphi(x, at); gu = tdphi(x, au_next);
		}
		else
		{
			if (tbphi(x, at, u) > tbphi(x, al, u)) { al_next = al; au_next = at; }
			else if ((tdbphi(x, at, u)*(al - at)) >= 0) { al_next = at; au_next = au; }
			else { al_next = at; au_next = al; }
			gl = tdbphi(x, al_next, u); gt = tdbphi(x, at, u); gu = tdbphi(x, au_next, u);
		}
		fl = tphi(x, al_next); ft = tphi(x, at); fu = tphi(x, au_next);

		double ac = cubeminpoint(al, at, fl, ft, gl, gt);
		double aq = intpolfgf(al, at, fl, ft, gl);
		double as = intpolfgf(al, at, fl, gl, gt);

		if (ft > fl)
		{
			if (abs(ac - al) < abs(aq - al)) { at_next = ac; }
			else { at_next = (aq + ac) / 2; }
		}
		else if (gt*gl < 0)
		{
			if (abs(ac - at) >= abs(as - at)) { at_next = ac; }
			else { at_next = as; }
		}
		else if (abs(gt)<abs(gl))
		{
			if (abs(ac - at) < abs(as - at)) { at_next = ac; }
			else { at_next = as; }

			atemp = at_next + 0.66*(au - at);
			if (at > al) { at_next = atemp < at_next ? atemp : at_next; }
			else { at_next = atemp > at_next ? atemp : at_next; }
		}
		else
		{
			at_next = cubeminpoint(au, at, fu, ft, gu, gt);
		}
		al = al_next;
		at = at_next;
		au = au_next;
	}
	at = tgetminp(x, al, au);
	x_next = x + at;
	return x_next;

}

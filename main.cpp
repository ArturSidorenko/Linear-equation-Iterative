#include"matrices.h"

using namespace std;

void plain_test();
void straightforward_tests();
void stopping_tests();
double max_norm_val(double a, double b, const method_tau &t);
matrix perebor(double a, double b, const method_tau &t);


//to estimate norms
double sf(double x, double tau);
double cumsf(double x, const method_tau &t);

int main() {
	//plain_test();
	try {
		straightforward_tests();
		stopping_tests();
	}
	catch (const std::out_of_range &r) {
		cerr << "RANGEERR: " << r.what() << "\n";
		system("pause");
		exit(-1);
	}
	catch (const std::invalid_argument &a) {
		cerr << "VALERR: " << a.what() << "\n";
		system("pause");
		exit(-1);
	}
	catch (...) {
		cerr << "SUPER UBER CRASH\n";
		exit(-1);
	}
	system("pause");
	return 0;
}

void plain_test() {
	ofstream s("text.txt");
	matrix A(2, 2);
	A.set(0, 0, 1); A.set(1, 1, 1);
	A.set(1, 0, 1);

	matrix B = 2 * A;

	matrix C = A * B;
	C = C;
	C = C + C;
	cout << C.get(0, 0) << "\n" << C.get(1, 1) << "\n";
	C.print(s);
	s.close();
	system("pause");
}

void straightforward_tests()
{
	cout << "\n\nstraightforward_test\n\n";
	size_t n = 200;
	matrix A = lagrange(n + 1);
	double h = 1. / (n + 1);
	matrix x(n, 1); //proposed solution

	for (size_t i = 0; i < n; i++) x.set(i, 0, 0.5*i + 1);
	matrix b = A * x; //proposed right hand of the equation

	permut clever = clever_met(4);
	size_t ord = clever.size(); //2^4
	method_tau met(ord);

	//naive order
	for (size_t i = 0; i < ord; i++) met[i] = inv_cheb_grid(ord, i, 4, 4. / h / h);

	matrix x0(n, 1); //zero vector

	cout << "Naive approach\n";
	cout << "Difference before: " << (x0 - x).linf_norm() << "\n";

	//norm research
	perebor(4, 4. / h / h, met).print("norms_naive.txt");

	matrix ans(n, 1);
	try {
		ans = test_n_steps_method(A, x0, b, met, x);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	double diff = (ans - x).linf_norm();
	cout << "Differecnce after = " << diff << "\n";

	//less naive order
	cout << "Non-naive processing\n";
	if (!(ord % 2)) for (size_t i = 0; i < (ord / 2); i++) {

		met[2 * i] = inv_cheb_grid(ord, ord - i - 1, 4, 4. / h / h);
		met[2 * i + 1] = inv_cheb_grid(ord, i, 4, 4. / h / h);
	}

	try {
		ans = test_n_steps_method(A, x0, b, met, x);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	diff = (ans - x).linf_norm();
	cout << "Differecnce in non-naive= " << diff << "\n";
	ofstream f("ans1.txt");
	ans.print(f);

	//norm research
	perebor(4, 4. / h / h, met).print("norms_cool.txt");

	//very clever processing

	cout << "Clever case processing\n";

	for (size_t i = 0; i < ord; i++) met[clever[i]] = inv_cheb_grid(ord, n - i - 1, 4, 4. / h / h);

	try {
		ans = test_n_steps_method(A, x0, b, met, x);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	diff = (ans - x).linf_norm();
	cout << "Differecnce in clever case= " << diff << "\n";
	f = ofstream("ans2.txt");
	ans.print(f);

	//norm research
	perebor(4, 4. / h / h, met).print("norms_optimal.txt");
}

void stopping_tests()
{
	cout << "\n\nstopping_time__test\n\n";
	size_t n = 100;
	matrix A = lagrange(n + 1);
	double h = 1. / (n + 1);
	matrix x(n, 1); //proposed solution

	for (size_t i = 0; i < n; i++) x.set(i, 0, 0.5*i + 1);
	matrix b = A * x; //proposed right hand of the equation

	size_t ord = 64;
	permut clever = clever_met(6);
	method_tau met(ord);

	//naive order
	for (size_t i = 0; i < ord; i++) met[i] = inv_cheb_grid(ord, i, 4, 4. / h / h);

	matrix x0(n, 1); //zero vector

	cout << "Naive approach\n";
	cout << "Difference before: " << (x0 - x).linf_norm() << "\n";

	matrix ans(n, 1);
	try {
		ans = n_steps_method(A, x0, b, met);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	double diff = (ans - x).linf_norm();
	cout << "Differecnce after = " << diff << "\n";

	//less naive order
	cout << "Non-naive processing\n";
	if (!(ord % 2)) for (size_t i = 0; i < (ord / 2); i++) {

		met[2 * i] = inv_cheb_grid(ord, ord - i - 1, 4, 4. / h / h);
		met[2 * i + 1] = inv_cheb_grid(ord, i, 4, 4. / h / h);
	}

	try {
		ans = n_steps_method(A, x0, b, met);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	diff = (ans - x).linf_norm();
	cout << "Differecnce in non-naive= " << diff << "\n";
	ofstream f("ans1.txt");
	ans.print(f);

	//very clever processing

	cout << "Clever case processing\n";

	for (size_t i = 0; i < ord; i++) met[clever[i]] = inv_cheb_grid(ord, i, 4, 4. / h / h);

	try {
		ans = n_steps_method(A, x0, b, met);
	}
	catch (const std::invalid_argument &ups) {
		cerr << "ERROR while solving: " << ups.what() << "\n";
		exit(-1);
	}
	diff = (ans - x).linf_norm();
	cout << "Differecnce in clever case= " << diff << "\n";
	f = ofstream("ans2.txt");
	ans.print(f);
}

double sf(double x, double tau) {
	return 1 - x * tau;
}

double cumsf(double x, const method_tau &t) {
	size_t n = t.size();
	double a = 1;
	for (size_t i = 0; i < n; i++) a *= sf(x, t[i]);
	return fabs(a);
}

double max_norm_val(double a, double b, const method_tau &t) {
	double ans = 0;
	double r;
	for (double c = a; c < b; c += 0.05) {
		r = cumsf(c, t);
		if (r > ans) ans = r;
	}
	return ans;
}

matrix perebor(double a, double b, const method_tau &t) {
	method_tau g;
	size_t n = t.size();
	matrix m(n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i; j < n; j++) {
			g.push_back(t[i]);
			m.set(i, j, max_norm_val(a, b, g));
		}
		g.clear();
	}
	return m;
}


#include"matrices.h"



//=======================constructors and assginments=================
matrix::matrix(size_t rows, size_t cols): 
	rows_(rows), cols_(cols) 
{
	data_ = new double*[cols];
	for (size_t i = 0; i < cols; i++) data_[i] = new double[rows];

	for (size_t j = 0; j < cols; j++) {
		for (size_t i = 0; i < rows; i++) data_[j][i] = 0;
	}
}

matrix::matrix(const matrix &other) :
	rows_(other.rows_), cols_(other.cols_)
{
	data_ = new double*[cols_];
	for (size_t i = 0; i < cols_; i++) data_[i] = new double[rows_];

	for (size_t j = 0; j < cols_; j++) {
		for (size_t i = 0; i < rows_; i++) data_[j][i] = other.data_[j][i];
	}
}

matrix::matrix(matrix &&other):
	rows_(other.rows_), cols_(other.cols_)
{
	data_ = other.data_;
	other.data_ = nullptr;
}

matrix::~matrix() {
	if (data_ != nullptr) {
		for (size_t i = 0; i < cols_; i++) delete[] data_[i];
		delete[] data_;
	}
}

matrix matrix::operator=(const matrix &other) {
	if (this != &other) {

		if (data_ != nullptr) {
			for (size_t i = 0; i < cols_; i++) delete[] data_[i];
			delete[] data_;
		}

		cols_ = other.cols_;
		rows_ = other.rows_;

		data_ = new double*[cols_];
		for (size_t i = 0; i < cols_; i++) data_[i] = new double[rows_];

		for (size_t j = 0; j < cols_; j++) {
			for (size_t i = 0; i < rows_; i++) data_[j][i] = other.data_[j][i];
		}

	}
	return *this;
}

matrix matrix::operator=(matrix &&other) {
	if (data_ != nullptr) {
		for (size_t i = 0; i < cols_; i++) delete[] data_[i];
		delete[] data_;
	}

	cols_ = other.cols_;
	rows_ = other.rows_;

	data_ = other.data_;

	other.data_ = nullptr;
	return *this;
}

//==============matrix operations================

matrix matrix::operator+(const matrix &M) const {
	if ((M.cols_ != cols_) || (M.rows_ != rows_))
		throw std::invalid_argument("Operator+: the sizes do not coincide");

	matrix ans(rows_, cols_);
	for (size_t j = 0; j < cols_; j++) {
		for (size_t i = 0; i < rows_; i++) ans.data_[j][i] = data_[j][i] + M.data_[j][i];
	}
	return ans;
}

matrix matrix::operator-(const matrix &M) const {
	if ((M.cols_ != cols_) || (M.rows_ != rows_))
		throw std::invalid_argument("Operator+: the sizes do not coincide");

	matrix ans(rows_, cols_);

	double *curcur, *curm, *cur_ans; //temps

	for (size_t j = 0; j < cols_; j++) {
		curcur = data_[j];
		curm = M.data_[j];
		cur_ans = ans.data_[j];
		for (size_t i = 0; i < rows_; i++) cur_ans[i] = curcur[i] - curm[i];
	}
	return ans;
}

matrix matrix::operator*(const matrix &M) const {
	if ((cols_) != (M.rows_)) throw std::invalid_argument("Operator*: wrong sizes");

	double *cur; //temp
	matrix ans(rows_, M.cols_);
	for (size_t j = 0; j < ans.cols_; j++) {
		cur = M.data_[j];
		for (size_t i = 0; i < ans.rows_; i++) {
			ans.data_[j][i] = 0;
			for (size_t k = 0; k < cols_; k++) {
				ans.data_[j][i] += data_[k][i] * cur[k];
			}
		}
	}
	return ans;
}

matrix matrix::operator*(double x) const {
	matrix ans = *this;

	for (size_t j = 0; j < cols_; j++) {
		for (size_t i = 0; i < rows_; i++) {
			ans.data_[j][i] *= x;
		}
	}
	return ans;
}

matrix operator*(double x, const matrix & A)
{
	matrix ans = A;
	for (size_t j = 0; j < ans.cols_; j++) {
		for (size_t i = 0; i < ans.rows_; i++) {
			ans.data_[j][i] *= x;
		}
	}
	return ans;
}

double matrix::get(size_t i, size_t j) const {
	if ((i >= rows_) || (j >= cols_)) throw std::out_of_range("(GET) The requested index is out of range");

	return data_[j][i];
}

void matrix::set(size_t i, size_t j, double x) {
	if ((i >= rows_) || (j >= cols_)) throw std::out_of_range("(SET) The requested index is out of range");

	data_[j][i] = x;
}

void matrix::print(std::ofstream &s) const {
	s << "\n============================\n";
	for (size_t i = 0; i < rows_; i++) {
		for(size_t j =0; j < cols_; j++)
			s << std::setw(10) << data_[j][i] << " "; 
		s << "\n";
	}
	s << "============================\n";
}

double matrix::l1_norm() const
{
	//col-wise sums
	double *maxes = new double[cols_];
	for (size_t k = 0; k < cols_; k++) maxes[k] = 0;

	for (size_t j = 0; j < cols_; j++) {
		for (size_t i = 0; i < rows_; i++) {
			maxes[j] += fabs(data_[j][i]);
		}
	}
	//finding the msximal value
	double c = maxes[0];
	for (size_t j = 1; j < cols_; j++)
		if (c < maxes[j]) c = maxes[j];

	return c;
}

double matrix::linf_norm() const
{
	//row-wise sums
	double *maxes = new double[rows_];
	for (size_t k = 0; k < rows_; k++) maxes[k] = 0;

	for (size_t j = 0; j < cols_; j++) {
		for (size_t i = 0; i < rows_; i++) {
			maxes[i] += fabs(data_[j][i]);
		}
	}
	//finding the msximal value
	double c = maxes[0];
	for (size_t i = 1; i < rows_; i++)
		if (c < maxes[i]) c = maxes[i];

	return c;
}

matrix id(size_t n)
{
	matrix m(n, n);
	for(size_t i = 0; i < n; i++) m.set(i, i, 1);
	return m;
}

matrix lagrange(size_t N)
{
	if (N <= 0) throw std::invalid_argument("The suze of the Lagrange matrix is (N-1)x(N-1)");
	matrix M(N - 1, N - 1);

	double h = 1. / N;
	
	for (size_t i = 0; i < N - 1; i++) {
		M.set(i, i, 2. / (h*h));
		if(i != 0)     M.set(i, i-1, -1. / (h*h));
		if(i != (N-2)) M.set(i, i+1, -1. / (h*h));
	}
	return M;
}

double cheb_grid(unsigned ord, unsigned k, double a, double b)
{
	if (ord == 0) throw std::invalid_argument("The zeroth Chebyshev polynomial does not have roots");
	if (k >= ord) throw std::invalid_argument("Cheb_grid: k >=ord");
	if (b < a)    throw std::invalid_argument("cheb_grid: right end < left end");
	k = ord - k - 1; //reverse

	double grid_van = cos(PI * (k + 0.5) / ord);

	return (a + b) * 0.5 + (b - a) * 0.5 * grid_van;
}

double inv_cheb_grid(unsigned ord, unsigned k, double a, double b) {
	return 1. / cheb_grid(ord, k, a, b);
}

matrix n_steps_method(const matrix & A, const matrix & X, const matrix &B, const method_tau & t)
{
	matrix ans = X;
	size_t ord = t.size();
	sol_logs logs;
	double res;

	std::cout << "    Size of the array of taus = " << ord << "\n";

	for (int g = 1; ; g++) {
		for (size_t i = 0; i < ord; i++) {
			ans = ans - t[i] * A * ans + t[i] * B; //(I - t*A)*x + t*b
		}

		if (!(g % 100)) {
			res = (A * ans - B).linf_norm();
			std::cout << "    Iteration #" << g << ",res = " << res << "\n";
			logs.push_back(res);
			if (break_loop(logs)) break;
			if (g > 2000) break;
		}
	}
	return ans;
}

matrix test_n_steps_method(const matrix & A, const matrix & X, const matrix &B, const method_tau & t, const matrix &true_ans)
{
	matrix ans = X;
	size_t ord = t.size();
	sol_logs logs;
	double diff;

	std::cout << "    Size of the array of taus = " << ord << "\n";

	for (int g = 1; ; g++) {
		for (size_t i = 0; i < ord; i++) {
			ans = ans - t[i] * A * ans + t[i] * B; //(I - t*A)*x + t*b
		}

		if (!(g % 100)) {
			diff = (ans - true_ans).linf_norm();
			std::cout << "    Iteration #" << g << ",L_inf diff = " << diff << "\n";
			logs.push_back(diff);
	
			if ((g > 2000) || (diff < 10e-12) || (diff > 10e+16)) break;
		}
	}
	return ans;
}

bool break_loop(const sol_logs &g) {
	
	size_t s = g.size();
	double cur = g[s-1];
	if (cur > 1e+15) return true;

	if (s < 6) return false; //to avoid multiple errors caused by to lack of data to alanyse
	double a = (g[s - 2] + g[s - 3] + g[s - 4]) / 3; //rolling average
	return (s > 0.95 * a);

}

permut clever_met(size_t t) {
	if (t == 0) {
		permut m(1);
		m[0] = 0;
		return m;
	}
	else {
		permut beg = clever_met(t - 1);
		permut adv(beg.size() * 2);
		for (int i = 0; i < beg.size(); i++) {
			adv[2 * i] = beg[i];
			adv[2 * i + 1] = adv.size() - beg[i] - 1;
		}
		return adv;
	}
	return permut();
}

#include"function_names.h"


vector<double> sweepMethod(int n,double a, double b, double mu1, double mu2, vector<double>f) {
	vector<double> result(n + 1); // вектор для записи результата
	result[0] = mu1; // левое граничное условие
	result[n] = mu2; //правое граничное условие
	vector<double> alpha(n), beta(n); //векторы для хранения коэффициентов альфа и бета
	alpha[0] = 0;
	beta[0] = 0;
	double h = (b-a) / n;

	//прямой ход прогонки
	for (int i = 1; i < n - 1; i++) {//у альфа и бета индексы сдвинуты на 1
		double Ai = h; // коэффициент Ai из метода прогонки
		double Bi = h; // коэффициент Вi из метода прогонки
		double Ci = -2*h; // коэффициент Сi из метода прогонки
		double fii = -6*((f[i+1]-f[i]) /h -(f[i]-f[i-1])/h); // коэффициент fi_i из разностной схемы
		alpha[i] = Bi / (Ci - Ai * alpha[i - 1]);
		beta[i] = (fii + Ai * beta[i - 1]) / (Ci - Ai * alpha[i - 1]);
	}

	//обратный ход прогонки
	for (int i = n - 1; i > 0; i--) {
		result[i] = alpha[i] * result[i + 1] + beta[i];
	}
	return result;
}
//pair<double, double>(*f)(double, double, double, double, double)
vector<vector<double>>get_coeffs(int n, double a_prm, double b_prm, double mu1, double mu2, double(*F)(double)) {
	int N = 2 * n;
	double h = (b_prm - a_prm) / N;
	vector<double> func(N+1);
	vector<double>a(N);
	vector<double>xs(N + 1);
	for (int i = 0; i < N; i++) {
		xs[i] = a_prm + i * h;
		func[i] = F(a_prm + i * h);
		a[i] = func[i];
	}
	func[N] = F(b_prm);
	xs[N] = b_prm;
	vector<double>c(N);
	c = sweepMethod(N, a_prm, b_prm, mu1, mu2, func);
	vector<double>b(N);
	for (int i = 1; i < N+1; i++) {
		b[i - 1] = (func[i] - func[i - 1]) / h + c[i] * h / 3. + c[i - 1] * h / 6.;
	}
	vector<double>d(N);
	for (int i = 1; i < N + 1; i++) {
		d[i - 1] = (c[i]-c[i-1])/h;
	}
	vector<vector<double>>res;
	res.push_back(xs);
	res.push_back(func);
	res.push_back(a);
	res.push_back(b);
	res.push_back(c);
	res.push_back(d);
	return res;
}
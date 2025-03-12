#include"function_names.h"


vector<double> sweepMethod(int n,double a, double b, double mu1, double mu2, vector<double>f) {
	vector<double> result(n + 1); // ������ ��� ������ ����������
	result[0] = mu1; // ����� ��������� �������
	result[n] = mu2; //������ ��������� �������
	vector<double> alpha(n), beta(n); //������� ��� �������� ������������� ����� � ����
	alpha[0] = 0;
	beta[0] = 0;
	double h = (b-a) / n;

	//������ ��� ��������
	for (int i = 1; i < n - 1; i++) {//� ����� � ���� ������� �������� �� 1
		double Ai = h; // ����������� Ai �� ������ ��������
		double Bi = h; // ����������� �i �� ������ ��������
		double Ci = -2*h; // ����������� �i �� ������ ��������
		double fii = -6*((f[i+1]-f[i]) /h -(f[i]-f[i-1])/h); // ����������� fi_i �� ���������� �����
		alpha[i] = Bi / (Ci - Ai * alpha[i - 1]);
		beta[i] = (fii + Ai * beta[i - 1]) / (Ci - Ai * alpha[i - 1]);
	}

	//�������� ��� ��������
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
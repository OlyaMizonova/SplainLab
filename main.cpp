#include"function_names.h"

int main() {
	int n = 1000;
	double a_prm = -1, b_prm = 1;
	vector<vector<double>> coeffs = get_coeffs(n, a_prm, b_prm, getFTest);
	//cout << "i    xi    fi      ai      bi       ci      di\n";
	//for (int i = 0; i < n; i++) {
	//	cout << i << "   ";
	//	for (int j = 0; j < 6; j++) {
	//		cout << coeffs[j][i]<<"   ";
	//	}
	//	cout << "\n";
	//}
	double maxDiff, maxDifffirstDerivate, maxDiffsecondDerivate;
	vector<vector<double>>errs = getSplainError(n * 4, &maxDiff, &maxDifffirstDerivate, &maxDiffsecondDerivate,
		getFTest, getFTestDerivate, getFTestSecondDerivate, coeffs[0], coeffs[2], coeffs[3], coeffs[4], coeffs[5]);
	//cout << "i    fi    si     fi-si     f'i      s'i     f'i-s'i      f''i      s''i     f''i-s''i\n";
	//for (int i = 0; i < n*4; i++) {
	//	cout << i << "   ";
	//	for (int j = 0; j < 9; j++) {
	//		cout << errs[j][i]<<"   ";
	//	}
	//	cout << "\n";
	//}
	cout << "max|F-S| = " << maxDiff << "    max|F'-S'| = " << maxDifffirstDerivate << "    max|F''-S''| = " << maxDiffsecondDerivate;
}
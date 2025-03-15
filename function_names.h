#pragma once
#include<vector>
#include <functional>
#include<iostream>
#include<cmath>
using namespace std;

vector<double> sweepMethod(int n, double a, double b, double mu1, double mu2, vector<double>f);
/*
����� �������� �� ������� ���� � ���������� ��������� �������������
*/


vector<vector<double>>get_coeffs(int n, double a_prm, double b_prm, double(*F)(double));
/*
������� ������������� ������� �� ����� ����������� n, 
double a_prm, double b_prm - ����� � ������ ������� �������
double(*F)(double) - ��������� �� �������, ������� �� ������������� ��������
���������� � ����� ������� - {xs, fi, a, b, c, d} - ������� ����, �������� �-�� � �����, ������������ �������
*/

vector<vector<double>> getSplainError(int N, double* maxDiff, double* maxDifffirstDerivate, double* maxDiffsecondDerivate,
	double(*F)(double), double(*F_firstDerivate)(double), double(*F_secondDerivate)(double),
	vector<double>xs, vector<double> a, vector<double>b, vector<double>c, vector<double>d);
/*
���������� �������� 2-3 ��� ������� �), �)
int N - ����� �����������, ������� n,
double* maxDiff, double* maxDifffirstDerivate, double* maxDiffsecondDerivate - ����������, 
					� ������� ������ ������������ �������� ��������� �-�� � ������� � �� �����������
	double(*F)(double), double(*F_firstDerivate)(double), double(*F_secondDerivate)(double) - �-�� � �����������, 
																	��� ������� �������� ������ 
	vector<double>xs - ���� ���������,
	vector<double> a, vector<double>b, vector<double>c, vector<double>d - ������������ �������
*/

vector<vector<double>> getSplainErrorForLastCase(int N, double* maxDiff, double* maxDifffirstDerivate, double* maxDiffsecondDerivate,
	double(*F)(double, double(*)(double)), double(*F_firstDerivate)(double, double(*)(double)), double(*F_secondDerivate)(double, double(*)(double)), double(*f)(double),
	vector<double>xs, vector<double> a, vector<double>b, vector<double>c, vector<double>d);
/*
����� �� �� �����, ��� � ��� ���������� �������, �� ����������� �������� ������� f(x)
*/

double getFTest(double x);
/*
�������� ������� - ����� �)
*/

double getFTestDerivate(double x);
/*
����������� �������� ������� - ����� �)
*/

double getFTestSecondDerivate(double x);
/*
������ ����������� �������� ������� - ����� �)
*/

double getFBfunc4(double x);
/*
������� � ����� � �� ������ �) ��� ������� 4 (�� ���� � �������)
*/

double getFBFirstDerivatefunc4(double x);
/*
����������� ������� � ����� � �� ������ �) ��� ������� 4 (�� ���� � �������)
*/

double getFBSecondDerivatefunc4(double x);
/*
������ ����������� ������� � ����� � �� ������ �) ��� ������� 4 (�� ���� � �������)
*/

double getFBfunc5(double x);
/*
������� � ����� � �� ������ �) ��� ������� 5 (�� ���� � �������)
*/

double getFBFirstDerivatefunc5(double x);
/*
����������� ������� � ����� � �� ������ �) ��� ������� 5 (�� ���� � �������)
*/

double getFBSecondDerivatefunc5(double x);
/*
������ ����������� ������� � ����� � �� ������ �) ��� ������� 5 (�� ���� � �������)
*/

double getFBfunc6(double x);
/*
������� � ����� � �� ������ �) ��� ������� 6 (�� ���� � �������)

*/

double getFBFirstDerivatefunc6(double x);
/*
����������� ������� � ����� � �� ������ �) ��� ������� 6 (�� ���� � �������)
*/

double getFBSecondDerivatefunc6(double x);
/*
������ ����������� ������� � ����� � �� ������ �) ��� ������� 6 (�� ���� � �������)
*/

double getFCfunc1(double x, double(*f)(double));
/*
������� � ����� � �� ������ c) f+cos(10x) (�� ���� � �������)
*/

double getFCFirstDerivatefunc1(double x, double(*f)(double));
/*
����������� ������� � ����� � �� ������ c) f+cos(10x) (�� ���� � �������)
*/

double getFCSecondDerivatefunc1(double x, double(*f)(double));
/*
������ ����������� ������� � ����� � �� ������ c) f+cos(10x) (�� ���� � �������)
*/

double getFCfunc2(double x, double(*f)(double));
/*
������� � ����� � �� ������ c) f+cos(100x) (�� ���� � �������)
*/

double getFCFirstDerivatefunc2(double x, double(*f)(double));
/*
����������� ������� � ����� � �� ������ c) f+cos(100x) (�� ���� � �������)
*/

double getFCSecondDerivatefunc2(double x, double(*f)(double));
/*
������ ����������� ������� � ����� � �� ������ c) f+cos(100x) (�� ���� � �������)
*/

double getSValue(double x, vector<double>xs, vector<double> a, vector<double>b, vector<double>c, vector<double>d);
/*
���������� �������� ������� � ����� �
*/

double getSDerivateValue(double x, vector<double>xs, vector<double>b, vector<double>c, vector<double>d);
/*
���������� �������� ����������� ������� � ����� �
*/

double getSSecondDerivateValue(double x, vector<double>xs, vector<double>c, vector<double>d);
/*
���������� �������� ������ ����������� ������� � ����� �
*/

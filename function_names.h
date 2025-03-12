#pragma once
#include<vector>
#include <functional>
#include<iostream>
#include<cmath>
using namespace std;

vector<double> sweepMethod(int n, double a, double b, double mu1, double mu2, vector<double>f);

vector<vector<double>>get_coeffs(int n, double a_prm, double b_prm, double mu1, double mu2, double(*F)(double));
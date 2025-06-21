#include <iostream>
#include <cmath>
#include <vector>
#include <list>
// author: Orlovskiy Gleb
// datе: 6.02.24
// name: Решение задачи Дирихле для двумерного уравнения Пауссона, а так же вспомогательные задачи.
//       Подробное условие в README.md
using namespace std;

const double pi = 3.14159265358979323846;


double findMaxAbsolute(const vector<vector<double>>& arr) {
    double maxAbs = fabs(arr[0][0]);
    for (int i = 0; i < arr.size(); ++i) {
        for (int j = 0; j < arr.size(); ++j) {
            double absValue = fabs(arr[i][j]);
            if (absValue > maxAbs) {
                maxAbs = absValue;
            }
        }
    }
    return maxAbs;
}


double operatorA(vector<vector<double>>& u, double h, double a, double b, int i, int j)
{
    return ((-a) / (h * h)) * (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j])
        + ((-b) / (h * h)) * (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]);
}



vector<double> powerMethod(const double a, const double b, const double eps, int n, void(*init)(vector<vector<double>>&), list<vector<int>>& query)
{
    vector<vector<double>> x0(n + 1, vector<double>(n + 1));
	init(x0);
    vector<vector<double>> x1(n + 1, vector<double>(n + 1));
	x1 = x0;
    double h = 1.0 / n;

    //степенной метод:
    double alpha;
    do
    {
        double d = findMaxAbsolute(x0);
        for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it) 
        {
			int i = (*it)[0];
			int j = (*it)[1];
            x1[i][j] = operatorA(x0, h, a, b, i, j) / d;
        }
        alpha = fabs((findMaxAbsolute(x1) / d) - 1.0);
        x0 = x1;
    } while (alpha >= eps);

    const double maxEValue = findMaxAbsolute(x1);

	init(x0);
    x1 = x0;
    do
    {
        double d = findMaxAbsolute(x0);
        for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
        {
            int i = (*it)[0];
            int j = (*it)[1];
            x1[i][j] = ((maxEValue * x0[i][j]) - operatorA(x0, h, a, b, i, j)) / d;
        }
        alpha = fabs((findMaxAbsolute(x1) / d) - 1.0);
        x0 = x1;
    } while (alpha >= eps);

    const double minEValue = maxEValue - findMaxAbsolute(x1);
	vector<double> result = { minEValue, maxEValue };
    return result;
}

void startTestEV(vector<vector<double>>& x0)
{
    int n = x0.size();
    //заполним x0 нулями на границе
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == 0 || j == 0 || j == n - 1 || i == n - 1)
            {
                x0[i][j] = 0.0;
            }
            else
            {
                x0[i][j] = 1.0 / (i + j + cos(i + j));
            }
        }
    }
}


void testEV(int n)
{
	double h = 1.0 / n;
    vector<vector<double>> x0(n + 1, vector<double>(n + 1));
    list<vector<int>> query;
    for (int i = 1; i < n; i++)
    {
		for (int j = 1; j < n; j++)
		{
			vector<int> v = { i, j };
			query.push_back(v);
		}
    }
	vector<double> result = powerMethod(1.0, 1.0, 0.000000001, n, startTestEV, query);
	double minEValue = result[0];
	double maxEValue = result[1];
    cout << "for n = " << n << " max eigenvalue = " << maxEValue << endl;
    cout << "for n = " << n << " min eigenvalue = " << minEValue << endl;
    cout << "for n = " << n << " |minEV - left| = " << fabs(minEValue - (8.0 / (h * h)) * pow(sin((pi * h) / 2.0), 2)) << endl;
    cout << "for n = " << n << " |maxEV - right| = " << fabs(maxEValue - (8.0 / (h * h)) * pow(cos((pi * h) / 2.0), 2)) << endl;
}

double u(double x, double y)
{
	return 1.0 / (1 + x * x + y * y);
}

void startEV(vector<vector<double>>& x0)
{
    int n = x0.size();
    double h = 1.0 / (n - 1);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x0[i][j] = 0;
        }
    }
    for (int i = (n - 1)/2; i < n; i++)
    {
		x0[i][i - (n - 1)/2] = u(i*h, (i - (n - 1)/2)*h);
		x0[i - (n - 1)/2][i] = u((i - (n - 1)/2)*h, i*h);
    }
    for (int i = 0; i < n / 2; i++)
    {
		x0[i][0] = u(i*h, 0);
        x0[0][i] = u(0, i*h);
    }
    for (int i = (n + 1) / 2; i < n; i++)
    {
        x0[i][n - 1] = u(i * h, (n - 1) * h);
		x0[n - 1][i] = u((n - 1) * h, i * h);
    }
}

bool checkInList(list<vector<int>>& l, const vector<int>& v)
{
    for (auto it = l.begin(); it != l.end(); it++)
    {
        if ((*it)[0] == v[0] && (*it)[1] == v[1])
        {
            return true;
        }
    }
    return false;
}

void show(const vector<vector<double>>& vec, list<vector<int>>& query)
{  
	int n = vec.size();
    for (int j = n - 1; j >= 0; j--) 
    {
            for (int i = 0; i < n; i++) 
            {
                if (checkInList(query, { i, j }))
                {
                    cout << "q";
                    if (vec[i][j] != 0.0)
                    {
                        cout << "ERROR(" << i << ", " << j << ")" << endl;
                    }
                }
                else
                {
                    if (vec[i][j] != 0.0)
                    {
                        cout << "1";
                    }
                    else
                    {
                        cout << "0";
                    }
                }
                cout << " ";
            }
            cout << endl;
    }
    for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it) {
		cout << "q(" << (*it)[0] << ", " << (*it)[1] <<")" << endl;
    }

}

//создание очереди обхода узлов
list<vector<int>> queryEV(const int n)
{
    vector<vector<double>> x0(n + 1, vector<double>(n + 1));
    list<vector<int>> query;
    startEV(x0);
    for (int j = 0; j < n + 1; j++)
    {
        int i = 0;
        //пробиваем границу
        while (x0[i][j] == 0)
        {
            i++;
            if (x0[i][j] != 0.0)
            {
                i++;
                break;
            }
        }
        if (x0[i][j] != 0) i++;
        //пробили, пошли дальше
        for (; i < n + 1 && j < n + 1; i++)
        {
            if (x0[i][j] == 0.0)
            {
                vector<int> v = { i, j };
                query.push_back(v);
            }
            else
            {
                if (i + 1 < n + 1)
                {
                    if (x0[i + 1][j] == 0)
                    {
                        break;
                    }
                }
            }
        }
    }
    return query;
}

//вычисление минимального и максимального собственных значений
vector<double> EV(const int n, const double a, const double b, list<vector<int>>& query)
{
    double h = 1.0 / n;
    //show(x0, query);//debug
    vector<double> result = powerMethod(a, b, 0.000000001, n, startEV, query);
    return result;
}

double f(const double a, const double b, const double x, const double y)
{
    return a * ((2.0 * pow(1.0 + x * x + y * y, 2)) - 8.0 * x * x * (1.0 + x * x + y * y)) / (pow(1.0 + x * x + y * y, 4)) +
        b * ((2.0 * pow(1.0 + x * x + y * y, 2)) - 8.0 * y * y * (1.0 + x * x + y * y)) / (pow(1.0 + x * x + y * y, 4));
}

double findMaxAbsoluteInQuery(list<vector<int>>& query, vector<vector<double>> matrix)
{
	double max = 0;
	for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
	{
		if (fabs(matrix[(*it)[0]][(*it)[1]]) > max)
		{
			max = fabs(matrix[(*it)[0]][(*it)[1]]);
		}
	}
	return max;
}

vector<vector<double>> upperRelaxationMethod(const int n, const double a, const double b, const double w, const double alpha, const double eps, void(*init)(vector<vector<double>>&), list<vector<int>> query)
{
    const double h = 1.0 / n;
    vector<vector<double>> x0(n + 1, vector<double>(n + 1));
    init(x0);
    vector <vector<double>> x1(n + 1, vector<double>(n + 1));
    vector <vector<double>> errorV(n + 1, vector<double>(n + 1));
	errorV = x0;
    double errorD = 0;
    show(x0, query);
    //забьем мусором внутренние точки
	for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
	{
		int i = (*it)[0];
		int j = (*it)[1];
        x0[i][j] = i + j + cos(i + j);;
	}
    x1 = x0;
    do
    {
		//метод релаксации:
        for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
        {
            int i = (*it)[0];
            int j = (*it)[1];
            x1[i][j] = ((1 - w) * 2.0 * (a + b) * x0[i][j] + w * (a * x1[i - 1][j] + b * x1[i][j - 1] + a * x0[i + 1][j] + b * x0[i][j + 1] + h * h * f(a, b, h * i, h * j))) / (2.0 * (a + b));
            //x0[i][j] = x1[i][j];
        }
        //вычисление ошибки:
		for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
		{
			int i = (*it)[0];
			int j = (*it)[1];
            errorV[i][j] = operatorA(x1, h, a, b, i, j) - f(a, b, i * h, j * h);
		}
		x0 = x1;
		errorD = findMaxAbsoluteInQuery(query, errorV);
        //cout << errorD << endl;
	} while (errorD >= alpha*eps);
    return x1;
}

int main()
{
    //testEV(8);
    const int n = 8;
    const double a = 0.9;
    const double b = 1.1;

    list<vector<int>> query = queryEV(n);
    const vector<double> ev = EV(n, a, b, query);
	cout << "for n = " << n << " min_ev(A_h) = " << ev[0] << " max_ev(A_h) = " << ev[1] << endl;
	const vector<vector<double>> solution = upperRelaxationMethod(n, a, b, 1.6, 0.00001, ev[0], startEV, query);
    
    double maxError = -1;
    const double h = 1.0 / n;
    for (list<vector<int>>::iterator it = query.begin(); it != query.end(); ++it)
    {
		int i = (*it)[0];
		int j = (*it)[1];
		maxError = max(maxError, fabs(solution[i][j] - u(i*h, j*h)));
	}
   for (int i = 0; i < n + 1; i++)
    {
        cout << solution[i][i] << " " << u(i * h, i * h) << " " << fabs(solution[i][i] - u(i * h, i * h)) <<  endl;
    }
	cout << "max error = " << maxError << endl;

    return 0;
}

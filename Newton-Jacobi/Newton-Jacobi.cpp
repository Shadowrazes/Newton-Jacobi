#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include<iomanip>
#include <cmath>

using namespace std;

const double E = 0.0001;

struct Var {
    double var = 1.0;
    double degree = 0.0;
};

struct eqMember {
    Var x;
    Var y;
};

void printVector(vector<eqMember>& equation) {
    cout.precision(1);
    for (const auto& it : equation) {
        cout << fixed << it.x.var << "*X^" << it.x.degree << "*";
        cout << fixed << it.y.var << "*Y^" << it.y.degree << "*";
        cout << " + ";
    }
        
    cout << endl << endl;
}

void printMatrix(vector<vector<eqMember>>& matrix) {
    cout.precision(1);
    for (auto& stroke : matrix) {
        for (auto& column : stroke) {
            cout << fixed << column.x.var << "*X^" << column.x.degree << "*";
            cout << fixed << column.y.var << "*Y^" << column.y.degree << "*";
            if(column.x.var != stroke[stroke.size() - 1].x.var)
                cout << " + ";
        }
        cout << endl << endl;
    }
}

void strokeWrap(string& stroke, vector<eqMember>& equation) {
    string buff = "";
    for (auto& it : stroke) {
        if (it == ' ') {
            if (buff.size() != 0) {
                if (buff.find("x") != -1 && buff.find("y") != -1) {
                    int xDegreeI = buff.find("^");
                    int yI = buff.find("y");
                    equation.push_back(eqMember{ {stod(buff.substr(0, xDegreeI - 1)), stod(buff.substr(xDegreeI + 1, yI - (xDegreeI + 1)))},
                        {1.0, stod(buff.substr(yI + 2))}});
                }
                else if (buff.find("x") != -1) {
                    int xDegreeI = buff.find("^");
                    equation.push_back(eqMember{ {stod(buff.substr(0, xDegreeI - 1)), stod(buff.substr(xDegreeI + 1))}, {} });
                }
                else if (buff.find("y") != -1) {
                    int yDegreeI = buff.find("^");
                    equation.push_back(eqMember{ {}, {stod(buff.substr(0, yDegreeI - 1)), stod(buff.substr(yDegreeI + 1))}});
                }
                else
                    equation.push_back(eqMember{ {stod(buff), 0}, {}});
                buff = "";
            }
            continue;
        }
        else
            buff += it;
    }
}

double func(vector<eqMember>& equation, vector<double>& X) {
    double sum = 0;
    for (auto &it : equation) {
        sum += it.x.var * pow(X[0], it.x.degree) * pow(X[1], it.y.degree);
    }
    return sum;
}

double determinant(vector<vector<double>> matrix) {
    return matrix[0][0] * matrix[1][1] - (matrix[0][1] * matrix[1][0]);
}

vector<double> matrixMultiply(vector<vector<double>>& matrix1, vector<double>& matrix2) {
    vector<double> multiply(matrix2.size(), 0);

    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix1.size(); j++) {
            multiply[i] += matrix1[i][j] * matrix2[j];
        }
    }
    return multiply;
}

vector<double> matrixSub(vector<double>& matrix1, vector<double> matrix2) {
    vector<double> sub(matrix2.size(), 0);

    for (int i = 0; i < matrix1.size(); i++) {
        sub[i] = matrix1[i] - matrix2[i];
    }
    return sub;
}

vector<vector<eqMember>> derivative(vector<vector<eqMember>> matrix, bool isX) {
    vector<vector<eqMember>> W = matrix;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (isX) {
                if (W[i][j].x.degree == 0) {
                    W[i][j].x = { 0, 0 };
                    W[i][j].y = { 0, 0 };
                    continue;
                }
                W[i][j].x = { matrix[i][j].x.var * matrix[i][j].x.degree , matrix[i][j].x.degree - 1 };
            }
            else
            {
                if (W[i][j].y.degree == 0) {
                    W[i][j].y = { 0, 0 };
                    W[i][j].x = { 0, 0 };
                    continue;
                }
                W[i][j].y = { W[i][j].y.var, matrix[i][j].y.degree - 1 };
                W[i][j].x.var *= matrix[i][j].y.degree;
            }
        }
    }
    return W;
}

vector<vector<double>> reverseMatrix(vector<vector<vector<eqMember>>> W, vector<double>& X) {
    vector<vector<double>> matrix(2, vector<double>(2));
    cout << "W(X): " << endl;
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            matrix[i][j] = func(W[i][j], X);
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    swap(matrix[0][0], matrix[1][1]);
    matrix[0][1] *= -1; matrix[1][0] *= -1;
    double detM1 = 1.0 / determinant(matrix);

    cout << "Wrev(X): " << endl;
    for (auto& stroke : matrix) {
        for (auto& column : stroke) {
            column *= detM1;
            cout << column << " ";
        }
        cout << endl;
    }
    cout << endl;
    return matrix;
}

bool checkE(vector<double>& x, vector<double>& xPrev) {
    for (int i = 0; i < x.size(); i++) {
        if (abs(x[i] - xPrev[i]) > E) {
            return true;
        }
    }
    return false;
}

int main()
{
    fstream matrixData("equation.txt");
    vector<vector<eqMember>> matrix;

    while (!matrixData.eof()) {
        vector<eqMember> tmp;
        matrix.push_back(tmp);
        string buff;
        getline(matrixData, buff);
        strokeWrap(buff, matrix[matrix.size() - 1]);
    }

    cout << "F: " << endl;
    printMatrix(matrix);
    cout << "---------------------------------------" << endl;
    cout.precision(10);
    vector<vector <double>> X;
    X.push_back({ 2, 1 });
    vector<vector<vector<eqMember>>> W(matrix.size());
    vector<vector<eqMember>> tmp = derivative(matrix, true);
    
    W[0].push_back(tmp[0]);
    W[1].push_back(tmp[1]);

    tmp = derivative(matrix, false);

    W[0].push_back(tmp[0]);
    W[1].push_back(tmp[1]);

    int i = 0;

    cout << "X" << i << ": ";
    for (auto it : X[i])
        cout << it << " ";
    cout << endl;

    vector<double> F;
    for (int j = 0; j < matrix.size(); j++)
        F.push_back(func(matrix[j], X[i]));

    auto Wrev = reverseMatrix(W, X[i]);
    X.push_back(matrixSub(X[i], matrixMultiply(Wrev, F)));
    cout << "---------------------------------------" << endl;
    i++;

    while (checkE(X[i], X[i - 1])) {
        cout << "X" << i << ": ";
        for (auto it : X[i])
            cout << it << " ";
        cout << endl;

        vector<double> F;
        for (int j = 0; j < matrix.size(); j++)
            F.push_back(func(matrix[j], X[i]));

        auto Wrev = reverseMatrix(W, X[i]);
        X.push_back(matrixSub(X[i], matrixMultiply(Wrev, F)));
        cout << "---------------------------------------" << endl;
        i++;
    }
    cout << "X" << i << ": ";
    for (auto it : X[i])
        cout << it << " ";
    cout << endl;
}
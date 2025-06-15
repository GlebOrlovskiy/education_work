
// author: Орловский Глеб Геннадьевич
// problem: Пересечение отрезков
// v6 лучшее решение
// data: 17.05.25

#include<iostream>
#include <cmath> 
#include <vector>
#include <iomanip>

using namespace std;


class Vector3D
{
public:
    double X;
    double Y;
    double Z;
    
    Vector3D()
    {
        this->X = 0.0;
        this->Y = 0.0;
        this->Z = 0.0;
    }
    
    Vector3D(double X, double Y, double Z)
    {
        this->X = X;
        this->Y = Y;
        this->Z = Z;
    }
    
    double norm() const
    {
        return sqrt(X*X + Y*Y + Z*Z);
    }
    
     // Скалярное произведение
    double dot(const Vector3D& other) const
    {
        return X * other.X + Y * other.Y + Z * other.Z;
    }
    
    Vector3D operator*(double scalar) const
    {
        return Vector3D(X * scalar, Y * scalar, Z * scalar);
    }
    
    Vector3D operator+(const Vector3D& other) const
    {
        return Vector3D(X + other.X, Y + other.Y, Z + other.Z);
    }

    
    Vector3D operator-(const Vector3D& other) const
    {
        return Vector3D(X - other.X, Y - other.Y, Z - other.Z);
    }
    
    // векторное произведение
    Vector3D cross(const Vector3D& other) const 
    {
        return Vector3D(
            Y * other.Z - Z * other.Y,
            Z * other.X - X * other.Z,
            X * other.Y - Y * other.X
        );
    }
};

ostream& operator<<(ostream& os, const Vector3D& v) 
{
        os << "(" << v.X << ", " << v.Y << ", " << v.Z << ")";
        return os;
}

class Matrix
{
private:
    vector<vector<double>> M;
public:
    Matrix(vector<vector<double>> matrix)
    {
        this->M = matrix;
    }
    Matrix(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3)
    {
        this->M = vector<vector<double>>({{v1.X, v1.Y, v1.Z},
                                         {v2.X, v2.Y, v2.Z},
                                         {v3.X, v3.Y, v3.Z}});
    }
    
    double det() const
    {
        return M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) +
              -M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) +
               M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
    }
};


class Plane
{
private:
    Vector3D start;
    Vector3D v1;
    Vector3D v2;
public:
    Plane(Vector3D start, Vector3D v1, Vector3D v2)
    {
        this->start = start;
        this->v1 = v1;
        this->v2 = v2;
    }
    double f(const Vector3D& point) const
    {
        return Matrix(point - start, v1, v2).det();
    }
};

class Segment3D
{
public:
    Vector3D start;
    Vector3D end;
    
    Segment3D(Vector3D start, Vector3D end)
    {
        this->start = start;
        this->end = end;
    }
    
    Vector3D getVector3D() const
    {
        return (this->end - this->start);
    }
    //Скрещиваются ли прямые сегментов
    bool isSkew(const Segment3D& seg, double eps) const
    {
        Vector3D v1 = seg.start - this->start;
        if(v1.norm() < eps) return false;
        Vector3D v2 = this->getVector3D();
        Plane plane(this->start, v1, v2);
        if(abs(plane.f(seg.end)) > eps) return true;
        return false;
    }
};


int sigm(double x)
{
    if(x < 0) return -1;
    if(x > 0) return 1;
    return 0;
}

bool isIntersect(const Segment3D& seg1, const Segment3D& seg2, const Vector3D& cross)
{
    Vector3D v1 = seg1.getVector3D();
    Vector3D v2 = seg2.getVector3D();

    Plane plane1(seg1.start, v1, cross);
    Plane plane2(seg2.start, v2, cross);
    
    //беру только знак для того чтобы избежать переполнения
    if(sigm(plane1.f(seg2.start))*sigm(plane1.f(seg2.end)) <= 0 && sigm(plane2.f(seg1.start))*sigm(plane2.f(seg1.end)) <= 0)
    {
        return true;
    }

    return false;
}

Vector3D Intersect(const Segment3D& seg1, const Segment3D& seg2, double eps_zero)
{
    /*
    IN: На входе сегменты не нулевой длины, eps_zero порог нуля
    OUT: на выходе точка пересечения отрезков,
         если пересечений больше одной точки или пересечений нет - на выходе NAN вектор 
    */
    
    Vector3D vectorNAN(NAN, NAN, NAN);
    
    if(seg1.isSkew(seg2, eps_zero)) return vectorNAN;
    
    Vector3D e1 = seg1.getVector3D();
    e1 = e1*(1.0/e1.norm());
    Vector3D e2 = seg2.getVector3D();
    e2 = e2*(1.0/e2.norm());
    Vector3D vcross = e1.cross(e2);
    
    if(vcross.norm() < eps_zero)
    {
        if(e1.dot(e2) > 0)
        {
            if((seg1.start - seg2.end).norm() < eps_zero)
            {
                return seg1.start;
            }
            else if((seg2.start - seg1.end).norm() < eps_zero)
            {
                return seg2.start;
            }
        }
        else
        {
            if((seg1.start - seg2.start).norm() < eps_zero)
            {
                return seg1.start;
            }
            else if((seg1.end - seg2.end).norm() < eps_zero)
            {
                return seg1.end;
            }
        }
        
        return vectorNAN;
    }
    
    if (!isIntersect(seg1, seg2, vcross))
    {
        return vectorNAN;
    }
    else
    {
        
        Vector3D b = seg2.start - seg1.start;
       
        /*
        l1 := seg1.start + t1*e1; l2 := seg2.start + t2*e2; где ti из [0, 1]
        l1 = l2 <=> t1*e1 - t2*e2 = seg2.start - seg1.start =: b
        
        уже знаем, что система ниже имеет одно решение
        |e1x -e2x| |t1| = |bx| <=> A*t = b
        |e1y -e2y| |t2|   |by|
        |e1z -e2z|        |bz|
        
        rank A >= 2
        Ранг матрицы — наивысший из порядков всевозможных ненулевых миноров этой матрицы.
        
        Значит некоторые из систем ниже имеют единственное решение
        
        |e1x -e2x||t1| = |bx|
        |e1y -e2y||t2|   |by|
        
        |e1x -e2x||t1| = |bx|
        |e1z -e2z||t2|   |bz|
        
        |e1y -e2y||t1| = |by|
        |e1z -e2z||t2|   |bz|

        */
        
        double det[3] = {(e1.X*(-e2.Y) - (-e2.X)*e1.Y),
                         (e1.X*(-e2.Z) - (-e2.X)*e1.Z),
                         (e1.Y*(-e2.Z) - (-e2.Y)*e1.Z)};
        int max = 0;
        
        for(int i = 1; i < 3; i++)
        {
            if(abs(det[max]) < abs(det[i]))
            {
                max = i;
            }
        }
        
        double t1 = 0.0;
        switch (max) 
        {
            case 0:   
                t1 = (b.X*(-e2.Y) - (-e2.X)*b.Y)/det[max];
                break;
            case 1:
                t1 = (b.X*(-e2.Z) - (-e2.X)*b.Z)/det[max];
                break;
            case 2:
                t1 = (b.Y*(-e2.Z) - (-e2.Y)*b.Z)/det[max];
                break;
        }
        
        return (seg1.start + e1*t1);
    }
}


int main()
{
    const double eps = 1e-6; // Настраиваемое
    const double eps_zero = 1e-6; // Настраиваемое
    
    cout << fixed << setprecision(3);
    
    cout << "Тест1 Точка пересечения, правильный ответ (0.5, 0.5, 0.0)" << endl;
    Segment3D seg1(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 0.0));
    Segment3D seg2(Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 1.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест2 Разнонаправленные на одной прямой, правильный ответ (0.0, 0.0, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    seg2 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(-1.0, 0.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест3 Разнонаправленные на одной прямой, правильный ответ (1.0, 0.0, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    seg2 = Segment3D(Vector3D(2.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест4 Точка персечения, правильный ответ (0.0, 0.0, 5.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0));
    seg2 = Segment3D(Vector3D(1.0, 0.0, 0.0), Vector3D(-1.0, 0.0, 1.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест5 Однонаправленные, правильный ответ (1.0, 1.0, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 0.0));
    seg2 = Segment3D(Vector3D(1.0, 1.0, 0.0), Vector3D(2.0, 2.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест6 Однонаправленные, правильный ответ (0.0, 0.0, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    seg2 = Segment3D(Vector3D(-1.0, 0.0, 0.0), Vector3D(0.0, 0.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест7 Точка пересечения, правильный ответ (0.5, 0.5, 0.5)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 1.0));
    seg2 = Segment3D(Vector3D(1.0, 1.0, 0.0), Vector3D(0.0, 0.0, 1.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест8 Вершина угла, правильный ответ (0.0, 0.0, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    seg2 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 1.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест9 Буква T, правильный ответ (0.5, 0.5, 0.0)" << endl;
    seg1 = Segment3D(Vector3D(0.5, 0.5, 0.0), Vector3D(5.0, 5.0, 1.0));
    seg2 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест10 Нет пересечения" << endl;
    seg1 = Segment3D(Vector3D(0.5, 0.5, 0.0), Vector3D(5.0, 5.0, 1.0));
    seg2 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(-1.0, -1.0, -1.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест11 Cкрещивающиеся" << endl;
    seg1 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 0.0));
    seg2 = Segment3D(Vector3D(0.1, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    cout << "Тест 12 Параллельные отрезки " << endl;
    seg1 = Segment3D(Vector3D(0.0, 1.0, 0.0), Vector3D(1.0, 1.0, 0.0));
    seg2 = Segment3D(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    cout << Intersect(seg1, seg2, eps_zero) << endl;
    
    
    return 0;
}

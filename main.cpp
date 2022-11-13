#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iomanip>
#include <numeric>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

class Vector
{
    vector<double> vektor;
public:
    explicit Vector(int n)
    {
        if (n <= 0) throw range_error("Bad dimension");
        vektor = vector<double> (n,0);
    }
    Vector(std::initializer_list<double> l)
    {
        if (l.size() == 0) throw range_error("Bad dimension");
        for (auto it = l.begin(); it != l.end(); it++) vektor.push_back(*it);
    }
    int NElems() const
    {
        return vektor.size();
    }
    double &operator[](int i)
    {
        return vektor[i];
    }
    double operator[](int i) const
    {
        return vektor[i];
    }
    double &operator()(int i)
    {
        if (i < 1 || i > vektor.size()) throw range_error("Invalid index");
        return vektor[i-1];
    }
    double operator()(int i) const
    {
        if (i < 1 || i > vektor.size()) throw range_error("Invalid index");
        return vektor[i-1];
    }
    double Norm() const
    {
        double s(0);
        for (int i = 0; i < vektor.size(); i++) {
            s += vektor[i] * vektor[i];
        }
        return sqrt(s);
    }
    friend double VectorNorm(const Vector &v)
    {
        double s(0);
        for (int i = 0; i < v.NElems(); i++) {
            s += v[i] * v[i];
        }
        return sqrt(s);
    }
    double GetEpsilon() const
    {
        return numeric_limits<double>::epsilon() * 10 * (this->Norm());
    }
    void Print(char separator = '\n', double eps = -1) const
    {
        for (int i = 0; i < vektor.size()-1; i++) {
            if (abs(vektor[i]) < eps) cout << 0 << separator;
            else cout << vektor[i] << separator;
        }
        if (abs(vektor[vektor.size()-1]) < eps) cout << 0;
        else cout << vektor[vektor.size()-1];
    }
    friend void PrintVector(const Vector &v, char separator = '\n', double eps = -1)
    {
        for (int i = 0; i < v.NElems()-1; i++) {
            if (abs(v[i]) < eps) cout << 0 << separator;
            else cout << v[i] << separator;
        }
        if (abs(v[v.NElems()-1]) < eps) cout << 0;
        else cout << v[v.NElems()-1];
    }
    friend Vector operator +(const Vector &v1, const Vector &v2)
    {
        if (v1.NElems() != v2.NElems()) throw domain_error("Incompatible formats");
        Vector v(v1.NElems());
        for (int i = 0; i<v.NElems(); i++) v[i] = v1[i]+v2[i];
        return v;
    }
    Vector &operator +=(const Vector &v)
    {
        if (this->NElems() != v.NElems()) throw domain_error("Incompatible formats");
        for (int i = 0; i < NElems(); i++) (*this)[i] += v[i];
        return *this;
    }
    friend Vector operator -(const Vector &v1, const Vector &v2)
    {
        if (v1.NElems() != v2.NElems()) throw domain_error("Incompatible formats");
        Vector v(v1.NElems());
        for (int i = 0; i<v.NElems(); i++) v[i] = v1[i]-v2[i];
        return v;
    }
    Vector &operator -=(const Vector &v)
    {
        if (this->NElems() != v.NElems()) throw domain_error("Incompatible formats");
        for (int i = 0; i < NElems(); i++) (*this)[i] -= v[i];
        return *this;
    }
    friend Vector operator *(double s, const Vector &v)
    {
        Vector povrat (v.NElems());
        for (int i = 0; i < povrat.NElems(); i++) povrat[i] = v[i]*s;
        return povrat;
    }
    friend Vector operator *(const Vector &v, double s)
    {
        Vector povrat (v.NElems());
        for (int i = 0; i < povrat.NElems(); i++) povrat[i] = v[i]*s;
        return povrat;
    }
    Vector &operator *=(double s)
    {
        for (int i = 0; i < NElems(); i++) (*this)[i]*=s;
        return *this;
    }
    friend double operator *(const Vector &v1, const Vector &v2)
    {
        if (v1.NElems() != v2.NElems()) throw domain_error("Incompatible formats");
        return inner_product(begin(v1.vektor),end(v1.vektor),begin(v2.vektor),0);

    }
    friend Vector operator /(const Vector &v, double s)
    {
        if (abs(s) < 0.00001) throw domain_error("Division by zero");
        Vector a(v.NElems());
        for (int i = 0; i < v.NElems(); i++) a[i] = (double)v[i]/s;
        return a;
    }
    Vector &operator /=(double s)
    {
        if (abs(s)<0.00001) throw domain_error("Division by zero");
        for (int i = 0; i<NElems(); i++) (*this)[i]/=s;
        return *this;
    }
    void Chop(double eps = -1)
    {
        double prag(eps);
        if (eps<0) prag = GetEpsilon();
        for (auto &x : vektor) if (x<prag) x = 0;
    }
    bool EqualTo(const Vector &v, double eps = -1) const
    {
        double tolerancija(eps);
        if (eps<0) tolerancija=GetEpsilon();
        if (NElems()!=v.NElems()) return false;
        for (int i = 0; i < NElems(); i++)
            if (fabs(vektor.at(i)-v.vektor.at(i))>tolerancija) return false;
        return true;
    }
};

class Matrix
{
    vector<vector<double>> matrica;
public:
    Matrix(int m, int n)
    {
        if (n <= 0 || m <= 0) throw range_error("Bad dimension");
        matrica = vector<vector<double>>(m,vector<double>(n,0));

    }
    Matrix(const Vector &v)
    {
        matrica = vector<vector<double>>(1,vector<double>(v.NElems(),0));
        for (int i = 0; i<v.NElems(); i++) matrica[0][i] = v[i];
    }
    Matrix(std::initializer_list<std::vector<double>> l)
    {
        if (l.size() == 0) throw range_error("Bad dimension");
        for (auto it2 = l.begin(); it2!=l.end(); it2++) {
            if (it2->size() == 0) throw range_error("Bad dimension");
        }
        for (auto it = l.begin() + 1; it != l.end(); it++)
            if (it->size()!=(it-1)->size()) throw logic_error("Bad matrix");
        matrica.resize(l.size());
        int i = 0;
        for (auto it = l.begin(); it!=l.end(); it++) {
            matrica[i] = *it;
            i++;
        }
    }
    Matrix(vector<vector<double>> mat) {
        matrica = mat;
    }

    int NRows() const
    {
        return matrica.size();
    }
    int NCols() const
    {
        return matrica.at(0).size();
    }
    double *operator[](int i)
    {
        return &matrica.at(i).at(0);
    }
    const double *operator[](int i) const
    {
        return &matrica.at(i).at(0);
    }
    double &operator()(int i, int j)
    {
        if (i<0 || j<0 || i>matrica.size() || j>matrica.at(0).size())
            throw range_error("Invalid index");
        return matrica.at(i-1).at(j-1);
    }
    double operator()(int i, int j) const
    {
        if (i<0 || j<0 || i>matrica.size() || j>matrica.at(0).size())
            throw range_error("Invalid index");
        return matrica.at(i-1).at(j-1);
    }
    double Norm() const
    {
        double s = 0;
        for (int i = 0; i < matrica.size(); i++) {
            for (int j = 0; j < matrica.at(i).size(); j++) {
                s += matrica.at(i).at(j) * matrica.at(i).at(j);
            }
        }
        return sqrt(s);
    }
    friend double MatrixNorm(const Matrix &m)
    {
        double s = 0;
        for (int i = 0; i < m.matrica.size(); i++) {
            for (int j = 0; j < m.matrica.at(i).size(); j++) {
                s += m.matrica.at(i).at(j) * m.matrica.at(i).at(j);
            }
        }
        return sqrt(s);
    }
    double GetEpsilon() const
    {
        return numeric_limits<double>::epsilon() * 10 * (this->Norm());
    }
    void Print(int width = 10, double eps = -1) const
    {
        for (int i = 0; i < matrica.size(); i++) {
            for (int j = 0; j < matrica.at(i).size(); j++) {
                if (abs(matrica.at(i).at(j)) < eps) cout << setw(width) << 0;
                else if (matrica.at(i).at(j)<0) cout << setw(width+1) << matrica.at(i).at(j);
                else cout<<setw(width)<<matrica.at(i).at(j);
            }
            if (i!=matrica.size()-1) cout << endl;
        }
    }
    friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1)
    {
        for (int i = 0; i < m.matrica.size(); i++) {
            for (int j = 0; j < m.matrica.at(i).size(); j++) {
                if (abs(m.matrica.at(i).at(j)) < eps) cout << setw(width) << 0;
                else cout << setw(width) << m.matrica.at(i).at(j);
            }
            if (i!=m.matrica.size()-1) cout << endl;
        }
    }
    friend Matrix operator +(const Matrix &m1, const Matrix &m2)
    {
        if (m1.NRows() != m2.NRows() || m1.NCols() != m2.NCols())
            throw domain_error("Incompatible formats");
        Matrix matrica(m1.NRows(),m1.NCols());
        for (int i = 0; i < matrica.NRows(); i++) {
            for (int j = 0; j < matrica.NCols(); j++) {
                matrica[i][j] = m1[i][j] + m2[i][j];
            }
        }
        return matrica;
    }
    Matrix &operator +=(const Matrix &m)
    {
        if (m.NRows()!=matrica.size() || m.NCols() != this->NCols())
            throw domain_error("Incompatible formats");
        for (int i = 0; i < matrica.size(); i++) {
            for (int j=0; j < matrica.at(0).size(); j++) {
                matrica[i][j] += m[i][j];
            }
        }
        return *this;
    }
    friend Matrix operator -(const Matrix &m1, const Matrix &m2)
    {
        if (m1.NRows() != m2.NRows() || m1.NCols() != m2.NCols())
            throw domain_error("Incompatible formats");
        Matrix matrica(m1.NRows(),m1.NCols());
        for (int i = 0; i < matrica.NRows(); i++) {
            for (int j = 0; j < matrica.NCols(); j++) {
                matrica[i][j] = m1[i][j] - m2[i][j];
            }
        }
        return matrica;
    }
    Matrix &operator -=(const Matrix &m)
    {
        if (m.NRows()!=matrica.size() || m.NCols() != this->NCols())
            throw domain_error("Incompatible formats");
        for (int i = 0; i < matrica.size(); i++) {
            for (int j=0; j < matrica.at(0).size(); j++) {
                matrica[i][j] -= m[i][j];
            }
        }
        return *this;
    }
    friend Matrix operator *(double s, const Matrix &m)
    {
        Matrix matrica(m.NRows(),m.NCols());
        for (int i = 0; i < matrica.NRows(); i++) {
            for (int j = 0; j < matrica.NCols(); j++) {
                matrica[i][j] = m[i][j]*s;
            }
        }
        return matrica;
    }
    friend Matrix operator *(const Matrix &m, double s)
    {
        Matrix matrica(m.NRows(),m.NCols());
        for (int i = 0; i < matrica.NRows(); i++) {
            for (int j = 0; j < matrica.NCols(); j++) {
                matrica[i][j] = m[i][j]*s;
            }
        }
        return matrica;
    }
    Matrix &operator *=(double s)
    {
        for (int i = 0; i < matrica.size(); i++) {
            for (int j = 0; j < matrica.at(i).size(); j++) {
                matrica[i][j]*=s;
            }
        }
        return *this;
    }
    friend Matrix operator *(const Matrix &m1, const Matrix &m2)
    {
        if (m1.NCols() != m2.NRows()) throw domain_error("Incompatible formats");
        Matrix matrica (m1.NRows(),m2.NCols());
        for (int i = 0; i < m1.NRows(); i++) {
            for (int j = 0; j < m2.NCols(); j++) {
                for (int k = 0; k < m1.NCols(); k++) {
                    matrica[i][j] += m1[i][k]*m2[k][j];
                }
            }
        }
        return matrica;
    }
    Matrix &operator *=(const Matrix &m)
    {
        if (NCols() != m.NRows()) throw domain_error("Incompatible formats");
        Matrix matrica2 (NRows(),m.NCols());
        for (int i = 0; i < NRows(); i++) {
            for (int j = 0; j < m.NCols(); j++) {
                for (int k = 0; k < NCols(); k++) {
                    matrica2[i][j] += matrica[i][k]*m[k][j];
                }
            }
        }
        *this = matrica2;
        return *this;
    }
    friend Vector operator *(const Matrix &m, const Vector &v)
    {
        if (m.NCols()!=v.NElems()) throw domain_error("Incompatible formats");
        Vector vektor (m.NRows());
        for (int i = 0; i < m.NRows(); i++) {
            for (int j = 0; j < m.NCols(); j++) {
                vektor[i] += m[i][j]*v[j];
            }
        }
        return vektor;
    }
    friend Matrix Transpose(const Matrix &m)
    {
        Matrix nova (m.NCols(),m.NRows());
        for (int i = 0; i < m.NRows(); i++) {
            for (int j= 0; j < m.NCols(); j++) {
                nova[j][i] = m[i][j];
            }
        }
        return nova;
    }
    void Transpose()
    {
        Matrix nova (NCols(),NRows());
        for (int i = 0; i < NRows(); i++) {
            for (int j = 0; j < NCols(); j++) {
                nova[j][i] = matrica[i][j];
            }
        }
        matrica=nova.matrica;
    }
    void Chop(double eps = -1)
    {
        double prag(eps);
        if (eps<0) prag=GetEpsilon();
        for (int i = 0; i < NRows(); i++)
            for (int j = 0; j < NCols(); j++)
                if (fabs(matrica.at(i).at(j))<prag) matrica.at(i).at(j)=0;
    }
    bool EqualTo(const Matrix &m, double eps = -1) const
    {
        double tolerancija(eps);
        if (eps<0) tolerancija=GetEpsilon();
        if (NRows()!=m.NRows() || NCols()!=m.NCols()) return false;
        for (int i = 0; i < NRows(); i++)
            for (int j = 0; j < NCols(); j++)
                if (fabs(matrica.at(i).at(j)-m.matrica.at(i).at(j))>tolerancija) return false;
        return true;
    }
    friend Matrix LeftDiv(Matrix m1, Matrix m2)
    {
        if (m1.NRows()!=m1.NCols()) throw domain_error("Divisor matrix is not square");
        if (m1.NRows()!=m2.NRows()) throw domain_error("Incompatible formats");
        int n = m1.NRows();
        for (int k = 0; k < n; k++) {
            int p = k;
            for (int i = k+1; i < n; i++) {
                if (fabs(m1.matrica.at(i).at(k))>fabs(m1.matrica.at(p).at(k)))
                    p=i;
            }
            if (fabs(m1.matrica.at(p).at(k))<m1.GetEpsilon()) throw domain_error("Divisor matrix is singular");
            if (p!=k) {
                m1.matrica[k].swap(m1.matrica[p]);
                m2.matrica[k].swap(m2.matrica[p]);
            }
            for (int i = k+1; i < n; i++) {
                double mi = m1.matrica[i][k] / m1.matrica[k][k];
                for (int j = k+1; j < n; j++)
                    m1.matrica[i][j]=m1.matrica[i][j]-mi*m1.matrica[k][j];
                for (int j = 0; j < m2.NCols(); j++)
                    m2.matrica[i][j]-=mi*m2.matrica[k][j];
            }
        }
        Matrix x(m1.NRows(),m2.NCols());
        for (int k = 0; k<m2.NCols(); k++) {
            for (int i = m1.NRows()-1; i>=0; i--) {
                double s = m2.matrica[i][k];
                for (int j = i+1; j < n; j++)
                    s-=m1.matrica[i][j]*x.matrica[j][k];
                x[i][k]=s/m1.matrica[i][i];
            }
        }
        return x;
    }
    friend Vector LeftDiv(Matrix m, Vector v)
    {
        if (m.NRows()!=m.NCols()) throw domain_error("Divisor matrix is not square");
        if (m.NRows()!=v.NElems()) throw domain_error("Incompatible formats");
        int n = m.NRows();
        for (int k = 0; k < n; k++) {
            int p = k;
            for (int i = k+1; i < n; i++) {
                if (fabs(m.matrica.at(i).at(k))>fabs(m.matrica.at(p).at(k)))
                    p=i;
            }
            if (fabs(m.matrica.at(p).at(k))<m.GetEpsilon()) throw domain_error("Divisor matrix is singular");
            if (p!=k) {
                m.matrica[k].swap(m.matrica[p]);
                swap(v[k],v[p]);
            }
            for (int i = k+1; i < n; i++) {
                double mi = m.matrica[i][k] / m.matrica[k][k];
                for (int j = k+1; j < n; j++)
                    m.matrica[i][j]=m.matrica[i][j]-mi*m.matrica[k][j];
                v[i]-=mi*v[k];
            }
        }
        Vector x(n);
        for (int i = n-1; i>=0; i--) {
            double s = v[i];
            for (int j=i+1; j<n; j++) {
                s-=m.matrica[i][j]*x[j];
            }
            x[i]=s/m.matrica.at(i).at(i);
        }
        return x;

    }
    friend Matrix operator /(const Matrix &m, double s)
    {
        if (s < 0.00001) throw domain_error("Division by zero");
        Matrix nova(m);
        nova/=s;
        return nova;
    }
    Matrix &operator /=(double s)
    {
        if (s < 0.00001) throw domain_error("Division by zero");
        for (int i = 0; i < NRows(); i++)
            for (int j = 0; j < NCols(); j++)
                matrica.at(i).at(j)/=s;
        return *this;
    }
    friend Matrix operator /(Matrix m1, Matrix m2)
    {
        return m1/=m2;
    }
    Matrix &operator /=(Matrix m)
    {
        if (m.NRows()!=m.NCols()) throw domain_error("Divisor matrix is not square");
        if (NCols()!=m.NRows()) throw domain_error("Incompatible formats");

        for (int k = 0; k < m.NRows(); k++) {
            int p = k;
            for (int i = k+1; i < m.NRows(); i++) {
                if (fabs(m.matrica.at(k).at(i))>fabs(m.matrica.at(k).at(p)))
                    p=i;
            }
            if (fabs(m.matrica.at(k).at(p))<m.GetEpsilon()) throw domain_error("Divisor matrix is singular");
            if (p!=k) {
                for (int i = 0; i < m.NRows(); i++) swap(m[i][k],m[i][p]);
                for (int i = 0; i < NRows(); i++) swap(matrica[i][k],matrica[i][p]);
            }
            for (int i = k+1; i < m.NRows(); i++) {
                double mi = m.matrica[k][i] / m.matrica[k][k];
                for (int j = k+1; j < m.NRows(); j++)
                    m.matrica[j][i]-=mi*m.matrica[j][k];
                for (int j = 0; j < NRows(); j++)
                    matrica[j][i]-=mi*matrica[j][k];
            }
        }
        Matrix x(NRows(),m.NRows());
        for (int k = 0; k<NRows(); k++) {
            for (int i = m.NRows()-1; i>=0; i--) {
                double s = matrica[k][i];
                for (int j = i+1; j < m.NRows(); j++)
                    s-=m.matrica[j][i]*x.matrica[k][j];
                x[k][i]=s/m.matrica[i][i];
            }
        }
        matrica=x.matrica;
        return *this;
    }
    friend double Det(Matrix m)
    {
        return m.Det();
    }
    double Det() const
    {
        Matrix n(*this);
        if (n.NRows()!=n.NCols()) throw domain_error("Matrix is not square");
        double d = 1;
        for (int k = 0; k < n.NRows(); k++) {
            int p = k;
            for (int i = k+1; i < n.NRows(); i++) {
                if (fabs(n.matrica.at(i).at(k))>fabs(n.matrica.at(p).at(k)))
                    p=i;
            }
            if (fabs(n.matrica.at(p).at(k))<n.GetEpsilon()) return 0;
            if (p!=k) {
                swap(n.matrica.at(k),n.matrica.at(p));
                d*=-1;
            }
            d*=n.matrica.at(k).at(k);
            for (int i = k+1; i < n.NRows(); i++) {
                double mi = n.matrica.at(i).at(k)/n.matrica.at(k).at(k);
                for (int j = k+1; j < n.NRows(); j++) n.matrica.at(i).at(j)-=mi*n.matrica.at(k).at(j);
            }
        }
        return d;
    }
    void Invert()
    {
        if (NCols()!=NRows()) throw domain_error("Matrix is not square");
        Vector w(NRows());
        for (int k = 0; k < NRows(); k++) {
            int p = k;
            for (int i = k+1; i < NRows(); i++)
                if (fabs(matrica.at(i).at(k))>matrica.at(p).at(k)) p=i;
            if (fabs(matrica.at(p).at(k))<GetEpsilon()) throw domain_error("Matrix is singular");
            if (p!=k) swap(matrica.at(p),matrica.at(k));
            w[k]=p;
            double mi = matrica.at(k).at(k);
            matrica.at(k).at(k)=1;
            for (int j = 0; j < NRows(); j++) matrica.at(k).at(j)/=mi;
            for (int i = 0; i < NRows(); i++) {
                if (i!=k) {
                    mi = matrica.at(i).at(k);
                    matrica.at(i).at(k)=0;
                    for (int j = 0; j < NRows(); j++)
                        matrica.at(i).at(j)-=mi*matrica.at(k).at(j);
                }
            }
        }
        for (int j = NRows()-1; j>=0; j--) {
            int p = w[j];
            if (p!=j) {
                for (int i=0; i<NRows(); i++) swap(matrica.at(i).at(p),matrica.at(i).at(j));
            }
        }
    }
    friend Matrix Inverse(Matrix m)
    {
        Matrix n(m);
        n.Invert();
        return n;
    }
    void ReduceToRREF()
    {
        int k = 0, l = 0;
        double eps = GetEpsilon();
        int m = NRows(), n = NCols();
        int p;
        vector<bool> w(NCols());
        for (int j = 0; j<n; j++)
            w.at(j)=false;
        while (k<m && l<n) {


            double v = 0;
            while (v<eps && l<n) {
                p=k;
                for (int i = k; i<m; i++) {
                    if (fabs(matrica.at(i).at(l))>v) {
                        v=fabs(matrica.at(i).at(l));
                        p=i;
                    }
                }
                if (v<eps)
                    l+=1;
            }
            if (l<n) {
                w.at(l)=true;
                if (p!=k)
                    swap(matrica.at(k),matrica.at(p));
                double mi = matrica.at(k).at(l);
                for (int j = l; j<n; j++) matrica.at(k).at(j)/=mi;
                for (int i = 0; i < m; i++) {
                    if (i!=k) {
                        mi = matrica.at(i).at(l);
                        for (int j = l; j <n; j++) matrica.at(i).at(j)-=mi*matrica.at(k).at(j);
                    }
                }
            }
            l+=1;
            k+=1;
        }
    }
    friend Matrix RREF(Matrix m)
    {
        m.ReduceToRREF();
        return m;
    }
    int Rank() const
    {
        Matrix pom(*this);
        int k = -1, l = -1;
        double eps = pom.GetEpsilon();
        int m = pom.NRows(), n = pom.NCols();
        int p;
        vector<bool> w(pom.NCols());
        for (int j = 0; j<n; j++)
            w.at(j)=false;
        while (k<m && l<n) {
            l+=1;
            k+=1;
            double v = 0;
            while (v<eps && l<n) {
                p=k;
                for (int i = k; i<m; i++) {
                    if (fabs(pom.matrica.at(i).at(l))>v) {
                        v=fabs(pom.matrica.at(i).at(l));
                        p=i;
                    }
                }
                if (v<eps)
                    l+=1;
            }
            if (l<n) {
                w.at(l)=true;
                if (p!=k)
                    swap(pom.matrica.at(k),pom.matrica.at(p));
                double mi = pom.matrica.at(k).at(l);
                for (int j = l; j<n; j++) pom.matrica.at(k).at(j)/=mi;
                for (int i = 0; i < m; i++) {
                    if (i!=k) {
                        mi = pom.matrica.at(i).at(l);
                        for (int j = l; j <n; j++) pom.matrica.at(i).at(j)-=mi*pom.matrica.at(k).at(j);
                    }
                }
            }

        }
        return k;
    }
    friend int Rank(Matrix m)
    {
        return m.Rank();
    }
};

class LUDecomposer
{
    Matrix mat;
    vector<int> v;
public:
    LUDecomposer(Matrix m) : v(m.NRows()), mat(m.NRows(),m.NCols())
    {
        if (m.NRows()!=m.NCols()) throw domain_error("Matrix is not square");
        int n = m.NRows();
        double e = m.GetEpsilon();
        double s;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < j+1; i++) {
                s=m[i][j];
                for (int k = 0; k < i; k++)
                    s-=m[i][k]*m[k][j];
                m[i][j]=s;
            }
            int p = j;
            for (int i = j+1; i < n; i++) {
                s=m[i][j];
                for (int k = 0; k < j; k++) {
                    s-=m[i][k]*m[k][j];
                }
                m[i][j]=s;
                if (fabs(s)>fabs(m[p][j])) p=i;
            }
            if (fabs(m[p][j])<e) throw domain_error("Matrix is singular");
            if (p!=j)
                for (int temp = 0; temp < n; temp++) swap(m[j][temp],m[p][temp]);
            v[j]=p;
            double mi = m[j][j];
            for (int i = j+1; i < n; i++) m[i][j]/=mi;
        }
        mat=m;
    }
    void Solve(const Vector &b, Vector &x) const
    {
        Vector a(b);
        if (mat.NRows()!=b.NElems() || mat.NRows()!=x.NElems()) throw domain_error("Incompatible formats");
        int p;
        double s;
        for (int i = 0; i < mat.NRows(); i++) {
            p = v[i];
            s = a[p];
            a[p] = a[i];
            for (int j = 0; j < i; j++) s-=mat[i][j]*x[j];
            x[i]=s;
        }
        for (int i = mat.NRows()-1; i>=0; i--) {
            s=x[i];
            for (int j=i+1; j<mat.NRows(); j++)
                s-=mat[i][j]*x[j];
            x[i]=s/mat[i][i];
        }
    }
    Vector Solve(Vector b) const
    {
        Vector x = b;
        this->Solve(x,b);
        return b;
    }
    void Solve(Matrix b, Matrix &x) const
    {
        if (x.NRows()!=mat.NRows() || x.NCols()!=b.NCols() || mat.NCols()!=b.NRows())
            throw domain_error("Incompatible formats");
        double s;
        //int prve = mat.NCols(), druge = b.NCols();
        for (int k = 0; k < b.NCols(); k++) {
            for (int i = 0; i < mat.NCols(); i++) {
                int p = v[i];
                s = b[p][k];
                b[p][k]=b[i][k];
                for (int j = 0; j < i; j++)
                    s-=mat[i][j]*x[j][k];
                x[i][k]=s;
            }
            for (int i = mat.NCols()-1; i>=0; i--) {
                s=x[i][k];
                for (int j = i+1; j < mat.NCols(); j++)
                    s-=mat[i][j]*x[j][k];
                x[i][k]=s/mat[i][i];
            }
        }
    }
    Matrix Solve(Matrix b) const
    {
        Matrix nova(b);
        this->Solve(nova,b);
        return b;
    }
    Matrix GetCompactLU() const
    {
        return mat;
    }
    Matrix GetL() const
    {
        Matrix L(mat.NRows(),mat.NCols());
        for (int i = 0; i < mat.NRows(); i++)
            for (int j = 0; j <= i; j++)
                L[i][j]=mat[i][j];
        for (int i = 0; i < mat.NRows(); i++) L[i][i]=1;
        return L;
    }
    Matrix GetU() const
    {
        Matrix U(mat.NRows(),mat.NCols());
        for (int i = 0; i < mat.NRows(); i++)
            for (int j = i; j < mat.NCols(); j++)
                U[i][j]=mat[i][j];
        return U;
    }
    Vector GetPermuation() const
    {
        Vector povrat(v.size());
        for (int i = 0; i < povrat.NElems(); i++) povrat[i]=v[i]+1;
        return povrat;
    }
};

int main ()
{
    try {
        vector<vector<double>> matrix;
        auto f = []() -> int { return rand() % 10000; };
        for (int i = 0; i < 100; i++) {
            vector<double> values(100);
            generate(values.begin(), values.end(), f);
            matrix.emplace_back(values);
        }
        Matrix M1(matrix);

        auto start = high_resolution_clock::now();
        LUDecomposer lu(M1);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << endl << "Time taken by function: " << duration.count() << " microseconds" << endl;
    } catch (std::domain_error e) {
        std::cout << "'" << e.what() << "'";
    }
    return 0;
}
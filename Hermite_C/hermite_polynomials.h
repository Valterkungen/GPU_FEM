#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;

class Polynomial {
private:
    vector<double> coeff;
public:
    Polynomial(const vector<double>& coeff) : coeff(coeff) {}

    double operator()(double x) {
        double result = 0;
        for (size_t n = 0; n < coeff.size(); ++n) {
            result += coeff[n] * pow(x, n);
        }
        return result;
    }

    Polynomial operator*(const Polynomial& other) const {
        vector<double> result_coeff(coeff.size() + other.coeff.size() - 1, 0);
        for (size_t i = 0; i < coeff.size(); ++i) {
            for (size_t j = 0; j < other.coeff.size(); ++j) {
                result_coeff[i + j] += coeff[i] * other.coeff[j];
            }
        }
        return Polynomial(result_coeff);
    }

    Polynomial derive() const {
        vector<double> derived_coeff(coeff.size() - 1, 0);
        for (size_t i = 1; i < coeff.size(); i++) {
            derived_coeff[i - 1] = i * coeff[i];
        }
        return Polynomial(derived_coeff);
    }

    double integrate(double low_bound, double high_bound) const {
        vector<double> integral_coeff(coeff.size() + 1, 0);
        for (size_t i = 0; i < coeff.size(); i++) {
            integral_coeff[i + 1] = coeff[i] / (i + 1);
        }
        Polynomial integral(integral_coeff);
        return integral(high_bound) - integral(low_bound);
    }

    friend ostream& operator<<(ostream& os, const Polynomial& poly) {
        stringstream ss;
        for (int i = poly.coeff.size() - 1; i >= 0; i--) {
            if (poly.coeff[i] != 0) {
                ss << poly.coeff[i] << "x^" << i;
                if (i != 0) {
                    ss << " + ";
                }
            }
        }
        string str = ss.str();
        size_t pos = str.find(" + -");
        if (pos != string::npos) {
            str.replace(pos, 4, "- ");
        }
        str.replace(str.find("x^0"), 3, "");
        os << str;
        return os;
    }
};



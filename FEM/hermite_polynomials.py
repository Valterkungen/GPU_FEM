class HermitePolynomial:
    def __init__(self, n, coeff=None):
        if coeff is not None:
            self.coeff = coeff
        elif n == 0:
            self.coeff = [1]
        elif n == 1:
            self.coeff = [0, 1]
        else:
            coeff1 = HermitePolynomial(n-1).coeff
            coeff2 = HermitePolynomial(n-2).coeff
            self.coeff = [0] + coeff1
            self.coeff = [
                a - (n-1)*b
                for a, b
                in zip(self.coeff, coeff2 + [0, 0])
            ]
    
    def __call__(self, x):
        return sum([
            k*x**n
            for n, k
            in enumerate(self.coeff)
        ])

    def __str__(self):
        return ' + '.join(reversed([
            f'{x}x^{n}'
            for n, x
            in enumerate(self.coeff)
            if x != 0
        ])).replace('+ -', '- ').replace('x^0', '')
    
    def __mul__(self, other):
        coeff = [0] * (len(self.coeff) + len(other.coeff) -1)

        for i, coeff_a in enumerate(self.coeff):
            for j, coeff_b in enumerate(other.coeff):
                coeff[i + j] += coeff_a * coeff_b
            
        return HermitePolynomial(0, coeff)

    def derive(self):
        return HermitePolynomial(0, [n*koeff for n, koeff in enumerate(self.coeff)][1:])

    def integrate(self, low_bound=0, high_bound=1):
        p = HermitePolynomial(0, [0] + [
            coeff/(n+1)
            for n, coeff
            in enumerate(self.coeff)
        ])

        return p(high_bound) - p(low_bound)

def main():
    print(HermitePolynomial(5)) # Print hermite polynomial of order 5

    p1 = HermitePolynomial(3) # Create hermite polynomial of order 3
    p2 = HermitePolynomial(3) # Create hermite polynomial of order 3

    print(p1(4)) # Evaluate polynomial p1 at x = 4

    p3 = p1 * p2# Calcualate the polynomial p1 * p2'
    print(p3.integrate()) # Calculate the integral of p3 from 0 to 1

if __name__ == '__main__':
    main()

class Polynomial:
    def __init__(self,coeff):
        self.coeff = coeff
    
    def __call__(self,x):
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

        return Polynomial(coeff=coeff)

    def derive(self):
        return Polynomial(coeff=[n*koeff for n, koeff in enumerate(self.coeff)][1:])

    def integrate(self,low_bound=0, high_bound=1):
        p = Polynomial(coeff=[0] + [
            coeff/(n+1)
            for n, coeff
            in enumerate(self.coeff)
        ])
        return p(high_bound) - p(low_bound)

def main():
    
    p1 = Polynomial([1,1,1,1]) # Create hermite polynomial of order 3
    #p2 = Polynomial(3) # Create hermite polynomial of order 3
    #print(p1(0)) # Evaluate polynomial p1 at x = 4

    p3 = p1 
    #print(p3.integrate()) # Calculate the integral of p3 from 0 to 1

if __name__ == '__main__':
    main()
using System.Linq;
using System.Numerics;
using static System.Math;
using static System.Console;

namespace DirichletSeries
{
    static class Program
    {
        /// <summary>
        /// Calculates the converged point for a Dirichlet series expansion.
        /// </summary>
        /// <param name="t">imaginary part of s. The first zero is at 14.134725</param>
        /// <param name="numberOfTerms">Use a higher number to find more accurate convergence.</param>
        /// <returns></returns>
        /// 
        static Complex CalcZetaZero(this double t, int numberOfTerms)
        {
            var range = Enumerable.Range(1, numberOfTerms);
            var zetaZero = Complex.Zero;

            foreach (int n in range)
            {
                var direction = n % 2 == 0 ? PI : 0;
                var newTerm = Complex.Exp(new Complex(-Log(n) * .5, -Log(n) * t + direction));
                zetaZero += newTerm;
            }

            return zetaZero;
        }

        static void Main(string[] args)
        {
            for (int terms = 100000; terms <= 1000000000; terms *= 10)
                WriteLine("{0}", CalcZetaZero(14.134725141734, terms));
            
            //for (double t=14.134725; t<14.134726; t+=1E-8)
            //    Console.WriteLine("{0:0.0000000000} : {1}", t, CalcZetaZero(t, 10000000).Magnitude);

            Write("Press Enter: ");
            ReadLine();
        }
    }
}

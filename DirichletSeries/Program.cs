using System;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using static System.Console;
using static System.Math;

namespace DirichletSeries
{
    static class Program
    {
        static private uint[] primes;

        static void prime_sieve(uint n)
        {
            uint p;
            primes.Initialize();

            primes[0] = 2;
            for (p = 0; primes[p] < n;)
            {
                for (uint i = primes[p] * primes[p]; i < n; i += primes[p])
                    primes[i] = 1;
                primes[++p] = primes[p - 1] + 1;
                for (; primes[p] < n && primes[primes[p]] == 1; primes[p]++) ;     //find next prime (where s[p]==0)
            }
            Array.Resize(ref primes, (int)p);
        }

        static void more_primes(uint n)
        {
            Random r = new Random();

            for (uint p = primes[primes.Length - 1] + 2; p < n; p += 2)
            {
                if (p % 5 == 0)
                    p += 2;

                bool isPrime = BigInteger.ModPow(r.Next((int)p), p - 1, p) == 1;

                if (isPrime)
                {
                    Array.Resize(ref primes, primes.Length + 1);
                    primes[primes.Length - 1] = p;
                    Write($"{p}\r");
                }
            }
            WriteLine($"\n\nprimes.Length: {primes.Length}\n");
        }

        // Lanczos approximation of complex gamma function
        // https://en.wikipedia.org/wiki/Lanczos_approximation
        static Complex gamma(Complex z)
        {
            double[] p = {
                676.5203681218851,
                -1259.1392167224028,
                771.32342877765313,
                -176.61502916214059,
                12.507343278686905,
                -0.13857109526572012,
                9.9843695780195716e-6,
                1.5056327351493116e-7
            };
            Complex y;

            if (z.Real < 0.5)
                y = Math.PI / (Complex.Sin(Math.PI * z) * gamma(1 - z));  // Reflection formula
            else
            {
                z -= 1;
                Complex x = 0.99999999999980993d;
                int i = 0;
                foreach (double pval in p)
                {
                    x += pval / (z + i + 1.0d);
                    i++;
                }
                Complex t = z + p.Length - 0.5d;
                y = Math.Sqrt(2 * Math.PI) * Complex.Pow(t, (z + 0.5)) * Complex.Exp(-t) * x;
            }
            if (Math.Abs(y.Imaginary) <= double.Epsilon)
                return y.Real;
            else
                return y;
        }

        // Stirling's formula for the gamma function
        // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
        static Complex loggamma(Complex z)
        {
            Complex log_z, z_sinh_z, z_pow6;

            log_z = Complex.Log(z);
            z_sinh_z = z * Complex.Sinh(Complex.Reciprocal(z));
            z_pow6 = Complex.Pow(z, -6.0) / 810.0;
            var log_sinh_z_minus_z_pow6 = Complex.Log(z_sinh_z + z_pow6);
            var _2_log_z = 2.0 * log_z;
            var last_term = z * (_2_log_z + log_sinh_z_minus_z_pow6 - 2.0);
            var _2_loggamma = Log(2.0 * PI) - log_z + last_term;

            return _2_loggamma / 2.0;
        }

        static Complex rsiegeltheta(double t)
        {
            Complex a, b, theta;
            a = loggamma(new Complex(0.25, t / 2.0));
            b = loggamma(new Complex(0.25, -t / 2.0));

            theta = (a - b) * new Complex(0, -0.5);
            theta -= Log(PI) / 2.0 * t;

            return theta;
        }

        static double rsiegelz(double t, int n)
        {
            double ZZ = 0.0;
            double R = 0.0;
            double j = 1.0;
            int k = 0;
            int N = (int)Floor(Sqrt(t * 0.5 / PI));
            double p = Sqrt(t * 0.5 / PI) - N;

            while (j <= N)
            {
                ZZ = ZZ + 1.0 / Sqrt(j) * Cos((rsiegeltheta(t).Real - t * Log(j)) % (2.0 * PI));
                j += 1.0;
            }
            ZZ = 2.0 * ZZ;

            while (k <= n)
            {
                R = R + C(2.0 * p - 1.0, k) * Pow(2.0 * PI / t, (double)k * 0.5);
                ++k;
            }
            R = Pow(2.0 * PI / t, 0.25) * R;
            if (N % 2 == 0)
                R = -R;
            return ZZ + R;
        }

        static double C(double z, int n)
        {
            switch (n)
            {
                case 0:
                    return (.38268343236508977173d * Pow(z, 0.0)
                + .43724046807752044936d * Pow(z, 2.0)
                + .13237657548034352332d * Pow(z, 4.0)
                - .01360502604767418865d * Pow(z, 6.0)
                - .01356762197010358089d * Pow(z, 8.0)
                - .00162372532314446528d * Pow(z, 10.0)
                + .00029705353733379691d * Pow(z, 12.0)
                + .00007943300879521470d * Pow(z, 14.0)
                + .00000046556124614505d * Pow(z, 16.0)
                - .00000143272516309551d * Pow(z, 18.0)
                - .00000010354847112313d * Pow(z, 20.0)
                + .00000001235792708386d * Pow(z, 22.0)
                + .00000000178810838580d * Pow(z, 24.0)
                - .00000000003391414390d * Pow(z, 26.0)
                - .00000000001632663390d * Pow(z, 28.0)
                - .00000000000037851093d * Pow(z, 30.0)
                + .00000000000009327423d * Pow(z, 32.0)
                + .00000000000000522184d * Pow(z, 34.0)
                - .00000000000000033507d * Pow(z, 36.0)
                - .00000000000000003412d * Pow(z, 38.0)
                + .00000000000000000058d * Pow(z, 40.0)
                + .00000000000000000015d * Pow(z, 42.0));
                case 1:
                    return (-.02682510262837534703d * Pow(z, 1.0)
                + .01378477342635185305d * Pow(z, 3.0)
                + .03849125048223508223d * Pow(z, 5.0)
                + .00987106629906207647d * Pow(z, 7.0)
                - .00331075976085840433d * Pow(z, 9.0)
                - .00146478085779541508d * Pow(z, 11.0)
                - .00001320794062487696d * Pow(z, 13.0)
                + .00005922748701847141d * Pow(z, 15.0)
                + .00000598024258537345d * Pow(z, 17.0)
                - .00000096413224561698d * Pow(z, 19.0)
                - .00000018334733722714d * Pow(z, 21.0)
                + .00000000446708756272d * Pow(z, 23.0)
                + .00000000270963508218d * Pow(z, 25.0)
                + .00000000007785288654d * Pow(z, 27.0)
                - .00000000002343762601d * Pow(z, 29.0)
                - .00000000000158301728d * Pow(z, 31.0)
                + .00000000000012119942d * Pow(z, 33.0)
                + .00000000000001458378d * Pow(z, 35.0)
                - .00000000000000028786d * Pow(z, 37.0)
                - .00000000000000008663d * Pow(z, 39.0)
                - .00000000000000000084d * Pow(z, 41.0)
                + .00000000000000000036d * Pow(z, 43.0)
                + .00000000000000000001d * Pow(z, 45.0));
                case 2:
                    return (+.00518854283029316849d * Pow(z, 0.0)
                   + .00030946583880634746d * Pow(z, 2.0)
                   - .01133594107822937338d * Pow(z, 4.0)
                   + .00223304574195814477d * Pow(z, 6.0)
                   + .00519663740886233021d * Pow(z, 8.0)
                   + .00034399144076208337d * Pow(z, 10.0)
                   - .00059106484274705828d * Pow(z, 12.0)
                   - .00010229972547935857d * Pow(z, 14.0)
                   + .00002088839221699276d * Pow(z, 16.0)
                   + .00000592766549309654d * Pow(z, 18.0)
                   - .00000016423838362436d * Pow(z, 20.0)
                   - .00000015161199700941d * Pow(z, 22.0)
                   - .00000000590780369821d * Pow(z, 24.0)
                   + .00000000209115148595d * Pow(z, 26.0)
                   + .00000000017815649583d * Pow(z, 28.0)
                   - .00000000001616407246d * Pow(z, 30.0)
                   - .00000000000238069625d * Pow(z, 32.0)
                   + .00000000000005398265d * Pow(z, 34.0)
                   + .00000000000001975014d * Pow(z, 36.0)
                   + .00000000000000023333d * Pow(z, 38.0)
                   - .00000000000000011188d * Pow(z, 40.0)
                   - .00000000000000000416d * Pow(z, 42.0)
                   + .00000000000000000044d * Pow(z, 44.0)
                   + .00000000000000000003d * Pow(z, 46.0));
                case 3:
                    return (-.00133971609071945690d * Pow(z, 1.0)
               + .00374421513637939370d * Pow(z, 3.0)
               - .00133031789193214681d * Pow(z, 5.0)
               - .00226546607654717871d * Pow(z, 7.0)
               + .00095484999985067304d * Pow(z, 9.0)
               + .00060100384589636039d * Pow(z, 11.0)
               - .00010128858286776622d * Pow(z, 13.0)
               - .00006865733449299826d * Pow(z, 15.0)
               + .00000059853667915386d * Pow(z, 17.0)
               + .00000333165985123995d * Pow(z, 19.0)
               + .00000021919289102435d * Pow(z, 21.0)
               - .00000007890884245681d * Pow(z, 23.0)
               - .00000000941468508130d * Pow(z, 25.0)
               + .00000000095701162109d * Pow(z, 27.0)
               + .00000000018763137453d * Pow(z, 29.0)
               - .00000000000443783768d * Pow(z, 31.0)
               - .00000000000224267385d * Pow(z, 33.0)
               - .00000000000003627687d * Pow(z, 35.0)
               + .00000000000001763981d * Pow(z, 37.0)
               + .00000000000000079608d * Pow(z, 39.0)
               - .00000000000000009420d * Pow(z, 41.0)
               - .00000000000000000713d * Pow(z, 43.0)
               + .00000000000000000033d * Pow(z, 45.0)
               + .00000000000000000004d * Pow(z, 47.0));
                default:
                    return (+.00046483389361763382d * Pow(z, 0.0)
              - .00100566073653404708d * Pow(z, 2.0)
              + .00024044856573725793d * Pow(z, 4.0)
              + .00102830861497023219d * Pow(z, 6.0)
              - .00076578610717556442d * Pow(z, 8.0)
              - .00020365286803084818d * Pow(z, 10.0)
              + .00023212290491068728d * Pow(z, 12.0)
              + .00003260214424386520d * Pow(z, 14.0)
              - .00002557906251794953d * Pow(z, 16.0)
              - .00000410746443891574d * Pow(z, 18.0)
              + .00000117811136403713d * Pow(z, 20.0)
              + .00000024456561422485d * Pow(z, 22.0)
              - .00000002391582476734d * Pow(z, 24.0)
              - .00000000750521420704d * Pow(z, 26.0)
              + .00000000013312279416d * Pow(z, 28.0)
              + .00000000013440626754d * Pow(z, 30.0)
              + .00000000000351377004d * Pow(z, 32.0)
              - .00000000000151915445d * Pow(z, 34.0)
              - .00000000000008915418d * Pow(z, 36.0)
              + .00000000000001119589d * Pow(z, 38.0)
              + .00000000000000105160d * Pow(z, 40.0)
              - .00000000000000005179d * Pow(z, 42.0)
              - .00000000000000000807d * Pow(z, 44.0)
              + .00000000000000000011d * Pow(z, 46.0)
              + .00000000000000000004d * Pow(z, 48.0));
            }
        }

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
            double sum = 0;
            object monitor = new object();
            object monitor_zeta = new object();
            //Parallel.ForEach(range, () => 0.0, (n, state, local) =>
            foreach (int n in range)
            {
                var direction = n % 2 == 0 ? PI : 0;
                var newTerm = Complex.Exp(new Complex(-Log(n) * 0.5d, -Log(n) * t + direction));
                lock (monitor_zeta)
                    zetaZero += newTerm;
                //return local++;
            }
            //, local => { lock (monitor) sum += local; });

            return zetaZero;
        }

        static Complex Zeta_Function(Complex s)
        {
            Complex product = 1.0d;
            foreach (int p in primes)
            {
                product *= (1.0d - Complex.Pow(p, -s));
            }
            return 1.0d / product;
        }

        static void Zeta_Zeros(double t, double step)
        {
            for (double i = t; i < t + step * 10; i += step)
            {
                Complex s = new Complex(0.5, i);
                Complex z = Zeta_Function(s);
                if (z.Magnitude < 0.01)
                    WriteLine("{0}\n{1,30:F10}\t{2:F10}\t{3:F10}", s, z, z.Magnitude, z.Phase);
            }
        }

        static Complex Dirichlet2(Complex s, int terms)
        {
            Complex zeta = 0;
            double sum = 0;
            var range = Enumerable.Range(1, terms);
            object monitor_zeta = new object();
            object monitor = new object();

            //for (int n = 1; n < terms; n++)
            Parallel.ForEach(range, () => 0.0, (n, state, local) =>
            {
                double d = (double)n;
                Complex t1 = d / Complex.Pow(d + 1.0, s);
                Complex t2 = (d - s) / Complex.Pow(d, s);
                lock (monitor_zeta)
                    zeta += t1 - t2;
                return local++;
            }, local => { lock (monitor) sum += local; });
            zeta /= (s - 1.0);
            return zeta;

            /*
             * >>> zeta(0.5-1j)
             * mpc(real='0.14393642707718907', imag='0.72209974353167306')
             */
        }

        static void Main(string[] args)
        {
            //for (int terms = 100000; terms <= 1000000000; terms *= 10)
            //    WriteLine("{0}", CalcZetaZero(14.134725141734, terms));

            const int TERMS = 1000000;      // 525000;    // 100000000;
            const double initial = 14.1d;
            double step = 0.01d;
            double lastMag = 0.5d;
            int lastSign = Sign(CalcZetaZero(initial, TERMS).Phase);
            int lastSignMag = Sign(CalcZetaZero(initial, TERMS).Magnitude);
            double t = initial;

            /*
            primes = new uint[TERMS];
            prime_sieve(TERMS);
            //foreach (int p in primes)
            //    Write("{0,8}", p);
            more_primes(6000000);
            WriteLine($"{primes.Length}");
            
            Zeta_Zeros(initial, step);
            //WriteLine("{0}", Zeta_Function(new Complex(0.5d, 14.134d)));
            */

            for (t = initial; t < initial + step * 100.0d; t += step)
            {
                Complex c = Dirichlet2(new Complex(0.5, t), TERMS);
                double magnitude = c.Magnitude;
                double phase = c.Phase;
                int sign = Sign(lastMag - magnitude);
                double z = Sin(Complex.Abs(rsiegeltheta(t))) * rsiegelz(t, 4);
                WriteLine("{0:0.0000000000} : {1}\t{2}\t{3}\t{4}", t, sign, magnitude, phase, z);
                lastMag = magnitude;
            }
            /*
            //WriteLine("{0}", Dirichlet2(new Complex(0.5, 14.134725d), TERMS));
            //WriteLine("{0}", Dirichlet2(new Complex(0.5, -1.0), TERMS));
            Complex s = new Complex(1.0, 14.134725141734);
            Complex s1 = new Complex(0.5, 236.52422966581619);
            WriteLine("{0}", gamma(s));
            WriteLine("{0}", gamma(s1));
            WriteLine("{0}", loggamma(s1));
            WriteLine("{0}", rsiegeltheta(s1.Imaginary));
            WriteLine("{0}", rsiegelz(s1.Imaginary, 4));
            WriteLine("{0}", Sin(Complex.Abs(rsiegeltheta(s1.Imaginary))) * rsiegelz(s1.Imaginary, 4));
            */

            Write("Press Enter: ");
            ReadLine();
            return;

            while (Abs(step) > 1e-3)
            {
                {
                    var Z = Sin(Complex.Abs(rsiegeltheta(t))) * rsiegelz(t, 4);    //Dirichlet2(new Complex(0.5d, t), TERMS);    //CalcZetaZero(t, TERMS);     //
                    var c = Complex.Exp(new Complex(0, rsiegeltheta(t).Imaginary)) * Z;
                    double magnitude = c.Magnitude;
                    int sign = Sign(Z);
                    int signMag = Sign(magnitude);
                    WriteLine("{0:0.0000000000} : {1}\t{2}\t|c|: {3}\targ(c): {4}\n{5}", t, sign, signMag, magnitude, c.Phase, Z);

                    if (sign != lastSign)
                    {
                        t -= step;
                        step /= 10.0d;
                        lastSign = -sign;
                        WriteLine("-------------------");
                    }
                    else if (signMag != lastSignMag)
                    {
                        t -= step;
                        step /= 10.0d;
                        lastSignMag = -signMag;
                        WriteLine("-------------------");
                    }
                    else
                    {
                        lastMag = magnitude;
                        t += step;
                    }
                }
            }


            Write("Press Enter: ");
            ReadLine();
        }
    }
}

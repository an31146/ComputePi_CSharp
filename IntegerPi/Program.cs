using System;
using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Threading.Tasks;

using static System.Console;

namespace IntegerPi
{
#pragma warning disable IDE1006
    class Program
    {
        const uint STEPS = 10000000;
        const uint DIGITS = 2000;

        static public BigInteger SquareRoot(BigInteger n)
        {
            BigInteger d = BigInteger.One, q, _d;
            bool bBreakExpr = true; 
            Stopwatch sw = new Stopwatch();

            //Newton's Method
            sw.Start();
            do
            {
                q = n / d + d;
                _d = q;
                q >>= 1;                                // _q = q / 2;
                d = q;                                  // _d = 2*d

                bBreakExpr = q * q > n;                 // Can't figure out a test expression without using multiplication here :-(
                //bBreakExpr = 2 * q < _d; // && q * 2 != d;
            } while (bBreakExpr);
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine("\nSquareRoot({0})\nElapsed time: {1}\n", n, strElapsed);
#endif
            return q;
        }

        public static BigInteger Sqrt(BigInteger x)
        {
            BigInteger div = BigInteger.One;
            div = BigInteger.Add(div, BigInteger.One);
            BigInteger TWO = new BigInteger(2);
            BigInteger div2 = div, y;
            Stopwatch sw = new Stopwatch();

            // Loop until we hit the same value twice in a row, or wind
            // up alternating.
            sw.Start();
            while (true)
            {
                y = BigInteger.Add(div, BigInteger.Divide(x, div));
                y >>= 1;
                if (y.Equals(div) || y.Equals(div2))
                    break;
                div2 = div;
                div = y;
            }
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine("\nSquareRoot({0})\nElapsed time: {1}\n", y, strElapsed);
#endif
            return y;
        }

        static double zeta_of_two_double()
        {
            double sum;
            const double unity = 1.0;

            sum = unity;
            object monitor = new object();
            //for (uint i = 2; i < STEPS; i++)
            Parallel.For<double>(2, STEPS, () => 2, (i, state, local) =>
            {
                double b = (double)i * i;           // b = i²
                local = unity / b;                  // local = 1 / i²
                sum += local;                       // sum = Ʃ (1 / i²)
                //Write("{0} \r", local.ToString("F18", CultureInfo.InvariantCulture));
                return (local);
            }, local => {
                            lock (monitor) { sum += local; }
                        });
            return sum;
        }

        static double zeta_of_four_double()
        {
            double sum;
            const double unity = 1.0;

            sum = unity;
            object monitor = new object();
            Parallel.For<double>(2, STEPS, () => 2, (i, state, local) =>
            //for (uint i = 2; i < STEPS; i++)
            {
                double b = (double)i * i;           // b = i²
                b *= b;                             // b = i⁴
                local = unity / b;                  // local = 1 / i⁴
                sum += local;                       // sum = Ʃ (1 / i⁴) 
                //Write("{0} \r", local.ToString("F18", CultureInfo.InvariantCulture));
                return (local);
            }, local => {
                            lock (monitor) { sum += local; }
                        });
            return sum;
        }


        static BigInteger normalize(uint n, uint digits)
        {
            BigInteger ret_val = new BigInteger(n);
            BigInteger exponent = 10;

            uint i = digits;
            while (i > 0)
            {
                if ((i & 1) == 1)
                    exponent *= 10;

                exponent *= exponent;
                i /= 2;
            }
            ret_val *= exponent;

            return ret_val;
        }

        static BigInteger zeta_of_two_bigint()
        {
            BigInteger unity = normalize(1, DIGITS);
            BigInteger exponent = unity;
            unity *= unity;
            BigInteger sum = unity;
            unity *= unity;

            for (uint i = 2; i < STEPS; i++)
            {
                BigInteger b = i * exponent;
                b *= b;
                sum += unity / b;
                //Write("{0}\r", i);
            }

            WriteLine();
            //WriteLine("{0}", unity);
            return sum;
        }

        static void Main(string[] args)
        {
            double pi_squared_over_six = zeta_of_two_double();
            WriteLine("pi²/6: {0}\n\n√(pi²/6): {1}\n\n", pi_squared_over_six, Math.Sqrt(pi_squared_over_six * 6.0d));

            double pi_to_the_fourth_over_ninety = zeta_of_four_double();
            WriteLine("pi⁴/90: {0}\n\n(pi⁴/90)^¼: {1}\n\n", pi_to_the_fourth_over_ninety, Math.Pow(pi_to_the_fourth_over_ninety * 90.0d, 0.25d));

            BigInteger TWO = normalize(2, DIGITS);
            WriteLine("SquareRoot({0}) =\n{1}\n", TWO, SquareRoot(TWO));
            WriteLine("Sqrt({0}) =\n{1}\n", TWO, Sqrt(TWO));

            BigInteger BigInt_pi_squared_over_six = zeta_of_two_bigint();
            BigInteger BigInt_pi = SquareRoot(BigInt_pi_squared_over_six * 6);
            WriteLine("BigInt_pi²/6: {0}\n\n√(BigInt_pi²/6): {1}\n\n", BigInt_pi_squared_over_six.ToString(), BigInt_pi.ToString());

            Write("Press Enter: ");
            ReadLine();
        }
    }
}

using System;
using System.Diagnostics;
using System.Numerics;
using System.Threading.Tasks;

using static System.Console;

#pragma warning disable IDE1006,IDE0049,CS0168,IDE0039,CS0219
/*
 * IDE1006 Suppress Naming Rule Violation IDE1006
 * IDE0049 Name Can Be Simplified
 * IDE0039 Use local function
 * CS0168  The variable 'var' is declared but never used.
 * CS0219  The variable 'b' is assigned but its value is never used
 */

namespace IntegerPi
{
    class Program
    {
        const uint STEPS = 80000000;
        const uint DIGITS = 150;

        public static BigInteger SquareRoot(BigInteger n)
        {
            BigInteger div = n >> 1, _div, quotient, _doubleOfQuotient, _remainder;
            bool bBreakExpr = false; 
            Stopwatch sw = new Stopwatch();

            //Newton's Method
            sw.Start();
            do
            {
                quotient = n / div + div;
                _div = div;                                     // preserve div for comparison later
                //_doubleOfQuotient = quotient;                 // _doubleOfquotient;       preserve double of quotient before halving
                quotient >>= 1;                                 // right-shift quotient;
                div = quotient;                                 // assign next divisor
                bBreakExpr = _div > quotient;                   // test if _div is still greater than the quotient - this is how to do it!

                //bBreakExpr = quotient * quotient > n;           // Can't figure out a test expression without using multiplication here :-(
                //bBreakExpr = quotient < _div;
                //_remainder = _doubleOfQuotient - div - div;     // 
                //bBreakExpr = quotient.Equals(div) || quotient.Equals(_doubleOfQuotient);

            } while (bBreakExpr);
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRoot:\n{quotient}\nElapsed time: {strElapsed}\n");
#endif
            return quotient;
        }   // SquareRoot

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

            WriteLine($"\nSqrt:\n{y}\nElapsed time: {strElapsed}\n");
#endif
            return y;
        }   // Sqrt

        public static BigInteger SquareRootFloor(BigInteger x)
        {
            if (x.CompareTo(BigInteger.Zero) < 0) {
                throw new ArgumentException("Negative argument.");
            }
            
            // square roots of 0 and 1 are trivial and
            // y == 0 will cause a divide-by-zero exception
            if (x.Equals(BigInteger.Zero) || x.Equals(BigInteger.One)) {
                return x;
            }

            BigInteger TWO = new BigInteger(2L);
            BigInteger y;
            Stopwatch sw = new Stopwatch();

            // starting with y = x / 2 avoids magnitude issues with x squared
            sw.Start();
            for (y = BigInteger.Divide(x, TWO);
                 y.CompareTo(BigInteger.Divide(x, y)) > 0;
                 y = BigInteger.Divide(BigInteger.Add(BigInteger.Divide(x, y), y), TWO));
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootFloor:\n{y})\nElapsed time: {strElapsed}\n");
#endif
            return y;
        }   // SquareRootFloor

        public static BigInteger SquareRootCeil(BigInteger x)
        {
            if (x.CompareTo(BigInteger.Zero) < 0) {
                throw new ArgumentException("Negative argument.");
            }
            
            // square roots of 0 and 1 are trivial and
            // y == 0 will cause a divide-by-zero exception
            if (x == BigInteger.Zero || x == BigInteger.One) {
                return x;
            }

            BigInteger TWO = new BigInteger(2L);
            BigInteger y;
            Stopwatch sw = new Stopwatch();

            // starting with y = x / 2 avoids magnitude issues with x squared
            sw.Start();
            for (y = BigInteger.Divide(x, TWO);
                 y.CompareTo(BigInteger.Divide(x, y)) > 0;
                 y = BigInteger.Divide(BigInteger.Add(BigInteger.Divide(x, y), y), TWO));
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootCeil:\n{y})\nElapsed time: {strElapsed}\n");
#endif

            if (x.CompareTo(BigInteger.Multiply(y, y)) == 0)
            {
                return y;
            }
            else
            {
                return BigInteger.Add(y, BigInteger.One);
            }
        }   // SquareRootCeil

        static double zeta_of_two_double()
        {
            double sum, term;
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
                lock (monitor)
                {
                    sum += local;
                }
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
                i >>= 1;
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
            unity *= unity;                     // normalize to become one

            object monitor = new object();
            Stopwatch sw = new Stopwatch();

            sw.Start();

            Parallel.For<BigInteger>(2, STEPS, () => 2, (i, loop, term) =>
            //for (uint i = 2; i < STEPS; i++)
            {
                BigInteger b = i * exponent;
                b *= b;
                term = unity / b;
                sum += term;
                //Write("{0}\r", i);
                return term;
            }, term => {
                lock (monitor)
                {
                    sum += term;
                }
            });
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nzeta_of_two_bigint:\n{sum})\nElapsed time: {strElapsed}\n");
#endif
            WriteLine();
            //WriteLine("{0}", unity);
            return sum;
        }

        /// <summary>Times the execution of a function and outputs both the elapsed time and the function's result.</summary>
        static dynamic TimeThis<T>(String strFuncName, Func<T> work)
        {
            if (work == null)
            {
                throw new ArgumentNullException(nameof(work));
            }

            Stopwatch sw = Stopwatch.StartNew();
                T result = work();
            sw.Stop();

            Write("--- {0} ---\t", strFuncName);
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", (float)sw.Elapsed.Milliseconds / 1000);
            //WriteLine("{0} ms\n", sw.ElapsedMilliseconds);
            WriteLine(strElapsed);

            return result;
        }

        static dynamic TimeThis<T, BigInteger>(String strFuncName, Func<T, BigInteger> work)
        {
            if (work == null)
            {
                throw new ArgumentNullException(nameof(work));
            }

            var sw = Stopwatch.StartNew();
            ArgIterator b = new ArgIterator();
            Func<T, BigInteger> selector = arg => work(arg);
            var result = selector;
            sw.Stop();
            Write("--- {0} ---\t", strFuncName);
            WriteLine("{0} ms\n", sw.ElapsedMilliseconds);
            return result;
        }

        static void Main(string[] args)
        {
#if DEBUG
            double pi_squared_over_six = zeta_of_two_double();
#else
            double pi_squared_over_six = TimeThis("zeta_of_two_double()", () => zeta_of_two_double());
#endif
            WriteLine("pi²/6: {0}\n\n√(pi²/6): {1}\n\n", pi_squared_over_six, Math.Sqrt(pi_squared_over_six * 6.0d));

#if DEBUG
            double pi_to_the_fourth_over_ninety = zeta_of_four_double();
#else
            double pi_to_the_fourth_over_ninety = TimeThis("zeta_of_four_double()", () => zeta_of_four_double());
#endif
            WriteLine("pi⁴/90: {0}\n\n(pi⁴/90)^¼: {1}\n\n", pi_to_the_fourth_over_ninety, Math.Pow(pi_to_the_fourth_over_ninety * 90.0d, 0.25d));


            BigInteger TWO = normalize(2, DIGITS);

#if DEBUG
            WriteLine("SquareRoot({0}) =\n{1}\n", TWO, SquareRoot(TWO));
            WriteLine("Sqrt({0}) =\n{1}\n", TWO, Sqrt(TWO));
            WriteLine("SquareRootFloor({0}) =\n{1}\n", TWO, SquareRootFloor(TWO));
            WriteLine("SquareRootCeil({0}) =\n{1}\n", TWO, SquareRootCeil(TWO));

            BigInteger BigInt_pi_squared_over_six = zeta_of_two_bigint();
#else
            WriteLine("SquareRoot({0}) =\n{1}\n", TWO, TimeThis("SquareRoot(TWO)", () => SquareRoot(TWO)));
            WriteLine("Sqrt({0}) =\n{1}\n", TWO, TimeThis("Sqrt(TWO)", () => Sqrt(TWO)));
            WriteLine("SquareRootFloor({0}) =\n{1}\n", TWO, TimeThis("SquareRootFloor(TWO)", () => SquareRootFloor(TWO)));
            WriteLine("SquareRootCeil({0}) =\n{1}\n", TWO, TimeThis("SquareRootCeil(TWO)", () => SquareRootCeil(TWO)));

            BigInteger BigInt_pi_squared_over_six = TimeThis("zeta_of_two_bigint()", () => zeta_of_two_bigint());
            BigInteger BigInt_pi = TimeThis("SquareRoot(BigInt_pi_squared_over_six * 6)", () => SquareRoot(BigInt_pi_squared_over_six * 6));
#endif
            WriteLine("BigInt_pi²/6: {0}\n\n√(BigInt_pi²/6): {1}\n\n", BigInt_pi_squared_over_six, BigInt_pi);

            Write("Press Enter: ");
            ReadLine();
        }
    }
}

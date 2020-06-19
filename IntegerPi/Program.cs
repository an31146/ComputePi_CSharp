﻿using System;
using System.Diagnostics;
using System.Numerics;
using System.Threading.Tasks;

using static System.Console;

#pragma warning disable IDE1006,IDE0049,CS0168,IDE0032,IDE0039,CS0219
/*
 * IDE1006 Suppress Naming Rule Violation IDE1006
 * IDE0049 Name Can Be Simplified
 * IDE0032 Use auto property
 * IDE0039 Use local function
 * CS0168  The variable 'var' is declared but never used.
 * CS0219  The variable 'b' is assigned but its value is never used
 */

namespace IntegerPi
{
    class PiFunctions
    {
        uint m_STEPS;
        uint m_DIGITS;
        public BigInteger ONE => normalize(1, m_DIGITS);
        public BigInteger TWO => normalize(2, m_DIGITS);

        public uint STEPS
        {
            get { return m_STEPS; }
            set { m_STEPS = value; }
        }

        public uint DIGITS
        {
            get { return m_DIGITS; }
            set { m_DIGITS = value; }
        }

        public PiFunctions(uint STEPS = 10000000, uint DIGITS = 53)
        {
            m_STEPS = STEPS;
            m_DIGITS = DIGITS;
        }

        /*
         * Sum series for exp(1)-1 between a, b, returning the result
         *  as an exact fraction(p, q).
         */
        private Tuple<BigInteger, BigInteger> BinarySplit(BigInteger a, BigInteger b)
        {
            if (b - a == 1)
                return new Tuple<BigInteger, BigInteger>(BigInteger.One, b);
            BigInteger m = (a + b) >> 1;
            Tuple<BigInteger, BigInteger> p = BinarySplit(a, m);
            Tuple< BigInteger, BigInteger> q = BinarySplit(m, b);
            
            return new Tuple<BigInteger, BigInteger>(p.Item1 * q.Item2 + q.Item1, p.Item2 * q.Item2);
        }

        // following adapted from https://bit.ly/3hC5fTr
        public BigInteger ExpFixed(uint digits)
        {
            uint prec = (uint)(digits / Math.Log(10, 2));
            int N = (int)(2.95 * digits / Math.Log(digits) + 35);
            Tuple<BigInteger, BigInteger> p = BinarySplit(0, N);

            //int len = Math.Max(p.Item1.ToString().Length, p.Item2.ToString().Length);
            return (p.Item1 + p.Item2) * normalize(1, digits >> 1) / p.Item2; 
        }

        public BigInteger LN2Fixed(uint digits)
        {
            uint N = (uint)(digits * 1.1); 
            BigInteger one_third = normalize(1, digits + 4 >> 1) / 3;
            BigInteger sum = one_third;

            for (uint n = 1; n < N; n++)    // How many iterations do we need to obtain desired precision?
            {
                one_third /= 9;
                sum += one_third / (2 * n + 1);
            }
            sum /= 1000;                    // truncate last few digits due to rounding
            return sum << 1;
        }

        // Assumes parameters are already normalized fixed-point to {m_DIGITS} places
        public BigInteger BigIntAGM(BigInteger a, BigInteger g)
        {
            BigInteger a1, g1;
            while (a != g)
            {
                a1 = (a + g) >> 1;
                g1 = Sqrt(a * g);
                a = a1; g = g1;
            }
            return a;
        }

        public BigInteger BigIntLogN(BigInteger n, uint digits)
        {
            // https://docs.microsoft.com/en-us/dotnet/csharp/language-reference/operators/
            uint N = (uint)(digits * 1.1);
            BigInteger ONE = normalize(1, digits + 4 >> 1);
            BigInteger N1 = ONE / 10240;
            int nLength = n.ToString().Length;

            BigInteger next = BigIntAGM(normalize(1, N1), N1);
            return next;
        }

        public BigInteger SquareRoot(BigInteger n)
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

        public BigInteger Sqrt(BigInteger x)
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

        public BigInteger SquareRootFloor(BigInteger x)
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
            //for (y = BigInteger.Divide(x, TWO);
            //     y.CompareTo(BigInteger.Divide(x, y)) > 0;
            //     y = BigInteger.Divide(BigInteger.Add(BigInteger.Divide(x, y), y), TWO)) ;
            for (y = x >> 1;
                 y.CompareTo(BigInteger.Divide(x, y)) > 0;
                 y = BigInteger.Add(BigInteger.Divide(x, y), y) >> 1) ;
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootFloor:\n{y}\nElapsed time: {strElapsed}\n");
#endif
            return y;
        }   // SquareRootFloor

        public BigInteger SquareRootCeil(BigInteger x)
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
                 y = BigInteger.Divide(BigInteger.Add(BigInteger.Divide(x, y), y), TWO)) ;
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootCeil:\n{y}\nElapsed time: {strElapsed}\n");
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

        // http://rosettacode.org/wiki/Integer_roots#C.23
        public BigInteger NthRoot(BigInteger @base, int n)
        {
            if (@base < 0 || n <= 0)
            {
                throw new ArgumentException();
            }

            int n1 = n - 1;
            BigInteger n2 = n;
            BigInteger n3 = n1;
            BigInteger c = 1;
            BigInteger d = (n3 + @base) / n2;
            BigInteger e = ((n3 * d) + (@base / BigInteger.Pow(d, n1))) / n2;
            while (c != d && c != e)
            {
                c = d;
                d = e;
                e = (n3 * e + @base / BigInteger.Pow(e, n1)) / n2;
            }
            if (d < e)
            {
                return d;
            }
            return e;
        }
        
        public BigInteger Factorial(uint n)
        {
            BigInteger ONE = new BigInteger(1L);
            BigInteger factorial = ONE;
            Stopwatch sw = new Stopwatch();

            sw.Start();
            if (n <= 1)
                return ONE;
            for (uint i = 2; i <= n; i++)
                factorial *= i;
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            //WriteLine($"\nFactorial({n}):\n{factorial}\nElapsed time: {strElapsed}\n");
#endif
            return factorial;
        }

        public BigInteger Ramanujan(uint n)
        {
            BigInteger ONE = normalize(1, DIGITS);
            BigInteger TWO = normalize(2, DIGITS);
            BigInteger sqrtTwo = SquareRoot(TWO);
            BigInteger sum = 0;

            // Can this loop be parallelized?
            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (uint i = 0; i <= n; i++)                           // each term adds ~10 d.p.
            {
                BigInteger term = ONE * (1103 + (26390 * i));
                term *= Factorial(4 * i);
                term /= BigInteger.Pow(Factorial(i), 4);
                term /= BigInteger.Pow(396, 4 * (int)i);
                sum += term;
            }
            BigInteger mult = TWO * sqrtTwo / normalize(9801, DIGITS);
            sum *= mult;
            //BigInteger reciprocal = normalize(1, (uint)sum.ToString().Length + 1) / sum;
            BigInteger reciprocal = ONE * ONE /  sum;
            sw.Stop();

            //WriteLine("ONE * ONE: {0}\n", ONE * ONE);
            WriteLine($"Sum of terms:\n{sum}\n");
            WriteLine($"1 / π:\n{reciprocal}\n");
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nRamanujan Multiplier:\n{mult}\nElapsed time: {strElapsed}\n");
#endif
            return reciprocal;
        }

        public double zeta_of_two_double()
        {
            double sum, term;
            const double unity = 1.0;

            sum = unity;
            //object monitor = new object();

            for (uint i = 2; i < STEPS; i++)
            //Parallel.For<double>(2, STEPS, () => 2, (i, state, local) =>
            {
                double b = (double)i * i;           // b = i²
                double local = unity / b;                  // local = 1 / i²
                sum += local;                       // sum = Ʃ (1 / i²)
                //Write("{0} \r", local.ToString("F18", CultureInfo.InvariantCulture));
            }

            //return (local);
            //}, local => {
            //    lock (monitor)
            //    {
            //        sum += local;
            //    }
            //});
            return sum;
        }

        public double zeta_of_four_double()
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

        public BigInteger normalize_old(uint n, uint digits)
        {
            BigInteger ret_val = new BigInteger(n);
            BigInteger exponent = 100;

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

        // 17/06/20 - new normalize function to make # of decimal places to {digits}
        public BigInteger normalize(uint n, uint digits)
        {
            _ = new BigInteger(n);
            string str_n = n.ToString() + new String('0', (int)digits * 2);
            return BigInteger.Parse(str_n);
        }

        // return normalized argument A as the same as B
        public BigInteger normalize(BigInteger A, BigInteger B)
        {
            int strLenA = A.ToString().Length;
            int strLenB = B.ToString().Length;
            int newlen = Math.Min(strLenA, strLenB);
            string strA = A.ToString();

            if (strLenA > newlen)
                return BigInteger.Parse(strA.Substring(0, newlen));
            else
                return BigInteger.Parse(strA + new String('0', strLenB - strLenA));
        }

        public BigInteger zeta_of_two_bigint()
        {
            BigInteger unity = normalize(1, DIGITS);
            BigInteger exponent = unity;
            unity *= unity;
            BigInteger sum = unity;
            unity *= unity;                     // normalize to become one

            object monitor = new object();
            Stopwatch sw = new Stopwatch();

            WriteLine("Calculating zeta_of_two_bigint() with {0} iterations and {1} digits.", STEPS, DIGITS);

            sw.Start();

            //Parallel.For<BigInteger>(2, STEPS, () => 2, (i, loop, term) =>
            for (uint i = 2; i < STEPS; i++)
            {
                BigInteger b = i * exponent;
                b *= b;
                BigInteger term = unity / b;
                sum += term;
                //Write("{0:F1}%\r", 100.0d * i / STEPS);
                //return term;
            //}, term => {
            //    lock (monitor) { sum += term; }
            } //);
            WriteLine();
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nzeta_of_two_bigint:\n{sum}\nElapsed time: {strElapsed}\n");
#endif
            WriteLine();
            //WriteLine("{0}", unity);
            return sum;
        }

        public BigInteger zeta_of_four_bigint()
        {
            BigInteger unity = normalize(1, DIGITS);
            BigInteger exponent = unity;
            unity *= unity;
            unity *= unity;
            BigInteger sum = unity;
            unity *= unity;                     // normalize to become one

            object monitor = new object();
            Stopwatch sw = new Stopwatch();

            WriteLine("Calculating zeta_of_four_bigint() with {0} iterations and {1} digits.", STEPS, DIGITS);

            sw.Start();

            Parallel.For<BigInteger>(2, STEPS, () => 2, (i, loop, term) =>
            //for (uint i = 2; i < STEPS; i++)
            {
                BigInteger b = i * exponent;
                b *= b;     // b^2
                b *= b;     // b^4
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

            WriteLine($"\nzeta_of_four_bigint:\n{sum}\nElapsed time: {strElapsed}\n");
#endif
            WriteLine();
            //WriteLine("{0}", unity);
            return sum;
        }

        /// <summary>Times the execution of a function and outputs both the elapsed time and the function's result.</summary>
        public dynamic TimeThis<T>(String strFuncName, Func<T> work)
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
                strElapsed = String.Format("{0:F1} s", (float)sw.Elapsed.Seconds);
            WriteLine(strElapsed);

            return result;
        }

        public dynamic TimeThis<T, BigInteger>(String strFuncName, Func<T, BigInteger> work)
        {
            if (work == null)
            {
                throw new ArgumentNullException(nameof(work));
            }

            var sw = Stopwatch.StartNew();
            ArgIterator b = new ArgIterator();                          // what is this line for?
            Func<T, BigInteger> selector = arg => work(arg);
            var result = selector;
            sw.Stop();
            Write("--- {0} ---\t", strFuncName);
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", (float)sw.Elapsed.Seconds);
            WriteLine(strElapsed);
            return result;
        }
    }   // class PiFunctions

    class Program
    { 
        static void Main(string[] args)
        {
            PiFunctions pf = new PiFunctions();

            if (args.Length >= 1)
            {
                pf.STEPS = UInt32.Parse(args[0]);
            }
            if (args.Length == 2)
            {
                pf.DIGITS = UInt32.Parse(args[1]);
            }
#if FLOAT
#if DEBUG
            double pi_squared_over_six = pf.zeta_of_two_double();
#else
            double pi_squared_over_six = pf.TimeThis("zeta_of_two_double()", () => pf.zeta_of_two_double());
#endif
            WriteLine("π²/6: {0}\n\n√(π²): {1}\n\n", pi_squared_over_six, Math.Sqrt(pi_squared_over_six * 6.0d));

#if DEBUG
            double pi_to_the_fourth_over_ninety = pf.zeta_of_four_double();
#else
            double pi_to_the_fourth_over_ninety = pf.TimeThis("zeta_of_four_double()", () => pf.zeta_of_four_double());
#endif
            WriteLine("π⁴/90: {0}\n\n(π⁴)^¼: {1}\n\n", pi_to_the_fourth_over_ninety, Math.Pow(pi_to_the_fourth_over_ninety * 90.0d, 0.25d));
#endif

#if DEBUG
#if SQRT
            WriteLine("SquareRoot({0}) =\n{1}\n", pf.TWO, pf.SquareRoot(pf.TWO));
            WriteLine("Sqrt({0}) =\n{1}\n", pf.TWO, pf.Sqrt(pf.TWO));
            WriteLine("SquareRootFloor({0}) =\n{1}\n", pf.TWO, pf.SquareRootFloor(pf.TWO));
            WriteLine("SquareRootCeil({0}) =\n{1}\n", pf.TWO, pf.SquareRootCeil(pf.TWO));
#endif
#if EXP
            WriteLine("ExpFixed =\n{0}\n", pf.ExpFixed(pf.DIGITS));
#endif
#if LN
            WriteLine("LN2Fixed =\n{0}\n", pf.LN2Fixed(pf.DIGITS));
            WriteLine("BigIntLogN(2) =\n{0}\n", pf.BigIntLogN(1024, pf.DIGITS));
            BigInteger sqrt2 = pf.SquareRoot(pf.TWO);
            BigInteger biAGM = pf.BigIntAGM(pf.normalize(pf.ONE, sqrt2), sqrt2);
            WriteLine("BigIntAGM(1, √2) =\n{0}\n", biAGM);
            WriteLine("Gauss's constant =\n{0}\n", pf.ONE / biAGM);
#endif
#if FACT
            WriteLine("Factorial({0}) = \n{1}\n", 100, pf.Factorial(100));
#endif
#if RAMANUJAN
            // 12 terms accurate to 103 d.p.
            // 6 req. for default prec to 53 d.p.
            WriteLine("Ramanujan series π:\n {0}\n\n", pf.Ramanujan(pf.STEPS));         // 44 terms accurate to 308 d.p. + rounding error
                                                                                        // 49 terms accurate to 397 d.p.
                                                                                        // 125 terms  accurate to 1002 d.p.
#endif
#if ZETA
            BigInteger BigInt_pi_squared_over_six = pf.zeta_of_two_bigint();
            BigInteger BigInt_pi = pf.SquareRoot(BigInt_pi_squared_over_six * 6);
            WriteLine("BigInt_pi²/6: {0}\n\n⁴√(BigInt_pi²*6): {1}\n\n", BigInt_pi_squared_over_six, BigInt_pi);

            BigInteger BigInt_pi_to_fourth_over_ninety = pf.zeta_of_four_bigint();
            BigInt_pi = pf.SquareRoot(pf.SquareRoot(BigInt_pi_to_fourth_over_ninety * 90));
            WriteLine("BigInt_pi⁴/90: {0}\n\n⁴√(BigInt_pi⁴*90): {1}\n\n", BigInt_pi_to_fourth_over_ninety, BigInt_pi);
#endif
#else
#if SQRT
            WriteLine("SquareRoot({0}) =\n{1}\n", pf.TWO, pf.TimeThis("SquareRoot(pf.TWO)", () => pf.SquareRoot(pf.TWO)));
            WriteLine("Sqrt({0}) =\n{1}\n", pf.TWO, pf.TimeThis("Sqrt(pf.TWO)", () => pf.Sqrt(pf.TWO)));
            WriteLine("SquareRootFloor({0}) =\n{1}\n", pf.TWO, pf.TimeThis("SquareRootFloor(pf.TWO)", () => pf.SquareRootFloor(pf.TWO)));
            WriteLine("SquareRootCeil({0}) =\n{1}\n", pf.TWO, pf.TimeThis("SquareRootCeil(pf.TWO)", () => pf.SquareRootCeil(pf.TWO)));
#endif
#if FACT
            WriteLine("Factorial({0}) = \n{1}\n", 100, pf.TimeThis("Factorial(100)", () => pf.Factorial(100)));
#endif

#if RAMANUJAN
            WriteLine("Ramanujan series π:\n {0}\n\n", pf.TimeThis("Ramanujan()", () => pf.Ramanujan(pf.STEPS)));
#endif

#endif
#if !DEBUG && ZETA
            BigInteger BigInt_pi_squared_over_six = pf.TimeThis("zeta_of_two_bigint()", () => pf.zeta_of_two_bigint());
            BigInteger BigInt_pi = pf.TimeThis("SquareRoot(BigInt_pi_squared_over_six * 6)", () => pf.SquareRoot(BigInt_pi_squared_over_six * 6));
            WriteLine("BigInt_pi²/6: {0}\n\n√(BigInt_pi²): {1}\n\n", BigInt_pi_squared_over_six, BigInt_pi);

            BigInteger BigInt_pi_to_fourth_over_ninety = pf.TimeThis("zeta_of_four_bigint()", () => pf.zeta_of_four_bigint());
            BigInt_pi = pf.TimeThis( "Root4(BigInt_pi_to_fourth_over_ninety)", () =>
                                     pf.SquareRoot(pf.SquareRoot(BigInt_pi_to_fourth_over_ninety * 90)) );
#endif
            Write("Press Enter: ");
            ReadLine();
        }
    }   // class Program
}   // namespace

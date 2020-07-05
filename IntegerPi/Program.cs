using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Threading;
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
        uint m_STEPS, m_DIGITS;
        BigInteger m_ONE, m_TWO;
        public BigInteger ONE => m_ONE;         // double digits precision
        public BigInteger _ONE;
        public BigInteger TWO => m_TWO;

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

        public PiFunctions(uint STEPS = 1000000, uint DIGITS = 53)
        {
            m_STEPS = STEPS;
            m_DIGITS = DIGITS;
            m_ONE = normalize(1, m_DIGITS);
            m_TWO = normalize(2, m_DIGITS);
            _ONE = BigInteger.Pow(10, (int)m_DIGITS);
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
        public BigInteger ExpFixed()
        {
            uint digits = m_DIGITS;
            uint prec = (uint)(digits / Math.Log(10, 2));
            int N = (int)(2.95 * digits / Math.Log(digits) + 35);
            Tuple<BigInteger, BigInteger> p = BinarySplit(0, N);

            //int len = Math.Max(p.Item1.ToString().Length, p.Item2.ToString().Length);
            return (p.Item1 + p.Item2) * _ONE / p.Item2; 
        }

        public BigInteger LN2Fixed()
        {
            //uint iterations = (uint)(m_DIGITS * 1.1); 
            BigInteger one_third = _ONE / 3;
            BigInteger sum = one_third, term = sum;

            for (uint n = 1; !term.IsZero; n++)    // How many iterations do we need to obtain desired precision?
            {
                one_third /= 9;
                term = one_third / (2 * n + 1);
                if (term.IsZero)
                    break;
                sum += term;
            }
            //sum = normalize(sum, _ONE / 10);         // truncate last few digits due to rounding
            return (sum << 1);                       // multiply by 2
        }

        // Assumes parameters are already normalized fixed-point to {m_DIGITS} places
        public BigInteger BigIntAGM(BigInteger a, BigInteger g)
        {
            BigInteger a1, g1;
            while (a != g)
            {
                a1 = (a + g) >> 1;
                g1 = SquareRootCeil(a * g);
                a = a1; g = g1;
            }
            return a;
        }

        public BigInteger BigIntLogN(BigInteger n)
        {
            // https://docs.microsoft.com/en-us/dotnet/csharp/language-reference/operators/
            int N = (int)(m_DIGITS * Math.Log(10, 2)) >> 1;     // if using desired number of digits by not shifting, N1 will = 0
            BigInteger N1 = 4 * _ONE / (n << N);                // but if using double the precision of ONE, PI division / agm loses precision
            BigInteger agm = BigIntAGM(_ONE, N1) << 1;

            BigInteger PI = BBPFixed();             // substitute with higher precision to calculate without halving m_DIGITS
                //RamanujanFixed();
            PI *= _ONE;                             // compensate precision for division by agm

            BigInteger LogN = PI / agm;             // this division causes loss of precision by half ~{m_DIGITS/2)
            BigInteger ln2 = LN2Fixed() * N;
            LogN -= ln2;
            LogN /= BigInteger.Pow(10, (int)m_DIGITS >> 1);    // truncate trailing ? digits due to rounding error
            return LogN;
        }

        public BigInteger SquareRoot(BigInteger n)
        {
            int LogBase2 = (int)BigInteger.Log(n, 2);
            BigInteger div = n >> (int)(LogBase2 >> 1);
            BigInteger _div, quotient, _doubleOfQuotient, _remainder;
            bool bBreakExpr = false;
            int iterations = 0;
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
                iterations++;
            } while (bBreakExpr);
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRoot() {iterations} iterations took: {strElapsed}\n");
#endif
            return quotient;
        }   // SquareRoot

        public BigInteger Sqrt(BigInteger x)
        {
            int LogBase2 = (int)BigInteger.Log(x, 2);
            BigInteger div = x >> (int)(LogBase2 >> 1);     // Guess by halving initial value's number of bits
            //BigInteger div = BigInteger.One << 1;         // Old intial guess value
            BigInteger y = x;
            int iterations = 0;
            Stopwatch sw = new Stopwatch();

            // Loop until we hit the same value twice in a row, or wind
            // up alternating.
            sw.Start();
            bool bBreak = false;
            while (!bBreak)
            {
                y = BigInteger.Add(div, BigInteger.Divide(x, div));
                y >>= 1;
                bBreak = y.Equals(div);
                //    break;
                div = y;
                iterations++;
            } 
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSqrt() {iterations} iterations took: {strElapsed}\n");
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

            BigInteger y;
            int iterations = 0;
            Stopwatch sw = new Stopwatch();

            // starting with y = x / 2 avoids magnitude issues with x squared
            sw.Start();
            for (y = x >> ((int)BigInteger.Log(x, 2) >> 1);
                 BigInteger.Multiply(y, y).CompareTo(x) > 0;
                 iterations++)
            {
                y = BigInteger.Add(BigInteger.Divide(x, y), y) >> 1;
            }
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds < 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootFloor() {iterations} iterations took: {strElapsed}\n");
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

            BigInteger y;
            int iterations = 0;
            Stopwatch sw = new Stopwatch();
            var _TWO = BigInteger.One * 2;

            // starting with y = x / 2 avoids magnitude issues with x squared
            sw.Start();
            for (y = BigInteger.Divide(x, x >> ((int)BigInteger.Log(x, 2) >> 1) + 1);
                 BigInteger.Multiply(y, y).CompareTo(x) > 0;
                 iterations++)
            {
                y = BigInteger.Divide(BigInteger.Add(BigInteger.Divide(x, y), y), _TWO);
            }
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds < 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nSquareRootCeil() {iterations} iterations took: {strElapsed}\n");
#endif

            if (x.CompareTo(BigInteger.Multiply(y, y)) == 0)
                return y;       // perfect square
            else
                return BigInteger.Add(y, BigInteger.One);

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
            BigInteger fact = 1;
            Stopwatch sw = new Stopwatch();

            sw.Start();
            if (n <= 1)
                return 1;
            for (uint i = 2; i <= n; i++)
                fact *= i;
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            //WriteLine($"\nFactorial({n}):\n{factorial}\nElapsed time: {strElapsed}\n");
#endif
            return fact;
        }

        public BigInteger RamanujanFixed()
        {
            BigInteger term1 = _ONE * 1103, a2 = _ONE * 26390, a3 = a2;
            BigInteger sum = term1, term = term1, fact = 1, fact4 = 1;
            //int j = 1, k = 1103;

            // Can this loop be parallelized?  No, not if the successive term now depends on the current one
            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 1; !term.IsZero; i++)            // each term adds ~10 d.p. decreases with higher digits
            {
                for (int j = (i - 1) * 4 + 1; j <= 4 * i; j++)
                    fact4 *= j;                             // (4*i)!
                term = fact4 * (term1 + a3);                // (4*i)! * (1103 + 26390*i)
                a3 += a2;                                   // a3 = 26390 * i
                fact *= i;                                  // i!
                term /= BigInteger.Pow(fact, 4);            // (i!)^4
                term /= BigInteger.Pow(396, 4 * i);         // 394^(4*i)
                sum += term;
            }
            BigInteger one_over_pi = 2 * SquareRootCeil(TWO) / 9801 * sum / _ONE;
            sw.Stop();

#if DEBUG
            //WriteLine($"Sum of terms:\n{sum}\n");
            WriteLine($"1 / π:\n{one_over_pi}\n");
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nElapsed time: {strElapsed}\n");
#endif
            return ONE / one_over_pi;           // Truncate last few digits for accuracy?
        }

#if BBPv2
        public BigInteger BBPFixed()
        {
            BigInteger top = _ONE * 47, bottom = 15;
            BigInteger sum = top / bottom;
            Stopwatch sw = new Stopwatch();
            bool bQuit = false;
            sw.Start();
            long i = 1;
            for (; !bQuit; i++)      // each term adds ~10 d.p. repeat until nothing to add, i.e. zero
            {
                top = 47; bottom = 15;
                long k = i * i;
                top += 120 * k + 151 * i;
                top *= _ONE;
                bottom += 712 * k + 194 * i;
                bottom += new BigInteger(k << 9) * k + new BigInteger(k << 10) * i;     // <!-- overflows for 32-bit int and 64-bit long -->
                
                BigInteger term = top / bottom;
                term >>= (int)(i << 2);
                bQuit = term.IsZero;
                sum += term;
            }
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"Sum of terms:\n{sum}\n");
            WriteLine($"\n{i - 1} iterations took: {strElapsed}\n");
#endif
            return sum;
        }
#else
        public BigInteger BBPFixed()
        {
            BigInteger sum = 0, term = 1;
            BigInteger _FOUR = _ONE << 2;
            BigInteger _TWO = _ONE << 1;
            // Can this loop be parallelized?
            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; !term.IsZero; i++)      // each term adds ~10 d.p. repeat until nothing to add, i.e. zero
            {
                int j = 8 * i;
                term = _FOUR / (j + 1);
                term -= _TWO / (j + 4);
                term -= _ONE / (j + 5);
                term -= _ONE / (j + 6);

                term >>= (i << 2);
                sum += term;
            }
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            //WriteLine($"Sum of terms:\n{sum}\n");
            WriteLine($"\nElapsed time: {strElapsed}\n");
#endif
            return sum;
        }
#endif
        public double SlowHarmonicSeries(long N)
        {
            double sum = 0, term;
            // Can this loop be parallelized?  Yes, but it doesn't work very well.
            Stopwatch sw = new Stopwatch();
            object Monitor = new object();
            sw.Start();
            /*
            Parallel.For<double>(1, 1000000000L, () => 1, (i, state, local) =>
            {
                local += 1.0 / i;
                return local;
            }, local => { lock (Monitor) sum += local; });
            */
            for (long i = 1; i < N; i++)
                sum += 1.0 / i;
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nElapsed time: {strElapsed}\n");
            WriteLine($"Sum of terms:\n{sum}\n");
#endif
            return sum;
        }

        public double Euler()
        {
            const long n = 4000000000L;
            return SlowHarmonicSeries(n) - Math.Log(n) - 1.0 / (n >> 1) + 1.0 / n / (n * 12.0);
        }

        public BigInteger SlowHarmonicSeries()
        {
            BigInteger sum = 0;

            for (BigInteger n = 1; n <= BigInteger.Pow(10, (int)m_DIGITS); n++)
                sum += _ONE / n;

            return sum;
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

        private BigInteger normalize_old(uint n, uint digits)
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
        private BigInteger normalize(uint n, uint digits)
        {
            return BigInteger.Pow(10, (int)digits * 2) * n;
            //string str_n = n.ToString() + new String('0', (int)digits * 2);
            //return BigInteger.Parse(str_n);               // expensive operation!
        }

        // return normalized argument A as the same as B
        private BigInteger normalize(BigInteger A, BigInteger B)
        {
            int Log10_A = (int)BigInteger.Log10(A);
            int Log10_B = (int)BigInteger.Log10(B);
            int newlen = Math.Min(Log10_A, Log10_B);

            if (Log10_A > newlen)
                return A / BigInteger.Pow(10, Log10_A - newlen);
            else
                return A * BigInteger.Pow(10, Log10_B - Log10_A);
        }

        public BigInteger BigIntZeta(uint s)    // must be +ve integer > 1
        {
            BigInteger unity = BigInteger.Pow(ONE, (int)(s >> 1));             // normalize to become 1 in the numerator 
            double sum = 1.0d;
            BigInteger zeta = unity;

            object monitor = new object();
            Stopwatch sw = new Stopwatch();

            WriteLine("Calculating BigIntZeta({2}) with {0} iterations and {1} digits.", STEPS, DIGITS, s);

            sw.Start();
            Parallel.For<Tuple<double, BigInteger>>(2, (int)STEPS, () => new Tuple<double, BigInteger>(2.0d, 2), (i, loop, T) =>
            //for (uint i = 2; i < STEPS; i++)
            {
                BigInteger div = BigInteger.Pow(i, (int)s);
                //zeta += unity / j;
                sum += 1.0d / Math.Pow(i, s);
                ValueTuple<double, BigInteger> X = T.ToValueTuple();
                X.Item1 += i;
                X.Item2 += unity / div;
                return X.ToTuple();
            }, t => {
                lock (monitor)
                {
                    sum += t.Item1;
                    zeta += t.Item2;
                }
            });
            sw.Stop();
#if DEBUG
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", sw.Elapsed.TotalSeconds);

            WriteLine($"\nBigIntZeta({s}):\n{zeta}\nElapsed time: {strElapsed}\n");
#endif
            WriteLine();
            return zeta;
        }

        public BigInteger zeta_of_four_bigint()
        {
            BigInteger unity = ONE * ONE;
            BigInteger sum = unity;

            object monitor = new object();
            Stopwatch sw = new Stopwatch();

            WriteLine("Calculating zeta_of_four_bigint() with {0} iterations and {1} digits.", STEPS, DIGITS);

            sw.Start();
            Parallel.For<BigInteger>(2, STEPS, () => 2, (i, loop, term) =>
            {
                BigInteger b = i;
                b = BigInteger.Pow(b, 4);     // b^4
                term += unity / b;
                //sum += term;
                return term;
            }, term => {
                lock (monitor) { 
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
            return sum;
        }

        /// <summary>Times the execution of a function and outputs both the elapsed time and the function's result.</summary>
        public dynamic TimeThis<T>(String strFuncName, Func<T> work)
        {
            if (work == null)
                throw new ArgumentNullException(nameof(work));

            Stopwatch sw = Stopwatch.StartNew();
            T result = work();
            sw.Stop();

            Write("--- {0} ---\t", strFuncName);
            string strElapsed;
            if (sw.ElapsedMilliseconds <= 1000)
                strElapsed = String.Format("{0} ms", sw.ElapsedMilliseconds);
            else
                strElapsed = String.Format("{0:F1} s", (float)sw.Elapsed.TotalSeconds);
            WriteLine(strElapsed);

            return result;
        }

        public dynamic TimeThis<T, BigInteger>(String strFuncName, Func<T, BigInteger> work)
        {
            if (work == null)
                throw new ArgumentNullException(nameof(work));

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
                strElapsed = String.Format("{0:F1} s", (float)sw.Elapsed.TotalSeconds);
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

                if (args.Length == 2)
                {
                    pf.DIGITS = UInt32.Parse(args[1]);
                }
                pf = new PiFunctions(pf.STEPS, pf.DIGITS);
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
            WriteLine("SquareRoot({0}) =\n{1}\n", pf.TWO.GetHashCode(), pf.SquareRoot(pf.TWO));
            WriteLine("Sqrt({0}) =\n{1}\n", pf.TWO.GetHashCode(), pf.Sqrt(pf.TWO));
            WriteLine("SquareRootFloor({0}) =\n{1}\n", pf.TWO.GetHashCode(), pf.SquareRootFloor(pf.TWO));
            WriteLine("SquareRootCeil({0}) =\n{1}\n", pf.TWO.GetHashCode(), pf.SquareRootCeil(pf.TWO));

            var sw = new Stopwatch();
            sw.Start();
            BigInteger sqrt2 = 0;
            for (int i = 1; i <= 1000; i++)
                sqrt2 = pf.SquareRoot(pf.ONE * i);
            sw.Stop();
            WriteLine("1000 runs of SquareRoot() took: {0:F1} s", sw.Elapsed.TotalSeconds);

            Write("Press Enter: "); ReadLine();

            sw.Restart();
            for (int i = 1; i <= 1000; i++)
                sqrt2 = pf.Sqrt(pf.ONE * i);
            sw.Stop();
            WriteLine("1000 runs of Sqrt() took: {0:F1} s", sw.Elapsed.TotalSeconds);

#endif
#if EXP || AGM
            WriteLine("ExpFixed =\n{0}\n", pf.ExpFixed());
            WriteLine("LN2Fixed =\n{0}\n", pf.LN2Fixed());

            BigInteger sqrt2 = pf.Sqrt(pf.TWO);
            BigInteger biAGM = pf.BigIntAGM(pf._ONE, sqrt2);
            WriteLine("BigIntAGM(1, √2) =\n{0}\n", biAGM);
            WriteLine("Gauss's constant =\n{0}\n", pf.ONE / biAGM);
#endif

#if LN
            for (long pow_of_ten = 10; pow_of_ten <= 10000; pow_of_ten *= 10)
            {
                WriteLine("BigIntLogN({0}) =\n{1}\n", pow_of_ten, pf.BigIntLogN(pow_of_ten));
            }
#endif
#if EULER
            //WriteLine("SlowHarmonicSeries()", pf.SlowHarmonicSeries());
            WriteLine("Euler-Mascheroni constant =\n{0}\n", pf.Euler());
#endif
#if FACT
            WriteLine("Factorial({0}) = \n{1}\n", 100, pf.Factorial(100));
#endif
#if RAMANUJAN
            /*
             * Timings are WITHOUT attaching the debugger AND breakpoints disabled!
             * 12 terms accurate to 103 d.p.
             * 6 req. for default prec to 53 d.p.
             * 19728 digits ~69.5 s
             */
            string strPI = pf.RamanujanFixed().ToString().Insert(1, ".");
            WriteLine("Ramanujan series π:\n{0}\n\n", strPI);         // 43 terms accurate to 309 d.p. + rounding error
                                                                      // 49 terms accurate to 397 d.p.
                                                                      // 125 terms  accurate to 1002 d.p.
#endif
#if BBP || BBPv2
            // 256 terms accurate to 309 d.p. 
            // 19728 digits ~2.2s for BBPFixed() v1 / ~1.2s for BBPFixed() v2
            // 39467 digits ~8.0s for BBPFixed() v1 / ~4.7s for BBPFixed() v2
            // 315655 digits ~582.6 s for BBPFixed() v1 / ~357.6 s for BBPFixed() v2
            string strPI = pf.BBPFixed().ToString().Insert(1, ".");
            WriteLine("BBP series π:\n{0}\n\n", strPI);
#endif
#if ZETA
            BigInteger BigIntZeta4 = pf.BigIntZeta(4);
            BigInteger BigInt_pi = pf.NthRoot(BigIntZeta4 * 90, 4);
            WriteLine("BigInt_pi⁴/6: {0}\n\n⁴√(BigInt_pi⁴*6): {1}\n\n", BigIntZeta4, BigInt_pi);

            BigInteger BigIntZeta6= pf.BigIntZeta(6);
            BigInt_pi = pf.NthRoot(BigIntZeta6 * 945, 6);
            WriteLine("BigInt_pi⁶/945: {0}\n\n⁶√(BigInt_pi⁶*945): {1}\n\n", BigIntZeta6, BigInt_pi);
#endif
#else
#if SQRT
            WriteLine("SquareRoot(pf.TWO.GetHashCode() = {0:x}) =\n{1}\n", pf.TWO.GetHashCode(), 
                                                    pf.TimeThis("SquareRoot(pf.TWO)", () => pf.SquareRoot(pf.TWO)));
            WriteLine("Sqrt(pf.TWO.GetHashCode() = {0:x}) =\n{1}\n", 
                                                    pf.TWO.GetHashCode(), pf.TimeThis("Sqrt(pf.TWO)", () => pf.Sqrt(pf.TWO)));
            WriteLine("SquareRootFloor(pf.TWO.GetHashCode() = {0:x}) =\n{1}\n", pf.TWO.GetHashCode(), 
                                                    pf.TimeThis("SquareRootFloor(pf.TWO)", () => pf.SquareRootFloor(pf.TWO)));
            WriteLine("SquareRootCeil(pf.TWO.GetHashCode() = {0:x}) =\n{1}\n", 
                                                    pf.TWO.GetHashCode(), 
                                                    pf.TimeThis("SquareRootCeil(pf.TWO)", () => pf.SquareRootCeil(pf.TWO)));
#endif
#if FACT
            WriteLine("Factorial({0}) = \n{1}\n", 100, pf.TimeThis("Factorial(100)", () => pf.Factorial(100)));
#endif


#if ZETA
            BigInteger BigInt_pi_squared_over_six = pf.TimeThis("zeta_of_two_bigint()", () => pf.zeta_of_two_bigint());
            BigInteger BigInt_pi =  pf.TimeThis("SquareRoot(BigInt_pi_squared_over_six * 6)", () => 
                                    pf.SquareRoot(BigInt_pi_squared_over_six * 6));
            WriteLine("BigInt_pi²/6: {0}\n\n√(BigInt_pi²): {1}\n\n", BigInt_pi_squared_over_six, BigInt_pi);

            BigInteger BigInt_pi_to_fourth_over_ninety = pf.TimeThis("zeta_of_four_bigint()", () => pf.zeta_of_four_bigint());
            BigInt_pi = pf.TimeThis( "Root4(BigInt_pi_to_fourth_over_ninety)", () =>
                                     pf.SquareRoot(pf.SquareRoot(BigInt_pi_to_fourth_over_ninety * 90)) );
#endif
#if RAMANUJAN
            // 19728 digits ~75.7 s
            WriteLine("Ramanujan series π:\n {0}\n\n", 
                pf.TimeThis("RamanujanFixed()", () => 
                pf.RamanujanFixed().ToString().Insert(1, ".")));
#endif
#if BBP
            // 256 terms accurate to 309 d.p. 
            WriteLine("BBP series π:\n{0}\n\n", pf.TimeThis("BBPFixed()", () => pf.BBPFixed().ToString().Insert(1, ".")));
#endif
#endif
            Write("Press Enter: ");
            ReadLine();
        }
    }   // class Program
}   // namespace

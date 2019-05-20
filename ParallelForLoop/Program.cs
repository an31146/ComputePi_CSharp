using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

using static System.Console;

namespace ParallelForLoop
{
    class Test
    {
        const int A = -100000000;
        const int B = 100000001;

        static void Main(string[] args)
        {
            int[] nums = Enumerable.Range(A, B-A+1).ToArray();
            long total = 0;
            object lockObj = new object();
            Stopwatch sw = new Stopwatch();

            sw.Start();
            {
                // Use type parameter to make subtotal a long, not an int
                WriteLine("{0} {1}", nums.Length, B-A);
                Parallel.For<long>(0, nums.Length, () => total, (j, loop, subtotal) =>
                    {
                        subtotal += nums[j];
                        //WriteLine("{0}", subtotal);
                        //if (nums[j] > 0)
                        //    loop.Stop();
                        return subtotal;
                    }, (x) => Interlocked.Add(ref total, x)
                );
            }
            sw.Stop();

            WriteLine("Sum of integers from {0:N0} to {1:N0} = {2:N0}\n", A, B, total);
            WriteLine("Time: {0} ms", sw.ElapsedMilliseconds);

            total = 0;
            sw.Restart();
            {
                foreach (int n in nums)
                    total += n;
            }
            sw.Stop();

            WriteLine("Sum of integers from {0:N0} to {1:N0} = {2:N0}\n", A, B, total);
            WriteLine("Time: {0} ms", sw.ElapsedMilliseconds);

            Write("Press ANY key... ");
            ReadKey();
        }
    }
}

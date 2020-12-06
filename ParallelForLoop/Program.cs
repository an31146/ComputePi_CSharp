using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

using static System.Console;

namespace ParallelForLoop
{
    class Test
    {
        const int A = -20000000;
        const int B = 20000001;

        static void Main(string[] args)
        {
            int[] nums = Enumerable.Range(A, B - A + 1).ToArray();
            long total = 0;
            object lockObj = new object();
            Stopwatch sw = new Stopwatch();

            WriteLine("sizeof nums[]: {0:N0} bytes\n", sizeof(int) * nums.Length);
            sw.Start();
            {
                WriteLine($"nums.Length: {nums.Length:N0}\nB - A:       {B - A:N0}\n");
                // Use type parameter to make subtotal a long, not an int
                Parallel.For<long>(0, nums.Length, () => total, (j, loop, subtotal) =>
                    {
                        subtotal += (long)nums[j];
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
                foreach (var n in nums)
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

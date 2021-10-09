using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

using static System.Console;

namespace ParallelForLoop
{
    class Test
    {
        const int A = -41000000;
        const int B = 41000001;

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
                ParallelOptions options = new ParallelOptions { MaxDegreeOfParallelism = 4 };
                Parallel.ForEach(nums, options, () => total, (j, loop, subtotal) =>
                    {
                        lock (lockObj)
                        {
                            //subtotal += nums[j];
                            total += j;
                            //System.Threading.Interlocked.Add(ref total, j);
                            //Write("{0}; ", nums[j]);
                        }
                        return subtotal;
                    }, (x) => { } //_ = total     // variable is corrupted if modified here!
                    //System.Threading.Interlocked.Add(ref total, x)
                );
            }
            sw.Stop();

            WriteLine("Parallel.ForEach: Sum of integers from {0:N0} to {1:N0} = {2:N0}\n", A, B, total);
            WriteLine("Time: {0} ms", sw.ElapsedMilliseconds);

            total = 0;
            sw.Restart();
            {
                foreach (var n in nums)
                    total += n;
            }
            sw.Stop();

            WriteLine("foreach: Sum of integers from {0:N0} to {1:N0} = {2:N0}\n", A, B, total);
            WriteLine("Time: {0} ms", sw.ElapsedMilliseconds);
        }
    }
}

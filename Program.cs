//--------------------------------------------------------------------------
// 
//  Copyright (c) Microsoft Corporation.  All rights reserved. 
// 
//  File: Program.cs
//
//--------------------------------------------------------------------------

using System;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

namespace ComputePi
{
    class Program
    {
        const int NUM_STEPS = 1000000000;
        const int REPS = 10;

        /// <summary>Main method to time various implementations of computing PI.</summary>
        static void Main(string[] args)
        {
            //ParallelForCancellation p = new ParallelForCancellation();
            //p.CancelDemo();

            Time("SerialLinqPi()", () => SerialLinqPi(), REPS);
            Time("ParallelLinqPi()", () => ParallelLinqPi(), REPS);
            Time("SerialPi()", () => SerialPi(), REPS);
            Time("ParallelPi()", () => ParallelPi(), REPS);
            Time("ParallelPartitionerPi()", () => ParallelPartitionerPi(), REPS);

            Console.WriteLine("---- Press Enter ----");
            Console.ReadLine();
        }

        /// <summary>Times the execution of a function and outputs both the elapsed time and the function's result.</summary>
        static void Time<T>(String strFuncName, Func<T> work, int reps)
        {
            dynamic result = 0.0;
            var sw = Stopwatch.StartNew();
            var prev_time = 0L;
            Console.WriteLine("\n--- {0} ---", strFuncName);
            foreach (var i in Enumerable.Range(0, reps))
            {
                result = work();
                var lap_time = sw.ElapsedMilliseconds;
                Console.WriteLine("#{0}: {1,5} ms", i, lap_time - prev_time);
                prev_time = lap_time;
            }
            sw.Stop();
            Console.WriteLine("\nAverage of {2} runs: {0,5} ms\nResult = {1}\n", sw.ElapsedMilliseconds / reps, result, reps);
        }

        /// <summary>Estimates the value of PI using a LINQ-based implementation.</summary>
        static double SerialLinqPi()
        {
            double step = 1.0 / NUM_STEPS;
            return 4.0d *
                (from i in Enumerable.Range(0, NUM_STEPS)
                 let x = (i + 0.5) * step
                 select 1.0 / (1.0 + x * x)).Sum() * step;
        }

        /// <summary>Estimates the value of PI using a PLINQ-based implementation.</summary>
        static double ParallelLinqPi()
        {
            double step = 1.0 / NUM_STEPS;
            return 4.0d * (from i in ParallelEnumerable.Range(0, NUM_STEPS)
                           let x = (i + 0.5) * step
                           select 1.0 / (1.0 + x * x)).Sum() * step;
        }

        /// <summary>Estimates the value of PI using a for loop.</summary>
        static double SerialPi()
        {
            double sum = 0.0;
            double step = 1.0 / NUM_STEPS;
            for (int i = 0; i < NUM_STEPS; i++)
            {
                double x = (i + 0.5) * step;
                sum = sum + 1.0 / (1.0 + x * x);
            }
            return 4.0d * step * sum;
        }

        /// <summary>Estimates the value of PI using a Parallel.For.</summary>
        static double ParallelPi()
        {
            double sum = 0.0;
            double step = 1.0 / NUM_STEPS;
            object monitor = new object();
            Parallel.For(0, NUM_STEPS, () => 0.0, (i, state, local) =>
            {
                double x = (i + 0.5) * step;
                return local + 1.0 / (1.0 + x * x);
            }, local => { lock (monitor) sum += local; });
            return 4.0d * step * sum;
        }

        /// <summary>Estimates the value of PI using a Parallel.ForEach and a range partitioner.</summary>
        static double ParallelPartitionerPi()
        {
            double sum = 0.0;
            double step = 1.0 / NUM_STEPS;
            object monitor = new object();
            //object monitor2 = new object();
            Parallel.ForEach(Partitioner.Create(0, NUM_STEPS), () => 0.0, (range, state, local) =>
            {
                //lock (monitor2) Console.WriteLine($"{range.Item1} {range.Item2}");
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    double x = (i + 0.5) * step;
                    local += 1.0 / (1.0 + x * x);
                }
                return local;
            }, local => { lock (monitor) sum += local; });
            return 4.0d * step * sum;
        }
    }
}
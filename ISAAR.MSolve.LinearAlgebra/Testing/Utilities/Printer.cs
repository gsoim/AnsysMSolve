﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.LinearAlgebra.Testing.Utilities
{
    // TODO: optional parameters that control formating, such as decimals, right alignment
    public class Printer
    {
        public string SectionSeparator { get; set; }

        public Printer(string rowSeparator = "\n", string colSeparator = " ")
        {
            this.SectionSeparator = "************************************************************************************";
        }

        public void Print(int[] array)
        {
            for (int i = 0; i < array.Length; ++i) Console.Write(array[i] + " ");
            Console.WriteLine();
        }

        public void Print(double[] vector)
        {
            (new Array1DWriter(vector)).WriteToConsole();
        }

        public void Print(double[,] matrix)
        {
            (new Array2DWriter(matrix)).WriteToConsole();
        }

        public void PrintFactorizationCholesky(bool isCorrect, double[,] matrix, double[,] uExpected, double[,] uComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine(SectionSeparator);
            Console.WriteLine("The following Cholesky factorization is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("U (expected) = ");
            Print(uExpected);
            Console.WriteLine();
            Console.Write("U (computed) = ");
            Print(uComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintFactorizationLU(bool isCorrect, double[,] matrix, double[,] lExpected, double[,] uExpected,
            double[,] lComputed, double[,] uComputed, bool isSingular)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine(SectionSeparator);
            if (isSingular) Console.Write("The matrix A is singular! ");
            Console.WriteLine("The following LU factorization is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("L (expected) = ");
            Print(lExpected);
            Console.WriteLine();
            Console.Write("L (computed) = ");
            Print(lComputed);
            Console.WriteLine();
            Console.Write("U (expected) = ");
            Print(uExpected);
            Console.WriteLine();
            Console.Write("U (computed) = ");
            Print(uComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintIndefiniteMatrix(double[,] matrix)
        {
            Console.WriteLine(SectionSeparator);
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine("is indefinite. Cannot apply Cholesky factorization.");
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintMatrixEquality(bool isCorrect, double[,] matrixExpected, double[,] matrixComputed)
        {
            string result = isCorrect ? "EQUAL" : "NOT EQUAL";
            Console.WriteLine(SectionSeparator);
            Console.WriteLine("The following matrices are " + result + ":");
            Console.Write("A (expected) = ");
            Print(matrixExpected);
            Console.WriteLine();
            Console.Write("A (computed) = ");
            Print(matrixComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintMatrixVectorMult(bool isCorrect, double[,] matrix, double[] x, double[] bExpected, double[] bComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine(SectionSeparator);
            Console.WriteLine("The following matrix multiplication is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("x = ");
            Print(x);
            Console.WriteLine();
            Console.Write("b (expected) = ");
            Print(bExpected);
            Console.WriteLine();
            Console.Write("b (computed) = ");
            Print(bComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintSingularMatrix(double[,] matrix, bool systemSolveUsecase = true)
        {
            Console.WriteLine(SectionSeparator);
            Console.Write("A = ");
            Print(matrix);
            Console.Write("is singular.");
            if (systemSolveUsecase) Console.Write(" Cannot solve the linear system");
            Console.WriteLine();
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintSystemSolution(bool isCorrect, double[,] matrix, double[] b, double[] xExpected, double[] xComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine(SectionSeparator);
            Console.WriteLine("The following linear system solution is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("b = ");
            Print(b);
            Console.WriteLine();
            Console.Write("x (expected) = ");
            Print(xExpected);
            Console.WriteLine();
            Console.Write("x (computed) = ");
            Print(xComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }

        public void PrintVectorEquality(bool isCorrect, double[] vectorExpected, double[] vectorComputed)
        {
            string result = isCorrect ? "EQUAL" : "NOT EQUAL";
            Console.WriteLine(SectionSeparator);
            Console.WriteLine("The following vectors are " + result + ":");
            Console.Write("v (expected) = ");
            Print(vectorExpected);
            Console.WriteLine();
            Console.Write("v (computed) = ");
            Print(vectorComputed);
            Console.WriteLine(SectionSeparator);
            Console.WriteLine();
        }
    }
}

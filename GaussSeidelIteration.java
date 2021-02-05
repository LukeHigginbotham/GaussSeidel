/*
 * Use The Gauss-Seidel Iteration method to solve a system of linear equations of any size n > 1. 
 *
 * CSC 2262 Gauss-Seidel Method Program Findings
 *
 * @author Luke Higginbotham 899568628
 * @since 11/19/2020 
 */
package gaussseideliteration;

import java.util.Arrays;

public class GaussSeidelIteration 
{

    static double EPSILON = 0.0001;
    static double[][] matrixA = 
        {
            {4, -1, 0, -1, 0, 0},
            {-1, 4, -1, 0, -1, 0},
            {0, -1, 4, 0, 0, -1},
            {-1, 0, 0, 4, -1, 0},
            {0, -1, 0, -1, 4, -1},
            {0, 0, -1, 0, -1, 4}
        };

    static double[] matrixb = {0, 5, 0, 6, -2, 6};

    public static void main(String[] args) 
    {
        GaussSeidelIteration.gaussSeidel();
        System.out.print("\nA = ");
        for (int k = 0; k < matrixA.length; k++) 
        {
            System.out.print("\n");
            for (int i = 0; i < matrixA.length; i++) 
            {
                System.out.print(" | " + matrixA[k][i]);
            }
            System.out.print(" | ");
        }
        System.out.println("\n\nx = | 1.0000 | 2.0000 | 1.0000 | 2.0000 | 1.0000 | 2.0000 |");
        GaussSeidelIteration.getNorm(matrixb);
        System.out.println("\nb = ");
        for (int k = 0; k < matrixA.length; k++) 
        {
            System.out.print(" | " + matrixb[k]);
        }
        System.out.print(" | ");
    }

    static double[] gaussSeidel() 
    {
        System.out.println("   k\t x0\tx1     x2     x3     x4     x5\t  Diff");
        int iteration = 0;
        double norm = 0;
        double diff = 0;
        double[] x = new double[matrixA.length];
        Arrays.fill(x, 0.0);
        
        System.out.printf("%4d %.4f %.4f %.4f %.4f %.4f %.4f\n", iteration, x[0], x[1], x[2], x[3], x[4], x[5], diff);

        while (true) 
        {
            for (int i = 0; i < matrixA.length; i++) 
            {
                double x0 = 0.0;
                for (int k = 0; k < matrixA.length; k++) 
                {
                    if (i != k) 
                    {
                        x0 += matrixA[i][k] * x[k];
                    }
                }
                x[i] = (matrixb[i] - x0) / matrixA[i][i];
            }
            iteration++;
            
            diff = norm - getNorm(x);
            norm = getNorm(x);
            System.out.printf("%4d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n", iteration, x[0], x[1], x[2], x[3], x[4], x[5], diff);
            if (Math.abs(diff) < EPSILON) 
            {
                break;
            }
        }
        return x;
    }

    static double getNorm(double[] x) 
    {
        double norm = 0;
        for (double i : x) 
        {
            norm += Math.pow(i, 2);
        }
        return Math.sqrt(norm);
    }
}

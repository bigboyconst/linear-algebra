using static System.Math;

namespace LinearAlgebra
{
    /// <summary>
    /// Represents a type of cubic spline.
    /// </summary>
    public enum SplineType
    {
        /// <summary>
        /// A closed spline; that is, one that forms a loop.
        /// </summary>
        Closed,
        /// <summary>
        /// A natural spline; that is, one that starts and ends smoothly.
        /// </summary>
        Natural
    }
    /// <summary>
    /// Class representing a cubic spline interpolating a set of points
    /// </summary>
    public static class Spline
    {
        const double _eps = 1e-12;

        /// <summary>
        /// Creates the matrices involved in the system of linear equations associated with a cubic closed spline interpolating the 
        /// values in <paramref name="pValues"/> at the respective value in <paramref name="tValues"/>
        /// </summary>
        public static List<double[,]> ClosedSplineMatrixSystem(List<int> tValues, List<double> pValues)
        {
            int nMoments = tValues.Count;
            int nPoints = pValues.Count;
            if (nMoments != nPoints)
                throw new InvalidLengthException();
            if (Abs(pValues.First() - pValues.Last()) >= _eps)
                throw new InvalidClosedSplineException();

            double[,] splineMatrix = new double[4 * (nMoments - 1), 4 * (nMoments - 1)];
            // left and right values
            for (int j = 0, i = 0; j < nMoments - 1 && i < nMoments - 1; j++, i++)
            {
                splineMatrix[2 * i, 4 * j] = 1;
                splineMatrix[2 * i, 4 * j + 1] = tValues[i];
                splineMatrix[2 * i, 4 * j + 2] = Pow(tValues[i], 2);
                splineMatrix[2 * i, 4 * j + 3] = Pow(tValues[i], 3);
                splineMatrix[2 * i + 1, 4 * j] = 1;
                splineMatrix[2 * i + 1, 4 * j + 1] = tValues[i + 1];
                splineMatrix[2 * i + 1, 4 * j + 2] = Pow(tValues[i + 1], 2);
                splineMatrix[2 * i + 1, 4 * j + 3] = Pow(tValues[i + 1], 3);
            }
            // derivatives at intermediary points
            for (int j = 0, i = 0; j < nMoments - 2 && i < nMoments - 2; i++, j++)
            {
                int k = 2 * (nMoments - 1) + i;
                splineMatrix[k, 4 * j] = 0;
                splineMatrix[k, 4 * j + 1] = 1;
                splineMatrix[k, 4 * j + 2] = 2 * tValues[i + 1];
                splineMatrix[k, 4 * j + 3] = 3 * Pow(tValues[i + 1], 2);
                splineMatrix[k, 4 * j + 4] = 0;
                splineMatrix[k, 4 * j + 5] = -1;
                splineMatrix[k, 4 * j + 6] = -2 * tValues[i + 1];
                splineMatrix[k, 4 * j + 7] = -3 * Pow(tValues[i + 1], 2);
            }
            // second derivatives at intermediary points
            for (int j = 0, i = 0; i < nMoments - 2 && j < nMoments - 2; i++, j++)
            {
                int k = 2 * (nMoments - 1) + (nMoments - 2) + i;
                splineMatrix[k, 4 * j] = 0;
                splineMatrix[k, 4 * j + 1] = 0;
                splineMatrix[k, 4 * j + 2] = 2;
                splineMatrix[k, 4 * j + 3] = 6 * tValues[i + 1];
                splineMatrix[k, 4 * j + 4] = 0;
                splineMatrix[k, 4 * j + 5] = 0;
                splineMatrix[k, 4 * j + 6] = -2;
                splineMatrix[k, 4 * j + 7] = -6 * tValues[i + 1];
            }
            // derivatives at endpoints
            splineMatrix[4 * (nMoments - 1) - 2, 0] = 0;
            splineMatrix[4 * (nMoments - 1) - 2, 1] = 1;
            splineMatrix[4 * (nMoments - 1) - 2, 2] = 2 * tValues[0];
            splineMatrix[4 * (nMoments - 1) - 2, 3] = 3 * Pow(tValues[0], 2);
            splineMatrix[4 * (nMoments - 1) - 2, 4 * (nMoments - 1) - 4] = 0;
            splineMatrix[4 * (nMoments - 1) - 2, 4 * (nMoments - 1) - 3] = -1;
            splineMatrix[4 * (nMoments - 1) - 2, 4 * (nMoments - 1) - 2] = -2 * tValues[nMoments - 1];
            splineMatrix[4 * (nMoments - 1) - 2, 4 * (nMoments - 1) - 1] = -3 * Pow(tValues[nMoments - 1], 2);

            // second derivatives at endpoints
            splineMatrix[4 * (nMoments - 1) - 1, 0] = 0;
            splineMatrix[4 * (nMoments - 1) - 1, 1] = 0;
            splineMatrix[4 * (nMoments - 1) - 1, 2] = 2;
            splineMatrix[4 * (nMoments - 1) - 1, 3] = 6 * tValues[0];
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 4] = 0;
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 3] = 0;
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 2] = -2;
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 1] = -6 * tValues[nMoments - 1];

            // matrix of constants
            double[,] constants = new double[4 * (nMoments - 1), 1];
            constants[0, 0] = pValues[0];
            constants[2 * (nMoments - 1) - 1, 0] = pValues[nMoments - 1];
            for (int i = 0; i < nMoments - 2; i++)
            {
                constants[2 * i + 1, 0] = pValues[i + 1];
                constants[2 * i + 2, 0] = pValues[i + 1];
            }
            return new List<double[,]> { splineMatrix, constants };
        }

        /// <summary>
        /// Creates the matrices involved in the system of linear equations associated with a cubic natural spline interpolating the 
        /// values in <paramref name="pValues"/> at the respective values in <paramref name="tValues"/>
        /// </summary>
        public static List<double[,]> NaturalSplineMatrixSystem(List<int> tValues, List<double> pValues)
        {
            int nMoments = tValues.Count;
            int nPoints = pValues.Count;
            if (nMoments != nPoints)
                throw new InvalidLengthException();

            double[,] splineMatrix = new double[4 * (nMoments - 1), 4 * (nMoments - 1)];
            // left and right values
            for (int j = 0, i = 0; j < nMoments - 1 && i < nMoments - 1; j++, i++)
            {
                splineMatrix[2 * i, 4 * j] = 1;
                splineMatrix[2 * i, 4 * j + 1] = tValues[i];
                splineMatrix[2 * i, 4 * j + 2] = Pow(tValues[i], 2);
                splineMatrix[2 * i, 4 * j + 3] = Pow(tValues[i], 3);
                splineMatrix[2 * i + 1, 4 * j] = 1;
                splineMatrix[2 * i + 1, 4 * j + 1] = tValues[i + 1];
                splineMatrix[2 * i + 1, 4 * j + 2] = Pow(tValues[i + 1], 2);
                splineMatrix[2 * i + 1, 4 * j + 3] = Pow(tValues[i + 1], 3);
            }
            // derivatives at intermediary points
            for (int j = 0, i = 0; j < nMoments - 2 && i < nMoments - 2; i++, j++)
            {
                int k = 2 * (nMoments - 1) + i;
                splineMatrix[k, 4 * j] = 0;
                splineMatrix[k, 4 * j + 1] = 1;
                splineMatrix[k, 4 * j + 2] = 2 * tValues[i + 1];
                splineMatrix[k, 4 * j + 3] = 3 * Pow(tValues[i + 1], 2);
                splineMatrix[k, 4 * j + 4] = 0;
                splineMatrix[k, 4 * j + 5] = -1;
                splineMatrix[k, 4 * j + 6] = -2 * tValues[i + 1];
                splineMatrix[k, 4 * j + 7] = -3 * Pow(tValues[i + 1], 2);
            }
            // second derivatives at intermediary points
            for (int j = 0, i = 0; i < nMoments - 2 && j < nMoments - 2; i++, j++)
            {
                int k = 2 * (nMoments - 1) + (nMoments - 2) + i;
                splineMatrix[k, 4 * j] = 0;
                splineMatrix[k, 4 * j + 1] = 0;
                splineMatrix[k, 4 * j + 2] = 2;
                splineMatrix[k, 4 * j + 3] = 6 * tValues[i + 1];
                splineMatrix[k, 4 * j + 4] = 0;
                splineMatrix[k, 4 * j + 5] = 0;
                splineMatrix[k, 4 * j + 6] = -2;
                splineMatrix[k, 4 * j + 7] = -6 * tValues[i + 1];
            }
            // second derivative at start
            splineMatrix[4 * (nMoments - 1) - 2, 2] = 2;
            splineMatrix[4 * (nMoments - 1) - 2, 3] = 6 * tValues[0];

            // second derivative at end
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 2] = 2;
            splineMatrix[4 * (nMoments - 1) - 1, 4 * (nMoments - 1) - 1] = 6 * tValues[nMoments - 1];

            // matrix of constants
            double[,] constants = new double[4 * (nMoments - 1), 1];
            constants[0, 0] = pValues[0];
            constants[2 * (nMoments - 1) - 1, 0] = pValues[nMoments - 1];
            for (int i = 0; i < nMoments - 2; i++)
            {
                constants[2 * i + 1, 0] = pValues[i + 1];
                constants[2 * i + 2, 0] = pValues[i + 1];
            }
            return new List<double[,]> { splineMatrix, constants };
        }

        /// <summary>
        /// Formats spline coefficients to live in an (n-1)x4 matrix instead of a 4(n-1)x1 matrix, so that they can easily be 
        /// read off as the coefficients for each constituent polynomial
        /// </summary>
        static double[,] FormatSplineCoefficients(double[,] coeffs)
        {
            double[,] result = new double[coeffs.GetLength(0) / 4, 4];
            for (int i = 0; i < coeffs.GetLength(0); i++)
            {
                result[i / 4, i % 4] = coeffs[i, 0];
            }
            return result;
        }

        /// <summary>
        /// Returns a matrix containing the coefficients for a cubic spline interpolating the list <paramref name="pValues"/> at 
        /// <paramref name="tValues"/> of type <paramref name="type"/>.
        /// </summary>
        /// <returns>
        /// An (n-1)x4 matrix where the column i represents the coefficient for the term of degree i and 
        /// each row represents a different cubic polynomial.
        /// </returns>
        public static Matrix CubicSplineInterpolation(List<int> tValues, List<double> pValues, SplineType type)
        {
            var system = type switch
            {
                SplineType.Closed => ClosedSplineMatrixSystem(tValues, pValues),
                SplineType.Natural => NaturalSplineMatrixSystem(tValues, pValues),
                _ => throw new ArgumentException("Something really bad went wrong!")
            };
            Matrix A = system[0];

            Vector b = Vector.ToVector(system[1]);

            return FormatSplineCoefficients(Solver.Solve(A, b));
        }
    }
}

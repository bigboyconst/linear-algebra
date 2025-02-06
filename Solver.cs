namespace LinearAlgebra
{
    /// <summary>
    /// Class containing many utilities related to solving systems of equations
    /// </summary>
    public static class Solver
    {
        /// <summary>
        /// Performs backward substitutions on an upper triangular matrix and a vector of constants.
        /// </summary>
        public static Vector BackwardSubstitute(Matrix upperTriangular, Vector consts)
        {
            int rows = upperTriangular.Rows;
            int cols = upperTriangular.Columns;
            if (rows != consts.Count)
            {
                throw new InvalidRowCountException();
            }
            Vector x = new(cols);
            double s = 0;
            x[cols - 1] = consts[cols - 1] / upperTriangular[rows - 1, cols - 1];
            for (int i = 1; i < rows; i++)
            {
                s = 0;
                for (int k = rows - i; k < cols; k++)
                {
                    s += upperTriangular[rows - i - 1, k] * x[k];
                }
                x[rows - i - 1, 0] = (consts[rows - i - 1] - s) / upperTriangular[rows - i - 1, cols - i - 1];
            }
            return x;
        }
        /// <summary>
        /// Performs forward substitutions on a lower triangular matrix and a vector of constants.
        /// </summary>
        public static Vector ForwardSubstitute(Matrix lowerTriangular, Vector consts)
        {
            if (lowerTriangular.Rows != consts.Count)
                throw new InvalidRowCountException();
            int rows = lowerTriangular.Rows;
            int cols = lowerTriangular.Columns;
            Vector x = new(cols, 1);
            double s = 0;
            x[0, 0] = consts[0] / lowerTriangular[0, 0];
            for (int i = 1; i < rows; i++)
            {
                s = 0;
                for (int k = 0; k < i; k++)
                    s += lowerTriangular[i, k] * x[k];
                x[i, 0] = (consts[i] - s) / lowerTriangular[i, i];
            }
            return x;
        }
        /// <summary>
        /// Solves the matrix equation Ax = b using the method of 
        /// LU decomposition
        /// </summary>
        public static Vector SolveLU(Matrix A, Vector b)
        {
            var PLU_A = A.LUDecomposition();
            var P = PLU_A[0];
            var L = PLU_A[1];
            var U = PLU_A[2];
            var Pinv = PLU_A[3];
            Vector Y, X;
            // If A has an LU decomposition
            if (P == Matrix.Id(A.Rows))
            {
                // LY = B
                Y = ForwardSubstitute(L, b);
                // UX = Y
                X = BackwardSubstitute(U, Y);
            }
            else
            {
                // B -> P^(-1)B
                b = Pinv * b;
                // LY = B
                Y = ForwardSubstitute(L, b);
                // UX = Y
                X = BackwardSubstitute(U, Y);
            }
            return X;
        }
        /// <summary>
        /// Solves the system Ax = b using the null space of A
        /// </summary>
        /// <returns>
        /// A matrix in which the first column is a particular solution 
        /// to the equation and the remaining columns (if the solution 
        /// is not unique) are the basis vectors for the null space.
        /// </returns>
        public static Matrix Solve(Matrix A, Vector b)
        {
            int rankA = A.Rank;
            int rankAb = A.Augment(b).Rank;
            if (rankA < rankAb)
                return new(0, 0); // there is no solution
            var rref = A.ExtendedRREF();
            var bStar = rref[1] * b;
            var link = A.LinkedVariables();
            Vector xStar = new(A.Columns);
            for (int i = 0; i < A.Columns; i++)
            {
                if (link.Contains(i))
                {
                    xStar[i] = bStar[link.IndexOf(i)];
                    continue;
                }
                xStar[i] = 0;
            }
            if (rankA == rankAb && rankA == A.Columns)
                return xStar;
            var ns = A.NullSpace();
            List<Vector> cols = new() { xStar };
            cols.AddRange(ns);
            return Matrix.AugmentList(cols);
        }

        /// <summary>
        /// Takes in a string representing a system of equations, with each equation separated by a semicolon (";") and no multiplication sign between the coefficients
        /// and their respective variables along with an array of strings stating what those variables are.  This function returns the coefficient matrix and constant 
        /// vector for the system and returns them as a tuple.
        /// </summary>
        /// <remarks>
        /// Every variable should be present in the represented equation, even if their coefficient is zero.  The equations should also be in the same order as the variables 
        /// and the first term in each equation has to have an explicitly written coefficient in front (even if it is 1).
        /// </remarks>
        public static (Matrix A, Vector b) ParseSystem(string system, string[] variables)
        {
            // equations are separated by a semicolon
            // no multiplication symbol between coefficients and variables
            // number of equations
            int nEqns = system.Count(c => c == '=');
            // nummber of variables
            int nVars = variables.Length;
            string[] splitStrings = new string[nVars + 1];
            splitStrings[0] = ";";
            for (int i = 0; i < nVars; i++)
                splitStrings[i + 1] = variables[i];

            List<string> coeffArray = system.Split(splitStrings, StringSplitOptions.None).ToList();

            for (int i = 0; i != coeffArray.Count; ++i)
            {
                if (string.IsNullOrWhiteSpace(coeffArray[i]))
                    coeffArray[i] = "0";
                if (coeffArray[i].RemoveSpaces() == "+" || coeffArray[i].RemoveSpaces() == "-")
                    coeffArray[i] += "1";
            }

            string[] constants = coeffArray.Where(s => s.Contains('=')).Select(s => s.Remove(s.IndexOf('='), 1)).ToArray();

            coeffArray.RemoveAll(s => s.Contains('='));

            coeffArray = coeffArray.Select(s => s.RemoveSpaces()).ToList();

            double[,] coefficients = new double[nEqns, nVars];
            for (int i = 0; i < coeffArray.Count; i++)
            {
                if (i / 3 >= nEqns || i % 3 >= nVars)
                    break;
                coefficients[i / 3, i % 3] = double.Parse(coeffArray[i]);
            }
            Vector consts = new(nEqns);
            for (int i = 0; i < constants.Length; i++)
                consts[i] = double.Parse(constants[i].RemoveSpaces());
            return (new(coefficients), consts);
        }
    }

    /// <summary>
    /// Different methods to solve a linear system
    /// </summary>
    public enum SolutionAlgorithm
    {
        /// <summary>
        /// Solving algorithm involving the null space of the coefficient
        /// matrix
        /// </summary>
        NullSpace,
        /// <summary>
        /// Solving algorithm using LU decomposition
        /// </summary>
        LU
    }
}

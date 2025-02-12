using static System.Math;

namespace LinearAlgebra
{
    /// <summary>
    /// A matrix with real entries
    /// </summary>
    public class Matrix : IEquatable<Matrix>
    {
        static int ValidateDimension(int candidate) =>
                                                    candidate < 0
                                                    ? throw new InvalidDimensionException()
                                                    : candidate;
        protected const double _eps = 1e-12;
        double[,] mat;
        int r, c;
        /// <summary>
        /// The number of rows in the matrix
        /// </summary>
        public int Rows
        {
            get => r;
            set
            {
                r = ValidateDimension(value);
            }
        }
        /// <summary>
        /// The number of columns in the matrix
        /// </summary>
        public int Columns
        {
            get => c;
            set
            {
                c = ValidateDimension(value);
            }
        }
        public Matrix(int rows, int columns)
        {
            Rows = rows;
            Columns = columns;
            mat = new double[Rows, Columns];
        }
        public Matrix(double[,] mat)
        {
            this.mat = mat;
            Rows = mat.GetLength(0);
            Columns = mat.GetLength(1);
        }
        /// <summary>
        /// Checks whether two matrices have the same dimensions
        /// </summary>
        public static bool HaveSameDimensions(Matrix m1, Matrix m2) => m1.Rows == m2.Rows && m1.Columns == m2.Columns;
        /// <summary>
        /// Whether or not the current instance of <see cref="Matrix"/> is square
        /// </summary>
        public bool IsSquare => Rows == Columns;
        public double this[int row, int col]
        {
            get => mat[row, col];
            set
            {
                mat[row, col] = value;
            }
        }
        public static implicit operator Matrix(double[,] mat) => new(mat);
        public static implicit operator double[,](Matrix mat) => mat.mat;

        #region Basic Operations
        // checks the equality of matrices
        static bool AreEqual(Matrix m1, Matrix m2)
        {
            if (!HaveSameDimensions(m1, m2))
                return false;
            for (int i = 0; i < m1.Rows; i++)
                for (int j = 0; j < m1.Columns; j++)
                    if (Abs(m1[i, j] - m2[i, j]) >= _eps)
                        return false;
            return true;
        }
        // adds two matrices
        static Matrix Add(Matrix m1, Matrix m2)
        {
            if (!HaveSameDimensions(m1, m2))
                throw new DifferentDimensionException();
            Matrix sum = new(m1);
            for (int i = 0; i < m1.Rows; i++)
                for (int j = 0; j < m1.Columns; j++)
                    sum[i, j] += m2[i, j];
            return sum;
        }
        // multiplies by a scalar
        static Matrix ProductByScalar(double k, Matrix m)
        {
            Matrix prod = new(m);
            for (int i = 0; i < m.Rows; i++)
                for (int j = 0; j < m.Columns; j++)
                    prod[i, j] *= k;
            return prod;
        }
        // multiplies two matrices
        static Matrix Product(Matrix m1, Matrix m2)
        {
            if (m1.Columns != m2.Rows)
                throw new InvalidMultiplicationException();
            Matrix prod = new(m1.Rows, m2.Columns);
            for (int i = 0; i < m1.Rows; i++)
                for (int j = 0; j < m2.Columns; j++)
                    for (int k = 0; k < m1.Columns; k++)
                        prod[i, j] += m1[i, k] * m2[k, j];
            return prod;
        }
        // forms the transpose
        static Matrix Transpose(Matrix m)
        {
            Matrix transpose = new(m.Columns, m.Rows);
            for (int i = 0; i < m.Rows; i++)
                for (int j = 0; j < m.Columns; j++)
                    transpose[j, i] = m[i, j];
            return transpose;
        }
        static double Trace(Matrix m)
        {
            if (!m.IsSquare)
                throw new SquareMatrixRequiredException();
            double s = 0;
            for (int i = 0; i < m.Rows; i++)
                s += m[i, i];
            return s;
        }
        /// <summary>
        /// Returns the <paramref name="n"/> x <paramref name="n"/> identity matrix.
        /// </summary>
        /// <param name="n">The order of the matrix</param>
        public static Matrix Id(int n)
        {
            Matrix id = new(n, n);
            for (int i = 0; i < n; i++)
                id[i, i] = 1;
            return id;
        }
        // matrix power
        static Matrix Pow(Matrix m, int p)
        {
            if (!m.IsSquare)
                throw new SquareMatrixRequiredException();

            ArgumentOutOfRangeException.ThrowIfNegative(p);
            if (p == 0)
                return Id(m.Rows);

            Matrix power = new(m);
            for (int i = 1; i < p; i++)
                power = Product(m, power);
            return power;
        }
        // Computes the matrix exponential (simply for fun)
        static Matrix Exp(Matrix m, int cutoff)
        {
            if (!m.IsSquare)
                throw new SquareMatrixRequiredException();
            var id = Id(m.Rows);
            if (cutoff == 0)
                return id;
            /*
                e^M = I + M + M^2 / 2! + M^3 / 3! + M^4 / 4! + ... + M^cutoff / cutoff!
            */
            int nFactorial = 1;
            Matrix exp = id;
            for (int i = 1; i <= cutoff; i++)
            {
                nFactorial *= i;
                var thisTerm = ProductByScalar(1.0 / nFactorial, Pow(m, i));
                exp = Add(exp, thisTerm);
            }
            return exp;
        }
        #endregion
        /// <summary>
        /// Calculates the matrix exponential up to a certain cutoff point (NUMERICALLY UNSTABLE)
        /// </summary>
        public Matrix Exp(int cutoff) => Exp(this, cutoff);
        /// <summary>
        /// Calculates the <paramref name="exponent"/>th power of the matrix.
        /// </summary>
        public Matrix Pow(int exponent) => Pow(this, exponent);
        /// <summary>
        /// Forms the transpose of the matrix.
        /// </summary>
        public Matrix Transpose() => Transpose(this);
        /// <summary>
        /// Shorthand for the transpose of the matrix
        /// </summary>
        public Matrix T => Transpose();
        /// <summary>
        /// Returns the trace of the matrix
        /// </summary>
        public double Trace() => Trace(this);

        public bool Equals(Matrix? other)
        {
            if (other as object == null)
                return false;
            return AreEqual(this, other);
        }

        #region Operations

        public static Matrix operator +(Matrix m) => m;
        public static Matrix operator +(Matrix m1, Matrix m2) => Add(m1, m2);
        public static Matrix operator *(double k, Matrix m) => ProductByScalar(k, m);
        public static Matrix operator *(Matrix m, double k) => k * m;
        public static Matrix operator *(Matrix m1, Matrix m2) => Product(m1, m2);
        public static Matrix operator -(Matrix m) => -1 * m;
        public static Matrix operator -(Matrix m1, Matrix m2) => m1 + (-m2);
        public static Matrix operator ^(Matrix m1, int p) => Pow(m1, p);
        public static bool operator ==(Matrix m1, Matrix m2) => m1.Equals(m2);
        public static bool operator !=(Matrix m1, Matrix m2) => m1.GetHashCode() != m2.GetHashCode();

        #endregion

        public override bool Equals(object? obj)
        {
            return Equals(obj as Matrix);
        }

        /// <summary>
        /// Returns the column vector of the column <paramref name="col"/>
        /// </summary>
        public Vector ColumnVector(int col)
        {
            Vector v = new(Rows);
            for (int i = 0; i < Rows; i++)
                v[i] = this[i, col];
            return v;
        }
        /// <summary>
        /// Returns the <paramref name="i"/>'th <paramref name="n"/>-dimensional 
        /// elementary vector (standard basis vector)
        /// </summary>
        public static Vector ElementaryVector(int n, int i)
        {
            var id = Id(n);
            return id.ColumnVector(i);
        }

        public override int GetHashCode() => mat.GetHashCode();

        #region Advanced Matrix Functionalities

        /*
        LU Decomposition, solving, determinants, inverses
        ref, rref
        solving linear system
        null space
        rank
        QR Decomposition (after creating a vector class)
        */

        // Elementary matrix functions
        /// <summary>
        /// Creates the elementary matrix for Gauss row reduction to obtain a pivot in 
        /// the column <paramref name="col"/>.
        /// </summary>
        public Matrix GaussElementaryMatrix(int col)
        {
            if (Abs(this[col, col]) < _eps)
                return Id(Rows);
            Matrix E = new(Rows, Rows);
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Rows; j++)
                {
                    if (i == j)
                        E[i, j] = 1;
                    if (j == col)
                        for (int s = 0; s < Rows; s++)
                            E[s, j] = -this[s, j] / this[col, j];
                }
            }
            return E;
        }
        /// <summary>
        /// Inverts a given Gauss elementary matrix
        /// </summary>
        public static Matrix GaussInverseElementaryMatrix(Matrix gaussElemMat)
        {
            Matrix Einv = new(gaussElemMat);
            for (int i = 0; i < gaussElemMat.Rows; i++)
            {
                for (int j = 0; j < gaussElemMat.Columns; j++)
                {
                    if (i != j && Abs(Einv[i, j]) >= _eps)
                        Einv[i, j] *= -1;
                }
            }
            return Einv;
        }
        /// <summary>
        /// Creates a <paramref name="dimension"/> x <paramref name="dimension"/> matrix that 
        /// permutes the rows <paramref name="row1"/> and <paramref name="row2"/>.
        /// </summary>
        public static Matrix PermutationMatrix(int dimension, int row1, int row2)
        {
            if (row1 == row2)
                return Id(dimension);
            Matrix P = new(dimension, dimension);
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    if (i == j && i != row1 && i != row2)
                        P[i, j] = 1;
                    if (i == row1)
                    {
                        P[i, row1] = 0;
                        P[i, row2] = 1;
                    }
                    if (i == row2)
                    {
                        P[i, row2] = 0;
                        P[i, row1] = 1;
                    }
                }
            }
            return P;
        }
        /// <summary>
        /// Creates the elementary matrix of order <paramref name="dimension"/> that multiplies row <paramref name="row"/>
        /// by the scalar <paramref name="k"/>.
        /// </summary>
        public static Matrix MultiplyRowMatrix(int dimension, double k, int row)
        {
            Matrix M = new(dimension, dimension);
            for (int i = 0; i < dimension; i++)
            {
                M[i, i] = (i == row) ? k : 1;
            }
            return M;
        }
        /// <summary>
        /// Computes the ordered product of the elements of the list <paramref name="lst"/>.
        /// </summary>
        public static Matrix ListProduct(List<Matrix> lst)
        {
            Matrix res = lst[0];
            for (int i = 1; i < lst.Count; i++)
                res *= lst[i];
            return res;
        }
        /// <summary>
        /// Computes the sign of the permutation given by the permutation matrix <paramref name="permutation"/>
        /// </summary>
        public static int Sign(Matrix permutation)
        {
            int[,] M = new int[2, permutation.Rows]; // 2-dimensional array representing all the permutations
            for (int i = 0; i < permutation.Rows; i++)
            {
                M[0, i] = i;
                for (int k = 0; k < permutation.Columns; k++)
                {
                    if (permutation[i, k] == 1)
                    {
                        M[1, i] = k;
                        break;
                    }
                }
            }
            List<int[]> pairs = new(); // List of all pairs of indices (i,j) such that j > i
            for (int i = 0; i < permutation.Rows; i++)
                for (int j = 0; j < permutation.Columns; j++)
                    if (j > i)
                        pairs.Add(new int[] { i, j });
            double prod = 1;
            foreach (var pair in pairs)
                prod *= (M[1, pair[1]] - M[1, pair[0]]) / (pair[1] - pair[0]);
            return (int)prod;
        }
        /// <summary>
        /// Searches for the row of an admissible and maximal pivot in the column <paramref name="col"/> 
        /// if <paramref name="nPivots"/> pivots have been found so far.  Returns -1 if the pivot
        /// entry has value 0.
        /// </summary>
        public int SearchForPivot(int nPivots, int col)
        {
            double m = Abs(this[nPivots, col]);
            int pivotRow = nPivots;
            for (int i = nPivots + 1; i < Rows; i++)
            {
                if (Abs(this[i, col]) > m)
                {
                    pivotRow = i;
                    m = Abs(this[i, col]);
                }
            }
            return (m < _eps) ? -1 : pivotRow;
        }

        /// <summary>
        /// Forms the elementary matrix for Gauss-Jordan row-reduction along the column <paramref name="col"/> 
        /// if <paramref name="nPivots"/> pivots have been found so far.
        /// </summary>
        public Matrix GaussJordanElementaryMatrix(int nPivots, int col)
        {
            Matrix E = Id(Rows);
            for (int i = 0; i < Rows; i++)
                if (i != nPivots)
                    E[i, nPivots] = -this[i, col] / this[nPivots, col];
            return E;
        }

        /// <summary>
        /// Searches for the row of a pivot for Gauss row-reduction in 
        /// the column <paramref name="col"/> after having found 
        /// <paramref name="nPivots"/> pivots
        /// </summary>
        public int SearchForPivotGauss(int nPivots, int col)
        {
            int row = nPivots;
            if (Abs(this[row, col]) >= _eps)
                return row;
            for (int i = row; i < Columns; i++)
                if (Abs(this[i, col]) >= _eps)
                    return i;
            return -1;
        }
        /// <summary>
        /// Computes the LU/PLU decomposition of this instance.
        /// </summary>
        /// <returns>
        /// A list containing, in order, 
        /// <list type="number">
        /// <item>The total necessary permutation</item>
        /// <item>The lower triangular matrix</item>
        /// <item>The upper triangular matrix (row-echelon form)</item>
        /// <item>The inverse of the permutation</item>
        /// </list>
        /// </returns>
        public List<Matrix> LUDecomposition()
        {
            if (!IsSquare)
                throw new SquareMatrixRequiredException();
            int nPivots = 0;
            var id = Id(Rows);
            Matrix reduced = this;
            Matrix P = id;
            Matrix E = new(Rows, Rows);
            List<Matrix> permutations = new();
            List<Matrix> elems = new();
            for (int i = 0; i < Columns; i++)
            {
                Matrix perm = id;
                if (nPivots == Rows)
                    break;
                int row = SearchForPivotGauss(nPivots, i);
                if (row != nPivots)
                    perm = PermutationMatrix(Rows, nPivots, row);
                permutations.Add(perm);
                reduced = perm * reduced;
                P = perm * P;
                E = GaussElementaryMatrix(i);
                reduced = E * reduced;
                elems.Add(E);
                nPivots++;
            }
            List<Matrix> inverses = new();
            for (int i = 0; i < elems.Count; i++)
                inverses.Add(GaussInverseElementaryMatrix(elems[i]));
            return new() { P, ListProduct(inverses), reduced, ListProduct(permutations) };
        }
        /// <summary>
        /// Find the row echelon form of the matrix
        /// </summary>
        public Matrix RowEchelonForm() => LUDecomposition()[2];
        /// <summary>
        /// Takes the product of all the elements down the diagonal of 
        /// the matrix
        /// </summary>
        public double DiagonalProduct()
        {
            if (!IsSquare)
                throw new SquareMatrixRequiredException();
            double prod = 1;
            for (int i = 0; i < Rows; i++)
                prod *= this[i, i];
            return prod;
        }
        /// <summary>
        /// Computes the determinant of the matrix
        /// </summary>
        public double Det()
        {
            if (!IsSquare)
                throw new SquareMatrixRequiredException();
            var PLU = LUDecomposition();
            var perm = PLU[0];
            var upper = PLU[2];
            double detPerm = Sign(perm);
            double detUpper = upper.DiagonalProduct();
            return detPerm * detUpper;
        }
        /// <summary>
        /// Forms the augmented matrix created by gluing <paramref name="rhs"/> 
        /// to the right of the other matrix
        /// </summary>
        public Matrix Augment(Matrix rhs)
        {
            if (Rows != rhs.Rows)
                throw new InvalidRowCountException();
            Matrix res = new(Rows, Columns + rhs.Columns);
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns + rhs.Columns; j++)
                {
                    if (j < Columns)
                        res[i, j] = this[i, j];
                    else if (j >= Columns)
                        res[i, j] = rhs[i, j - Columns];
                }
            }
            return res;
        }

        /// <summary>
        /// Creates the matrix formed by gluing together all matrices in <paramref name="lst"/>
        /// </summary>
        public static Matrix AugmentList(List<Matrix> lst)
        {
            Matrix mat = lst[0];
            for (int i = 1; i < lst.Count; i++)
                mat = mat.Augment(lst[i]);
            return mat;
        }
        public static Matrix AugmentList(List<Vector> lst)
        {
            Matrix mat = lst[0];
            for (int i = 1; i < lst.Count; i++)
                mat = mat.Augment(lst[i]);
            return mat;
        }
        /// <summary>
        /// Calculates the inverse of the matrix
        /// </summary>
        public Matrix Inverse()
        {
            if (!IsSquare)
                throw new SquareMatrixRequiredException();
            List<Matrix> Xis = new();
            for (int i = 0; i < Rows; i++)
                Xis.Add(Solver.SolveLU(this, ElementaryVector(Rows, i)));

            return AugmentList(Xis);
        }

        /// <summary>
        /// Extension of the simple reduced row echelon form.
        /// </summary>
        /// <returns>
        /// A list containing the reduced row echelon form and the pseudo-inverse
        /// of this matrix
        /// </returns>
        public List<Matrix> ExtendedRREF()
        {
            int nPivots = 0;
            var id = Id(Rows);
            var pseudoInverse = id;
            Matrix E = new(Rows, Rows);
            Matrix reduced = this;
            Matrix perm;
            for (int i = 0; i < Columns; i++)
            {
                if (nPivots == Rows)
                    break;
                int row = reduced.SearchForPivot(nPivots, i);
                perm = id;
                if (row != nPivots)
                    perm = PermutationMatrix(Rows, nPivots, row);
                reduced = perm * reduced;
                pseudoInverse = perm * pseudoInverse;
                E = reduced.GaussJordanElementaryMatrix(nPivots, i);
                reduced = E * reduced;
                pseudoInverse = E * pseudoInverse;
                nPivots++;
            }

            // normalizing pivots
            for (int i = 0; i < Rows; i++)
            {
                double nonZero = reduced.FirstNonzeroElement(i);
                if (Abs(nonZero) >= _eps)
                {
                    E = MultiplyRowMatrix(Rows, 1 / nonZero, i);
                    reduced = E * reduced;
                    pseudoInverse = E * pseudoInverse;
                }
            }
            return new() { reduced, pseudoInverse };
        }

        /// <summary>
        /// Computes the reduced row echelon form of this matrix
        /// </summary>
        public Matrix RREF() => ExtendedRREF()[0];

        /// <summary>
        /// Calculates the pseudo-inverse of this matrix; that is, the 
        /// product of all elementary matrices that put this matrix into
        /// reduced row echelon form.
        /// </summary>
        public Matrix PseudoInverse() => ExtendedRREF()[1];

        /// <summary>
        /// Returns the value of the first nonzero element in the row <paramref name="row"/>
        /// </summary>
        public double FirstNonzeroElement(int row)
        {
            for (int i = 0; i < Columns; i++)
                if (Abs(this[row, i]) >= _eps)
                    return this[row, i];
            return 0;
        }

        /// <summary>
        /// Returns a list of tuples representing all pairs of indices 
        /// corresponding to pivots in rref(A).
        /// </summary>
        public List<(int row, int col)> PivotPositions()
        {
            Matrix rref = RREF();
            List<(int row, int col)> pivotPositions = new();
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    if (Abs(rref[i, j]) >= _eps)
                    {
                        pivotPositions.Add((i, j));
                        break;
                    }
                }
            }
            return pivotPositions;
        }

        /// <summary>
        /// The rank of the matrix
        /// </summary>
        public int Rank => PivotPositions().Count;

        /// <summary>
        /// Returns a list of indices representing the linked variables in 
        /// the system
        /// </summary>
        public List<int> LinkedVariables()
        {
            var pivotPositions = PivotPositions();
            List<int> pivotColumns = pivotPositions.Select((i, j) => j).ToList();
            return pivotColumns;
        }
        /// <summary>
        /// Returns a list of indices representing the free variables in 
        /// the system
        /// </summary>
        public List<int> FreeVariables() => Enumerable.Range(0, Columns)
                                            .Where(
                                            n => !LinkedVariables().Contains(n)
                                            ).ToList();
        /// <summary>
        /// Returns the null space of the matrix, represented as a 
        /// list of vectors forming a basis of the null space.
        /// </summary>
        public List<Vector> NullSpace()
        {
            List<Vector> ns = new();
            var rref = RREF();
            var freeVars = FreeVariables();
            double[,] variableValues = new double[Columns, freeVars.Count];
            int treatedLinkedVariables = 0;
            for (int i = 0; i < Columns; i++)
            {
                if (freeVars.Contains(i))
                    variableValues[i, freeVars.IndexOf(i)] = 1;
                else
                {
                    for (int j = 0; j < freeVars.Count; j++)
                        variableValues[i, j] = -rref[treatedLinkedVariables, freeVars[j]];
                    treatedLinkedVariables++;
                }
            }
            for (int i = 0; i < freeVars.Count; i++)
                ns.Add(((Matrix)variableValues).ColumnVector(i));
            return ns;
        }

        public List<Vector> ToColumnVectorList()
        {
            List<Vector> cols = new();
            for (int i = 0; i < Columns; i++)
                cols.Add(ColumnVector(i));
            return cols;
        }

        /// <summary>
        /// Produces the orthogonal matrix in the QR decomposition of 
        /// this matrix
        /// </summary>
        public Matrix OrthogonalQ()
        {
            var gramSchmidt = Vector.GramSchmidt(ToColumnVectorList());
            return AugmentList(gramSchmidt);
        }
        /// <summary>
        /// Produces the upper triangular matrix R in the QR decomposition 
        /// of this matrix
        /// </summary>
        public Matrix OrthogonalR() => OrthogonalQ().T * this;

        /// <summary>
        /// Finds the QR decomposition of this matrix
        /// </summary>
        public (Matrix Q, Matrix R) QRDecomposition()
        {
            var Q = OrthogonalQ();
            return (Q, Q.T * this);
        }
        #endregion

        #region Eigenvectors, eigenvalues, diagonalization, ...
        /// <summary>
        /// Returns the list of the diagonal entries in the matrix <paramref name="mat"/> (does not check if <paramref name="mat"/> is square)
        /// </summary>
        public static List<double> DiagonalEntries(Matrix mat)
        {
            List<double> diag = new();
            for (int i = 0; i < mat.Rows; i++)
                diag.Add(mat[i, i]);
            return diag;
        }
        /// <summary>
        /// Checks if two matrices are within <paramref name="tolerance"/> distance from eachother (meaning that no entry in 
        /// (<paramref name="A"/> - <paramref name="B"/>) (in absolute value) is greater than or equal to <paramref name="tolerance"/>.
        /// </summary>
        public static bool AllClose(Matrix A, Matrix B, double tolerance)
        {
            if (!HaveSameDimensions(A, B))
                return false;
            for (int i = 0; i < A.Rows; i++)
                for (int j = 0; j < A.Columns; j++)
                    if (Abs(A[i, j] - B[i, j]) >= tolerance)
                        return false;
            return true;
        }
        /// <summary>
        /// Checks if the current instance is a symmetric matrix.
        /// </summary>
        public bool IsSymmetric()
        {
            if (!IsSquare)
                return false;
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    if (Abs(this[i, j] - this[j, i]) >= _eps)
                        return false;
            return true;
        }
        /// <summary>
        /// Turns the matrix <paramref name="A"/> into an upper triangular matrix by 
        /// setting all entries in the lower triangle to zero.
        /// </summary>
        public static Matrix ToUpTri(Matrix A)
        {
            if (!A.IsSquare)
                throw new SquareMatrixRequiredException();
            Matrix U = new(A.Rows, A.Rows);
            for (int i = 0; i < A.Rows; i++)
                for (int j = i; j < A.Columns; j++)
                    U[i, j] = A[i, j];
            return U;
        }
        /// <summary>
        /// Approximates the eigenvalues of a symmetric matrix using the QR algorithm
        /// </summary>
        /// <remarks>
        /// It is important to note that numerical imprecision is bound to happen given the 
        /// sheer magnitude of calculations necessary for this algorithm. Answers risk only 
        /// being precise up to the second decimal place.
        /// </remarks>
        public List<double> SymmetricEigenvalues()
        {
            if (!IsSymmetric())
                throw new ArgumentException("Matrix must be symmetric to apply this algorithm.");

            const double tolerance = 1e-12;

            var QR = QRDecomposition();
            Matrix A_i = QR.R * QR.Q;
            do
            {
                QR = A_i.QRDecomposition();
                A_i = QR.R * QR.Q;
            }
            while (!AllClose(A_i, ToUpTri(A_i), tolerance));
            return DiagonalEntries(A_i);
        }

        #endregion
    }
}
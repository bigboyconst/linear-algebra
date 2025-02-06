using static System.Math;

namespace LinearAlgebra
{
    /// <summary>
    /// A column vector
    /// </summary>
    public class Vector : Matrix, IEquatable<Vector>
    {
        public bool Equals(Vector? other) => Equals(other as Matrix);
        public Vector(int count) : base(count, 1) { }
        public int Count => Rows;
        static double[,] To2DimensionalArray(double[] oneD)
        {
            double[,] twoDee = new double[oneD.Length, 1];
            for (int i = 0; i < oneD.Length; i++)
                twoDee[i, 0] = oneD[i];
            return twoDee;
        }
        public static Vector ToVector(Matrix mat)
        {
            if (mat.Columns != 1)
                throw new InvalidVectorException();
            double[] v = new double[mat.Rows];
            for (int i = 0; i < mat.Rows; i++)
                v[i] = mat[i, 0];
            return new(v);
        }
        public Vector(params double[] vec) : base(To2DimensionalArray(vec)) { }
        public override bool Equals(object? obj)
        {
            return Equals(obj as Vector);
        }

        public double this[int i]
        {
            get => base[i, 0];
            set
            {
                base[i, 0] = value;
            }
        }

        public override int GetHashCode() => base.GetHashCode();

        /// <summary>
        /// Computes the dot product of the vectors <paramref name="a"/> and <paramref name="b"/>.
        /// </summary>
        public static double Dot(Vector a, Vector b)
        {
            if (!HaveSameDimensions(a, b))
                throw new DifferentDimensionException();
            double s = 0;
            for (int i = 0; i < a.Count; i++)
                s += a[i] * b[i];
            return s;
        }
        Vector MultiplyByScalar( double k)
        {
            Vector v = new(Count);
            for (int i = 0; i < Count; i++)
                v[i] = k * this[i];
            return v;
        }

        /// <summary>
        /// Calculates the cross product of 2- or 3-dimensional vectors.
        /// </summary>
        public static Vector Cross(Vector a, Vector b)
        {
            if (a.Count != 3 || b.Count != 3)
                throw new InvalidCrossProductException();
            if (a.Count == 2)
                a = new(a[0], a[1], 0);
            if (b.Count == 2)
                b = new(b[0], b[1], 0);
            Vector cross = new(3);
            cross[0] = a[1] * b[2] - a[2] * b[1];
            cross[1] = -(a[0] * b[2] - a[2] * b[0]);
            cross[2] = a[0] * b[1] - a[1] * b[0];
            return cross;
        }
        public static Vector operator +(Vector a) => a;
        public static Vector operator -(Vector a) => ToVector(-(a as Matrix));
        public static Vector operator +(Vector a, Vector b) => ToVector((a as Matrix) + b);
        public static Vector operator -(Vector a, Vector b) => ToVector((a as Matrix) - b);
        public static Vector operator *(double k, Vector a) => a.MultiplyByScalar(k);
        public static Vector operator /(Vector a, double k) => 1 / k * a;
        public static Vector operator *(Vector a, double k) => k * a;
        public static Vector operator *(Matrix M, Vector a) => ToVector(M * (a as Matrix));
        
        /// <summary>
        /// The magnitude of the vector
        /// </summary>
        public double Norm => Sqrt(NormSquared);
        /// <summary>
        /// The squared magnitude (norm) of the vector
        /// </summary>
        public double NormSquared => Dot(this, this);
        /// <summary>
        /// The normalized version of the vector.  Note that this does not check if the vector's 
        /// magnitude is zero.
        /// </summary>
        public Vector Normalized => 1 / Norm * this;
        /// <summary>
        /// Computes the orthogonal projection of this instance onto the vector <paramref name="onto"/>.
        /// </summary>
        public Vector ProjectOnto(Vector onto) => Dot(this, onto) / onto.NormSquared * onto;
        /// <summary>
        /// Projects the vector <paramref name="projected"/> onto the vector <paramref name="onto"/>.
        /// </summary>
        public static Vector Proj(Vector onto, Vector projected) => projected.ProjectOnto(onto);
        /// <summary>
        /// Takes in the basis <paramref name="basis"/> and performs the Gram-Schmidt process to return 
        /// an orthonormal basis
        /// </summary>
        public static List<Vector> GramSchmidt(List<Vector> basis)
        {
            List<Vector> orthonormalBasis = new();
            if (basis[0].Norm < _eps)
                return orthonormalBasis;
            orthonormalBasis.Add(basis[0].Normalized);
            for (int i = 1; i < basis.Count; i++)
            {
                var u_i = basis[i];
                for (int k = 0; k < i; k++)
                {
                    var e_k = orthonormalBasis[k];
                    // projects v_i onto the k'th normalized vector
                    u_i -= Proj(e_k, basis[i]);
                }
                if (u_i.Norm < _eps)
                    break;
                orthonormalBasis.Add(u_i.Normalized);
            }
            return orthonormalBasis;
        }
    }
}

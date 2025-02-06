namespace LinearAlgebra
{
    /// <summary>
    /// Exception that is thrown when an invalid value is passed as the dimension of a matrix and/or vector
    /// </summary>
    public class InvalidDimensionException : Exception
    {
        public static readonly string message = "The value given to the dimension is not valid.";
        public InvalidDimensionException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when two arguments must have the same dimensions
    /// </summary>
    public class DifferentDimensionException : Exception
    {
        public static readonly string message = "The arguments must have the same dimensions.";
        public DifferentDimensionException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when the multiplication between two given matrices is not defined
    /// </summary>
    public class InvalidMultiplicationException : Exception
    {
        public static readonly string message = "The left hand side of the matrix product must have as " +
            "many columns as the right hand side has rows.";
        public InvalidMultiplicationException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when a square matrix is required for the operation
    /// </summary>
    public class SquareMatrixRequiredException : Exception
    {
        public static readonly string message = "A square matrix is required to perform this operation.";
        public SquareMatrixRequiredException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when the number of columns of two matrices do not match
    /// </summary>
    public class InvalidColumnCountException : Exception
    {
        public static readonly string message = "Matrices must have the same number of columns";
        public InvalidColumnCountException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when the number of rows in two matrices do not match
    /// </summary>
    public class InvalidRowCountException : Exception
    {
        public static readonly string message = "Matrices must have the same number of rows";
        public InvalidRowCountException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when trying to perform the cross product on vectors that are not 3 dimensional
    /// </summary>
    public class InvalidCrossProductException : Exception
    {
        public static readonly string message = "Cross products are only defined in 3 dimensions";
        public InvalidCrossProductException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when we try to create a vector with multiple columns
    /// </summary>
    public class InvalidVectorException : Exception
    {
        public static readonly string message = "A vector must be a matrix with one column";
        public InvalidVectorException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when two lists have different lengths
    /// </summary>
    public class InvalidLengthException : Exception
    {
        public static readonly string message = "The lists must have the same length";
        public InvalidLengthException() : base(message) { }
    }
    /// <summary>
    /// Exception that is thrown when given data points do not reflect the behaviour of a closed spline
    /// </summary>
    public class InvalidClosedSplineException : Exception
    {
        public static readonly string message = "The data does not represent a closed spline";
        public InvalidClosedSplineException() : base(message) { }
    }
}

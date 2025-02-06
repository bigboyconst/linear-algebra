#nullable disable
using System.Text;
using LinearAlgebra;
using static System.Math;


static void GenerateRandomLatexMatrix()
{
    Random rand = new();

    bool keepGenerating = true;

    do
    {
        Console.Write("Dimensions (format: [rows]x[cols]): ");

        var dims = ParseDims(Console.ReadLine());

        Console.WriteLine($"{dims.m} rows and {dims.n} cols");

        Matrix generated = RandomIntMatrix(dims.m, dims.n, rand);

        Console.WriteLine(ToLatexString(generated));

        Console.Write("Generate another? (n to exit): ");

        if (Console.ReadLine().Equals("N", StringComparison.CurrentCultureIgnoreCase))
            keepGenerating = false;
    }
    while (keepGenerating);
}
static (int m, int n) ParseDims(string dims)
{
    string[] dimsElems = dims.Split('x');
    return (int.Parse(dimsElems[0]), int.Parse(dimsElems[1]));
}
static Matrix RandomIntMatrix(int m, int n, Random rnd)
{
    Matrix mat = new(m, n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            mat[i, j] = rnd.Next(-10, 10);
    return mat;
}
static string ToLatexString(Matrix m)
{
    StringBuilder sb = new();
    sb.AppendLine("\\begin{bmatrix}");
    for (int i = 0; i < m.Rows; i++)
    {
        for (int j = 0; j < m.Columns; j++)
            sb.Append(m[i, j] + ((j == m.Columns - 1) ? "\\\\ \n" : " & "));
    }
    sb.AppendLine("\\end{bmatrix}");
    return sb.ToString();
}

static void Write(Matrix m, int nDigits = 4)
{
    for (int i = 0; i < m.Rows; i++)
    {
        for (int j = 0; j < m.Columns; j++)
            Console.Write($"{Math.Round(m[i, j], nDigits)}    ");
        Console.WriteLine();
    }
    Console.WriteLine();
}
static void WritePolynomials(double[,] formattedCoefficients, List<int> Lt)
{
    for (int i = 0; i < formattedCoefficients.GetLength(0); i++)
    {
        Console.WriteLine($"{Round(formattedCoefficients[i, 0], 4)} " +
            $"{SignString(formattedCoefficients[i, 1])}t " +
            $"{SignString(formattedCoefficients[i, 2])}t^2 " +
            $"{SignString(formattedCoefficients[i, 3])}t^3" +
            $"    t < {Lt[i] + 1}");
    }
}
static string SignString(double x)
{
    x = Round(x, 4);
    return x < 0 ? "- " + Abs(x).ToString() : "+ " + x.ToString();
}
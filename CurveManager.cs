using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LinearAlgebra
{
    public static class CurveManager
    {
        public static Vector[] SplineInterpolation(Vector[] controlPoints, int numSubdivisions, SplineType type = SplineType.Natural)
        {
            float step = 1f / numSubdivisions;
            List<double> x = new();
            List<double> y = new();
            List<double> z = new();
            for (int i = 0; i < controlPoints.Length; i++)
            {
                x.Add(controlPoints[i][0]);
                y.Add(controlPoints[i][1]);
                z.Add(controlPoints[i][2]);
            }
            List<int> tValues = Enumerable.Range(0, controlPoints.Length).ToList();

            var splineX = Spline.CubicSplineInterpolation(tValues, x, type);
            var splineY = Spline.CubicSplineInterpolation(tValues, y, type);
            var splineZ = Spline.CubicSplineInterpolation(tValues, z, type);

            int nPts = numSubdivisions * (controlPoints.Length - 1) + 1;

            double[] xSamples = new double[nPts];
            double[] ySamples = new double[nPts];
            double[] zSamples = new double[nPts];
            int iter = 0;
            for (double t = 0; t < controlPoints.Length - 1; t += step)
            {
                int row = iter / numSubdivisions;
                xSamples[iter] = splineX[row, 0]
                    + splineX[row, 1] * t
                    + splineX[row, 2] * Math.Pow(t, 2)
                    + splineX[row, 3] * Math.Pow(t, 3);
                ySamples[iter] = splineY[row, 0]
                    + splineY[row, 1] * t
                    + splineY[row, 2] * Math.Pow(t, 2)
                    + splineY[row, 3] * Math.Pow(t, 3);
                zSamples[iter] = splineZ[row, 0]
                    + splineZ[row, 1] * t
                    + splineZ[row, 2] * Math.Pow(t, 2)
                    + splineZ[row, 3] * Math.Pow(t, 3);
                iter++;
            }
            xSamples[iter] = x.Last();
            ySamples[iter] = y.Last();
            zSamples[iter] = z.Last();
            return ToVectorArray(xSamples, ySamples, zSamples);
        }
        static Vector[] ToVectorArray(double[] x, double[] y, double[] z)
        {
            Vector[] vecs = new Vector[x.Length];
            for (int i = 0; i < x.Length; i++)
                vecs[i] = new([x[i], y[i], z[i]]);
            return vecs;
        }
    }
}

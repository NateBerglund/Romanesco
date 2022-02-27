using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace Romanesco
{
    class Program
    {
        const double A = 0.5;
        const double B = 0.01;
        const double C = 1.0 / 6.0;
        const double ScaleFactor = 0.1;

        static void Main(string[] args)
        {
            // Note: fractal equations are modified version(s) of: https://medium.com/@rodrigosetti/the-broccolis-equation-2ed8faaa3f3a
            Matrix<double> S1 = Matrix<double>.Build.Dense(3, 3, (i, j) => i == j ? ScaleFactor : 0.0);

            // Broccoli = Cone Union S1 * Union_{theta = -100*pi by pi/4 to 100*pi} (S2 * R * T * Broccoli)

            //Equation of logarithmic spiral:
            //cylindrical coords:
            //r = A * e ^ (B * theta)
            //z = -slope * r

            //x = A * e ^ (B * theta) * cos(theta)
            //y = A * e ^ (B * theta) * sin(theta)
            //z = A * (e ^ (B * 20 * pi) - 2 * e ^ (B * theta))


        }

        /// <summary>
        /// Matrix that scales by C * e^(B * theta)
        /// </summary>
        /// <param name="theta">Angle</param>
        /// <returns>Scaling matrix</returns>
        static Matrix<double> S2(double theta)
        {
            return Matrix<double>.Build.Dense(3, 3, (i, j) => i == j ? C * Math.Exp(B * theta) : 0.0);
        }

        /// <summary>
        /// Matrix that rotates about ...
        /// </summary>
        /// <param name="theta">Angle</param>
        /// <returns>Scaling matrix</returns>
        static Matrix<double> R(double theta)
        {
            return Matrix<double>.Build.Dense(3, 3, (i, j) => i == j ? C * Math.Exp(B * theta) : 0.0);  // todo: put in correct forumula
        }

        /// <summary>
        /// Matrix that translates about ...
        /// </summary>
        /// <param name="theta">Angle</param>
        /// <returns>Scaling matrix</returns>
        static Vector<double> T(double theta)
        {
            return Vector<double>.Build.Dense(3, i => 0.0); // todo: put in correct forumula
        }
    }
}

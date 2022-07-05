using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace Romanesco
{
    class Program
    {
        const string outputFile = "romanesco.stl";

        const double sqrtHalf = 0.70710678118654757; // sqrt(0.5)
        const double latticeScaling = 1.0001; // amount to scale from lattice-index coordinates to world coordinates
        const double offsetX = -10.001; // Offset to apply to lattice-origin in world coordinates
        const double offsetY = -10.001;
        const double offsetZ = -6.499;

        const double A = 0.5;
        const double B = 0.01;
        const double C = 1.0 / 6.0;
        const double ScaleFactor = 0.1;

        const int nX = 20;
        const int nY = 20;
        const int nZ = 10;
        static double[,,,] latticePoints;

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

            InitializeLattice();
            ExportLatticeToMesh();
        }

        static void InitializeLattice()
        {
            latticePoints = new double[nX, nY, nZ, 3];
            for (int xIdx = 0; xIdx < nX; xIdx++)
            {
                for (int yIdx = 0; yIdx < nY; yIdx++)
                {
                    for (int zIdx = 0; zIdx < nZ; zIdx++)
                    {
                        double x = xIdx;
                        double y = yIdx;
                        double z = zIdx * sqrtHalf;
                        if (zIdx % 2 == 1) // odd z index
                        {
                            x += 0.5;
                            y += 0.5;
                        }
                        // (x,y,z) are now in "lattice coordinates"

                        x *= latticeScaling;
                        y *= latticeScaling;
                        z *= latticeScaling;
                        x += offsetX;
                        y += offsetY;
                        z += offsetZ;
                        // (x,y,z) are now in "world coordinates"

                        // Test if the point (x,y,z) is inside the starting cone
                        if (z <= 0 && (z * z) >= (x * x) + (y * y))
                        {
                            latticePoints[xIdx, yIdx, zIdx, 0] = x;
                            latticePoints[xIdx, yIdx, zIdx, 1] = y;
                            latticePoints[xIdx, yIdx, zIdx, 2] = z;
                        }
                    }
                }
            }
        }

        static void ExportLatticeToMesh()
        {
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.AppendLine("solid romanesco");

            for (int xIdx = -2; xIdx < nX + 2; xIdx++)
            {
                for (int yIdx = -2; yIdx < nY + 2; yIdx++)
                {
                    for (int zIdx = -1; zIdx < nZ + 1; zIdx++)
                    {
                        double x = xIdx;
                        double y = yIdx;
                        double z = zIdx * sqrtHalf;
                        if ((zIdx + 2) % 2 == 1) // odd z index
                        {
                            x += 0.5;
                            y += 0.5;
                        }
                        // (x,y,z) are now in "lattice coordinates"

                        bool currentCellIsOccupied =
                            xIdx >= 0 && xIdx < nX
                            && yIdx >= 0 && yIdx < nY
                            && zIdx >= 0 && zIdx < nZ &&
                            (latticePoints[xIdx, yIdx, zIdx, 0] != 0 ||
                            latticePoints[xIdx, yIdx, zIdx, 1] != 0 ||
                            latticePoints[xIdx, yIdx, zIdx, 2] != 0); // whether current lattice cell is occupied

                        bool neighborCellIsOccupied = false;

                        // Check the forwards neighbors in each of the x and y directions
                        int xIdxNbr = xIdx;
                        int yIdxNbr = yIdx;
                        int zIdxNbr = zIdx;
                        for (int dim = 0; dim < 2; dim++)
                        {
                            if (dim == 0)
                            {
                                xIdxNbr = xIdx + 1;
                                yIdxNbr = yIdx;
                            }
                            else
                            {
                                xIdxNbr = xIdx;
                                yIdxNbr = yIdx + 1;
                            }
                            neighborCellIsOccupied =
                                xIdxNbr >= 0 && xIdxNbr < nX
                                && yIdxNbr >= 0 && yIdxNbr < nY
                                && zIdxNbr >= 0 && zIdxNbr < nZ &&
                                (latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 0] != 0 ||
                                latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 1] != 0 ||
                                latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 2] != 0); // whether neigboring lattice cell is occupied

                            if (neighborCellIsOccupied != currentCellIsOccupied)
                            {
                                double[] coords = new double[3];
                                double[] relLatCoords = new double[3]; // lattice coords relative to the current cell
                                relLatCoords[dim] = 0.5;
                                int orien = currentCellIsOccupied ? 1 : -1; // surface orientation
                                for (int sgn = -1; sgn <= 1; sgn += 2)
                                {
                                    // Facet normals will be recomputed by the software we use to open the stl file, so it's not critical that they are correct
                                    stringBuilder.AppendLine("facet normal 0.0 0.0 1.0");
                                    stringBuilder.AppendLine("    outer loop");

                                    relLatCoords[1 - dim] = (2 * dim - 1) * sgn * 0.5;
                                    relLatCoords[2] = 0;
                                    coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                    coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                    coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                    stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                    relLatCoords[1 - dim] = 0;
                                    relLatCoords[2] = -sgn * orien * 0.5 * sqrtHalf;
                                    if (zIdx == 0 && relLatCoords[2] < 0)
                                    {
                                        relLatCoords[2] = 0;
                                    }
                                    coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                    coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                    coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                    stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                    relLatCoords[1 - dim] = 0;
                                    relLatCoords[2] = sgn * orien * 0.5 * sqrtHalf;
                                    if (zIdx == 0 && relLatCoords[2] < 0)
                                    {
                                        relLatCoords[2] = 0;
                                    }
                                    coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                    coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                    coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                    stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                    stringBuilder.AppendLine("    endloop");
                                    stringBuilder.AppendLine("endfacet");
                                }
                            }
                        }
                        // Check the upwards neighbors in each direction
                        int xShift = 0;
                        int yShift = 0;
                        if ((zIdx + 2) % 2 == 1) // odd z index
                        {
                            xShift = 1;
                            yShift = 1;
                        }
                        zIdxNbr = zIdx + 1;
                        for (int xDir = -1; xDir <= 1; xDir += 2)
                        {
                            for (int yDir = -1; yDir <= 1; yDir += 2)
                            {
                                int xToAdd = (xDir - 1) / 2 + xShift;
                                int yToAdd = (yDir - 1) / 2 + yShift;
                                xIdxNbr = xIdx + xToAdd;
                                yIdxNbr = yIdx + yToAdd;

                                neighborCellIsOccupied =
                                xIdxNbr >= 0 && xIdxNbr < nX
                                && yIdxNbr >= 0 && yIdxNbr < nY
                                && zIdxNbr >= 0 && zIdxNbr < nZ &&
                                (latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 0] != 0 ||
                                latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 1] != 0 ||
                                latticePoints[xIdxNbr, yIdxNbr, zIdxNbr, 2] != 0); // whether neigboring lattice cell is occupied

                                if (neighborCellIsOccupied != currentCellIsOccupied)
                                {
                                    double[] coords = new double[3];
                                    double[] relLatCoords = new double[3]; // lattice coords relative to the current cell
                                    int orien = currentCellIsOccupied ? 1 : -1; // surface orientation
                                    for (int sgn = -1; sgn <= 1; sgn += 2)
                                    {
                                        // Facet normals will be recomputed by the software we use to open the stl file, so it's not critical that they are correct
                                        stringBuilder.AppendLine("facet normal 0.0 0.0 1.0");
                                        stringBuilder.AppendLine("    outer loop");

                                        relLatCoords[0] = sgn < 0 ? 0 : 0.5 * xDir;
                                        relLatCoords[1] = sgn < 0 ? 0 : 0.5 * yDir;
                                        relLatCoords[2] = sgn < 0 ? sqrtHalf : 0;
                                        if (zIdx == -1)
                                        {
                                            relLatCoords[2] = sqrtHalf;
                                        }
                                        coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                        coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                        coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                        stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                        relLatCoords[0] = ((1 - sgn * xDir * yDir * orien) / 2) * 0.5 * xDir;
                                        relLatCoords[1] = ((sgn * xDir * yDir * orien + 1) / 2) * 0.5 * yDir;
                                        relLatCoords[2] = 0.5 * sqrtHalf;
                                        if (zIdx == -1)
                                        {
                                            relLatCoords[2] = sqrtHalf;
                                        }
                                        coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                        coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                        coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                        stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                        relLatCoords[0] = ((sgn * xDir * yDir * orien + 1) / 2) * 0.5 * xDir;
                                        relLatCoords[1] = ((1 - sgn * xDir * yDir * orien) / 2) * 0.5 * yDir;
                                        relLatCoords[2] = 0.5 * sqrtHalf;
                                        if (zIdx == -1)
                                        {
                                            relLatCoords[2] = sqrtHalf;
                                        }
                                        coords[0] = latticeScaling * (x + relLatCoords[0]) + offsetX;
                                        coords[1] = latticeScaling * (y + relLatCoords[1]) + offsetY;
                                        coords[2] = latticeScaling * (z + relLatCoords[2]) + offsetZ;
                                        stringBuilder.AppendLine("        vertex " + coords[0].ToString() + " " + coords[1].ToString() + " " + coords[2].ToString());

                                        stringBuilder.AppendLine("    endloop");
                                        stringBuilder.AppendLine("endfacet");
                                    }
                                }
                            }
                        }
                    }
                }
            }

            stringBuilder.AppendLine("endsolid");

            File.WriteAllText(outputFile, stringBuilder.ToString());
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

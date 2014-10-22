using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using log4net;

namespace VectorMath
{

    public class Vector
    {

        public class EdgeFamily
        {
            public List<MemorySafe_CartCoord> startendpt { get; set; }
            public string sbdec { get; set; }
            public List<EdgeFamily> relatedEdges { get; set; }
        }

        private static readonly ILog log = LogManager.GetLogger(typeof(Vector));
        public enum WalkDirection
        {
            Clockwise,
            Counterclockwise,
        };

        public class CartCoord
        {
            public double X { get; set; }
            public double Y { get; set; }
            public double Z { get; set; }
        }

        public class MemorySafe_CartCoord
        {
            private readonly double _X, _Y, _Z;

            public MemorySafe_CartCoord(double X, double Y, double Z)
            {
                _X = X;
                _Y = Y;
                _Z = Z;
            }

            public double X { get { return _X; } }
            public double Y { get { return _Y; } }
            public double Z { get { return _Z; } }
        }

        public class PointVector
        {
            public MemorySafe_CartCoord p { get; set; }
            public MemorySafe_CartVect v { get; set; }
            public double s { get; set; }
            public PointVector(MemorySafe_CartCoord p1, MemorySafe_CartCoord p2)
            {
                p = p1;
                double vX = p2.X - p1.X;
                double vY = p2.Y - p1.Y;
                double vZ = p2.Z - p1.Z;
                v = new MemorySafe_CartVect(vX, vY, vZ);
                s = 1;
            }

        }

        public class MemorySafe_LongCoord
        {
            private readonly long _X, _Y, _Z;

            public MemorySafe_LongCoord(long X, long Y, long Z)
            {
                _X = X;
                _Y = Y;
                _Z = Z;
            }

            public long X { get { return _X; } }
            public long Y { get { return _Y; } }
            public long Z { get { return _Z; } }
        }

        public class CartVect
        {
            public CartVect()
            {

            }

            private double _X,_Y,_Z;
            public CartVect(double X, double Y, double Z)
            {
                _X = X;
                _Y = Y;
                _Z = Z;
            }
            public double X { get; set; }
            public double Y { get; set; }
            public double Z { get; set; }

        }

        public class MemorySafe_CartVect
        {
            private readonly double _X, _Y, _Z;

            public MemorySafe_CartVect(double X, double Y, double Z)
            {
                _X = X;
                _Y = Y;
                _Z = Z;
            }

            public double X { get { return _X; } }
            public double Y { get { return _Y; } }
            public double Z { get { return _Z; } }

        }

        public static CartVect CreateVector(CartCoord cd1, CartCoord cd2)
        {
            CartVect vector = new CartVect();
            vector.X = cd2.X - cd1.X;
            vector.Y = cd2.Y - cd1.Y;
            vector.Z = cd2.Z - cd1.Z;
            return vector;
        }

        public static CartVect CreateVector(MemorySafe_CartCoord cd1, MemorySafe_CartCoord cd2)
        {
            CartVect vector = new CartVect();
            vector.X = cd2.X - cd1.X;
            vector.Y = cd2.Y - cd1.Y;
            vector.Z = cd2.Z - cd1.Z;
            return vector;
        }

        public static MemorySafe_CartVect CreateMemorySafe_Vector(MemorySafe_CartCoord cd1, MemorySafe_CartCoord cd2)
        {
            double X = cd2.X - cd1.X;
            double Y = cd2.Y - cd1.Y;
            double Z = cd2.Z - cd1.Z;
            MemorySafe_CartVect vector = new MemorySafe_CartVect(X,Y,Z);
            return vector;
        }

        public static MemorySafe_CartVect CreateMemorySafe_Vector(CartCoord cd1, CartCoord cd2)
        {
            double X = cd2.X - cd1.X;
            double Y = cd2.Y - cd1.Y;
            double Z = cd2.Z - cd1.Z;
            MemorySafe_CartVect vector = new MemorySafe_CartVect(X, Y, Z);
            return vector;
        }

        public static CartVect convertToTempVector(MemorySafe_CartVect vector)
        {
            CartVect Vect = new CartVect();
            Vect.X = vector.X;
            Vect.Y = vector.Y;
            Vect.Z = vector.Z;
            return Vect;
        }

        public static MemorySafe_CartVect convertToMemorySafeVector(CartVect vector)
        {
            MemorySafe_CartVect memVect = new MemorySafe_CartVect(vector.X, vector.Y, vector.Z);
            return memVect;
        }

        public static MemorySafe_CartCoord convertToMemorySafeCoord(CartCoord coord)
        {
            MemorySafe_CartCoord memCoord = new MemorySafe_CartCoord(coord.X, coord.Y, coord.Z);
            return memCoord;
        }

        public static Double VectorMagnitude(CartVect vector)
        {
            double magnitude= 0.0;
            
            magnitude = Math.Sqrt(Math.Pow((vector.X),2) + Math.Pow((vector.Y),2) + Math.Pow((vector.Z),2));
            return magnitude;
        }

        public static Double VectorMagnitude(MemorySafe_CartVect vector)
        {
            double magnitude = 0.0;

            magnitude = Math.Sqrt(Math.Pow((vector.X), 2) + Math.Pow((vector.Y), 2) + Math.Pow((vector.Z), 2));
            return magnitude;
        }

        public static CartVect UnitVector(CartVect vector)
        {
            CartVect UV = new CartVect();
            double magnitude = VectorMagnitude(vector);

            UV.X = vector.X / magnitude;
            UV.Y = vector.Y / magnitude;
            UV.Z = vector.Z / magnitude;
            return UV;
        }

        public static MemorySafe_CartVect UnitVector(MemorySafe_CartVect vector)
        {
            double magnitude = VectorMagnitude(vector);
            double X = vector.X / magnitude;
            double Y = vector.Y / magnitude;
            double Z = vector.Z / magnitude;
            MemorySafe_CartVect UV = new MemorySafe_CartVect(X,Y,Z);
            return UV;
        }

        public static MemorySafe_CartVect VectorTimesScalar(MemorySafe_CartVect vect, double scalar)
        {
            double newx = vect.X * scalar;
            double newy = vect.Y * scalar;
            double newz = vect.Z * scalar;
            Vector.MemorySafe_CartVect newvect = new Vector.MemorySafe_CartVect(newx, newy, newz);
            return newvect;
        }

        public static MemorySafe_CartCoord SumPointAndLine(MemorySafe_CartCoord coord,MemorySafe_CartVect vect)
        {
            try
            {
                double retx = coord.X + vect.X;
                double rety = coord.Y + vect.Y;
                double retz = coord.Z + vect.Z;
                Vector.MemorySafe_CartCoord retcoord = new MemorySafe_CartCoord(retx,rety,retz);
                return retcoord;
            }
            catch
            {

            }
            //on the off chance something goes awry
            return new MemorySafe_CartCoord(-999,-999,-999);
        }

        public static CartVect CrossProduct(CartVect vector1, CartVect vector2)
        {
            CartVect xProd = new CartVect();

            xProd.X = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            xProd.Y = -1*(vector2.Z * vector1.X - vector1.Z * vector2.X);
            xProd.Z = vector2.Y * vector1.X - vector1.Y * vector2.X;
            return xProd;
        }

        public static CartVect CrossProductNVRetMSNV(MemorySafe_CartVect vector1, CartVect vector2)
        {
            double xProdX = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            double xProdY = -1 * (vector2.Z * vector1.X - vector1.Z * vector2.X);
            double xProdZ = vector2.Y * vector1.X - vector1.Y * vector2.X;
            CartVect xProd = new CartVect(xProdX, xProdY, xProdZ);
            xProd.X = xProdX;
            xProd.Y = xProdY;
            xProd.Z = xProdZ;
            return xProd;
        }

        public static CartVect CrossProductNVRetMSMS(MemorySafe_CartVect vector1, MemorySafe_CartVect vector2)
        {
            double xProdX = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            double xProdY = -1 * (vector2.Z * vector1.X - vector1.Z * vector2.X);
            double xProdZ = vector2.Y * vector1.X - vector1.Y * vector2.X;
            CartVect xProd = new CartVect(xProdX, xProdY, xProdZ);
            xProd.X = xProdX;
            xProd.Y = xProdY;
            xProd.Z = xProdZ;
            return xProd;
        }

        public static CartVect CrossProductNVRetNVMS(CartVect vector1, MemorySafe_CartVect vector2)
        {
            double xProdX = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            double xProdY = -1 * (vector2.Z * vector1.X - vector1.Z * vector2.X);
            double xProdZ = vector2.Y * vector1.X - vector1.Y * vector2.X;
            CartVect xProd = new CartVect(xProdX, xProdY, xProdZ);
            xProd.X = xProdX;
            xProd.Y = xProdY;
            xProd.Z = xProdZ;
            return xProd;
        }

        public static MemorySafe_CartVect CrossProduct(MemorySafe_CartVect vector1, MemorySafe_CartVect vector2)
        {
            double xProdX = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            double xProdY = -1 * (vector2.Z * vector1.X - vector1.Z * vector2.X);
            double xProdZ = vector2.Y * vector1.X - vector1.Y * vector2.X;
            MemorySafe_CartVect xProd = new MemorySafe_CartVect(xProdX, xProdY, xProdZ);
            return xProd;
        }

        public static MemorySafe_CartVect CrossProductMSRetMSNV(MemorySafe_CartVect vector1, CartVect vector2)
        {
            double xProdX = vector2.Z * vector1.Y - vector1.Z * vector2.Y;
            double xProdY = -1 * (vector2.Z * vector1.X - vector1.Z * vector2.X);
            double xProdZ = vector2.Y * vector1.X - vector1.Y * vector2.X;
            MemorySafe_CartVect xProd = new MemorySafe_CartVect(xProdX, xProdY, xProdZ);
            return xProd;
        }

        public static double DotProduct(MemorySafe_CartVect v1, MemorySafe_CartVect v2)
        {
            try
            {
                //dot product definition
                return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
            }
            catch (Exception e)
            {
                log.Debug("DotProduct Failed." +e.ToString());
                throw e;
            }
        }

        public static double DotProductMag(MemorySafe_CartVect v1, MemorySafe_CartVect v2)
        {
            try
            {
                //dot product definition
                double magv1 = VectorMagnitude(v1);
                double magv2 = VectorMagnitude(v2);
                return (v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z) / (magv1 * magv2);
            }
            catch (Exception e)
            {
                
            }
            return 0;
        }

        public static double GetAngle(MemorySafe_CartVect v1, MemorySafe_CartVect v2)
        {
            try
            {
                double dot = DotProduct(v1, v2);
                double angle = (180/Math.PI) * Math.Acos(dot / (Vector.VectorMagnitude(v1) * Vector.VectorMagnitude(v2)));
                return angle;
            }
            catch (Exception e)
            {
                log.Debug("Get Angle Failed: " + e.ToString());
                throw e;
            }
        }

        static public double FindTilt(Vector.MemorySafe_CartVect normalVector)
        {
            double calculatedTilt = -999;
            //may need to also take into account other factors that, at this stage, seem to not be important
            //building Direction of Relative North
            //zone Direction of Relative North
            //GlobalGeometryRules coordinate system
            //I may need to know this in the future then rotate the axis vectors I am making below

            //x-axis [1 0 0] points east, y-axis [0 1 0] points north, z-axis[0 0 1] points up to the sky
            //alignment with y axis means north pointing, alignment with z-axis means it is pointing up to the sky (like a flat roof)
            double nX = 0;
            double nY = 1;
            double nZ = 0;
            Vector.MemorySafe_CartVect northVector = new Vector.MemorySafe_CartVect(nX, nY, nZ);

            double uX = 0;
            double uY = 0;
            double uZ = 1;
            Vector.MemorySafe_CartVect upVector = new Vector.MemorySafe_CartVect(uX, uY, uZ);

            //rotate the axis vectors for the future

            //ensure the vector passed into the function is a unit vector
            normalVector = Vector.UnitVector(normalVector);
            //get tilt:  cross product of normal vector and upVector
            //since parallel and anti parallel vectors will return the same cross product [0,0,0] I need to filter out the antiparalll case
            if (normalVector.X == upVector.X * -1 && normalVector.Y == upVector.Y * -1 && normalVector.Z == upVector.Z * -1)
            {
                calculatedTilt = 180;
                return calculatedTilt;
            }
            else
            {
                Vector.MemorySafe_CartVect tiltVector = Vector.CrossProduct(normalVector, upVector);
                double tiltVectorMagnitude = Vector.VectorMagnitude(tiltVector);
                calculatedTilt = Math.Round(Math.Asin(tiltVectorMagnitude) * 180 / Math.PI, 2);
                return calculatedTilt;
            }
        }

        public static WalkDirection isCounterClockwise(List<MemorySafe_CartCoord> coords)
        {
            WalkDirection wd = new WalkDirection();
            wd = WalkDirection.Clockwise;

            try
            {
                log.Info("Starting CounterClockwise winding checks.");
                Vector.CartVect plRHRVect = new Vector.CartVect();
                //this list will store all of the rhr values returned by any arbitrary polyloop
                List<Vector.MemorySafe_CartVect> RHRs = new List<Vector.MemorySafe_CartVect>();

                int coordCount = coords.Count;
                for (int i = 0; i < coordCount - 2; i++)
                {
                    Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(coords[i], coords[i + 1]);
                    Vector.MemorySafe_CartVect v2 = Vector.CreateMemorySafe_Vector(coords[i + 1], coords[i + 2]);
                    Vector.MemorySafe_CartVect uv = Vector.UnitVector(Vector.CrossProduct(v1, v2));
                    RHRs.Add(uv);
                    log.Debug("Making vectors between points " + i.ToString() + ", " + (i + 1).ToString() + ", and " + (i+2).ToString());
                }
                int RHRVectorCount = RHRs.Count;
                int clockwiseCount = 0;
                //the Distinct().ToList() routine did not work because, we believe, the item in the list is not recognized by Distinct()
                //distinctRHRs = RHRs.Distinct().ToList();
                //so we took the following approach to try and find unique vectors and store them
                List<int> uniqueIndices = new List<int>();
                for(int j = 0; j < RHRVectorCount; j++)
                {

                    if (RHRs[j].X == 0 && RHRs[j].Y == 0 && RHRs[j].Z == 1)
                    {
                        //means that the vectors are not facing in the same direction
                        clockwiseCount++;
                    }
                    else
                    {
                        clockwiseCount--;
                    }

                }

                if (clockwiseCount > 0)
                {

                    return wd;
                }
                else
                {
                    wd = WalkDirection.Counterclockwise;
                    return wd;
                }
            }
            catch
            {

            }
            return wd;
        }

        public static bool isPlanar(List<MemorySafe_CartCoord> coords)
        {
            try
            {
                Dictionary<string, List<Vector.CartVect>> surfaceXProducts = new Dictionary<string, List<Vector.CartVect>>();
                List<Vector.CartVect> xProducts = new List<Vector.CartVect>();
                for (int i = 0; i < coords.Count - 2; i++)
                {
                    //Get the Cross Product
                    VectorMath.Vector.CartVect v1 = VectorMath.Vector.CreateVector(coords[i], coords[i + 1]);
                    VectorMath.Vector.CartVect v2 = VectorMath.Vector.CreateVector(coords[i + 1], coords[i + 2]);
                    Vector.CartVect xProd = Vector.CrossProduct(v1, v2);
                    if (xProd.X == 0 && xProd.Y == 0 && xProd.Z == 0)
                    {
                        log.Debug("Point "+i+" and "+(i+1)+" and "+(i+2)+" form parallel vectors.");
                        log.Debug("Cross product will not be included in the analysis.");
                        log.Info("Points form parallel vectors that are end to end.");
                    }
                    else
                    {
                        xProd = Vector.UnitVector(xProd);
                        xProducts.Add(xProd);
                    }
                    
                    
                }
                for (int j = 0; j < xProducts.Count - 1; j++)
                {
                    //parallel and anti parallel
                    if (xProducts[j].X == xProducts[j + 1].X && xProducts[j].Y == xProducts[j + 1].Y && xProducts[j].Z == xProducts[j + 1].Z)
                    {
                        continue;
                    }
                    //anti-parallel
                    else if (xProducts[j].X == -1 * xProducts[j + 1].X && xProducts[j].Y == -1 * xProducts[j + 1].Y && xProducts[j].Z == -1 * xProducts[j + 1].Z)
                    {
                        continue;
                    }
                    else if (Math.Abs(xProducts[j].X) - Math.Abs(xProducts[j + 1].X) < .0001 && Math.Abs(xProducts[j].Y) - Math.Abs(xProducts[j + 1].Y) < .0001 &&
                        Math.Abs(xProducts[j].Z) - Math.Abs(xProducts[j + 1].Z) < 0.0001)
                    {
                        continue;
                    }
                    else
                    {
                        return false;
                    }
                }
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
            return false;
        }

        public static List<MemorySafe_CartCoord> ScalePolygonByPerimDepth(List<MemorySafe_CartCoord> coordlist, double perimdepth, bool isRegular)
        {
            
            //list of all the coordinates that intersect
            List<MemorySafe_CartCoord> intlist = new List<MemorySafe_CartCoord>();
            List<Vector.MemorySafe_CartVect> origpointvec = new List<Vector.MemorySafe_CartVect>();
            List<Vector.MemorySafe_CartVect> centpointvec = new List<Vector.MemorySafe_CartVect>();
            try
            {
                //point vectors of the origninal coordinates
                origpointvec = MakePointVecList(coordlist);
                if (origpointvec.Count == 0)
                {
                    //throw an error
                }
                
                //find centroid, to create the proper vector framework
                MemorySafe_CartCoord C = GetCentroid(coordlist, isRegular);

                //create a new list of coordinates to make a vector dictionary based on the centroid
                List<MemorySafe_CartCoord> centroidcoords = new List<MemorySafe_CartCoord>();
                List<MemorySafe_CartCoord> shuffledcoords = new List<MemorySafe_CartCoord>();
                for (int i = 1; i < coordlist.Count(); i++)
                {
                    shuffledcoords.Add(coordlist[i]);
                    centroidcoords.Add(coordlist[i]);
                    centroidcoords.Add(C);
                }
                shuffledcoords.Add(coordlist[0]);
                centroidcoords.Add(coordlist[0]);
                centroidcoords.Add(C);
                centpointvec = MakeCentroidPointVecList(centroidcoords);
                if (centpointvec.Count == 0)
                {
                    //throw some sort of exception
                }

                //start the algorithm for vector translation and new point definitions
                for (int i = 0; i < origpointvec.Count(); i++)
                {
                    Vector.MemorySafe_CartVect L1 = origpointvec[i];
                    Vector.MemorySafe_CartCoord P11 = coordlist[i];
                    Vector.MemorySafe_CartCoord P12 = SumPointAndLine(P11, L1);

                    MemorySafe_CartCoord midpoint = GetMidpoint(P11, P12);
                    List<double> slope = GetSlopeComponents(P11, P12);
                    double[] normalm = new double[2];
                    normalm[0] = slope[1]*-1;
                    normalm[1] = slope[0];
                    double magn = Math.Sqrt(Math.Pow(normalm[0], 2) + Math.Pow(normalm[1], 2));
                    MemorySafe_CartVect normalunit = new MemorySafe_CartVect(normalm[0] / magn, normalm[1] / magn, 0);

                    MemorySafe_CartVect perimnormal = new MemorySafe_CartVect(normalunit.X * perimdepth, normalunit.Y * perimdepth, normalunit.Z * perimdepth);
                    Vector.MemorySafe_CartCoord transP11 = new MemorySafe_CartCoord(P11.X + perimnormal.X, P11.Y + perimnormal.Y, P11.Z + perimnormal.Z);
                    Vector.MemorySafe_CartCoord transP12 = new MemorySafe_CartCoord(P12.X + perimnormal.X, P12.Y + perimnormal.Y, P12.Z + perimnormal.Z);

                    //January 13 2014 - :
                    //

                    Vector.MemorySafe_CartVect L2 = centpointvec[i];
                    Vector.MemorySafe_CartCoord P21 = P12;
                    //alternative to above
                    //Vector.MemorySafe_CartCoord P21 = shuffledcoords[i];
                    Vector.MemorySafe_CartCoord P22 = C;

                    double[] arr1 = new double[3];
                    double A1 = transP11.Y - transP12.Y;
                    double B1 = transP12.X - transP11.X;
                    double C1 = A1 * transP11.X + B1 * transP11.Y;
                    arr1[0] = A1;
                    arr1[1] = B1;
                    arr1[2] = C1;

                    double[] arr2 = new double[3];
                    double A2 = P21.Y - P22.Y;
                    double B2 = P22.X - P21.X;
                    double C2 = A2 * P21.X + B2 * P21.Y;
                    arr2[0] = A2;
                    arr2[1] = B2;
                    arr2[2] = C2;

                    double[][] comparr = new double[2][];
                    comparr[0] = arr1;
                    comparr[1] = arr2;
                    comparr = ReducedRowEchelonForm(comparr);
                    double xmin1 = Math.Min(transP11.X, transP12.X);
                    double xmax1 = Math.Max(transP11.X, transP12.X);
                    double xmin2 = Math.Min(P21.X, P22.X);
                    double xmax2 = Math.Max(P21.X, P22.X);
                    double ymin1 = Math.Min(transP11.Y, transP12.Y);
                    double ymax1 = Math.Max(transP11.Y, transP12.Y);
                    double ymin2 = Math.Min(P21.Y, P22.Y);
                    double ymax2 = Math.Max(P21.Y, P22.Y);
                    //if it reduces to reduced row echelon form, then extract the x and y coordinates
                    if (comparr[0][0] == 1 && comparr[0][1] == 0 && comparr[1][0] == 0 && comparr[1][1] == 1)
                    {
                        //find the intersection point of this new translation with the centroidvectors
                        double xint = comparr[0][2];
                        double yint = comparr[1][2];
                        MemorySafe_CartCoord intersection = new MemorySafe_CartCoord(xint, yint, 0);
                        intlist.Add(intersection);
                    }
                    else
                    {
                        //this seems pretty bad
                    }
                }
                //reshuffle intlist
                MemorySafe_CartCoord firstcoord = intlist.Last();
                intlist.Remove(firstcoord);
                intlist.Insert(0, firstcoord);
                return intlist;

            }
            catch (Exception e)
            {
                e.ToString();
            }
            return intlist;
        }

        public static List<MemorySafe_CartCoord> ScalePolygonByFactor(List<MemorySafe_CartCoord> coordlist, double pct, bool isRegular)
        {
            List<MemorySafe_CartCoord> retlist = new List<MemorySafe_CartCoord>();
            //list of all the coordinates that intersect
            List<MemorySafe_CartCoord> intlist = new List<MemorySafe_CartCoord>();
            List<Vector.MemorySafe_CartVect> origpointvec = new List<Vector.MemorySafe_CartVect>();
            List<Vector.MemorySafe_CartVect> centpointvec = new List<Vector.MemorySafe_CartVect>();
            try
            {
                //point vectors of the origninal coordinates
                origpointvec = MakePointVecList(coordlist);
                if (origpointvec.Count == 0)
                {
                    //throw an error
                }
                //find centroid, to create the proper vector framework
                MemorySafe_CartCoord C = GetCentroid(coordlist, isRegular);

                //create a new list of coordinates to make a vector dictionary based on the centroid
                List<MemorySafe_CartCoord> centroidcoords = new List<MemorySafe_CartCoord>();
                List<MemorySafe_CartCoord> shuffledcoords = new List<MemorySafe_CartCoord>();
                for (int i = 1; i < coordlist.Count(); i++)
                {
                    shuffledcoords.Add(coordlist[i]);
                    centroidcoords.Add(coordlist[i]);
                    centroidcoords.Add(C);
                }
                shuffledcoords.Add(coordlist[0]);
                centroidcoords.Add(coordlist[0]);
                centroidcoords.Add(C);
                centpointvec = MakeCentroidPointVecList(centroidcoords);
                if (centpointvec.Count == 0)
                {
                    //throw some sort of exception
                }
                //shuffle vectorlist
                MemorySafe_CartVect firstvect = centpointvec.Last();
                centpointvec.Remove(firstvect);
                centpointvec.Insert(0, firstvect);

                //start the algorithm for new scaled polygon
                for (int i = 0; i < centpointvec.Count(); i++)
                {
                    double coreX = coordlist[i].X + centpointvec[i].X * (1 - pct);
                    double coreY = coordlist[i].Y + centpointvec[i].Y * (1 - pct);
                    double coreZ = coordlist[i].Z + centpointvec[i].Z * (1 - pct);
                    MemorySafe_CartCoord corecoord = new MemorySafe_CartCoord(coreX, coreY, coreZ);
                    retlist.Add(corecoord);
                }
                return retlist;

            }
            catch (Exception e)
            {
                e.ToString();
                throw e;
            }
            
        }

        public static List<List<Vector.MemorySafe_CartCoord>> FinishCore(List<List<MemorySafe_CartCoord>> cds, double pct, bool isRegular)
        {
            //Jan 14 - simple scaling algorithm that does not really look to verify the result or ensure it meets any complex area requirements
            List<List<List<Vector.MemorySafe_CartCoord>>> returncoords = new List<List<List<MemorySafe_CartCoord>>>();
            try
            {
                
                double grossarea = GetAreaFrom2DPolyLoop(cds[0]);
                double corearea = GetAreaFrom2DPolyLoop(cds[1]);
                double servicecorearea = grossarea * pct;
                double coreratio = servicecorearea / corearea;
                //get a new depth based upon these calculated values, then scale polygon
                //to find new depth, we actually don't know of a mathematical way to do this, so we use an iterative approach
                List<MemorySafe_CartCoord> corepoints = ScalePolygonByFactor(cds[0], pct, isRegular);
                cds.Add(corepoints);

                return cds;
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public static Dictionary<MemorySafe_LongCoord, MemorySafe_CartVect> MakeLongPointVec(Dictionary<MemorySafe_CartCoord, MemorySafe_CartVect> centpointvec)
        {
            Dictionary<MemorySafe_LongCoord, MemorySafe_CartVect> longcentdict = new Dictionary<MemorySafe_LongCoord, MemorySafe_CartVect>();
            foreach (KeyValuePair<MemorySafe_CartCoord, MemorySafe_CartVect> kp in centpointvec)
            {
                //this number is arbitrary but necessary
                long X = Convert.ToInt64(kp.Key.X * 100000000);
                long Y = Convert.ToInt64(kp.Key.Y * 100000000);
                long Z = Convert.ToInt64(kp.Key.Z * 100000000);
                MemorySafe_LongCoord lc = new MemorySafe_LongCoord(X, Y, Z);
                longcentdict[lc] = kp.Value;
            }
            return longcentdict;
        }

        public static MemorySafe_LongCoord ConvertCoordToLong(MemorySafe_CartCoord P12)
        {
            long X = Convert.ToInt64(P12.X * 100000000);
            long Y = Convert.ToInt64(P12.Y * 100000000);
            long Z = Convert.ToInt64(P12.Z * 100000000);
            MemorySafe_LongCoord lc = new MemorySafe_LongCoord(X, Y, Z);
            return lc;
        }

        public static Dictionary<MemorySafe_CartCoord, MemorySafe_CartVect> MakeCentroidPointVecDict(List<MemorySafe_CartCoord> centroidcoords)
        {
            Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec = new Dictionary<MemorySafe_CartCoord, MemorySafe_CartVect>();
            try
            {
                //turn coordinates into a point vector dictionary
                //we can assume for now that the floor plane represented by the coordinates is in the X,Y Plane
                int coordCount = centroidcoords.Count;
                for (int i = 0; i < coordCount; i++)
                {
                    MemorySafe_CartVect vec = CreateMemorySafe_Vector(centroidcoords[i], centroidcoords[i + 1]);
                    pointvec[centroidcoords[i]] = vec;
                    i++;
                }
                return pointvec;
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public static List<MemorySafe_CartVect> MakeCentroidPointVecList(List<MemorySafe_CartCoord> centroidcoords)
        {
            List<Vector.MemorySafe_CartVect> pointvec = new List<MemorySafe_CartVect>();
            try
            {
                //turn coordinates into a point vector dictionary
                //we can assume for now that the floor plane represented by the coordinates is in the X,Y Plane
                int coordCount = centroidcoords.Count;
                for (int i = 0; i < coordCount; i++)
                {
                    MemorySafe_CartVect vec = CreateMemorySafe_Vector(centroidcoords[i], centroidcoords[i + 1]);
                    pointvec.Add(vec);
                    i++;
                }
                return pointvec;
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        private static List<double> GetSlopeComponents(MemorySafe_CartCoord P11, MemorySafe_CartCoord P12)
        {
            List<double> components = new List<double>();
            components.Add(P12.X - P11.X);
            components.Add(P12.Y - P11.Y);
            return components;
        }

        public static double GetSlope(MemorySafe_CartCoord P11, MemorySafe_CartCoord P12)
        {
            return (P12.Y - P11.Y) / (P12.Z - P11.X);
        }

        public static MemorySafe_CartCoord GetMidpoint(MemorySafe_CartCoord P11, MemorySafe_CartCoord P12)
        {
            double x = (P11.X + P12.X) / 2;
            double y = (P11.Y + P12.Y) / 2;
            double z = (P11.Z + P12.Z) / 2;
            return new MemorySafe_CartCoord(x,y,z);
        }

        //rseed
        public static bool BruteForceIntersectionTest(List<Vector.MemorySafe_CartCoord> coords)
        {
            Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec = new Dictionary<MemorySafe_CartCoord, Vector.MemorySafe_CartVect>();
            try
            {
                //turn coordinates into a point vector dictionary
                //we can assume for now that the floor plane represented by the coordinates is in the X,Y Plane
                pointvec = MakePointVecDict(coords);
                if (pointvec.Count == 0) return false;
                if (CheckIntersections(pointvec, coords))
                {
                    log.Info("Self Intersection Tests Passed.");
                    return true;
                }
                else
                {
                    //this is bad
                }
                
            }
            catch (Exception e)
            {
                //log something
            }
            return false;
        }

        public static void MakeBasicTemplate()
        {

        }

        public static Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> MakePointVecDict(List<Vector.MemorySafe_CartCoord> coords)
        {
            Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec = new Dictionary<MemorySafe_CartCoord, MemorySafe_CartVect>();
            try
            {
                //turn coordinates into a point vector dictionary
                //we can assume for now that the floor plane represented by the coordinates is in the X,Y Plane
                int coordCount = coords.Count;
                for (int i = 0; i < coordCount; i++)
                {

                    if (i < coordCount - 1)
                    {
                        Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(coords[i], coords[i + 1]);
                        pointvec[coords[i]] = v1;
                        log.Debug("Vector " + i + ": (" + coords[i].X + "," + coords[i].Y + "," + coords[i].Z + ");[" + v1.X + "," + v1.Y + "," + v1.Z + "]");
                    }
                    else
                    {
                        Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(coords[i], coords[0]);
                        pointvec[coords[i]] = v1;
                        log.Debug("Vector " + i + ": (" + coords[i].X + "," + coords[i].Y + "," + coords[i].Z + ");[" + v1.X + "," + v1.Y + "," + v1.Z + "]");
                    }

                }
            }
            catch (Exception e)
            {
                log.Error(e.ToString());
            }
            return pointvec;
        }

        public static List<Vector.MemorySafe_CartVect> MakePointVecList(List<Vector.MemorySafe_CartCoord> coords)
        {
            List<Vector.MemorySafe_CartVect> pointvec = new List<MemorySafe_CartVect>();
            try
            {
                //turn coordinates into a point vector dictionary
                //we can assume for now that the floor plane represented by the coordinates is in the X,Y Plane
                int coordCount = coords.Count;
                for (int i = 0; i < coordCount; i++)
                {

                    if (i < coordCount - 1)
                    {
                        Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(coords[i], coords[i + 1]);
                        pointvec.Add(v1);
                        log.Debug("Vector " + i + ": (" + coords[i].X + "," + coords[i].Y + "," + coords[i].Z + ");[" + v1.X + "," + v1.Y + "," + v1.Z + "]");
                    }
                    else
                    {
                        Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(coords[i], coords[0]);
                        pointvec.Add(v1);
                        log.Debug("Vector " + i + ": (" + coords[i].X + "," + coords[i].Y + "," + coords[i].Z + ");[" + v1.X + "," + v1.Y + "," + v1.Z + "]");
                    }

                }
            }
            catch (Exception e)
            {
                log.Error(e.ToString());
            }
            return pointvec;
        }

        public static bool CheckIntersections(Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec, List<Vector.MemorySafe_CartCoord> coords)
        {
            try
            {
                log.Debug("Start Intersection Checking of Coordinate List.");
                int coordCount = pointvec.Count;

                for (int j = 0; j < coordCount; j++)
                {
                    if (j == 0)
                    {
                        //expect the same nearest neighbor test to hold true
                        if (CheckAdjacencyIntersections(pointvec, coords, j)) continue;
                        else return false;
                    }
                    else if (j == coordCount - 1)
                    {
                        log.Debug("Do nothing.  This is the last vector in the set and has already been checked.");
                        return true;
                    }
                    else
                    {
                        if (CheckAdjacencyIntersections(pointvec, coords, j)) continue;
                        else return false;
                    }
                }
                
            }
            catch (Exception e)
            {
                log.Error("Exception caught." + e.ToString());
            }
            return false;
        }

        //private static bool CheckAdjacencyIntersections(Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec, List<Vector.MemorySafe_CartCoord> coords, int j)
        //{
        //    List<bool> IntTruth = new List<bool>();

        //    try
        //    {
                
        //        int coordCount = pointvec.Count();
        //        log.Debug("Check Adjacency Intersections for vector " + j.ToString());

        //        for (int k = 1; k < coordCount - j; k++)
        //        {
        //            Vector.MemorySafe_CartVect L1 = pointvec[coords[j]];
        //            Vector.MemorySafe_CartCoord P1 = coords[j];

        //            //log.Debug("Testing adjacent vector " + k.ToString());
        //            Vector.MemorySafe_CartVect L2 = pointvec[coords[j + k]];
        //            Vector.MemorySafe_CartCoord P2 = coords[j + k];

        //            Vector.MemorySafe_CartVect tempvec = CreateMemorySafe_Vector(P1, P2);
        //            Vector.MemorySafe_CartVect numcross = CrossProduct(tempvec, L2);
        //            double numcrmag = VectorMagnitude(numcross);
        //            Vector.MemorySafe_CartVect dencross = CrossProduct(L1, L2);
        //            double dencrmag = VectorMagnitude(dencross);
        //            //basic checks suggested:  are numcross and dencross parallel?
        //            //is numcross not the zero vector?
        //            if (numcrmag == 0 && dencrmag == 0)
        //            {
        //                //the two lines are clearly parallel
                        
        //                //not sure if this check is needed, but we ensure that these two parallel lines are unique and do not overlap
        //                //it gets a bit more complicated from here
                        
        //                double p12x = P1.X + L1.X;
        //                double p12y = P1.Y + L1.Y;
        //                double p12z = P1.Z + L1.Z;
        //                MemorySafe_CartCoord P12 = new Vector.MemorySafe_CartCoord(p12x,p12y,p12z);

        //                double p22x = P2.X + L2.X;
        //                double p22y = P2.Y + L2.Y;
        //                double p22z = P2.Z + L2.Z;
        //                MemorySafe_CartCoord P22 = new Vector.MemorySafe_CartCoord(p22x, p22y, p22z);
        //                double diffx = P12.X - P1.X;
        //                double diffy = P12.Y - P1.Y;
        //                double p1xmin = Math.Min(P1.X, P12.X);
        //                double p1xmax = Math.Max(P1.X, P12.X);
        //                double p1ymin = Math.Min(P1.Y, P12.Y);
        //                double p1ymax = Math.Max(P1.Y, P12.Y);
        //                if (diffy == 0)
        //                {
        //                    //at this point all we need to do is check to see if X or Y of one vector fall into the range of the other
        //                    if ((P2.X >= p1xmin && P2.X <= p1xmax) || (P22.X >= p1xmin && P22.X <= p1xmax))
        //                    {
        //                        log.Info("FAIL:  Lines intersect illegally.");
        //                        return false;
        //                    }
        //                    else
        //                    {
        //                        IntTruth.Add(true);
        //                        continue;
        //                    }
        //                }
        //                else if (diffx == 0)
        //                {
        //                    //at this point all we need to do is check to see if X or Y of one vector fall into the range of the other
        //                    if ((P2.Y >= p1ymin && P2.Y <= p1ymax) || (P22.Y >= p1ymin && P22.Y <= p1ymax))
        //                    {
        //                        log.Info("FAIL:  Lines intersect illegally.");
        //                        return false;
        //                    }
        //                    else
        //                    {
        //                        IntTruth.Add(true);
        //                        continue;
        //                    }
        //                }
        //                else
        //                {
        //                    //I just need to see if the two lines intersect somewhere using common cartesian math
        //                    double slope1 = diffy / diffx;
        //                    double slope2 = (P22.Y - P2.Y) / (P22.X - P2.X);
        //                    double b1 = slope1 * P1.X - P1.Y;
        //                    double b2 = slope2 * P2.X - P2.Y;

        //                    double xint = (b2 - b1) / (slope1 - slope2);
        //                    double yint = slope1 * P1.X + b1;

        //                    if ((xint >= p1xmin && xint <= p1xmax) || (yint >= p1ymin && yint <= p1ymax))
        //                    {
        //                        log.Info("FAIL:  Lines intersect illegally.");
        //                        return false;
        //                    }
        //                    else
        //                    {
        //                        IntTruth.Add(true);
        //                        continue;
        //                    }

        //                }
        //            }
        //            double a = numcrmag / dencrmag;

        //            if (k == 1)
        //            {
        //                log.Info("Vector "+(j+k).ToString()+" is a nearest neighbor.");
        //                //we expect a to be one and P2 to be equal to the resulting coord
        //                //for the first nearest neighbor
        //                if (a == 1.0)
        //                {
        //                    log.Info("PASS:  Nearest neighbor "+(j+k).ToString()+" is intersected at its startpoint by vector "+j.ToString()+".");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "PASS");
                            
        //                    Vector.MemorySafe_CartCoord rescoord = SumPointAndLine(P1, L1);
        //                    if (rescoord.X == P2.X && rescoord.Y == P2.Y && rescoord.Z == P2.Z)
        //                    {
        //                        //and b == 0.  This truth statement is the equivalent
        //                        IntTruth.Add(true);
        //                        continue;
        //                    }
        //                }
        //                else
        //                {
        //                    log.Info("FAIL:  Nearest neighbor "+(j+k).ToString()+" is not intersected at its endpoint by vector "+j.ToString()+".  This indicates an invalid enclosure.");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "FAIL");
        //                    return false;
        //                }
        //            }
        //            else
        //            {
        //                //the last vector in the set of vector 0
        //                if (j == 0 && a == 0 && k == coordCount - 1) { log.Info("Vector " + k.ToString() + " is a nearest neighbor."); }
        //                else { log.Info("Vector " + (j+k).ToString() + " is not a nearest neighbor."); }
        //                //parallel lines
        //                if (a == double.PositiveInfinity || a == double.NegativeInfinity)
        //                {
        //                    log.Info("PASS:  Vectors "+j.ToString()+ " and " + (j+k).ToString()+" are parallel vectors.");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "PASS");
        //                    IntTruth.Add(true);
        //                }
        //                //ok when j==0, but not otherwise
        //                else if (j==0 && a == 0 && k == coordCount - 1)
        //                {
        //                    //This may take the place of the test at the last vector in the list
        //                    //I believe this only happens when the line intersects at the starting point
        //                    //we may not need 
        //                    log.Info("The vector " + (j+k).ToString() + " intersects vector " + (j).ToString() + " at its starting point.");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "PASS");
        //                    return true;
        //                }
        //                else if (j != 0 && a == 0)
        //                {
        //                    //means vectors are perpendicular to each other
        //                    //this needs to be updated, to be similar to a<1.0
        //                    log.Error("The vector " + (j + k).ToString() + " intersects vector " + j.ToString() + " at its starting point.");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "FAIL");
        //                    return false;
        //                }
        //                //ok
        //                else if (a > 1.0)
        //                {
        //                    //this means the points intersect, but some point beyond the length of the vector
        //                    log.Info("Pass.  Vectors would intersect but beyond the bounds of the vectors.");
        //                    log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "PASS");
        //                    IntTruth.Add(true);
        //                }
        //                //bad
        //                else if (a < 1.0)
        //                {
        //                    //this means they may intersect, I need to check the magnitude of the proposed intersection
        //                    //maybe, solve for a and b and report
        //                    log.Debug("a is equal to " + a);
        //                    if (a == 0) log.Debug("Vector "+j+" and "+k+"are perpendicular.");

        //                    Vector.MemorySafe_CartVect F1 = VectorTimesScalar(L1, a);
        //                    Vector.MemorySafe_CartCoord rescoord = SumPointAndLine(P1, F1);
        //                    Vector.MemorySafe_CartVect newvect = CreateMemorySafe_Vector(P2, rescoord);
        //                    //Dec 28, 2013
        //                    //this new vector may not face in the direction of L2.  Since we just draw a line from P2 to the new resulting coordinate.
        //                    //if the new vector is not parallel to L2, this is another indication, similar to case 2.
        //                    Vector.MemorySafe_CartVect L2unitvec = Vector.UnitVector(L2);
        //                    Vector.MemorySafe_CartVect newunitvec = Vector.UnitVector(newvect);
        //                    double diffx = Math.Abs(L2unitvec.X - newunitvec.X);
        //                    double diffy = Math.Abs(L2unitvec.Y - newunitvec.Y);
        //                    double diffz = Math.Abs(L2unitvec.Z - newunitvec.Z);

        //                    if (diffx <= 0.01 || diffy <= 0.01 && diffz <= 0.01)
        //                    {
        //                        double newvectmag = Vector.VectorMagnitude(newvect);
        //                        double L2mag = Vector.VectorMagnitude(L2);
        //                        if (newvectmag <= L2mag)
        //                        {
        //                            //this clearly means that the vectors intersect
        //                            log.Info("An illegal intersection in the polygon has been detected.");
        //                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
        //                            return false;
        //                        }
        //                        else
        //                        {
        //                            log.Info("PASS:  Vectors could intersect but beyond the bounds of the vectors.");
        //                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "PASS");
        //                            IntTruth.Add(true);
        //                        }
        //                    }
        //                    //Dec 28 2013
        //                    //If we arrive here, this essentially means that the vectors don't cross
        //                    else
        //                    {
        //                        log.Info("PASS:  Vectors could intersect but beyond the bounds of the vectors.");
        //                        log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
        //                        IntTruth.Add(true);
        //                    }
                            
        //                }
        //            }
        //        }
        //    }
        //    catch (Exception e)
        //    {
        //        log.Debug(e.ToString());
        //    }

        //    if (IntTruth.Contains(false))
        //    {
        //        return false;
        //    }
        //    else
        //    {
        //        return true;
        //    }
        //}

        public static Vector.MemorySafe_CartCoord FindVolumetricCentroid(List<List<Vector.MemorySafe_CartCoord>> encdef)
        {
            //initialize centroid coordinates
            double x = 0;
            double y = 0;
            double z = 0;
            try
            {
                //we first do just a very simple implementation where we find the volume's centroid as the average of all coordinates
                
                int coordcount = 0;
                foreach (List<Vector.MemorySafe_CartCoord> surface in encdef)
                {
                    
                    foreach (MemorySafe_CartCoord c in surface)
                    {
                        x += c.X;
                        y += c.Y;
                        z += c.Z;
                        coordcount++;
                    }
                }
                x = x / coordcount;
                y /= coordcount;
                z = z / coordcount;
            }
            catch (Exception e)
            {

            }
            MemorySafe_CartCoord centroid = new MemorySafe_CartCoord(x, y, z);
            return centroid;
        }

        //April 14,2014
        //Added to Accommodate named surfaces.
        public static Vector.MemorySafe_CartCoord FindVolumetricCentroid(Dictionary<string,List<Vector.MemorySafe_CartCoord>> encdef)
        {
            //initialize centroid coordinates
            double x = 0;
            double y = 0;
            double z = 0;
            try
            {
                //we first do just a very simple implementation where we find the volume's centroid as the average of all coordinates

                int coordcount = 0;
                foreach (KeyValuePair<string,List<Vector.MemorySafe_CartCoord>> surface in encdef)
                {

                    foreach (MemorySafe_CartCoord c in surface.Value)
                    {
                        x += c.X;
                        y += c.Y;
                        z += c.Z;
                        coordcount++;
                    }
                }
                x = x / coordcount;
                y /= coordcount;
                z = z / coordcount;
            }
            catch (Exception e)
            {

            }
            MemorySafe_CartCoord centroid = new MemorySafe_CartCoord(x, y, z);
            return centroid;
        }

        public static Dictionary<string, bool> SurfacesCCWound(MemorySafe_CartCoord centroid, List<List<MemorySafe_CartCoord>> encdef, int spacecount)
        {
            Dictionary<string, bool> results = new Dictionary<string, bool>();
            try
            {
                //is each surface normal vector pointing away from the centroid?
                int surfcount = 1;
                foreach (List<MemorySafe_CartCoord> surface in encdef)
                {
                    string surfid = spacecount + "#" + "shellgeometry-"+ surfcount.ToString();
                    MemorySafe_CartVect RHRVector = GetMemRHR(surface);
                    MemorySafe_CartCoord scent = GetCentroid(surface, true);
                    //make vector from centroid to surface centroid
                    MemorySafe_CartVect cent2cent = UnitVector(new MemorySafe_CartVect(scent.X - centroid.X, scent.Y - centroid.Y, scent.Z - centroid.Z));

                    //take dot product
                    double dot = DotProduct(RHRVector, cent2cent);
                    //dot product should always indicate an angle < 90 between the RHR and the vector cent2cent
                    if (dot < 0)
                    {
                        //then there is something wrong
                        log.Error("A normal is facing in the wrong direction");
                        results.Add(surfid, false);
                    }
                    else
                    {
                        results.Add(surfid, true);
                    }
                    surfcount++;
                }
            }
            catch (Exception e)
            {

            }
            return results;
        }

        //April 14, 2014
        //Chien Harriman Carmel Software Corporation
        //added this to accommodate surfaces that are already uniquely named.  e.g. - Surface or Space Boundaries
        public static Dictionary<string, bool> NamedSurfacesCCWound(MemorySafe_CartCoord centroid, Dictionary<string,List<MemorySafe_CartCoord>> encdef, int spacecount, bool fromspace)
        {
            Dictionary<string, bool> results = new Dictionary<string, bool>();
            try
            {
                //is each surface normal vector pointing away from the centroid?
                int surfcount = 1;
                foreach (KeyValuePair<string,List<MemorySafe_CartCoord>> surf in encdef)
                {
                    string surfid = "";
                    if (fromspace) { surfid = spacecount+"#"+surf.Key; }
                    else { surfid = surf.Key; }
                    MemorySafe_CartVect RHRVector = GetMemRHR(surf.Value);
                    MemorySafe_CartCoord scent = GetCentroid(surf.Value, true);
                    //make vector from centroid to surface centroid
                    MemorySafe_CartVect cent2cent = UnitVector(new MemorySafe_CartVect(scent.X - centroid.X, scent.Y - centroid.Y, scent.Z - centroid.Z));

                    //take dot product
                    double dot = DotProduct(RHRVector, cent2cent);
                    //dot product should always indicate an angle < 90 between the RHR and the vector cent2cent
                    if (dot < 0)
                    {
                        //then there is something wrong
                        log.Error("A normal is facing in the wrong direction");
                        results.Add(surfid, false);
                    }
                    else
                    {
                        results.Add(surfid, true);
                    }
                    surfcount++;
                }
            }
            catch (Exception e)
            {

            }
            return results;
        }
        //what is the point of this method
        //It takes the ring of coordinates, and does a branch tree test to ensure that only nearest neighbors intersect.  It exhaustively searches the entire chain
        private static bool CheckAdjacencyIntersections(Dictionary<Vector.MemorySafe_CartCoord, Vector.MemorySafe_CartVect> pointvec, List<Vector.MemorySafe_CartCoord> coords, int j)
        {
            List<bool> IntTruth = new List<bool>();

            try
            {

                int coordCount = pointvec.Count();
                log.Debug("Check Adjacency Intersections for vector " + j.ToString());

                for (int k = 1; k < coordCount - j; k++)
                {
                    Vector.MemorySafe_CartVect L1 = pointvec[coords[j]];
                    Vector.MemorySafe_CartCoord P11 = coords[j];
                    Vector.MemorySafe_CartCoord P12 = SumPointAndLine(P11, L1);

                    Vector.MemorySafe_CartVect L2 = pointvec[coords[j + k]];
                    Vector.MemorySafe_CartCoord P21 = coords[j + k];
                    Vector.MemorySafe_CartCoord P22 = SumPointAndLine(P21,L2);

                    double[] arr1 = new double[3];
                    double A1 = P11.Y - P12.Y;
                    double B1 = P12.X - P11.X;
                    double C1 = A1 * P11.X + B1 * P11.Y;
                    arr1[0] = A1;
                    arr1[1] = B1;
                    arr1[2] = C1;

                    double[] arr2 = new double[3];
                    double A2 = P21.Y - P22.Y;
                    double B2 = P22.X - P21.X;
                    double C2 = A2 * P21.X + B2 * P21.Y;
                    arr2[0] = A2;
                    arr2[1] = B2;
                    arr2[2] = C2;

                    double[][] comparr = new double[2][];
                    comparr[0] = arr1;
                    comparr[1] = arr2;
                    comparr = ReducedRowEchelonForm(comparr);
                    double xmin1 = Math.Min(P11.X, P12.X);
                    double xmax1 = Math.Max(P11.X, P12.X);
                    double xmin2 = Math.Min(P21.X, P22.X);
                    double xmax2 = Math.Max(P21.X, P22.X);
                    double ymin1 = Math.Min(P11.Y, P12.Y);
                    double ymax1 = Math.Max(P11.Y, P12.Y);
                    double ymin2 = Math.Min(P21.Y, P22.Y);
                    double ymax2 = Math.Max(P21.Y, P22.Y);
                    //if it reduces to reduced row echelon form, then extract the x and y coordinates
                    if (comparr[0][0] == 1 && comparr[0][1] == 0 && comparr[1][0] == 0 && comparr[1][1] == 1)
                    {
                        if (k == 1)
                        {
                            //this is supposed to happen this is good, since it is a nearest neighbor
                            log.Info("Vector " + (j + k).ToString() + " is a nearest neighbor.");
                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNY" + ";" + "PASS");
                            IntTruth.Add(true);
                            continue;
                        }
                        else if (j == 0 && k == coordCount - 1)
                        {
                            //this is also supposed to happen for since this is a nearest neighbor check.
                            log.Info("Vector " + (j + k).ToString() + " is a nearest neighbor.");
                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNY" + ";" + "PASS");
                            IntTruth.Add(true);
                            continue;
                        }
                        else
                        {
                            //being interested only in the line that we currently are checking...
                            double slope = (P12.Y - P11.Y) / (P12.X - P11.X);
                            if (slope == 0)
                            {
                                log.Debug("Vector " + j + " is horizontal.");
                                if (isInBounds(comparr[0][2], xmin1, xmin2, xmax1, xmax2) && isInBounds(comparr[1][2], ymin1, ymin2, ymax1, ymax2))
                                {
                                    log.Info("FAIL:  Illegal intersection");
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
                                    return false;
                                }
                                else
                                {
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "PASS");
                                    continue;
                                }
                            }
                            else if (slope == double.PositiveInfinity || slope == double.NegativeInfinity)
                            {
                                if (isInBounds(comparr[0][2], xmin1, xmin2, xmax1, xmax2) && isInBounds(comparr[1][2], ymin1, ymin2, ymax1, ymax2))
                                {
                                    log.Info("FAIL: Illegal intersection");
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
                                    return false;
                                }
                                else
                                {
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "PASS");
                                    continue;
                                }
                            }
                            else
                            {
                                if (!isInBounds(comparr[0][2], xmin1, xmin2, xmax1, xmax2) && !isInBounds(comparr[1][2], ymin1, ymin2, ymax1, ymax2))
                                {
                                    //this is good
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "PASS");
                                    continue;
                                }
                                else
                                {
                                    log.Info("FAIL:  Illegal intersection");
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FALSE");
                                    return false;
                                }
                            }
                            //this should not happen and we need to return false
                        }
                    }
                    else
                    {
                        if (k == 1)
                        {
                            //this should never happen, it is bad, since it is a nearest neighbor
                            log.Info("Vector " + (j + k).ToString() + " is a nearest neighbor.");
                            log.Info("FAIL: Not enclosed.");
                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNY" + ";" + "FAIL");
                            return false;
                        }
                        else if (j == 0 && k == coordCount - 1)
                        {
                            log.Info("Vector " + (j + k).ToString() + " is a nearest neighbor.");
                            log.Info("FAIL: Not enclosed.");
                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNY" + ";" + "FAIL");
                            return false;
                        }
                        else
                        {
                            //this is good
                            log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "TRUE");
                            IntTruth.Add(true);
                            continue;
                        }
                    }
                }
            }
            catch (Exception e)
            {
                log.Debug(e.ToString());
            }

            if (IntTruth.Contains(false))
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        private static bool isInBounds(double coord, double min1, double min2, double max1, double max2)
        {
            if (coord >= min1 && coord >= min2 && coord <= max1 && coord <= max2)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public static double[][] ReducedRowEchelonForm(double[][] comparr)
        {
            
            int lead = 0;
            int rowcount = comparr.Length;
            int columncount = comparr[0].Length;

            for (int r = 0; r < rowcount; r++)
            {
                //reached the end
                if (columncount <= lead) break;
                int i = r;
                //checks to see if two rows should be switched
                while (comparr[i][lead] == 0)
                {
                    i++;
                    if (i == rowcount)
                    {
                        i = r;
                        lead++;
                        if (columncount == lead)
                        {
                            //reached the end
                            lead--;
                            break;
                        }
                    }
                }
                //switches the rows if needed
                for (int j = 0; j < columncount; j++)
                {
                    double temp = comparr[r][j];
                    comparr[r][j] = comparr[i][j];
                    comparr[i][j] = temp;

                }
                //do the math to reduce to echelon form
                double div = comparr[r][lead];
                for (int j = 0; j < columncount; j++)
                {
                    comparr[r][j] /= div;
                }
                for (int j = 0; j < rowcount; j++)
                {
                    if (j != r)
                    {
                        double sub = comparr[j][lead];
                        for (int k = 0; k < columncount; k++)
                        {
                            comparr[j][k] -= sub * comparr[r][k];
                        }
                    }
                }
                lead++;
            }

            return comparr;
        }
        public static MemorySafe_CartCoord GetCentroid(List<MemorySafe_CartCoord> coordinates, bool isRegular)
        {

            try
            {
                //from Mathworld, a simple average of coordinates
                double x = 0;
                double y = 0;
                double z = 0;
                int coordcount = 0;
                foreach (MemorySafe_CartCoord c in coordinates)
                {
                    x += c.X;
                    y += c.Y;
                    z += c.Z;
                    coordcount++;
                }
                x = x / coordcount;
                y = y / coordcount;
                z = z / coordcount;

                //old, does not properly scale to 3 dimensions, as of yet
                //CartCoord centroid = new CartCoord();
                //centroid.X = 0;
                //centroid.Y = 0;
                //centroid.Z = 0;
                
                //if (isRegular)
                //{
                //    //calculate the area of the polyLoop
                //    //from Wikipedia
                //    double area = GetAreaFrom2DPolyLoop(coordinates);
                //    for (int i = 0; i < coordinates.Count - 1; i++)
                //    {
                //        centroid.X += (coordinates[i].X + coordinates[i + 1].X) * (coordinates[i].X * coordinates[i + 1].Y - coordinates[i + 1].X * coordinates[i].Y);
                //        centroid.Y += (coordinates[i].Y + coordinates[i + 1].Y) * (coordinates[i].X * coordinates[i + 1].Y - coordinates[i + 1].X * coordinates[i].Y);
                //    }
                //    centroid.X /= 6 * area;
                //    centroid.Y /= 6 * area;
                //    MemorySafe_CartCoord memC = new MemorySafe_CartCoord(centroid.X, centroid.Y, centroid.Z);
                //    return memC;
                //}
                MemorySafe_CartCoord ret = new MemorySafe_CartCoord(x, y, z);
                return ret;
            }
            catch (Exception e)
            {
                log.Debug("Fail:" + e.ToString());
            }
            MemorySafe_CartCoord memF = new MemorySafe_CartCoord(-1000, -1000, -1000);
            return memF;
            
        }

        //February 20 2013
        //Created by Chien Si Harriman Senior Product Manager for the Carmel Software Corporation
        //Currently the tool assumes that the polyloop is a valid one (counterclockwise coordinates)  Previous checks ensure this is the case?
        //and the segments of the polygon are not self-intersecting  (there are no previous tests for this as of the date above)
        public static double GetAreaFrom2DPolyLoop(List<Vector.CartCoord> coordList)
        {
            int count = coordList.Count;
            double areaprod = 0;
            bool XisZero = true;
            bool YisZero = true;
            bool ZisZero = true;
            //the following calculates the area of any irregular polygon
            foreach (Vector.CartCoord coord in coordList)
            {
                if (coord.X != 0) XisZero = false;
                if (coord.Y != 0) YisZero = false;
                if (coord.Z != 0) ZisZero = false;
            }

            //since we can only get the area of a 2-D projection, if all three of the coordinates are not zero, then we would have to
            //return an error since this is not possible
            if (!XisZero && !YisZero && !ZisZero) return -999;

            //the rest uses greens theorem
            if (XisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].Y * coordList[i + 1].Z - coordList[i].Z * coordList[i + 1].Y);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].Y * coordList[0].Z - coordList[i].Z * coordList[0].Y);
                    }
                }
                areaprod /= 2;
            }
            else if (YisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[i + 1].Z - coordList[i].Z * coordList[i + 1].X);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[0].Z - coordList[i].Z * coordList[0].X);
                    }
                }
                areaprod /= 2;

            }
            else if (ZisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[i + 1].Y - coordList[i].Y * coordList[i + 1].X);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[0].Y - coordList[i].Y * coordList[0].X);
                    }
                }
                areaprod /= 2;
            }
            return Math.Abs(areaprod);
        }

        public static double GetAreaFrom2DPolyLoop(List<Vector.MemorySafe_CartCoord> coordList)
        {
            int count = coordList.Count;
            double areaprod = 0;
            bool XisZero = true;
            bool YisZero = true;
            bool ZisZero = true;
            //the following calculates the area of any irregular polygon
            foreach (Vector.MemorySafe_CartCoord coord in coordList)
            {
                if (coord.X != 0) XisZero = false;
                if (coord.Y != 0) YisZero = false;
                if (coord.Z != 0) ZisZero = false;
            }

            //since we can only get the area of a 2-D projection, if all three of the coordinates are not zero, then we would have to
            //return an error since this is not possible
            if (!XisZero && !YisZero && !ZisZero) return -999;

            //the rest uses greens theorem
            if (XisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].Y * coordList[i + 1].Z - coordList[i].Z * coordList[i + 1].Y);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].Y * coordList[0].Z - coordList[i].Z * coordList[0].Y);
                    }
                }
                areaprod /= 2;
            }
            else if (YisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[i + 1].Z - coordList[i].Z * coordList[i + 1].X);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[0].Z - coordList[i].Z * coordList[0].X);
                    }
                }
                areaprod /= 2;

            }
            else if (ZisZero)
            {
                for (int i = 0; i < count; i++)
                {
                    if (i < count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[i + 1].Y - coordList[i].Y * coordList[i + 1].X);
                    }
                    else if (i == count - 1)
                    {
                        areaprod += (coordList[i].X * coordList[0].Y - coordList[i].Y * coordList[0].X);
                    }
                }
                areaprod /= 2;
            }
            return Math.Abs(areaprod);
        }

        public static double GetAreaofSurface(ModelingUtilities.BuildingObjects.Surface surface)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            CartVect RHRVector = GetRHR(surface.SurfaceCoords);

            //now that I have this, I can move on

            //there are two basic cases for calculating the area that we cover here, 
            //one where we get the area using greens theorem when the surface is parallel to one of the axes of the project global reference frame
            //and the second where the surface is not parallel to one of the axes of the global reference frame

            //Surface normal Parallel to global reference frame X Axis
            if (Math.Abs(RHRVector.X) == 1 && RHRVector.Y == 0 && RHRVector.Z == 0)
            {
                List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
                foreach (Vector.CartCoord coord in surface.SurfaceCoords)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    coord.X = 0;
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;


            }
            //Surface normal Parallel to global reference frame y Axis
            else if (RHRVector.X == 0 && Math.Abs(RHRVector.Y) == 1 && RHRVector.Z == 0)
            {
                List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
                foreach (Vector.CartCoord coord in surface.SurfaceCoords)
                {
                    //only take the X and Z coordinates and throw out the Y because we can assume that they are all the same
                    coord.Y = 0;
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }
            else if (RHRVector.X == 0 && RHRVector.Y == 0 && Math.Abs(RHRVector.Z) == 1)
            {
                List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
                foreach (Vector.CartCoord coord in surface.SurfaceCoords)
                {
                    //only take the X and Y coordinates and throw out the Z because we can assume that they are all the same
                    coord.Z = 0;
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }

            //the surface is not aligned with one of the reference frame axes, which requires a bit more work to determine the right answer.
            else
            {

                //New Z Axis for this plane is the normal vector already calculated, does not need to be created
                //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please
                Vector.CartVect localY = new Vector.CartVect();
                Vector.CartVect globalReferenceX = new Vector.CartVect();
                globalReferenceX.X = 1;
                globalReferenceX.Y = 0;
                globalReferenceX.Z = 0;
                localY = Vector.CrossProduct(RHRVector, globalReferenceX);
                localY = Vector.UnitVector(localY);

                //new X axis is the localY cross the surface normal vector
                Vector.CartVect localX = new Vector.CartVect();
                localX = Vector.CrossProduct(localY, RHRVector);
                localX = Vector.UnitVector(localX);

                //convert the polyloop coordinates to a local 2-D reference frame
                //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
                List<Vector.CartCoord> translatedCoordinates = new List<Vector.CartCoord>();
                //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
                Vector.CartCoord newOrigin = new Vector.CartCoord();
                newOrigin.X = 0;
                newOrigin.Y = 0;
                newOrigin.Z = 0;
                translatedCoordinates.Add(newOrigin);
                for (int j = 1; j < surface.SurfaceCoords.Count; j++)
                {
                    //randomly assigns the first polyLoop coordinate as the origin
                    Vector.CartCoord origin = surface.SurfaceCoords[0];
                    //captures the components of a vector drawn from the new origin to the 
                    Vector.CartVect distance = new Vector.CartVect();
                    distance.X = surface.SurfaceCoords[j].X - origin.X;
                    distance.Y = surface.SurfaceCoords[j].Y - origin.Y;
                    distance.Z = surface.SurfaceCoords[j].Z - origin.Z;
                    Vector.CartCoord translatedPt = new Vector.CartCoord();
                    //x coordinate is distance vector dot the new local X axis
                    translatedPt.X = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
                    //y coordinate is distance vector dot the new local Y axis
                    translatedPt.Y = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
                    translatedPt.Z = 0;
                    translatedCoordinates.Add(translatedPt);

                }
                double area = GetAreaFrom2DPolyLoop(translatedCoordinates);
                return area;
            }
            
        }

        //CHarriman Feb 14 2014 <3<3<3
        //a revision of the get area of a surface.  Used in the Rhino Grasshopper library
        public static double GetAreaofMemSafeCoords(List<Vector.MemorySafe_CartCoord> coords)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            MemorySafe_CartVect RHRVector = GetMemRHR(coords);

            //now that I have this, I can move on

            //there are two basic cases for calculating the area that we cover here, 
            //one where we get the area using greens theorem when the surface is parallel to one of the axes of the project global reference frame
            //and the second where the surface is not parallel to one of the axes of the global reference frame

            //Surface normal Parallel to global reference frame X Axis
            if (Math.Abs(RHRVector.X) == 1 && RHRVector.Y == 0 && RHRVector.Z == 0)
            {
                List<CartCoord> coordList = new List<CartCoord>();
                foreach (MemorySafe_CartCoord coord in coords)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    CartCoord tempcoord = new CartCoord();
                    tempcoord.X = 0;
                    tempcoord.Y = coord.Y;
                    tempcoord.Z = coord.Z;
                    coordList.Add(tempcoord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;


            }
            //Surface normal Parallel to global reference frame y Axis
            else if (RHRVector.X == 0 && Math.Abs(RHRVector.Y) == 1 && RHRVector.Z == 0)
            {
                List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
                foreach (MemorySafe_CartCoord coord in coords)
                {
                    //only take the X and Z coordinates and throw out the Y because we can assume that they are all the same
                    CartCoord tempcoord = new CartCoord();
                    tempcoord.X = coord.X;
                    tempcoord.Y = 0;
                    tempcoord.Z = coord.Z;
                    coordList.Add(tempcoord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }
            else if (RHRVector.X == 0 && RHRVector.Y == 0 && Math.Abs(RHRVector.Z) == 1)
            {
                List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
                foreach (MemorySafe_CartCoord coord in coords)
                {
                    //only take the X and Y coordinates and throw out the Z because we can assume that they are all the same
                    CartCoord tempcoord = new CartCoord();
                    tempcoord.X = coord.X;
                    tempcoord.Y = coord.Y;
                    tempcoord.Z = 0;
                    coordList.Add(tempcoord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }

            //the surface is not aligned with one of the reference frame axes, which requires a bit more work to determine the right answer.
            else
            {

                //New Z Axis for this plane is the normal vector already calculated, does not need to be created
                //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please

                CartVect globalReferenceX = new CartVect();
                globalReferenceX.X = 1;
                globalReferenceX.Y = 0;
                globalReferenceX.Z = 0;

                CartVect localY = CrossProductNVRetMSNV(RHRVector, globalReferenceX);
                localY = UnitVector(localY);

                
                //new X axis is the localY cross the surface normal vector
                Vector.CartVect localX = new Vector.CartVect();

                localX = CrossProductNVRetNVMS(localY, RHRVector);
                localX = Vector.UnitVector(localX);

                //convert the polyloop coordinates to a local 2-D reference frame
                //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
                List<Vector.CartCoord> translatedCoordinates = new List<Vector.CartCoord>();
                //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
                Vector.CartCoord newOrigin = new Vector.CartCoord();
                newOrigin.X = 0;
                newOrigin.Y = 0;
                newOrigin.Z = 0;
                translatedCoordinates.Add(newOrigin);
                for (int j = 1; j < coords.Count; j++)
                {
                    //randomly assigns the first polyLoop coordinate as the origin
                    Vector.CartCoord origin = new CartCoord();
                    origin.X = coords[0].X;
                    origin.Y = coords[0].Y;
                    origin.Z = coords[0].Z;
                    //captures the components of a vector drawn from the new origin to the 
                    Vector.CartVect distance = new Vector.CartVect();
                    
                    distance.X = coords[j].X - origin.X;
                    distance.Y = coords[j].Y - origin.Y;
                    distance.Z = coords[j].Z - origin.Z;
                    Vector.CartCoord translatedPt = new Vector.CartCoord();
                    //x coordinate is distance vector dot the new local X axis
                    translatedPt.X = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
                    //y coordinate is distance vector dot the new local Y axis
                    translatedPt.Y = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
                    translatedPt.Z = 0;
                    translatedCoordinates.Add(translatedPt);

                }
                double area = GetAreaFrom2DPolyLoop(translatedCoordinates);
                return area;
            }

        }

        //public static double GetAreaofPolyList(List<Vector.MemorySafe_CartCoord> coordlist)
        //{
        //    //Used to figure out how best to calculate the area from a given surfacce. 
        //    //Get the coordinates that define the surface
        //    //get the area based on the coordinates

        //    //Get the RHRVector (the actual direction is not important
        //    MemorySafe_CartVect RHRVector = GetMemRHR(coordlist);

        //    //now that I have this, I can move on

        //    //there are two basic cases for calculating the area that we cover here, 
        //    //one where we get the area using greens theorem when the surface is parallel to one of the axes of the project global reference frame
        //    //and the second where the surface is not parallel to one of the axes of the global reference frame

        //    //Surface normal Parallel to global reference frame X Axis
        //    if (Math.Abs(RHRVector.X) == 1 && RHRVector.Y == 0 && RHRVector.Z == 0)
        //    {
        //        List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
        //        foreach (Vector.CartCoord coord in coordList)
        //        {
        //            //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
        //            coord.X = 0;
        //            double Y = surface.SurfaceCoords[i].Y;
        //            double Z = surface.SurfaceCoords[i].Z;
        //            Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
        //            coordList.Add(coord);

        //        }
        //        double area = GetAreaFrom2DPolyLoop(coordList);
        //        return area;


        //    }
        //    //Surface normal Parallel to global reference frame y Axis
        //    else if (RHRVector.X == 0 && Math.Abs(RHRVector.Y) == 1 && RHRVector.Z == 0)
        //    {
        //        List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
        //        foreach (Vector.CartCoord coord in surface.SurfaceCoords)
        //        {
        //            //only take the X and Z coordinates and throw out the Y because we can assume that they are all the same
        //            coord.Y = 0;
        //            coordList.Add(coord);

        //        }
        //        double area = GetAreaFrom2DPolyLoop(coordList);
        //        return area;
        //    }
        //    else if (RHRVector.X == 0 && RHRVector.Y == 0 && Math.Abs(RHRVector.Z) == 1)
        //    {
        //        List<Vector.CartCoord> coordList = new List<Vector.CartCoord>();
        //        foreach (Vector.CartCoord coord in surface.SurfaceCoords)
        //        {
        //            //only take the X and Y coordinates and throw out the Z because we can assume that they are all the same
        //            coord.Z = 0;
        //            coordList.Add(coord);

        //        }
        //        double area = GetAreaFrom2DPolyLoop(coordList);
        //        return area;
        //    }

        //    //the surface is not aligned with one of the reference frame axes, which requires a bit more work to determine the right answer.
        //    else
        //    {

        //        //New Z Axis for this plane is the normal vector already calculated, does not need to be created
        //        //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please
        //        Vector.CartVect localY = new Vector.CartVect();
        //        Vector.CartVect globalReferenceX = new Vector.CartVect();
        //        globalReferenceX.X = 1;
        //        globalReferenceX.Y = 0;
        //        globalReferenceX.Z = 0;
        //        localY = Vector.CrossProduct(RHRVector, globalReferenceX);
        //        localY = Vector.UnitVector(localY);

        //        //new X axis is the localY cross the surface normal vector
        //        Vector.CartVect localX = new Vector.CartVect();
        //        localX = Vector.CrossProduct(localY, RHRVector);
        //        localX = Vector.UnitVector(localX);

        //        //convert the polyloop coordinates to a local 2-D reference frame
        //        //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
        //        List<Vector.CartCoord> translatedCoordinates = new List<Vector.CartCoord>();
        //        //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
        //        Vector.CartCoord newOrigin = new Vector.CartCoord();
        //        newOrigin.X = 0;
        //        newOrigin.Y = 0;
        //        newOrigin.Z = 0;
        //        translatedCoordinates.Add(newOrigin);
        //        for (int j = 1; j < surface.SurfaceCoords.Count; j++)
        //        {
        //            //randomly assigns the first polyLoop coordinate as the origin
        //            Vector.CartCoord origin = surface.SurfaceCoords[0];
        //            //captures the components of a vector drawn from the new origin to the 
        //            Vector.CartVect distance = new Vector.CartVect();
        //            distance.X = surface.SurfaceCoords[j].X - origin.X;
        //            distance.Y = surface.SurfaceCoords[j].Y - origin.Y;
        //            distance.Z = surface.SurfaceCoords[j].Z - origin.Z;
        //            Vector.CartCoord translatedPt = new Vector.CartCoord();
        //            //x coordinate is distance vector dot the new local X axis
        //            translatedPt.X = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
        //            //y coordinate is distance vector dot the new local Y axis
        //            translatedPt.Y = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
        //            translatedPt.Z = 0;
        //            translatedCoordinates.Add(translatedPt);

        //        }
        //        double area = GetAreaFrom2DPolyLoop(translatedCoordinates);
        //        return area;
        //    }

        //}

        public static double GetAreaofSurface(ModelingUtilities.BuildingObjects.MemorySafe_Surface surface)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            MemorySafe_CartVect RHRVector = GetMemRHR(surface.SurfaceCoords);

            //now that I have this, I can move on

            //there are two basic cases for calculating the area that we cover here, 
            //one where we get the area using greens theorem when the surface is parallel to one of the axes of the project global reference frame
            //and the second where the surface is not parallel to one of the axes of the global reference frame

            //Surface normal Parallel to global reference frame X Axis
            if (Math.Abs(RHRVector.X) == 1 && RHRVector.Y == 0 && RHRVector.Z == 0)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for(int i = 0; i< surface.SurfaceCoords.Count; i++)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    double X = 0;
                    double Y = surface.SurfaceCoords[i].Y;
                    double Z = surface.SurfaceCoords[i].Z;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;


            }
            //Surface normal Parallel to global reference frame y Axis
            else if (RHRVector.X == 0 && Math.Abs(RHRVector.Y) == 1 && RHRVector.Z == 0)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for (int i = 0; i < surface.SurfaceCoords.Count; i++ )
                {
                    //only take the X and Z coordinates and throw out the Y because we can assume that they are all the same
                    double X = surface.SurfaceCoords[i].X;
                    double Y = 0;
                    double Z = surface.SurfaceCoords[i].Z;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }
            else if (RHRVector.X == 0 && RHRVector.Y == 0 && Math.Abs(RHRVector.Z) == 1)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for (int i = 0; i < surface.SurfaceCoords.Count; i++)
                {
                    //only take the X and Y coordinates and throw out the Z because we can assume that they are all the same
                    double X = surface.SurfaceCoords[i].X;
                    double Y = surface.SurfaceCoords[i].Y;
                    double Z = 0;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);
                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }

            //the surface is not aligned with one of the reference frame axes, which requires a bit more work to determine the right answer.
            else
            {

                //New Z Axis for this plane is the normal vector already calculated, does not need to be created
                //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please
                double X = 1;
                double Y = 0;
                double Z = 0;
                Vector.MemorySafe_CartVect globalReferenceX = new Vector.MemorySafe_CartVect(X,Y,Z);

                Vector.MemorySafe_CartVect localY = Vector.CrossProduct(RHRVector, globalReferenceX);
                localY = Vector.UnitVector(localY);

                //new X axis is the localY cross the surface normal vector
                Vector.MemorySafe_CartVect localX = Vector.CrossProduct(localY, RHRVector);
                localX = Vector.UnitVector(localX);

                //convert the polyloop coordinates to a local 2-D reference frame
                //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
                List<Vector.MemorySafe_CartCoord> translatedCoordinates = new List<Vector.MemorySafe_CartCoord>();
                //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
                double originX = 0;
                double originY = 0;
                double originZ = 0;
                Vector.MemorySafe_CartCoord newOrigin = new Vector.MemorySafe_CartCoord(originX, originY, originZ);
                translatedCoordinates.Add(newOrigin);
                for (int j = 1; j < surface.SurfaceCoords.Count; j++)
                {
                    //randomly assigns the first polyLoop coordinate as the origin
                    Vector.MemorySafe_CartCoord origin = surface.SurfaceCoords[0];
                    //captures the components of a vector drawn from the new origin to the 
                    double xDistance = surface.SurfaceCoords[j].X - origin.X;
                    double yDist = surface.SurfaceCoords[j].Y - origin.Y;
                    double zDist = surface.SurfaceCoords[j].Z - origin.Z;
                    Vector.MemorySafe_CartVect distance = new Vector.MemorySafe_CartVect(xDistance, yDist, zDist);
                    double translPtX = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
                    double translPtY = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
                    double translPtZ = 0;
                    Vector.MemorySafe_CartCoord translatedPt = new Vector.MemorySafe_CartCoord(translPtX, translPtY, translPtZ);
                    translatedCoordinates.Add(translatedPt);

                }
                double area = GetAreaFrom2DPolyLoop(translatedCoordinates);
                return area;
            }

        }

        public static List<Vector.MemorySafe_CartCoord> flattenPolyLoop(List<Vector.MemorySafe_CartCoord> coordlist)
        {

            MemorySafe_CartVect RHRVector = GetMemRHR(coordlist);
            //New Z Axis for this plane is the normal vector already calculated, does not need to be created
                //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please
                double X = 1;
                double Y = 0;
                double Z = 0;
                Vector.MemorySafe_CartVect globalReferenceX = new Vector.MemorySafe_CartVect(X,Y,Z);

                Vector.MemorySafe_CartVect localY = Vector.CrossProduct(RHRVector, globalReferenceX);
                localY = Vector.UnitVector(localY);

                //new X axis is the localY cross the surface normal vector
                Vector.MemorySafe_CartVect localX = Vector.CrossProduct(localY, RHRVector);
                localX = Vector.UnitVector(localX);

                //convert the polyloop coordinates to a local 2-D reference frame
                //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
                List<Vector.MemorySafe_CartCoord> translatedCoordinates = new List<Vector.MemorySafe_CartCoord>();
                //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
                double originX = 0;
                double originY = 0;
                double originZ = 0;
                Vector.MemorySafe_CartCoord newOrigin = new Vector.MemorySafe_CartCoord(originX, originY, originZ);
                translatedCoordinates.Add(newOrigin);
                for (int j = 1; j < coordlist.Count; j++)
                {
                    //randomly assigns the first polyLoop coordinate as the origin
                    Vector.MemorySafe_CartCoord origin = coordlist[0];
                    //captures the components of a vector drawn from the new origin to the 
                    double xDistance = coordlist[j].X - origin.X;
                    double yDist = coordlist[j].Y - origin.Y;
                    double zDist = coordlist[j].Z - origin.Z;
                    Vector.MemorySafe_CartVect distance = new Vector.MemorySafe_CartVect(xDistance, yDist, zDist);
                    double translPtX = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
                    double translPtY = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
                    double translPtZ = 0;
                    Vector.MemorySafe_CartCoord translatedPt = new Vector.MemorySafe_CartCoord(translPtX, translPtY, translPtZ);
                    translatedCoordinates.Add(translatedPt);
                }
                return translatedCoordinates;
        }

        public static double GetAreaofWindow(ModelingUtilities.BuildingObjects.MemorySafe_OpeningDefinitions opening)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            MemorySafe_CartVect RHRVector = GetMemRHR(opening.coordinateList);

            //now that I have this, I can move on

            //there are two basic cases for calculating the area that we cover here, 
            //one where we get the area using greens theorem when the surface is parallel to one of the axes of the project global reference frame
            //and the second where the surface is not parallel to one of the axes of the global reference frame

            //Surface normal Parallel to global reference frame X Axis
            if (Math.Abs(RHRVector.X) == 1 && RHRVector.Y == 0 && RHRVector.Z == 0)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for (int i = 0; i < opening.coordinateList.Count; i++)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    double X = 0;
                    double Y = opening.coordinateList[i].Y;
                    double Z = opening.coordinateList[i].Z;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;


            }
            //Surface normal Parallel to global reference frame y Axis
            else if (RHRVector.X == 0 && Math.Abs(RHRVector.Y) == 1 && RHRVector.Z == 0)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for (int i = 0; i < opening.coordinateList.Count; i++)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    double X = opening.coordinateList[i].X;
                    double Y = 0;
                    double Z = opening.coordinateList[i].Z;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }
            else if (RHRVector.X == 0 && RHRVector.Y == 0 && Math.Abs(RHRVector.Z) == 1)
            {
                List<Vector.MemorySafe_CartCoord> coordList = new List<Vector.MemorySafe_CartCoord>();
                for (int i = 0; i < opening.coordinateList.Count; i++)
                {
                    //only take the Y and Z coordinates and throw out the X because we can assume that they are all the same
                    double X = opening.coordinateList[i].X ;
                    double Y = opening.coordinateList[i].Y;
                    double Z = 0;
                    Vector.MemorySafe_CartCoord coord = new MemorySafe_CartCoord(X, Y, Z);
                    coordList.Add(coord);

                }
                double area = GetAreaFrom2DPolyLoop(coordList);
                return area;
            }

            //the surface is not aligned with one of the reference frame axes, which requires a bit more work to determine the right answer.
            else
            {

                //New Z Axis for this plane is the normal vector already calculated, does not need to be created
                //Get New Y Axis which is the surface Normal Vector cross the original global reference X unit vector (all unit vectors please
                double X = 1;
                double Y = 0;
                double Z = 0;
                Vector.MemorySafe_CartVect globalReferenceX = new Vector.MemorySafe_CartVect(X, Y, Z);

                Vector.MemorySafe_CartVect localY = Vector.CrossProduct(RHRVector, globalReferenceX);
                localY = Vector.UnitVector(localY);

                //new X axis is the localY cross the surface normal vector
                Vector.MemorySafe_CartVect localX = Vector.CrossProduct(localY, RHRVector);
                localX = Vector.UnitVector(localX);

                //convert the polyloop coordinates to a local 2-D reference frame
                //using a trick employed by video game programmers found here http://stackoverflow.com/questions/1023948/rotate-normal-vector-onto-axis-plane
                List<Vector.MemorySafe_CartCoord> translatedCoordinates = new List<Vector.MemorySafe_CartCoord>();
                //put the origin in place in these translated coordinates since our loop skips over this first arbitrary point
                double originX = 0;
                double originY = 0;
                double originZ = 0;
                Vector.MemorySafe_CartCoord newOrigin = new Vector.MemorySafe_CartCoord(originX, originY, originZ);
                translatedCoordinates.Add(newOrigin);
                for (int j = 1; j < opening.coordinateList.Count; j++)
                {
                    //randomly assigns the first polyLoop coordinate as the origin
                    Vector.MemorySafe_CartCoord origin = opening.coordinateList[0];
                    //captures the components of a vector drawn from the new origin to the 
                    double xDistance = opening.coordinateList[j].X - origin.X;
                    double yDist = opening.coordinateList[j].Y - origin.Y;
                    double zDist = opening.coordinateList[j].Z - origin.Z;
                    Vector.MemorySafe_CartVect distance = new Vector.MemorySafe_CartVect(xDistance, yDist, zDist);
                    double translPtX = distance.X * localX.X + distance.Y * localX.Y + distance.Z * localX.Z;
                    double translPtY = distance.X * localY.X + distance.Y * localY.Y + distance.Z * localY.Z;
                    double translPtZ = 0;
                    Vector.MemorySafe_CartCoord translatedPt = new Vector.MemorySafe_CartCoord(translPtX, translPtY, translPtZ);
                    translatedCoordinates.Add(translatedPt);

                }
                double area = GetAreaFrom2DPolyLoop(translatedCoordinates);
                return area;
            }

        }  
                    
        


        public static Vector.CartVect GetRHR(List<Vector.CartCoord> plCoords)
        {
            Vector.CartVect plRHRVect = new Vector.CartVect();
            //this list will store all of the rhr values returned by any arbitrary polyloop
            List<Vector.CartVect> RHRs = new List<Vector.CartVect>();

            int coordCount = plCoords.Count;
            for (int i = 0; i < coordCount - 2; i++)
            {
                Vector.CartVect v1 = Vector.CreateVector(plCoords[i], plCoords[i + 1]);
                Vector.CartVect v2 = Vector.CreateVector(plCoords[i + 1], plCoords[i + 2]);
                Vector.CartVect uv = Vector.UnitVector(Vector.CrossProduct(v1, v2));
                RHRs.Add(uv);
            }
            int RHRVectorCount = RHRs.Count;
            List<Vector.CartVect> distinctRHRs = new List<Vector.CartVect>();
            int parallelCount = 0;
            int antiParallelCount = 0;
            //the Distinct().ToList() routine did not work because, we believe, the item in the list is not recognized by Distinct()
            //distinctRHRs = RHRs.Distinct().ToList();
            //so we took the following approach to try and find unique vectors and store them
            distinctRHRs.Add(RHRs[0]);
            List<int> uniqueIndices = new List<int>();
            for (int j = 1; j < RHRVectorCount; j++)
            {

                if (RHRs[j].X == distinctRHRs[0].X*-1 && RHRs[j].Y == distinctRHRs[0].Y*-1 && RHRs[j].Z == distinctRHRs[0].Z*-1)
                {
                    //means that the vectors are not facing in the same direction
                    antiParallelCount++;
                }
                else
                {
                    parallelCount++;
                }
                
            }

            if (antiParallelCount > parallelCount)
            {
                Vector.CartVect antiParallel = new CartVect { };
                antiParallel.X = distinctRHRs[0].X * -1;
                antiParallel.Y = distinctRHRs[0].Y * -1;
                antiParallel.Z = distinctRHRs[0].Z* -1;

                return antiParallel;
            }
            else
            {
                return distinctRHRs[0];
            }
        }

        public static Vector.MemorySafe_CartVect GetMemRHR(List<Vector.MemorySafe_CartCoord> plCoords)
        {
            Vector.CartVect plRHRVect = new Vector.CartVect();
            //this list will store all of the rhr values returned by any arbitrary polyloop
            List<Vector.MemorySafe_CartVect> RHRs = new List<Vector.MemorySafe_CartVect>();

            int coordCount = plCoords.Count;
            for (int i = 0; i < coordCount - 2; i++)
            {
                Vector.MemorySafe_CartVect v1 = Vector.CreateMemorySafe_Vector(plCoords[i], plCoords[i + 1]);
                Vector.MemorySafe_CartVect v2 = Vector.CreateMemorySafe_Vector(plCoords[i + 1], plCoords[i + 2]);
                Vector.MemorySafe_CartVect uv = Vector.UnitVector(Vector.CrossProduct(v1, v2));
                RHRs.Add(uv);
            }
            int RHRVectorCount = RHRs.Count;
            List<Vector.MemorySafe_CartVect> distinctRHRs = new List<Vector.MemorySafe_CartVect>();
            int parallelCount = 0;
            int antiParallelCount = 0;
            //the Distinct().ToList() routine did not work because, we believe, the item in the list is not recognized by Distinct()
            //distinctRHRs = RHRs.Distinct().ToList();
            //so we took the following approach to try and find unique vectors and store them
            distinctRHRs.Add(RHRs[0]);
            List<int> uniqueIndices = new List<int>();
            for (int j = 1; j < RHRVectorCount; j++)
            {

                if (RHRs[j].X == distinctRHRs[0].X * -1 && RHRs[j].Y == distinctRHRs[0].Y * -1 && RHRs[j].Z == distinctRHRs[0].Z * -1)
                {
                    //means that the vectors are not facing in the same direction
                    antiParallelCount++;
                }
                else
                {
                    parallelCount++;
                }

            }

            if (antiParallelCount > parallelCount)
            {
                double X = distinctRHRs[0].X * -1;
                double Y = distinctRHRs[0].Y * -1;
                double Z = distinctRHRs[0].Z * -1;
                Vector.MemorySafe_CartVect antiParallel = new Vector.MemorySafe_CartVect(X,Y,Z);
                return antiParallel;
            }
            else
            {
                return distinctRHRs[0];
            }
        }

        private static bool IsSurfaceRegular(ModelingUtilities.BuildingObjects.Surface Surface)
        {
            //tests to see if all candidate surfaces and the standard surface are regular (rectangular polygons)

            bool isRegularPolygon = true;
            //see if the standard surface has four coordinates defining its polyloop (one marker of a rectangle)
            int standSurfaceCoordinateCount = Surface.SurfaceCoords.Count;
            if (standSurfaceCoordinateCount == 4)
            {
                //check the two potentially parallel sides, to ensure they are indeed parallel
                Vector.CartVect v1 = Vector.CreateVector(Surface.SurfaceCoords[0], Surface.SurfaceCoords[1]);
                Vector.CartVect v2 = Vector.CreateVector(Surface.SurfaceCoords[2], Surface.SurfaceCoords[3]);
                Vector.CartVect v1xv2 = Vector.CrossProduct(v1, v2);
                v1xv2 = Vector.UnitVector(v1xv2);
                double magnitudev1xv2 = Vector.VectorMagnitude(v1xv2);
                Vector.CartVect v3 = Vector.CreateVector(Surface.SurfaceCoords[1], Surface.SurfaceCoords[2]);
                Vector.CartVect v4 = Vector.CreateVector(Surface.SurfaceCoords[3], Surface.SurfaceCoords[0]);
                Vector.CartVect v3xv4 = Vector.CrossProduct(v3, v4);
                v3xv4 = Vector.UnitVector(v3xv4);
                double magnitudev3xv4 = Vector.VectorMagnitude(v3xv4);
                //the unit vector will not be a number NaN if the Cross product detects a zero vector (indicating parallel vectors)
                if (double.IsNaN(magnitudev1xv2) && double.IsNaN(magnitudev3xv4))
                {
                    isRegularPolygon = true;
                }
                else
                {
                    isRegularPolygon = false;
                }
            }
            else
            {
                //might as well stop here because 
                isRegularPolygon = false;
                return isRegularPolygon;
            }
            return isRegularPolygon;

        }

        //surfaceID is the name of the surface
        //dictionary holds all of the edges that are related
        //coords holds the coordinates of the given surface
        //coordtol is the allowable tolerance when comparing coordinates of surfaces
        //vectortol is the allowable tolerance when trying to determine alignment of edges (usually a cross product or dot product)
        public static Dictionary<int, EdgeFamily> GetEdgeFamilies(string surfaceId, Dictionary<int, EdgeFamily> uniqueedges, List<Vector.MemorySafe_CartCoord> coords, double coordtol, double dotprodtol)
        {
            //for any given surface group that is made up of a series of coordinates
            try
            {
                bool perfectMatch = false;
                int coordcount = coords.Count;

                for (int i = 0; i < coordcount; i++)
                {
                    int uniqueedgect = uniqueedges.Count();
                    //initialize the edge being tested, the test edge
                    EdgeFamily edge = new EdgeFamily();
                    edge.sbdec = surfaceId;
                    edge.relatedEdges = new List<EdgeFamily>();
                    edge.startendpt = new List<Vector.MemorySafe_CartCoord>();
                    if (uniqueedges.Count == 0)
                    {
                        uniqueedges[uniqueedgect] = edge;
                        edge.startendpt.Add(coords[i]);
                        edge.startendpt.Add(coords[i + 1]);
                        uniqueedgect++;
                        continue;

                    }
                    //the first one is easy becaues it will always be a uniqueedge
                    #region
                    if (i < coordcount - 1)
                    {
                        edge.startendpt.Add(coords[i]);
                        edge.startendpt.Add(coords[i + 1]);
                    }
                    else
                    {
                        edge.startendpt.Add(coords[i]);
                        edge.startendpt.Add(coords[0]);
                    }
                    //search through existing edges to try and find a perfect match
                    int edgeloopcounter = 0; //keeps track of how many guest edges in the unique edge dictionary I've searched through
                    foreach (KeyValuePair<int, EdgeFamily> kp in uniqueedges)
                    {
                        Vector.MemorySafe_CartCoord gueststartpt = kp.Value.startendpt[0];
                        Vector.MemorySafe_CartCoord guestendpt = kp.Value.startendpt[1];
                        double diffx1 = gueststartpt.X - edge.startendpt[0].X;
                        double diffy1 = gueststartpt.Y - edge.startendpt[0].Y;
                        double diffz1 = gueststartpt.Z - edge.startendpt[0].Z;
                        double diffx2 = gueststartpt.X - edge.startendpt[1].X;
                        double diffy2 = gueststartpt.Y - edge.startendpt[1].Y;
                        double diffz2 = gueststartpt.Z - edge.startendpt[1].Z;
                        double diffx3 = guestendpt.X - edge.startendpt[0].X;
                        double diffy3 = guestendpt.Y - edge.startendpt[0].Y;
                        double diffz3 = guestendpt.Z - edge.startendpt[0].Z;
                        double diffx4 = guestendpt.X - edge.startendpt[1].X;
                        double diffy4 = guestendpt.Y - edge.startendpt[1].Y;
                        double diffz4 = guestendpt.Z - edge.startendpt[1].Z;
                        
                        //In the unique edge dictionary, I have located at least one point that is similar to my test edge.
                        if (Math.Abs(diffx1) < coordtol && Math.Abs(diffy1) < coordtol && Math.Abs(diffz1) < coordtol)
                        {
                            //diffx1 = guestendpt.X - edge.startendpt[1].X;
                            //diffy1 = guestendpt.Y - edge.startendpt[1].Y;
                            //diffz1 = guestendpt.Z - edge.startendpt[1].Z;
                            if (Math.Abs(diffx4) < coordtol && Math.Abs(diffy4) < coordtol && Math.Abs(diffz4) < coordtol)
                            {
                                //both match, means the match is perfect, so the unique edge has found it complement related edge.  Great!
                                kp.Value.relatedEdges.Add(edge);
                                //I am done searching this test edge, and I can start over again with the next edge in question
                                perfectMatch = true;
                                break;

                            }

                            else
                            {
                                //so far, I have found only one thing in common, sharing of one point, even though second point did not match, the edges could still align
                                //draw vector A
                                double Ax = guestendpt.X - edge.startendpt[1].X;
                                double Ay = guestendpt.Y - edge.startendpt[1].Y;
                                double Az = guestendpt.Z - edge.startendpt[1].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                double Bx = gueststartpt.X - edge.startendpt[1].X;
                                double By = gueststartpt.Y - edge.startendpt[1].Y;
                                double Bz = gueststartpt.Z - edge.startendpt[1].Z;
                                Vector.MemorySafe_CartVect B = new MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);

                                //this is the easiest intersection because A and B prove to be inside the guest edge
                                //these vectors are antiparallel in this case
                                double dotAB = DotProductMag(A, B);
                                if (dotAB < 0 && Math.Abs(Math.Abs(dotAB) - 1) < dotprodtol)
                                {
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //in this case the point is outside of the guest edge, because both vectors prove to be parallel
                                else if (Math.Abs(dotAB-1) < dotprodtol)
                                {
                                    //the question is whether it completely overlaps the guest edge, or does not intersect it
                                    double ex = edge.startendpt[1].X - edge.startendpt[0].X;
                                    double ey = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                    double ez = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                    //take cross product to see if they are parallel
                                    double evX = guestendpt.X - gueststartpt.X;
                                    double evY = guestendpt.Y - gueststartpt.Y;
                                    double evZ = guestendpt.Z - gueststartpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    double dot = Vector.DotProductMag(A, ev);
                                    if (dot < 0 && Math.Abs(Math.Abs(dot) - 1) < coordtol)
                                    {
                                        //antiparallel is an indication that they overlap
                                        kp.Value.relatedEdges.Add(edge);
                                        //the edge is still unique but accumulates a neighbor
                                        edge.relatedEdges.Add(kp.Value);
                                        edgeloopcounter++;
                                        continue;
                                    }
                                    else
                                    {
                                        //they don't intersect
                                        edgeloopcounter++;
                                        continue;
                                    }

                                }

                            }
                        }
                        else if (Math.Abs(diffx2) < coordtol && Math.Abs(diffy2) < coordtol && Math.Abs(diffz2) < coordtol)
                        {
                            if (Math.Abs(diffx3) <= coordtol && Math.Abs(diffy3) <= coordtol && Math.Abs(diffz3) <= coordtol)
                            {
                                //both match, means the match is perfect, so add it to the related surfaces list
                                kp.Value.relatedEdges.Add(edge);
                                perfectMatch = true;
                                break;

                            }
                            else
                            {
                                //the edge may be unique, though it could still have neighboring relationships
                                double Ax = guestendpt.X - edge.startendpt[0].X;
                                double Ay = guestendpt.Y - edge.startendpt[0].Y;
                                double Az = guestendpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                double Bx = gueststartpt.X - edge.startendpt[0].X;
                                double By = gueststartpt.Y - edge.startendpt[0].Y;
                                double Bz = gueststartpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect B = new MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);

                                //this is the easiest intersection because A and B prove to be inside the guest edge
                                //these vectors are antiparallel in this case
                                double dotAB = DotProductMag(A, B);
                                if (dotAB<0 && Math.Abs(Math.Abs(dotAB) - 1) < dotprodtol)
                                {
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //in this case the point is outside of the guest edge, because both vectors prove to be parallel
                                else if (Math.Abs(dotAB-1) < dotprodtol)
                                {
                                    //the question is whether it completely overlaps the guest edge, or does not intersect it
                                    double ex = edge.startendpt[1].X - edge.startendpt[0].X;
                                    double ey = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                    double ez = edge.startendpt[1].Z - edge.startendpt[0].Z;


                                    //take cross product to see if they are antiparallel
                                    double evX = guestendpt.X - gueststartpt.X;
                                    double evY = guestendpt.Y - gueststartpt.Y;
                                    double evZ = guestendpt.Z - gueststartpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    double dot = Vector.DotProductMag(A, ev);
                                    //in this case, if they are anti-parallel, then they overlap
                                    if (dot < 0 && Math.Abs(Math.Abs(dot) - 1) < coordtol)
                                    {
                                        kp.Value.relatedEdges.Add(edge);
                                        //the edge is still unique but accumulates a neighbor
                                        edge.relatedEdges.Add(kp.Value);
                                        edgeloopcounter++;
                                        continue;
                                    }
                                    else
                                    {
                                        //they don't intersect
                                        edgeloopcounter++;
                                        continue;
                                    }
                                }

                            }
                        }
                        else if (Math.Abs(diffx3) < coordtol && Math.Abs(diffy3) < coordtol && Math.Abs(diffz3) < coordtol)
                        {
                            if (Math.Abs(diffx2) <= coordtol && Math.Abs(diffy2) <= coordtol && Math.Abs(diffz2) <= coordtol)
                            {
                                //both match, means the match is perfect, so add it to the related surfaces list
                                kp.Value.relatedEdges.Add(edge);
                                perfectMatch = true;
                                break;

                            }
                            else
                            {
                                //the edge may be unique, and we can prove it still have neighboring relationships
                                //same as diff2
                                double Ax = gueststartpt.X - edge.startendpt[1].X;
                                double Ay = gueststartpt.Y - edge.startendpt[1].Y;
                                double Az = gueststartpt.Z - edge.startendpt[1].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                //same as diff4
                                double Bx = guestendpt.X - edge.startendpt[1].X;
                                double By = guestendpt.Y - edge.startendpt[1].Y;
                                double Bz = guestendpt.Z - edge.startendpt[1].Z;
                                Vector.MemorySafe_CartVect B = new MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);

                                //this is the easiest intersection because A and B prove to be inside the guest edge
                                //these vectors are antiparallel in this case
                                double dotAB = DotProductMag(A, B);
                                if (dotAB < 0 && Math.Abs(Math.Abs(dotAB) - 1) < dotprodtol)
                                {
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //in this case the point is outside of the guest edge, because both vectors prove to be parallel
                                else if (Math.Abs(dotAB-1) < dotprodtol)
                                {
                                    //the question is whether it completely overlaps the guest edge, or does not intersect it
                                    double ex = edge.startendpt[1].X - edge.startendpt[0].X;
                                    double ey = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                    double ez = edge.startendpt[1].Z - edge.startendpt[0].Z;

                                    double evX = guestendpt.X - gueststartpt.X;
                                    double evY = guestendpt.Y - gueststartpt.Y;
                                    double evZ = guestendpt.Z - gueststartpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    double dot = Vector.DotProductMag(A, ev);
                                    //in this case, if they are parallel, then they overlap
                                    if (Math.Abs(dot - 1) < dotprodtol)
                                    {
                                        kp.Value.relatedEdges.Add(edge);
                                        //the edge is still unique but accumulates a neighbor
                                        edge.relatedEdges.Add(kp.Value);
                                        edgeloopcounter++;
                                        continue;
                                    }
                                    else
                                    {
                                        //they don't intersect
                                        edgeloopcounter++;
                                        continue;
                                    }
                                }

                            }
                        }
                        else if (Math.Abs(diffx4) < coordtol && Math.Abs(diffy4) < coordtol && Math.Abs(diffz4) < coordtol)
                        {
                            if (Math.Abs(diffx1) < coordtol && Math.Abs(diffy1) < coordtol && Math.Abs(diffz1) < coordtol)
                            {
                                //both match, means the match is perfect, so the unique edge has found it complement related edge.  Great!
                                kp.Value.relatedEdges.Add(edge);
                                //I am done searching this test edge, and I can start over again with the next edge in question
                                perfectMatch = true;
                                break;

                            }

                            else
                            {
                                //so far, I have found only one thing in common, sharing of one point, even though second point did not match, the edges could still align
                                //draw vector A
                                double Ax = guestendpt.X - edge.startendpt[0].X;
                                double Ay = guestendpt.Y - edge.startendpt[0].Y;
                                double Az = guestendpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                double Bx = gueststartpt.X - edge.startendpt[0].X;
                                double By = gueststartpt.Y - edge.startendpt[0].Y;
                                double Bz = gueststartpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect B = new MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);

                                //this is the easiest intersection because A and B prove to be inside the guest edge
                                //these vectors are antiparallel in this case
                                double dotAB = DotProductMag(A, B);
                                if (dotAB < 0 && Math.Abs(Math.Abs(dotAB) - 1) < dotprodtol)
                                {
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //in this case the point is outside of the guest edge, because both vectors prove to be parallel
                                else if (Math.Abs(dotAB-1) < coordtol)
                                {
                                    //the question is whether it completely overlaps the guest edge, or does not intersect it
                                    double ex = edge.startendpt[1].X - edge.startendpt[0].X;
                                    double ey = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                    double ez = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                    //take cross product to see if they are parallel
                                    double evX = guestendpt.X - gueststartpt.X;
                                    double evY = guestendpt.Y - gueststartpt.Y;
                                    double evZ = guestendpt.Z - gueststartpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    double dot = Vector.DotProductMag(A, ev);
                                    if (Math.Abs(dot-1) < dotprodtol)
                                    {
                                        //parallel is an indication that they overlap
                                        kp.Value.relatedEdges.Add(edge);
                                        //the edge is still unique but accumulates a neighbor
                                        edge.relatedEdges.Add(kp.Value);
                                        edgeloopcounter++;
                                        continue;
                                    }
                                    else
                                    {
                                        //they don't intersect
                                        edgeloopcounter++;
                                        continue;
                                    }

                                }

                            }
                        }
                        
                        //neither points perfectly coincide, so we do an exhaustive overlap check.
                        else
                        {
                            
                            //are the two vectors even parallel?  because if they are not, no need to get more complex
                            double evX = guestendpt.X - gueststartpt.X;
                            double evY = guestendpt.Y - gueststartpt.Y;
                            double evZ = guestendpt.Z - gueststartpt.Z;
                            Vector.MemorySafe_CartVect guestvec = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                            

                            double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                            double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                            double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                            Vector.MemorySafe_CartVect testedgev = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                            //tolerance?
                            double dot = DotProductMag(guestvec, testedgev);
                            if (Math.Abs(Math.Abs(dot) - 1) > dotprodtol)
                            {
                                //they are not even parallel so move on
                                edgeloopcounter++;
                                continue;
                            }
                            //
                            //try to determine how the two edges overlap, or if they overlap at all.  We only know that they are parallel at this point
                            //test edge point 1

                            double Ax = gueststartpt.X - edge.startendpt[0].X;
                            double Ay = gueststartpt.Y - edge.startendpt[0].Y;
                            double Az = gueststartpt.Z - edge.startendpt[0].Z;
                            Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                            double Amag = Vector.VectorMagnitude(A);


                            double Bx = gueststartpt.X - edge.startendpt[1].X;
                            double By = gueststartpt.Y - edge.startendpt[1].Y;
                            double Bz = gueststartpt.Z - edge.startendpt[1].Z;
                            Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                            double Bmag = Vector.VectorMagnitude(B);
                            double crossB = Vector.VectorMagnitude(Vector.CrossProduct(B, guestvec));

                            double Cx = guestendpt.X - edge.startendpt[0].X;
                            double Cy = guestendpt.Y - edge.startendpt[0].Y;
                            double Cz = guestendpt.Z - edge.startendpt[0].Z;
                            Vector.MemorySafe_CartVect C = new MemorySafe_CartVect(Cx, Cy, Cz);

                            double Dx = guestendpt.X - edge.startendpt[1].X;
                            double Dy = guestendpt.Y - edge.startendpt[1].Y;
                            double Dz = guestendpt.Z - edge.startendpt[1].Z;
                            Vector.MemorySafe_CartVect D = new MemorySafe_CartVect(Dx, Dy, Dz);

                            double dotCD = Vector.DotProductMag(C, D);

                            //take the dot product to see if the point falls inside of the test edge
                            //this is actually the harder case first......
                            double dotAB = Vector.DotProductMag(A, B);
                            //if dotAB comes back as 1, this means the two vectors are parallel and thus falls outside of the guest edge
                            if (Math.Abs(dotAB - 1) < dotprodtol)
                            {
                                //does the other point fall outside of the test edge?  
                                
                                if (Math.Abs(dotCD - 1) < dotprodtol)
                                {
                                    //this means that this point also falls outside.  Now the question is, does it envelop
                                    //or does it completely miss the mark?  Dot A or B with C or D
                                    double dotAC = DotProductMag(A, C);
                                    //if this comes back as negative 1, the guest edge envelops because A and C are antiparallel
                                    if (dotAC < 0 && Math.Abs(Math.Abs(dotAC) - 1) < dotprodtol)
                                    {
                                        //parallel is an indication that they overlap
                                        kp.Value.relatedEdges.Add(edge);
                                        //the edge is still unique but accumulates a neighbor
                                        edge.relatedEdges.Add(kp.Value);
                                        edgeloopcounter++;
                                        continue;
                                    }
                                    else
                                    {
                                        //falls totally outside
                                        edgeloopcounter++;
                                        continue;
                                    }

                                }
                                //alternatively, this means the point falls inside while the other falls outside
                                else if (dotCD < 0 && Math.Abs(Math.Abs(dotCD) - 1) < dotprodtol)
                                {
                                    //we already know that they are parallel from the ccross product test first conducted
                                    //so we can conclude that regardless of what happens from here on out, at least one
                                    //point is overlapping, so we call it quits here
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                else
                                {
                                    edgeloopcounter++;
                                    continue;
                                }

                            }
                            //easier case...anti-parallel this means the point falls inside of the test edge
                            else if (Math.Abs(Math.Abs(dotAB) - 1) < dotprodtol)
                            {
                                //we already know that they are parallel from the ccross product test first conducted
                                //so we can conclude that regardless of what happens from here on out, at least one
                                //point is overlapping, so we call it quits here
                                kp.Value.relatedEdges.Add(edge);
                                //the edge is still unique but accumulates a neighbor
                                edge.relatedEdges.Add(kp.Value);
                                edgeloopcounter++;
                                continue;

                            }
                            //easier case...this is the only other condition where an overlap could exists
                            else if (dotCD < 0 && Math.Abs(Math.Abs(dotCD) - 1) < dotprodtol)
                            {
                                //we already know that they are parallel from the ccross product test first conducted
                                //so we can conclude that regardless of what happens from here on out, at least one
                                //point is overlapping, so we call it quits here
                                kp.Value.relatedEdges.Add(edge);
                                //the edge is still unique but accumulates a neighbor
                                edge.relatedEdges.Add(kp.Value);
                                edgeloopcounter++;
                                continue;
                            }
                            else
                            {
                                edgeloopcounter++;
                                continue;
                            }
                            
                        }


                    } //for each edge in unique edges
                    #endregion
                    //if I made it here, it means I did not find a perfect match for the test edge, so I add it to unique edges, with the relationships it has thus
                    //far accumulated
                    if (!perfectMatch)
                    {
                        uniqueedges[uniqueedgect] = edge;
                        uniqueedgect++;
                    }
                    perfectMatch = false;
                } //for loop of coordinates
            }
            catch (Exception e)
            {
                log.Error(e.ToString());
                log.Error("Error found when creating edge families.");

            }
            return uniqueedges;

        }

        public double getPlanarSA(List<CartVect> polygonVect)
        {
            List<CartVect> normalizedPlane = new List<CartVect>();
            //the new plane's first coordinate is arbitrarily set to zero
            normalizedPlane[0].X = 0;
            normalizedPlane[0].Y = 0;
            normalizedPlane[0].Z = 0;
            double diffX = 0;
            double diffY = 0;
            double diffZ = 0;

            double surfaceArea = -1;
            int numPoints = polygonVect.Count;
            for(int i=0; i<numPoints; i++)
            {
                
                if (i > 0)
                {

                }
            }
            return surfaceArea;
        }

    }
}

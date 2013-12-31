using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using log4net;

namespace VectorMath
{

    public class Vector
    {
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
                    Vector.MemorySafe_CartCoord P1 = coords[j];

                    //log.Debug("Testing adjacent vector " + k.ToString());
                    Vector.MemorySafe_CartVect L2 = pointvec[coords[j + k]];
                    Vector.MemorySafe_CartCoord P2 = coords[j + k];

                    Vector.MemorySafe_CartVect tempvec = CreateMemorySafe_Vector(P1, P2);
                    Vector.MemorySafe_CartVect numcross = CrossProduct(tempvec, L2);
                    double numcrmag = VectorMagnitude(numcross);
                    Vector.MemorySafe_CartVect dencross = CrossProduct(L1, L2);
                    double dencrmag = VectorMagnitude(dencross);
                    //basic checks suggested:  are numcross and dencross parallel?
                    //is numcross not the zero vector?
                    if (numcrmag == 0 && dencrmag == 0)
                    {
                        //the two lines are clearly parallel
                        //not sure if this check is needed, but we ensure that these two parallel lines are unique and do not overlap
                        //it gets a bit more complicated from here

                    }
                    double a = numcrmag / dencrmag;

                    if (k == 1)
                    {
                        log.Info("Vector "+(j+k).ToString()+" is a nearest neighbor.");
                        //we expect a to be one and P2 to be equal to the resulting coord
                        //for the first nearest neighbor
                        if (a == 1.0)
                        {
                            log.Info("PASS:  Nearest neighbor "+(j+k).ToString()+" is intersected at its startpoint by vector "+j.ToString()+".");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "PASS");
                            
                            Vector.MemorySafe_CartCoord rescoord = SumPointAndLine(P1, L1);
                            if (rescoord.X == P2.X && rescoord.Y == P2.Y && rescoord.Z == P2.Z)
                            {
                                //and b == 0.  This truth statement is the equivalent
                                IntTruth.Add(true);
                                continue;
                            }
                        }
                        else
                        {
                            log.Info("FAIL:  Nearest neighbor "+(j+k).ToString()+" is not intersected at its endpoint by vector "+j.ToString()+".  This indicates an invalid enclosure.");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "FAIL");
                            return false;
                        }
                    }
                    else
                    {
                        //the last vector in the set of vector 0
                        if (j == 0 && a == 0 && k == coordCount - 1) { log.Info("Vector " + k.ToString() + " is a nearest neighbor."); }
                        else { log.Info("Vector " + (j+k).ToString() + " is not a nearest neighbor."); }
                        //parallel lines
                        if (a == double.PositiveInfinity || a == double.NegativeInfinity)
                        {
                            log.Info("PASS:  Vectors "+j.ToString()+ " and " + (j+k).ToString()+" are parallel vectors.");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "PASS");
                            IntTruth.Add(true);
                        }
                        //ok when j==0, but not otherwise
                        else if (j==0 && a == 0 && k == coordCount - 1)
                        {
                            //This may take the place of the test at the last vector in the list
                            //I believe this only happens when the line intersects at the starting point
                            //we may not need 
                            log.Info("The vector " + (j+k).ToString() + " intersects vector " + (j).ToString() + " at its starting point.");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNY" + ";" + "PASS");
                            return true;
                        }
                        else if (j != 0 && a == 0)
                        {
                            log.Error("The vector " + (j + k).ToString() + " intersects vector " + j.ToString() + " at its starting point.");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "FAIL");
                            return false;
                        }
                        //ok
                        else if (a > 1.0)
                        {
                            //this means the points intersect, but some point beyond the length of the vector
                            log.Info("Pass.  Vectors would intersect but beyond the bounds of the vectors.");
                            log.Debug(j.ToString() + ":" + (j+k).ToString() + ";" + "NNN" + ";" + "PASS");
                            IntTruth.Add(true);
                        }
                        //bad
                        else if (a < 1.0)
                        {
                            //this means they may intersect, I need to check the magnitude of the proposed intersection
                            //maybe, solve for a and b and report
                            log.Debug("a is less than 1.");

                            
                            Vector.MemorySafe_CartVect F1 = VectorTimesScalar(L1, a);
                            Vector.MemorySafe_CartCoord rescoord = SumPointAndLine(P1, F1);
                            Vector.MemorySafe_CartVect newvect = CreateMemorySafe_Vector(P2, rescoord);
                            //Dec 28, 2013
                            //this new vector may not face in the direction of L2.  Since we just draw a line from P2 to the new resulting coordinate.
                            //if the new vector is not parallel to L2, this is another indication, similar to case 2.
                            Vector.MemorySafe_CartVect L2unitvec = Vector.UnitVector(L2);
                            Vector.MemorySafe_CartVect newunitvec = Vector.UnitVector(newvect);
                            double diffx = Math.Abs(L2unitvec.X - newunitvec.X);
                            double diffy = Math.Abs(L2unitvec.Y - newunitvec.Y);
                            double diffz = Math.Abs(L2unitvec.Z - newunitvec.Z);

                            if (diffx <= 0.01 || diffy <= 0.01 && diffz <= 0.01)
                            {
                                double newvectmag = Vector.VectorMagnitude(newvect);
                                double L2mag = Vector.VectorMagnitude(L2);
                                if (newvectmag <= L2mag)
                                {
                                    //this clearly means that the vectors intersect
                                    log.Info("An illegal intersection in the polygon has been detected.");
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
                                    return false;
                                }
                                else
                                {
                                    log.Info("PASS:  Vectors could intersect but beyond the bounds of the vectors.");
                                    log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "PASS");
                                    IntTruth.Add(true);
                                }
                            }
                            //Dec 28 2013
                            //If we arrive here, this essentially means that the vectors don't cross
                            else
                            {
                                log.Info("PASS:  Vectors could intersect but beyond the bounds of the vectors.");
                                log.Debug(j.ToString() + ":" + (j + k).ToString() + ";" + "NNN" + ";" + "FAIL");
                                IntTruth.Add(true);
                            }
                            
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

        public static CartCoord GetCentroid(List<CartCoord> coordinates)
        {
            //based on a formula for a centroid of a planar surface defined with three coordinates
            //Cx = (1/6A)*sig(i=0->N-1) (x(i)+x(i+1))(x(i)y(i+1)-x(i+1)y(i))
            //Cy = (1/6A)*sig(i=0->N-1) (y(i)+y(i+1))(x(i)y(i+1)-x(i+1)y(i))
            //Cz = (1/6A)*sig(i=0->N-1) (z(i)+z(i+1))(x(i)z(i+1)-x(i+1)z(i))
            //we recognize this may not be entirely accurate but should be for the types of surfaces we are going to deal with we assume the centroid will fall somewhere inside the coordinates)

            CartCoord centroid = new CartCoord();
            //calculate the area of the polyLoop

            return centroid;
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

        public static double GetAreaofSurface(ModelingUtilities.BuildingObjects.MemorySafe_Surface surface)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            MemorySafe_CartVect RHRVector = GetRHR(surface.SurfaceCoords);

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

        public static double GetAreaofWindow(ModelingUtilities.BuildingObjects.MemorySafe_OpeningDefinitions opening)
        {
            //Used to figure out how best to calculate the area from a given surfacce. 
            //Get the coordinates that define the surface
            //get the area based on the coordinates

            //Get the RHRVector (the actual direction is not important
            MemorySafe_CartVect RHRVector = GetRHR(opening.coordinateList);

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

        public static Vector.MemorySafe_CartVect GetRHR(List<Vector.MemorySafe_CartCoord> plCoords)
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

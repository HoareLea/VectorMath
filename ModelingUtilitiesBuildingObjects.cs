using System;
using VectorMath.Vector;


namespace ModelingUtilities
{
    public class BuildingObjects
    {
        public class Surface
        {
            
            public string name;
            public SurfaceTypes surfaceType;
            public string constructionName;
            public OutsideBoundary outsideBoundary;
            public enum SurfaceTypes
            {
                Wall,
                Floor,
                Ceiling,
                Roof,
                Blank
            }
            public enum OutsideBoundary
            {
                Surface,
                Ground,
                Outdoors,
                Zone,
                OtherSideCoefficients,
                OtherSideConditionsModel,
                Blank

            }
            public string zoneName;
            public string outsideBoundaryCondition;
            public string sunExposureVar;
            public string windExposureVar;
            public double viewFactor;
            public int numVertices;
            public List<Vector.CartCoord> SurfaceCoords;
            public double tilt;
            public double azimuth;

            public void Clear()
            {
                name = "";
                surfaceType = SurfaceTypes.Blank;
                constructionName = "";
                outsideBoundary = OutsideBoundary.Blank;
                zoneName = "";
                sunExposureVar = "";
                windExposureVar = "";
                numVertices = 0;
                SurfaceCoords.Clear();

            }

        }
    }

}
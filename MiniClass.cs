using System;

public class Building
{
	public Building()
	{
	}

    public class MemorySafe_Spaces
    {
        private readonly string _name;
        private readonly int _multiplier;
        private readonly List<MemorySafe_Surface> _spaceSurfaces;

        public MemorySafe_Spaces(string name, int multiplier, List<MemorySafe_Surface> spaceSurfaces)
        {
            _name = name;
            _multiplier = multiplier;
            _spaceSurfaces = spaceSurfaces;
        }

        public string name { get { return _name; } }
        public int multiplier { get { return _multiplier; } }
        public List<MemorySafe_Surface> spaceSurfaces { get { return _spaceSurfaces; } }

    }

    public class Spaces
    {
        public string name;
        public int multiplier;
        public List<Surface> spaceSurfaces;

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
    public enum SurfaceTypes
    {
        Wall,
        Floor,
        Ceiling,
        Roof,
        Blank
    }

    public class MemorySafe_Surface
    {

        public string _name;
        public int _multiplier;
        public SurfaceTypes _surfaceType;
        public string _constructioName;
        public OutsideBoundary _outsideBoundary;
        public string _zoneName;
        public string _outsideBoundaryCondition;
        public string _sunExposureVar;
        public string _windExposureVar;
        public double _viewFactor;
        public int _numVertices;
        public List<Vector.MemorySafe_CartCoord> _SurfaceCoords;
        public double _tilt;
        public double _azimuth;

        public MemorySafe_Surface(string name, int multiplier, SurfaceTypes surfaceType, string constName, OutsideBoundary ob,
            string zoneName, string outsideBC, string sunExp, string windExp, double vF, int numVert, List<Vector.MemorySafe_CartCoord> sC,
            double tilt, double az)
        {
            _name = name;
            _multiplier = multiplier;
            _surfaceType = surfaceType;
            _constructioName = constName;
            _outsideBoundary = ob;
            _zoneName = zoneName;
            _outsideBoundaryCondition = outsideBC;
            _sunExposureVar = sunExp;
            _windExposureVar = windExp;
            _viewFactor = vF;
            _numVertices = numVert;
            _SurfaceCoords = sC;
            _tilt = tilt;
            _azimuth = az;

        }

        public string name { get { return _name; } }
        public int multiplier { get { return _multiplier; } }
        public SurfaceTypes surfaceType { get { return _surfaceType; } }
        public string constructionName { get { return _constructioName; } }
        public OutsideBoundary outsideBoundary { get { return _outsideBoundary; } }
        public string zoneName { get { return _zoneName; } }
        public string outsideBoundaryCondition { get { return _outsideBoundaryCondition; } }
        public string sunExposureVar { get { return _sunExposureVar; } }
        public string windExposureVar { get { return _windExposureVar; } }
        public double viewFactor { get { return _viewFactor; } }
        public int numVertices { get { return _numVertices; } }
        public List<Vector.MemorySafe_CartCoord> SurfaceCoords { get { return _SurfaceCoords; } }
        public double tilt { get { return _tilt; } }
        public double azimuth { get { return _azimuth; } }

    }
}

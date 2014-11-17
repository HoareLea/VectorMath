using System;
using System.Collections.Generic;
using VectorMath;



namespace ModelingUtilities
{
    public class BuildingObjects
    {
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
            public string name { get; set; }
            public double area { get; set; }
            public int multiplier { get; set; }
            public List<Surface> spaceSurfaces { get; set; }


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

            public MemorySafe_Surface(string name,int multiplier,SurfaceTypes surfaceType,string constName,OutsideBoundary ob,
                string zoneName,string outsideBC,string sunExp,string windExp,double vF, int numVert, List<Vector.MemorySafe_CartCoord> sC,
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
            public int multiplier { get {return _multiplier;}}
            public SurfaceTypes surfaceType { get { return _surfaceType; }}
            public string constructionName { get { return _constructioName; }}
            public OutsideBoundary outsideBoundary { get { return _outsideBoundary; }}
            public string zoneName { get { return _zoneName; }}
            public string outsideBoundaryCondition { get { return _outsideBoundaryCondition; }}
            public string sunExposureVar { get {return _sunExposureVar; }}
            public string windExposureVar { get {return _windExposureVar; }}
            public double viewFactor { get { return _viewFactor; }}
            public int numVertices { get { return _numVertices; }}
            public List<Vector.MemorySafe_CartCoord> SurfaceCoords { get { return _SurfaceCoords; }}
            public double tilt { get { return _tilt; }}
            public double azimuth { get { return _azimuth; }}

        }
        public class Surface
        {
            public string name;
            public int multiplier;
            public SurfaceTypes surfaceType;
            public string constructionName;
            public OutsideBoundary outsideBoundary;
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

        public class MemorySafe_OpeningDefinitions
        {
            private readonly string _nameId;
            private readonly string _openingType;
            private readonly string _parentSurfaceNameId;
            private readonly double _parentAzimuth;
            private readonly double _parentTilt;
            private readonly string _outsideBoundaryConditionObj;
            private readonly double _viewFactortoGround;
            private readonly string _shadeControlSch;
            private readonly List<Vector.MemorySafe_CartCoord> _coordinateList;
            private readonly double _Azimuth;
            private readonly double _Tilt;
            private readonly Vector.MemorySafe_CartVect _rHRVector;
            private readonly string _constructionName;

            private readonly string _frameAndDividerName;
            private readonly int _multiplier;
            private readonly int _numVertices;
            private readonly double _area;

            public MemorySafe_OpeningDefinitions(string nameId, string openType,string parentSurf,double parentAz,double parentTilt,
                string oBCond,double vF, string shadeCntrlSch, List<Vector.MemorySafe_CartCoord> coordList, double az,
                double tilt, Vector.MemorySafe_CartVect RHR, string constName, string frameandDivider, int mult, int numVert, double area)
            {
                _nameId = nameId;
                _openingType = openType;
                _parentSurfaceNameId = parentSurf;
                _parentAzimuth=parentAz;
                _parentTilt=parentTilt;
            
                 _outsideBoundaryConditionObj = oBCond;
                 _viewFactortoGround = vF;
                 _shadeControlSch = shadeCntrlSch;
                 _coordinateList = coordList;
                 _Azimuth = az;
                 _Tilt = tilt;
                  Vector.MemorySafe_CartVect _rHRVector = RHR;
                 _constructionName = constName;

                 _frameAndDividerName = frameandDivider;
                 _multiplier = mult;
                 _numVertices = numVert;
                 _area = area;
            }

            public string nameId { get { return _nameId; }}
            public string openingType { get { return _openingType; }}
            public string parentSurfaceNameId { get { return _parentSurfaceNameId; }}
            public double parentAzimuth { get { return _parentAzimuth; }}
            public double parentTilt { get { return _parentTilt; }}
            public string outsideBoundaryConditionObj { get { return _outsideBoundaryConditionObj; }}
            public double viewFactortoGround { get { return _viewFactortoGround; }}
            public string shadeControlSch { get { return _shadeControlSch; }}
            public List<Vector.MemorySafe_CartCoord> coordinateList { get {return _coordinateList; }}
            public double Azimuth { get {return _Azimuth; }}
            public double Tilt { get { return _Tilt; }}
            public Vector.MemorySafe_CartVect rHRVector { get { return _rHRVector; }}
            public string constructionName { get {return _constructionName; }}

            public string frameAndDividerName { get { return _frameAndDividerName; }}
            public int multiplier { get { return _multiplier; }}
            public int numVertices { get { return _numVertices; }}
            public double area { get { return _area; }}

            public MemorySafe_OpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening)
            {
                _nameId = previousOpening.nameId;
                _openingType = previousOpening.openingType;
                _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                _parentAzimuth = previousOpening.parentAzimuth;;
                _parentTilt = previousOpening.parentTilt;
                _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                _viewFactortoGround = previousOpening.viewFactortoGround;
                _shadeControlSch = previousOpening.shadeControlSch;
                _coordinateList = previousOpening.coordinateList;
                _Azimuth = previousOpening.Azimuth;
                _Tilt = previousOpening.Tilt;
                _rHRVector = previousOpening.rHRVector;
                _constructionName = previousOpening.constructionName;

                _frameAndDividerName = previousOpening.frameAndDividerName;
                _multiplier = previousOpening.multiplier;
                _numVertices = previousOpening.numVertices;
                _area = previousOpening.area;
            }

            public MemorySafe_OpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, int multiplier)
            {
                _nameId = previousOpening.nameId;
                _openingType = previousOpening.openingType;
                _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                _parentAzimuth = previousOpening.parentAzimuth; ;
                _parentTilt = previousOpening.parentTilt;
                _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                _viewFactortoGround = previousOpening.viewFactortoGround;
                _shadeControlSch = previousOpening.shadeControlSch;
                _coordinateList = previousOpening.coordinateList;
                _Azimuth = previousOpening.Azimuth;
                _Tilt = previousOpening.Tilt;
                _rHRVector = previousOpening.rHRVector;
                _constructionName = previousOpening.constructionName;

                _frameAndDividerName = previousOpening.frameAndDividerName;
                _multiplier = multiplier;
                _numVertices = previousOpening.numVertices;
                _area = previousOpening.area;
            }

            public MemorySafe_OpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, List<Vector.MemorySafe_CartCoord> desiredCoords)
            {
                _nameId = previousOpening.nameId;
                _openingType = previousOpening.openingType;
                _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                _parentAzimuth = previousOpening.parentAzimuth; ;
                _parentTilt = previousOpening.parentTilt;
                _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                _viewFactortoGround = previousOpening.viewFactortoGround;
                _shadeControlSch = previousOpening.shadeControlSch;
                //gives a new set of coordinates as provided
                _coordinateList = desiredCoords;
                _Azimuth = previousOpening.Azimuth;
                _Tilt = previousOpening.Tilt;
                _rHRVector = previousOpening.rHRVector;
                _constructionName = previousOpening.constructionName;

                _frameAndDividerName = previousOpening.frameAndDividerName;
                _multiplier = previousOpening.multiplier;
                _numVertices = previousOpening.numVertices;
                _area = previousOpening.area;
            }

            public MemorySafe_OpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, List<Vector.MemorySafe_CartCoord> desiredCoords, int newMultiplier)
            {
                _nameId = previousOpening.nameId;
                _openingType = previousOpening.openingType;
                _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                _parentAzimuth = previousOpening.parentAzimuth; ;
                _parentTilt = previousOpening.parentTilt;
                _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                _viewFactortoGround = previousOpening.viewFactortoGround;
                _shadeControlSch = previousOpening.shadeControlSch;
                //gives a new set of coordinates as provided
                _coordinateList = desiredCoords;
                _Azimuth = previousOpening.Azimuth;
                _Tilt = previousOpening.Tilt;
                _rHRVector = previousOpening.rHRVector;
                _constructionName = previousOpening.constructionName;

                _frameAndDividerName = previousOpening.frameAndDividerName;
                _multiplier = newMultiplier;
                _numVertices = previousOpening.numVertices;
                _area = previousOpening.area;
            }
        }

        public class MemorySafe_ADOpeningDefinitions
        {
            private readonly string _nameId;
            private readonly string _openingType;
            private readonly string _parentSurfaceNameId;
            private readonly double _parentAzimuth;
            private readonly double _parentTilt;
            private readonly string _outsideBoundaryConditionObj;
            private readonly double _viewFactortoGround;
            private readonly string _shadeControlSch;
            private readonly List<Vector.MemorySafe_CartCoord> _coordinateList;
            private readonly double _Azimuth;
            private readonly double _Tilt;
            private readonly Vector.MemorySafe_CartVect _rHRVector;
            private readonly string _constructionName;

            private readonly string _frameAndDividerName;
            private readonly int _multiplier;
            private readonly int _numVertices;
            private readonly double _area;
            private readonly double _x;
            private readonly double _z;
            private readonly double _height;
            private readonly double _length;
            

            public MemorySafe_ADOpeningDefinitions(string nameId, string openType, string parentSurf, double parentAz, double parentTilt,
                string oBCond, double vF, string shadeCntrlSch, List<Vector.MemorySafe_CartCoord> coordList, double az,
                double tilt, Vector.MemorySafe_CartVect RHR, string constName, string frameandDivider, int mult, int numVert, double area, 
                double x, double z, double height, double length)
            {
                _nameId = nameId;
                _openingType = openType;
                _parentSurfaceNameId = parentSurf;
                _parentAzimuth = parentAz;
                _parentTilt = parentTilt;

                _outsideBoundaryConditionObj = oBCond;
                _viewFactortoGround = vF;
                _shadeControlSch = shadeCntrlSch;
                _coordinateList = coordList;
                _Azimuth = az;
                _Tilt = tilt;
                Vector.MemorySafe_CartVect _rHRVector = RHR;
                _constructionName = constName;

                _frameAndDividerName = frameandDivider;
                _multiplier = mult;
                _numVertices = numVert;
                _area = area;
                _x = x;
                _z = z;
                _height = height;
                _length = length;
            }

            public string nameId { get { return _nameId; } }
            public string openingType { get { return _openingType; } }
            public string parentSurfaceNameId { get { return _parentSurfaceNameId; } }
            public double parentAzimuth { get { return _parentAzimuth; } }
            public double parentTilt { get { return _parentTilt; } }
            public string outsideBoundaryConditionObj { get { return _outsideBoundaryConditionObj; } }
            public double viewFactortoGround { get { return _viewFactortoGround; } }
            public string shadeControlSch { get { return _shadeControlSch; } }
            public List<Vector.MemorySafe_CartCoord> coordinateList { get { return _coordinateList; } }
            public double Azimuth { get { return _Azimuth; } }
            public double Tilt { get { return _Tilt; } }
            public Vector.MemorySafe_CartVect rHRVector { get { return _rHRVector; } }
            public string constructionName { get { return _constructionName; } }

            public string frameAndDividerName { get { return _frameAndDividerName; } }
            public int multiplier { get { return _multiplier; } }
            public int numVertices { get { return _numVertices; } }
            public double area { get { return _area; } }
            public double x { get { return _x; } }
            public double z  { get { return _z;} }
            public double height { get { return _height; } }
            public double length { get { return _length; } }

            //public MemorySafe_ADOpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening)
            //{
            //    _nameId = previousOpening.nameId;
            //    _openingType = previousOpening.openingType;
            //    _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
            //    _parentAzimuth = previousOpening.parentAzimuth; ;
            //    _parentTilt = previousOpening.parentTilt;
            //    _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
            //    _viewFactortoGround = previousOpening.viewFactortoGround;
            //    _shadeControlSch = previousOpening.shadeControlSch;
            //    _coordinateList = previousOpening.coordinateList;
            //    _Azimuth = previousOpening.Azimuth;
            //    _Tilt = previousOpening.Tilt;
            //    _rHRVector = previousOpening.rHRVector;
            //    _constructionName = previousOpening.constructionName;

            //    _frameAndDividerName = previousOpening.frameAndDividerName;
            //    _multiplier = previousOpening.multiplier;
            //    _numVertices = previousOpening.numVertices;
            //    _area = previousOpening.area;
            //}

            //public MemorySafe_ADOpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, int multiplier)
            //{
            //    _nameId = previousOpening.nameId;
            //    _openingType = previousOpening.openingType;
            //    _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
            //    _parentAzimuth = previousOpening.parentAzimuth; ;
            //    _parentTilt = previousOpening.parentTilt;
            //    _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
            //    _viewFactortoGround = previousOpening.viewFactortoGround;
            //    _shadeControlSch = previousOpening.shadeControlSch;
            //    _coordinateList = previousOpening.coordinateList;
            //    _Azimuth = previousOpening.Azimuth;
            //    _Tilt = previousOpening.Tilt;
            //    _rHRVector = previousOpening.rHRVector;
            //    _constructionName = previousOpening.constructionName;

            //    _frameAndDividerName = previousOpening.frameAndDividerName;
            //    _multiplier = multiplier;
            //    _numVertices = previousOpening.numVertices;
            //    _area = previousOpening.area;
            //}

            //public MemorySafe_ADOpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, List<Vector.MemorySafe_CartCoord> desiredCoords)
            //{
            //    _nameId = previousOpening.nameId;
            //    _openingType = previousOpening.openingType;
            //    _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
            //    _parentAzimuth = previousOpening.parentAzimuth; ;
            //    _parentTilt = previousOpening.parentTilt;
            //    _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
            //    _viewFactortoGround = previousOpening.viewFactortoGround;
            //    _shadeControlSch = previousOpening.shadeControlSch;
            //    //gives a new set of coordinates as provided
            //    _coordinateList = desiredCoords;
            //    _Azimuth = previousOpening.Azimuth;
            //    _Tilt = previousOpening.Tilt;
            //    _rHRVector = previousOpening.rHRVector;
            //    _constructionName = previousOpening.constructionName;

            //    _frameAndDividerName = previousOpening.frameAndDividerName;
            //    _multiplier = previousOpening.multiplier;
            //    _numVertices = previousOpening.numVertices;
            //    _area = previousOpening.area;
            //}

            //public MemorySafe_ADOpeningDefinitions(MemorySafe_OpeningDefinitions previousOpening, List<Vector.MemorySafe_CartCoord> desiredCoords, int newMultiplier)
            //{
            //    _nameId = previousOpening.nameId;
            //    _openingType = previousOpening.openingType;
            //    _parentSurfaceNameId = previousOpening.parentSurfaceNameId;
            //    _parentAzimuth = previousOpening.parentAzimuth; ;
            //    _parentTilt = previousOpening.parentTilt;
            //    _outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
            //    _viewFactortoGround = previousOpening.viewFactortoGround;
            //    _shadeControlSch = previousOpening.shadeControlSch;
            //    //gives a new set of coordinates as provided
            //    _coordinateList = desiredCoords;
            //    _Azimuth = previousOpening.Azimuth;
            //    _Tilt = previousOpening.Tilt;
            //    _rHRVector = previousOpening.rHRVector;
            //    _constructionName = previousOpening.constructionName;

            //    _frameAndDividerName = previousOpening.frameAndDividerName;
            //    _multiplier = newMultiplier;
            //    _numVertices = previousOpening.numVertices;
            //    _area = previousOpening.area;
            //}
        }

        public static MemorySafe_Surface convert2MemorySafeSurface(Surface surface)
        {
            List<Vector.MemorySafe_CartCoord> surfaceCoords = new List<Vector.MemorySafe_CartCoord>();
            foreach(Vector.CartCoord coord in surface.SurfaceCoords)
            {
                Vector.MemorySafe_CartCoord surfaceCoord = new Vector.MemorySafe_CartCoord(coord.X,coord.Y,coord.Z);
                surfaceCoords.Add(surfaceCoord);
            }

            MemorySafe_Surface memSurface = new MemorySafe_Surface(surface.name, surface.multiplier, surface.surfaceType, surface.constructionName,
                surface.outsideBoundary, surface.zoneName, surface.outsideBoundaryCondition, surface.sunExposureVar, surface.windExposureVar,
                surface.viewFactor, surface.numVertices, surfaceCoords, surface.tilt, surface.azimuth);

            return memSurface;

        }

        public static MemorySafe_OpeningDefinitions convert2MemorySafeOpening(OpeningDefinitions opening)
        {
            List<Vector.MemorySafe_CartCoord> surfaceCoords = new List<Vector.MemorySafe_CartCoord>();
            foreach (Vector.CartCoord coord in opening.coordinateList)
            {
                Vector.MemorySafe_CartCoord surfaceCoord = new Vector.MemorySafe_CartCoord(coord.X, coord.Y, coord.Z);
                surfaceCoords.Add(surfaceCoord);
            }
            Vector.MemorySafe_CartVect memSafeRHR = Vector.convertToMemorySafeVector(opening.rHRVector);

            MemorySafe_OpeningDefinitions memOpening = new MemorySafe_OpeningDefinitions(opening.nameId,
                opening.openingType, opening.parentSurfaceNameId,
                opening.parentAzimuth, opening.parentTilt, opening.outsideBoundaryConditionObj, opening.viewFactortoGround,
                opening.shadeControlSch, surfaceCoords, opening.Azimuth, opening.Tilt, memSafeRHR,
                opening.constructionName, opening.frameAndDividerName, opening.multiplier, opening.numVertices,
                opening.area);

            return memOpening;

        }

        public static MemorySafe_ADOpeningDefinitions convert2ADMemorySafeOpening(OpeningDefinitions opening)
        {
            List<Vector.MemorySafe_CartCoord> surfaceCoords = new List<Vector.MemorySafe_CartCoord>();

            Vector.CartVect dummy = new Vector.CartVect();
            dummy.X = -999;
            dummy.Y = -999;
            dummy.Z = -999;
            Vector.MemorySafe_CartVect memSafeRHR = Vector.convertToMemorySafeVector(dummy);
            string openingType = "N/A";
            double parentAzimuth = 0;
            double parentTilt = 0;
            string outsideBoundary = "Outdoors";
            double viewFactor = 0;
            double az = -999;
            double tilt = -999;
            int numvertices = 4;
            double area = opening.height * opening.length;
            MemorySafe_ADOpeningDefinitions memOpening = new MemorySafe_ADOpeningDefinitions(opening.nameId,
                openingType, opening.parentSurfaceNameId,
                parentAzimuth, parentTilt, outsideBoundary, viewFactor,
                opening.shadeControlSch, surfaceCoords, az, tilt, memSafeRHR,
                opening.constructionName, opening.frameAndDividerName, opening.multiplier, numvertices,
                opening.area, opening.X, opening.Z, opening.height, opening.length);

            return memOpening;

        }

        public static Surface convert2TempSurface(MemorySafe_Surface surface)
        {
            Surface tempSurface = new Surface();
            
            foreach (Vector.MemorySafe_CartCoord coord in surface.SurfaceCoords)
            {
                Vector.CartCoord surfaceCoord = new Vector.CartCoord();
                surfaceCoord.X = coord.X;
                surfaceCoord.Y = coord.Y;
                surfaceCoord.Z = coord.Z;
                tempSurface.SurfaceCoords.Add(surfaceCoord);
            }

            tempSurface.name = surface.name;
            tempSurface.multiplier = surface.multiplier;
            tempSurface.surfaceType = surface.surfaceType;
            tempSurface.constructionName = surface.constructionName;
            tempSurface.outsideBoundary = surface.outsideBoundary;
            tempSurface.zoneName = surface.zoneName;
            tempSurface.outsideBoundaryCondition = surface.outsideBoundaryCondition;
            tempSurface.sunExposureVar = surface.sunExposureVar;
            tempSurface.windExposureVar = surface.windExposureVar;
            tempSurface.viewFactor = surface.viewFactor;
            tempSurface.numVertices = surface.numVertices;
            tempSurface.tilt =  surface.tilt;
            tempSurface.azimuth = surface.azimuth;

            return tempSurface;

        }
        public static OpeningDefinitions convert2TempOpening(MemorySafe_OpeningDefinitions opening)
        {
            OpeningDefinitions Opening = new OpeningDefinitions();
            
            foreach (Vector.MemorySafe_CartCoord coord in opening.coordinateList)
            {
                Vector.CartCoord opCoord = new Vector.CartCoord();
                opCoord.X = coord.X;
                opCoord.Y = coord.Y;
                opCoord.Z = coord.Z;
                Opening.coordinateList.Add(opCoord);
            }
            Vector.CartVect RHR = Vector.convertToTempVector(opening.rHRVector);

            
            Opening.nameId = opening.nameId;
            Opening.openingType = opening.openingType;
            Opening.parentSurfaceNameId = opening.parentSurfaceNameId;
            Opening.parentAzimuth = opening.parentAzimuth;
            Opening.parentTilt = opening.parentTilt;
            Opening.outsideBoundaryConditionObj = opening.outsideBoundaryConditionObj;
            Opening.viewFactortoGround = opening.viewFactortoGround;
            Opening.shadeControlSch = opening.shadeControlSch;
            
            Opening.Azimuth = opening.Azimuth;
            Opening.Tilt = opening.Tilt;
            Opening.rHRVector = RHR;
            Opening.constructionName =  opening.constructionName;
            Opening.frameAndDividerName = opening.frameAndDividerName;
            Opening.multiplier = opening.multiplier;
            Opening.numVertices = opening.numVertices;
            Opening.area = opening.area;

            return Opening;

        }
        public class OpeningDefinitions
        {
            //an empty constructor
            public OpeningDefinitions()
            {
                
            }
            public OpeningDefinitions(OpeningDefinitions previousOpening)
            {
                nameId = previousOpening.nameId;
                openingType = previousOpening.openingType;
                parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                parentAzimuth = previousOpening.parentAzimuth;;
                parentTilt = previousOpening.parentTilt;
                outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                viewFactortoGround = previousOpening.viewFactortoGround;
                shadeControlSch = previousOpening.shadeControlSch;
                coordinateList = previousOpening.coordinateList;
                Azimuth = previousOpening.Azimuth;
                Tilt = previousOpening.Tilt;
                rHRVector = previousOpening.rHRVector;
                constructionName = previousOpening.constructionName;

                frameAndDividerName = previousOpening.frameAndDividerName;
                multiplier = previousOpening.multiplier;
                numVertices = previousOpening.numVertices;
                area = previousOpening.area;
            }

            public OpeningDefinitions(OpeningDefinitions previousOpening, List<Vector.CartCoord> desiredCoords)
            {
                nameId = previousOpening.nameId;
                openingType = previousOpening.openingType;
                parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                parentAzimuth = previousOpening.parentAzimuth; ;
                parentTilt = previousOpening.parentTilt;
                outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                viewFactortoGround = previousOpening.viewFactortoGround;
                shadeControlSch = previousOpening.shadeControlSch;
                //gives a new set of coordinates as provided
                coordinateList = desiredCoords;
                Azimuth = previousOpening.Azimuth;
                Tilt = previousOpening.Tilt;
                rHRVector = previousOpening.rHRVector;
                constructionName = previousOpening.constructionName;

                frameAndDividerName = previousOpening.frameAndDividerName;
                multiplier = previousOpening.multiplier;
                numVertices = previousOpening.numVertices;
                area = previousOpening.area;
            }

            public OpeningDefinitions(OpeningDefinitions previousOpening, List<Vector.CartCoord> desiredCoords, int newMultiplier)
            {
                nameId = previousOpening.nameId;
                openingType = previousOpening.openingType;
                parentSurfaceNameId = previousOpening.parentSurfaceNameId;
                parentAzimuth = previousOpening.parentAzimuth; ;
                parentTilt = previousOpening.parentTilt;
                outsideBoundaryConditionObj = previousOpening.outsideBoundaryConditionObj;
                viewFactortoGround = previousOpening.viewFactortoGround;
                shadeControlSch = previousOpening.shadeControlSch;
                //gives a new set of coordinates as provided
                coordinateList = desiredCoords;
                Azimuth = previousOpening.Azimuth;
                Tilt = previousOpening.Tilt;
                rHRVector = previousOpening.rHRVector;
                constructionName = previousOpening.constructionName;

                frameAndDividerName = previousOpening.frameAndDividerName;
                multiplier = newMultiplier;
                numVertices = previousOpening.numVertices;
                area = previousOpening.area;
            }
            public string nameId;
            public string openingType;
            public string parentSurfaceNameId;
            public double parentAzimuth;
            public double parentTilt;
            public string outsideBoundaryConditionObj;
            public double viewFactortoGround;
            public string shadeControlSch;
            public List<Vector.CartCoord> coordinateList;
            public double Azimuth;
            public double Tilt;
            public Vector.CartVect rHRVector;
            public string constructionName;

            public string frameAndDividerName;
            public int multiplier;
            public int numVertices;
            public double area;
            public double X;
            public double Z;
            public double length;
            public double height;
        }
    }

    public class ProjectFeatures
    {
        //the constructor
        public ProjectFeatures(string location, int year)
        {
            currentProjectWWR = 0.0;
            currentNorthWWR = 0.0;
            currentEastWWR = 0;
            currentSouthWWR = 0;
            currentWestWWR = 0;
            Dictionary<string, double> projectWWRs = new Dictionary<string, double>();
            targetedProjectWWR = 0.0;
            targetedNorthWWR = 0.0;
            targetedEastWWR = 0.0;
            targetedSouthWWR = 0.0;
            targetedWestWWR = 0.0;

            //these are some examples
            if (location == "California" && year > 2010 && year <= 2013)
            {
                energyCode = "Title24";
                energyCodeYear = "2010";
                energyCodeVersion = "1";
                ventilationCode = "ASH62+Title24";
                ventilationCodeYear = "2007";
                ventilationcodeVersion = "1";
            }
            else if (location == "LEED" && year > 2010 && year <= 2013)
            {
                energyCode = "ASHRAE90-1";
                energyCodeYear = "2007";
                energyCodeVersion = "1";
                ventilationCode = "ASHRAE62-1";
                ventilationCodeYear = "2007";
                ventilationcodeVersion = "1";
            }
        }

        double currentProjectWWR;
        double currentNorthWWR;
        double currentEastWWR;
        double currentSouthWWR;
        double currentWestWWR;
        Dictionary<string, double> projectWWRs = new Dictionary<string, double>();
        double targetedProjectWWR;
        double targetedNorthWWR;
        double targetedEastWWR;
        double targetedSouthWWR;
        double targetedWestWWR;

        //relevant codes
        string energyCode;
        string energyCodeYear;
        string energyCodeVersion;
        string ventilationCode;
        string ventilationCodeYear;
        string ventilationcodeVersion;
    }

    public class SupportFiles
    {
       
    }

    public class OrientationUtil
    {
        public static bool isNorth(BuildingObjects.OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isNorth = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (((opening.Azimuth + globalAzimuthRotation) >= 315 && (opening.Azimuth + globalAzimuthRotation) <=360) ||
                ((opening.Azimuth + globalAzimuthRotation) >=0 && (opening.Azimuth + globalAzimuthRotation) < 45))
            {
                isNorth = true;
            }
            return isNorth;
        }

        public static bool isEast(BuildingObjects.OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isEast = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 45 && opening.Azimuth + globalAzimuthRotation < 135)
            {
                isEast = true;
            }
            return isEast;
        }

        public static bool isSouth(BuildingObjects.OpeningDefinitions opening,double globalAzimuthRotation)
        {
            bool isSouth = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 135 && opening.Azimuth + globalAzimuthRotation < 225)
            {
                isSouth = true;
            }
            return isSouth;
        }

        public static bool isWest(BuildingObjects.OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isWest = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 215 && opening.Azimuth + globalAzimuthRotation < 315)
            {
                isWest = true;
            }
            return isWest;
        }

        public static bool isNorth(BuildingObjects.MemorySafe_OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isNorth = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (((opening.Azimuth + globalAzimuthRotation) >= 315 && (opening.Azimuth + globalAzimuthRotation) <= 360) ||
                ((opening.Azimuth + globalAzimuthRotation) >= 0 && (opening.Azimuth + globalAzimuthRotation) < 45))
            {
                isNorth = true;
            }
            return isNorth;
        }

        public static bool isEast(BuildingObjects.MemorySafe_OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isEast = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 45 && opening.Azimuth + globalAzimuthRotation < 135)
            {
                isEast = true;
            }
            return isEast;
        }

        public static bool isSouth(BuildingObjects.MemorySafe_OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isSouth = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 135 && opening.Azimuth + globalAzimuthRotation < 225)
            {
                isSouth = true;
            }
            return isSouth;
        }

        public static bool isWest(BuildingObjects.MemorySafe_OpeningDefinitions opening, double globalAzimuthRotation)
        {
            bool isWest = false;
            //assumes that globalAzimuthRotation is positive counter clockwise
            if (opening.Azimuth + globalAzimuthRotation >= 215 && opening.Azimuth + globalAzimuthRotation < 315)
            {
                isWest = true;
            }
            return isWest;
        }
    }
}
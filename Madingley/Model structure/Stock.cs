using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace Madingley
{
    /// <summary>
    /// Hold individual stocks
    /// </summary>
    public class Stock
    {
        /// <summary>
        /// The index of the functional group that the stock belongs to
        /// </summary>
        private byte _FunctionalGroupIndex;
        /// <summary>
        /// Get and set the functional group that the stock belongs to
        /// </summary>
        public byte FunctionalGroupIndex { get { return _FunctionalGroupIndex; } }

        /// <summary>
        /// The mean body mass of an individual in this stock
        /// </summary>
        private double _IndividualBodyMass;
        /// <summary>
        /// Get and set the mean body mass of an individual in this stock
        /// </summary>
        public double IndividualBodyMass
        {
            get { return _IndividualBodyMass; }
            set { _IndividualBodyMass = value; }
        }

        /// <summary>
        /// The total biomass of the stock
        /// </summary>
        private double _TotalBiomass;
        /// <summary>
        /// Get and set the total biomass of this stock
        /// </summary>
        public double TotalBiomass
        {
            get { return _TotalBiomass; }
            set { _TotalBiomass = value; }
        }

        // The name of the stock
        private string _StockName;

        public string StockName
        {
            get { return _StockName; }
            set { _StockName = value; }
        }


        private double _FractionalArea;

        public double FractionalArea
        {
            get { return _FractionalArea; }
            set { _FractionalArea = value; }
        }

        private double[] _SizeStructure;
        public double[] SizeStructure
        {
            get { return _SizeStructure; }
            set { _SizeStructure = value; }
        }
        private double[] _SizeBinEdges;
        public double[] SizeBinEdges
        {
            get { return _SizeBinEdges; }
            set { _SizeBinEdges = value; }
        }
        private double[] _SizeBinCentres;
        public double[] SizeBinCentres
        {
            get { return _SizeBinCentres; }
            set { _SizeBinCentres = value; }
        }
        /// <summary>
        /// Constructor for stock class. Assigns stock starting properties
        /// </summary>
        /// <param name="functionalGroupIndex">The functional group index of the stock being generated</param>
        /// <param name="individualMass">The individual mass of the stock</param>
        /// <param name="initialTotalBiomass">The initial total biomass of the stock</param>
        //public Stock(byte functionalGroupIndex, double individualMass, double initialTotalBiomass)
        //{
        //    _FunctionalGroupIndex = functionalGroupIndex;
        //    _IndividualBodyMass = individualMass;
        //    _TotalBiomass = initialTotalBiomass;
            
        //}

        /// <summary>
        /// Constructor for stock class. Assigns stock starting properties
        /// </summary>
        /// <param name="functionalGroupIndex">The functional group index of the stock being generated</param>
        /// <param name="individualMass">The individual mass of the stock</param>
        /// <param name="initialTotalBiomass">The initial total biomass of the stock</param>
        /// <param name="stockName">The name of the stock</param>
        public Stock(byte functionalGroupIndex, double individualMass, double initialTotalBiomass, double fractionalArea, string stockName)
        {
            _FunctionalGroupIndex = functionalGroupIndex;
            _IndividualBodyMass = individualMass;
            _TotalBiomass = initialTotalBiomass;
            _FractionalArea = fractionalArea;
            _StockName = stockName;

        }

        public Stock(Stock s)
        {
            _FunctionalGroupIndex = s._FunctionalGroupIndex;
            _IndividualBodyMass = s._IndividualBodyMass;
            _TotalBiomass = s._TotalBiomass;
            _FractionalArea = s.FractionalArea;
            _StockName = s.StockName;
        }


        public void SetSizeDistribution(double[] binCentres, double[] binEdges, double[] nsfCentres)
        {
            _SizeStructure = nsfCentres;
            _SizeBinEdges = binEdges;
            _SizeBinCentres = binCentres;
        }

    }
}

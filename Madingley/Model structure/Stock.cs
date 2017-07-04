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

        /// <summary>
        /// Set and get the biomass from the NPZ model
        /// </summary>
        private double[] _BiomassFromNPZ;
        public double[] BiomassFromNPZ { get { return _BiomassFromNPZ; } }
        

        /// <summary>
        /// Constructor for stock class. Assigns stock starting properties
        /// </summary>
        /// <param name="functionalGroupIndex">The functional group index of the stock being generated</param>
        /// <param name="individualMass">The individual mass of the stock</param>
        /// <param name="initialTotalBiomass">The initial total biomass of the stock</param>
        public Stock(byte functionalGroupIndex, double individualMass, double initialTotalBiomass)
        {
            _FunctionalGroupIndex = functionalGroupIndex;
            _IndividualBodyMass = individualMass;
            _TotalBiomass = initialTotalBiomass;
            
        }

        /// <summary>
        /// Constructor for stock class. Assigns stock starting properties
        /// </summary>
        /// <param name="functionalGroupIndex">The functional group index of the stock being generated</param>
        /// <param name="individualMass">The individual mass of the stock</param>
        /// <param name="initialTotalBiomass">The initial total biomass of the stock</param>
        /// <param name="stockName">The name of the stock</param>
        public Stock(byte functionalGroupIndex, double individualMass, double initialTotalBiomass, string stockName)
        {
            _FunctionalGroupIndex = functionalGroupIndex;
            _IndividualBodyMass = individualMass;
            _TotalBiomass = initialTotalBiomass;
            _StockName = stockName;
        }

        public Stock(Stock s)
        {
            _FunctionalGroupIndex = s._FunctionalGroupIndex;
            _IndividualBodyMass = s._IndividualBodyMass;
            _TotalBiomass = s._TotalBiomass;
        }
    }
}

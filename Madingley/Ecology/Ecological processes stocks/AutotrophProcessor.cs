using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using MathNet.Numerics;

namespace Madingley
{
    /// <summary>
    /// Class for converting primary productivity estimates to autotroph biomass
    /// </summary>
    public class AutotrophProcessor
    {
        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        private UtilityFunctions Utilities;

        // Conversion ratio for phytoplankton from grams carbon to grams wet weight
        /// <summary>
        /// Factor to convert phytoplankton biomass from grams carbon to grams wet weight
        /// </summary>
        /// <remarks>Currently derived from Ho et al. (2003) J. Phycol., Dalsgaard and Pauly (1997) and Strickland (1966)</remarks>
        private double _PhytoplanktonConversionRatio;
        /// <summary>
        /// Get the conversion ratio for phytoplankton from grams carbon to grams wet weight
        /// </summary>
        public double PhytoplanktonConversionRatio { get { return _PhytoplanktonConversionRatio; } }

        /// <summary>
        /// Factor to convert NPP from units per m^2 to units per km^2
        /// </summary>
        private const double _MsqToKmSqConversion = 1000000.0;
        /// <summary>
        /// Get the factor to convert NPP from units per m^2 to units per km^2
        /// </summary>
        public double MsqToKmSqConversion { get { return _MsqToKmSqConversion; } }

        /// <summary>
        /// Constructor for the autotroph processor: initialises necessary classes
        /// </summary>
        public AutotrophProcessor()
        {

            _PhytoplanktonConversionRatio = EcologicalParameters.Parameters["AutotrophProcessor.ConvertNPPtoAutotroph.PhytoplanktonConversionRatio"];

            // Initialise the utility functions
            Utilities = new UtilityFunctions();
        }

        /// <summary>
        /// Calculate coefficients for the abundance-size slope
        /// </summary>
        /// <param name="temperature"></param>
        /// <param name="no3"></param>
        /// <returns></returns>
        public Tuple<double, double> GetPhytoAbundSizeSlope(double temperature, double no3)
        {
            double[] phyto = new double[3];
            double[] nsfPhyto = new double[3];
            double[] phytoSize = new double[3] { Math.Log10(1E-7), Math.Log10(1E-5), Math.Log10(1E-3) };
            double[] regCoefs = new double[2];

            phyto[0] = 1.145 - 0.021 * no3 - 6.936E-6 * temperature;
            phyto[1] = 1.146 + 0.013 * no3 - 0.064 * temperature;
            phyto[2] = 0.804 - 0.002 * no3 - 0.077 * temperature;
            nsfPhyto[0] = Math.Log10(Math.Exp(phyto[0]) / (Math.Exp(phyto[0]) + Math.Exp(phyto[1]) + Math.Exp(phyto[2])));
            nsfPhyto[1] = Math.Log10(Math.Exp(phyto[1]) / (Math.Exp(phyto[0]) + Math.Exp(phyto[1]) + Math.Exp(phyto[2])));
            nsfPhyto[2] = Math.Log10(Math.Exp(phyto[2]) / (Math.Exp(phyto[0]) + Math.Exp(phyto[1]) + Math.Exp(phyto[2])));

            Tuple<double, double> linReg = Fit.Line(phytoSize, nsfPhyto);
            regCoefs[0] = linReg.Item1; // intercept
            regCoefs[1] = linReg.Item2; // regression coefficient
            return Tuple.Create(regCoefs[0], regCoefs[1]);
        }


        public double[] GetPhytoDistributionEnvironment(double temperature, double no3,double[] BinEdges,
            double[] BinCentres)
        {
            Tuple<double, double> regCoefs = GetPhytoAbundSizeSlope(temperature, no3);

            double[] NsfEdges = new double[BinEdges.Length];
            double[] NsfCentres = new double[BinCentres.Length];

            //int b = 0;
            //for (double i = -8; i < 0; i++, b++)
            //{
            //    NsfEdges[b] = Math.Pow(10, regCoefs.Item2 * BinEdges[b] + regCoefs.Item1);
            //}

            //Calculate the mean NSF value (in log space) because bins are uniform and 
            for (int b =0; b <BinCentres.Length; b++)
            {
                NsfCentres[b] = Math.Pow(10, regCoefs.Item2 * BinCentres[b] + regCoefs.Item1);//(NsfEdges[b] + NsfEdges[b + 1]) / 2.0;
            }

            //Normalise so that the sum of Nsf centres = 1
            double TotalNsf = NsfCentres.Sum();

            for (int b = 0; b < NsfCentres.Length; b++)
            {
                NsfCentres[b] = NsfCentres[b] / TotalNsf;
            }

            return NsfCentres;
        }


        
        /// <summary>
        /// Convert NPP estimate into biomass of an autotroph stock
        /// </summary>
        /// <param name="cellEnvironment">The environment of the current grid cell</param>
        /// <param name="gridCellStockHandler">The stock handler for the current stock</param>
        /// <param name="actingStock">The location of the stock to add biomass to</param>
        /// <param name="terrestrialNPPUnits">The units of the terrestrial NPP data</param>
        /// <param name="oceanicNPPUnits">The units of the oceanic NPP data</param>
        /// <param name="currentTimestep">The current model time step</param>
        /// <param name="GlobalModelTimeStepUnit">The time step unit used in the model</param>
        /// <param name="trackProcesses">Whether to output data describing the ecological processes</param>
        /// <param name="globalTracker">Whether to output data describing the global-scale environment</param>
        /// <param name="outputDetail">The level of output detail to use for the outputs</param>
        /// <param name="specificLocations">Whether the model is being run for specific locations</param>
        /// <param name="currentMonth">The current month in the model run</param>
        public double ConvertNPPToAutotroph(FunctionalGroupDefinitions cohortDefinitions, FunctionalGroupDefinitions stockDefinitions,
            SortedList<string,double[]> cellEnvironment, GridCellStockHandler gridCellStockHandler, int[] 
            actingStock, string terrestrialNPPUnits, string oceanicNPPUnits, uint currentTimestep, string GlobalModelTimeStepUnit,
            ProcessTracker trackProcesses, FunctionalGroupTracker functionalTracker, HighResFGTracker highResFGTracker, GlobalProcessTracker globalTracker, string outputDetail, bool specificLocations,uint currentMonth)
        {
            double NPP = new double();

            // Check that this is an ocean cell
            if(cellEnvironment["Realm"][0] == 2.0)
            {
                switch (gridCellStockHandler[actingStock].StockName)
                {
                    case "oceannpp":
                        // Get NPP from cell environment
                        NPP = cellEnvironment["NPP"][currentMonth];

                        // If NPP is a missing value then set to zero
                        if (NPP == cellEnvironment["Missing Value"][0]) NPP = 0.0;

                        break;

                    default:
                        Debug.Fail("Can't find phytoplankton name in stocks");
                        Console.WriteLine("Can't find phytoplankton name in stocks");
                        break;
                }

                // Check that the units of oceanic NPP are gC per m2 per day
                Debug.Assert(oceanicNPPUnits == "gC/m2/day", "Oceanic NPP data are not in the correct units for this formulation of the model");

                // Convert to g/cell/month
                NPP *= _MsqToKmSqConversion;

                // Multiply by cell area to get g/cell/day
                NPP *= cellEnvironment["Cell Area"][0];

                // Convert to g wet matter, assuming carbon content of phytoplankton is 10% of wet matter
                NPP *= _PhytoplanktonConversionRatio;

                // Finally convert to g/cell/month and add to the stock totalbiomass
                NPP *= Utilities.ConvertTimeUnits(GlobalModelTimeStepUnit, "day");
                
            }

            // Else if neither on land or in the ocean
            else
            {
                Debug.Fail("This is not a marine cell!");
                // Set the autotroph biomass to zero
                gridCellStockHandler[actingStock].TotalBiomass = 0.0;
            }
            Debug.Assert(gridCellStockHandler[actingStock].TotalBiomass >= 0.0, "stock negative");

            return NPP;
        }
    }
}

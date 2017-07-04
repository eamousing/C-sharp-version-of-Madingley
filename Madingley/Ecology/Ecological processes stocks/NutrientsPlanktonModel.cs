using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Timing;

namespace Madingley
{
    /// <summary>
    /// Class for running the nutrients-plankton model (NPZ model)
    /// </summary>
    class NutrientsPlanktonModel
    {
        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        private UtilityFunctions Utilities;

        /// <summary>
        /// An instance of StopWatch to time individual time steps
        /// </summary>
        private StopWatch NPZTimer;

        /// <summary>
        /// Set and get the number of size levels (stocks) in the NPZ model
        /// </summary>
        private int _nLevels;
        public int nLevels { get { return _nLevels; } }

        /// <summary>
        /// Set and get the internal timestep for the NPZ model
        /// </summary>
        private double _dt;
        public double dt { get { return _dt; } }

        /// <summary>
        /// Set and get the number of days to run the NPZ model. 
        /// This is the same as the global model time step (i.e. about 30 days if run monthly).
        /// </summary>
        private double _maxTime;
        public double maxTime { get { return _maxTime; } }

        /// <summary>
        /// Set and get the cell carbon quota factor (Lovdal et al.)
        /// </summary>
        private double _cellCarbonQuotaFactor;
        public double cellCarbonQuotaFactor { get { return _cellCarbonQuotaFactor; } }

        /// <summary>
        /// Set and get the cell carbon quota exponent (Lovdal et al.)
        /// </summary>
        private double _cellCarbonQuotaExponent;
        public double cellCarbonQuotaExponent { get { return _cellCarbonQuotaExponent; } }

        /// <summary>
        /// Set and get the constant mortality term (metabolism/respiration). 
        /// In the future this will be changed to be size dependent
        /// </summary>
        private double _constantMortalityTerm;
        public double constantMortalityTerm { get { return _constantMortalityTerm; } }

        /// <summary>
        /// Constructor for the NPZ model processor: initialises necessary variables and classes
        /// </summary>
        public NutrientsPlanktonModel()
        {
            _cellCarbonQuotaFactor = NutrientPlanktonModelParameters.NPZParameters["aCell"];
            _cellCarbonQuotaExponent = NutrientPlanktonModelParameters.NPZParameters["bCell"];
            _constantMortalityTerm = NutrientPlanktonModelParameters.NPZParameters["kResp"];
            _nLevels = 3; // Hard coded for now. Has to be the number of stocks defined in the stocksDefinitionFile.
            _dt = 0.01; // Corresponds to a time scale of approximately 7-8 hours.
            _maxTime = 365; // Hard coded for now. Has to be the number of the days in the Madingley model time step.

            // Initialize the utility functions
            Utilities = new UtilityFunctions();
            NPZTimer = new StopWatch();
        }


        public double[] RunNPZModel(SortedList<string, double[]> cellEnvironment, uint currentMonth)
        {
            // Get nutrient supply rate from the environment
            double supplyNitrate = cellEnvironment["microNPP"][currentMonth]; // nitrate supply layer missing. Using microNPP as a proxy.

            // If the nutrient supply rate is a missing value, then set to zero
            if (supplyNitrate == cellEnvironment["Missing Value"][0]) supplyNitrate = 0.0;

            // Calculate the maximum number of steps in the NPZ model
            int nStepMax = Convert.ToInt32(Math.Round(maxTime / dt));
            int nStepOut = nStepMax / 10;

            // Create working arrays
            double[] biomassN = new double[nLevels]; // Total biomass in nitrate units
            for(int i = 0; i < nLevels; i++)
            {
                biomassN[i] = 1.0;
            }
            double[] dBiomassNDt = new double[nLevels]; // Changes in nitrate biomass per time step
            double[] biomassC = new double[nLevels]; // Total biomass in carbon units.
            double[] biomassCPerCell = new double[nLevels]; // Carbon biomass per cell

            // Create arrays for outputting
            int nStepOutMax = nStepMax / nStepOut;
            double[,] biomassNOut = new double[nLevels, nStepOutMax];
            double[] biomassNOutFinal = new double[nLevels];

            // Set cell volumes for all size levels/stocks (import from stock definition file)
            double[] cellVol = new double[nLevels];
            //for(int i = 0; i < nLevels; i++)
            //{
            //    cellVol[i] = Convert.ToDouble(stockDefinitions.GetTraitValuesAllFunctionalGroups("individual mass")[i]);
            //}
            double[] sizeExponents = new double[6] { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };
            for(int i = 0; i < nLevels; i++)
            {
                cellVol[i] = Math.Pow(10.0, sizeExponents[i]);
            }

            // Calculate carbon cell quotas for each size group
            double[] quotaCarbon = new double[nLevels];
            for(int i = 0; i < quotaCarbon.Length; i++)
            {
                quotaCarbon[i] = cellCarbonQuotaFactor * Math.Pow(cellVol[i], cellCarbonQuotaExponent);
            }

            // Calculate nitrate cell quotas for each size group based on Redfield proportions
            double[] quotaNitrate = new double[nLevels];
            for(int i = 0; i < quotaNitrate.Length; i++)
            {
                quotaNitrate[i] = quotaCarbon[i] * (16.0 / 106.0);
            }

            // Convert nitrate quota to micromoles nitrate per cell
            // Ask Mick about where 1.0e15 comes from.
            for(int i = 0; i < quotaNitrate.Length; i++)
            {
                quotaNitrate[i] = quotaNitrate[i] * 1.0e6 / 1.0e15;
            }

            // Calculate nitrate maximum uptake rate (micromoles N cell-1 day-1; Parametes from Litchman et al., 2007)
            double[] vmaxN = new double[nLevels];
            for (int i = 0; i < vmaxN.Length; i++)
            {
                vmaxN[i] = 9.1e-9 * Math.Pow(cellVol[i], 0.67);
            }

            // Calculate half-saturation constant for nitrate uptake (micromoles N; parameters from Litchman et al., 2007)
            double[] knN = new double[nLevels];
            for(int i = 0; i < knN.Length; i++)
            {
                knN[i] = 0.17 * Math.Pow(cellVol[i], 0.27);
            }

            // Set mortality term/background respiration/metabolism to be the same for all organisms
            // In the future metabolism have to be size dependent.
            double[] kMort = new double[nLevels];
            for(int i = 0; i < kMort.Length; i++)
            {
                kMort[i] = constantMortalityTerm;
            }

            // Setup some more working arrays and variables
            double[] autotrophy = new double[nLevels];
            double[] respiration = new double[nLevels];
            int imax = nLevels;

            // Initiate time related variables within the loop.
            int nCount = 0;
            int nOut = 0;
            double time = 0.0;
            double dNitrateDt = 0.0;
            double[] prod = new double[nLevels];

            // Re-initialize biomass (micromoles N per liter)
            // Here we set the biomass to be equal the biomass at the end of the former step.
            for(int i = 0; i < biomassC.Length; i++)
            {
                //int[] actingStock = new int[2] { i, 1 };
                //biomassC[i] = gridCellStocks[actingStock].TotalBiomass;
                biomassN[i] = biomassN[i] * 1.0e-2; // Temporary set as a constant for testing purposes
            }

            // Re-initialize nitrate concentration
            // Here we set the nitrate concentration to be equal to the nitrate concentration
            // at the end of the former step
            double nitrateConcentration = 1.0; // Temporary set as a constant for testing purposes


            // Start the main time loop
            for(int nStep = 0; nStep < nStepMax; nStep++)
            {
                // Evaluate autotrophy
                for(int i = 0; i < autotrophy.Length; i++)
                {
                    autotrophy[i] = (vmaxN[i] / quotaNitrate[i]) * (nitrateConcentration / (nitrateConcentration + knN[i])) * biomassN[i];
                }

                // Maintain respiration/background loss
                for(int i = 0; i < respiration.Length; i++)
                {
                    respiration[i] = kMort[i] * biomassN[i];
                }

                // To avoid numerical problems, do not let biomass decline below a low threshold value
                for(int i = 0; i < nLevels; i++)
                {
                    if(biomassN[i] < 1.0e-25)
                    {
                        respiration[i] = 0.0;
                    }
                }

                // Evaluate rates of change
                dNitrateDt = -autotrophy.Sum() + supplyNitrate;
                for(int i = 0; i < dBiomassNDt.Length; i++)
                {
                    dBiomassNDt[i] = autotrophy[i] - respiration[i];
                }

                // Euler forward step
                double[] biomassNNew = new double[nLevels];
                double nitrateNew;
                for(int i = 0; i < biomassNNew.Length; i++)
                {
                    biomassNNew[i] = biomassN[i] + dBiomassNDt[i] * dt;
                    biomassN[i] = biomassNNew[i];
                }

                nitrateNew = nitrateConcentration + dNitrateDt * dt;
                nitrateConcentration = nitrateNew;

                // Increment time
                time = time + dt;
            }
            // End of main time loop

            // Save final biomasses
            for(int i = 0; i < nLevels; i++)
            {
                biomassNOutFinal[i] = biomassN[i];
            }

            // Calculate and return carbon biomass for each stock
            double[] nitrateToCarbon = new double[nLevels];
            for (int i = 0; i < nLevels; i++)
            {
                // Convert from micromoles nitrate per L to carbon
                nitrateToCarbon[i] = biomassNOutFinal[i] * (106.0 / 16.0);
                // Convert to grams carbon per cubic meter
                biomassC[i] = nitrateToCarbon[i] * 12.0107 * 1.0e3 * 1.0e-6;
                // Convert to grams carbon per square kilometer (assuming a one layer system where m3 = m2)
                biomassC[i] *= 1000000;
                // Calculate grams per cell
                biomassC[i] *= cellEnvironment["Cell Area"][0];
                // Increase for testing purposes
                biomassC[i] *= 1000;
            }

            Console.WriteLine("Plankton carbon biomass after 30 days");
            Console.WriteLine("[{0}]", string.Join(",", biomassC));

            return biomassC;

        }

        //public double BiomassConverter(string fromString, string toString, double fromDouble, double toDouble)
        //{
        //    switch (fromString)
        //    {
        //        case "nitrate":
        //            switch (toString)
        //            {
        //                case "carbon":
        //                    return fromDouble 
        //                    break;
        //                default:
        //                    break;
        //            }
        //            break;
        //        case "carbon":
        //        default:
        //            break;
        //    }
        //}

    }
}

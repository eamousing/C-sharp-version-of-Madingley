using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// Class for running the nutrients-phytoplankton-zooplankton (NPZ model)
    /// </summary>
    class NutrientsPlanktonModel
    {

        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        /// 
        private UtilityFunctions Utilities;

        // Define variables

        // Set number of trophic levels to model
        int nLevels = 3;

        // Set time step parameters
        double dt = 0.01;
        double maxTime = 30;
        
        
        // Set ecological parameters

        // Parameters for calculating nitrate and carbon cell quotas
        double aCell = 18.7;
        double bCell = 0.89;
        // Parameters for calculating grazing rates
        double aGraz = 0.5;
        double bGraz = -0.16;

        // Set nutrient supply rate (moles nitrate per square meter per day)
        // double and trible check the units!
        double sNitrate = 0.0;

        /// <summary>
        /// Factor to convert NPP from units per m^2 to units per km^2
        /// </summary>
        private const double _MsqToKmSqConversion = 1000000.0;
        /// <summary>
        /// Get the factor to convert NPP from units per m^2 to units per km^2
        /// </summary>
        public double MsqToKmSqConversion { get { return _MsqToKmSqConversion; } }


        /// <summary>
        /// Constructor for the NPZ model processor: initialises necessary classes
        /// </summary>
        public NutrientsPlanktonModel()
        {
            Utilities = new UtilityFunctions();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="cellEnvironment"></param>
        /// <param name="currentMonth"></param>
        /// <returns></returns>
        public double[] RunNPZModel(SortedList<string, double[]> cellEnvironment, uint currentMonth)
        {

            // Get nutrient supply rate from the cell environment
                sNitrate = cellEnvironment["microNPP"][currentMonth];

                // if nutrient supply rate is a missing value, then set to zero
                if (sNitrate == cellEnvironment["Missing Value"][0]) sNitrate = 0.0;


                int nStepMax = Convert.ToInt32(Math.Round(maxTime / dt));
                int nStepOut = nStepMax / 10;

                // Create working arrays
                double[] biomassN = new double[nLevels];
                for (int i = 0; i < nLevels; i++)
                {
                    biomassN[i] = 1.0;
                }
                double[] dBiomassNDt = new double[nLevels];
                double[] biomassC = new double[nLevels];
                double[] biomassCPerCell = new double[nLevels];

                // Create output array
                int nStepOutMax = nStepMax / nStepOut;
                double[,] biomassN_Out = new double[nLevels, nStepOutMax];
                double[] biomassNout = new double[nLevels];

                // Define cell volumes (for now they are hard coded)
                double[] cellVol = new double[nLevels];
                cellVol[0] = 4.2; // cell diameter = 2 um
                cellVol[1] = 113.1; // cell diameter = 6 um
                cellVol[2] = 1767.1; // cell diameter = 15 um
                //cellVol[3] = 14377.2; // cell diameter = 30 um
                //cellVol[4] = 381703.5; // cell diameter = 90 um
                //cellVol[5] = 904778.7; // cell diameter = 120 um

                // Print cell volumes to console
                Console.WriteLine("Setting cell volumes to:");
                Console.WriteLine("[{0}]", string.Join(", ", cellVol));

                // Create arrays for cell quotas
                double[] quotaCarbon = new double[nLevels];
                double[] quotaNitrate = new double[nLevels];

                // Calculate carbon cell quota for each size group
                for (var i = 0; i < quotaCarbon.Length; i++)
                {
                    quotaCarbon[i] = aCell * Math.Pow(cellVol[i], bCell);
                }

                // Calculate nitrate cell quota for each size group based on Redfield proportions
                for (var i = 0; i < quotaNitrate.Length; i++)
                {
                    quotaNitrate[i] = quotaCarbon[i] * (16.0 / 106.0);
                }

                // Convert nitrate quota to micromoles N cell-1
                // ask Mick about where 1.0e15 comes from.
                for (var i = 0; i < quotaNitrate.Length; i++)
                {
                    quotaNitrate[i] = quotaNitrate[i] * 1.0e6 / 1.0e15;
                }

                // Create arrays for autotrophic traits
                double[] vmaxN = new double[nLevels];
                double[] kn = new double[nLevels];
                double[] kResp = new double[nLevels];
                double[] specVmaxN = new double[nLevels];

                // Calculate growth rates from cell volumes
                for (int i = 0; i < vmaxN.Length; i++)
                {
                    vmaxN[i] = 9.1e-9 * Math.Pow(cellVol[i], 0.67);
                }

                // Calculate half-saturation constants from cell volumes
                for (int i = 0; i < kn.Length; i++)
                {
                    kn[i] = 0.17 * Math.Pow(cellVol[i], 0.27);
                }

                // Set mortality term to be the same for all types
                for (var i = 0; i < kResp.Length; i++)
                {
                    kResp[i] = 0.03;
                }

                // Initiate arrays for heterotrophic traits
                double[,] gmax = new double[nLevels, nLevels];
                double[,] gmax1 = new double[nLevels, nLevels];
                double[,] kbN = new double[nLevels, nLevels];
                double[,] gamma = new double[nLevels, nLevels];
                double[,] eye = new double[nLevels, nLevels];
                double[,] grazInteract = new double[nLevels, nLevels];

                // Populate grazing interaction matrix
                // All plankton are mixotrophic. One prey each.
                //grazInteract[0, 1] = 1.0;
                //grazInteract[0, 2] = 0.2;
                //grazInteract[1, 2] = 1.0;
                //grazInteract[1, 3] = 0.2;
                //grazInteract[2, 3] = 1.0;
                //grazInteract[2, 4] = 0.2;
                //grazInteract[3, 4] = 1.0;
                //grazInteract[3, 5] = 0.2;
                //grazInteract[4, 5] = 1.0;

                // Calculate maximum grazing rate for all types
                for (var i = 0; i < gmax.GetLength(0); i++)
                {
                    for (var j = 0; j < gmax.GetLength(1); j++)
                    {
                        gmax[i, j] = aGraz * Math.Pow(cellVol[j], bGraz) * grazInteract[i, j];
                        eye[i, j] = i;
                    }
                }

                // Set gmax1 = grazing eate for Holling I ((mol N l-1 day)-1)
                for (var i = 0; i < gmax1.GetLength(0); i++)
                {
                    for (var j = 0; j < gmax1.GetLength(1); j++)
                    {
                        // Ask Mick about this.. gmax and gmax1 are identical(?)
                        gmax1[i, j] = gmax[i, j] / 1.0;
                    }
                }


                // Calculate half-saturation constant for grazing (from Ward et al. (2013))
                for (var i = 0; i < kbN.GetLength(0); i++)
                {
                    for (var j = 0; j < kbN.GetLength(1); j++)
                    {
                        kbN[i, j] = 0.5 * (16.0 / 106.0) * grazInteract[i, j]; // micromol N biomass L-1 
                    }
                }

                // Calculate trophic transfer efficiency
                for (var i = 0; i < gamma.GetLength(0); i++)
                {
                    for (var j = 0; j < gamma.GetLength(1); j++)
                    {
                        gamma[i, j] = 0.6 * grazInteract[i, j];
                    }
                }

                // Setup some more working arrays and variables
                double[] heterotrophy = new double[nLevels];
                double[] autotrophy = new double[nLevels];
                double[] predation = new double[nLevels];
                double[] respiration = new double[nLevels];
                int imax = nLevels;
                int jmax = nLevels;
                int kmax = nLevels;

                // Initiate time related variables within the loop
                int nCount = 0;
                int nOut = 0;
                double time = 0.0;
                double dBmassNDtOld = 0.0;
                double dNitrateDtOld = 0.0;
                double dNitrateDt = 0.0;
                double[] prod = new double[nLevels];

                // Re-initialize biomass (convert to micromoles N l-1??)
                for (int i = 0; i < biomassN.Length; i++)
                {
                    biomassN[i] = biomassN[i] * 1.0e-2;
                }

                // Set nitrate concentration
                double nitrate = 1.0;

                // Start main time loop
                for (int nStep = 0; nStep < nStepMax; nStep++)
                {
                    // Evaluate autotrophy
                    for (int i = 0; i < autotrophy.Length; i++)
                    {
                        autotrophy[i] = (vmaxN[i] / quotaNitrate[i]) * (nitrate / (nitrate + kn[i])) * biomassN[i];
                    }

                    // Maintain respiration/background loss
                    for (int i = 0; i < respiration.Length; i++)
                    {
                        respiration[i] = kResp[i] * biomassN[i];
                    }

                    // Evaluate heterotrophy
                    for (int i = 0; i < imax; i++)
                    {
                        // 
                        double sum1 = 0.0;

                        for (int j = 0; j < jmax; j++)
                        {
                            sum1 = sum1 + (gamma[j, i] * gmax1[j, i] * biomassN[i] * biomassN[j]);
                            heterotrophy[i] = sum1;
                        }

                        // Predation by all others (ask Mick)
                        double sum2 = 0.0;
                        for (int k = 0; k < kmax; k++)
                        {
                            sum2 = sum2 + (gmax1[i, k] * biomassN[k] * biomassN[i]);
                            predation[i] = sum2;
                        }
                    }


                    // To avoid numerical problems, do not let biomass decline below a low threshold value
                    for (int i = 0; i < nLevels; i++)
                    {
                        if (biomassN[i] < 1.0e-25)
                        {
                            respiration[i] = 0.0;
                            predation[i] = 0.0;
                        }
                    }

                    // Evaluate rates of change
                    dNitrateDt = -autotrophy.Sum() + sNitrate;

                    for (var i = 0; i < dBiomassNDt.Length; i++)
                    {
                        dBiomassNDt[i] = autotrophy[i] + heterotrophy[i] - respiration[i] - predation[i];
                    }

                    // Store total growth
                    for (var i = 0; i < prod.Length; i++)
                    {
                        prod[i] = autotrophy[i] + heterotrophy[i];
                    }

                    // Euler forward step
                    var bmassNNew = new double[nLevels];
                    double nitrateNew;

                    for (var i = 0; i < bmassNNew.Length; i++)
                    {
                        bmassNNew[i] = biomassN[i] + dBiomassNDt[i] * dt;
                        biomassN[i] = bmassNNew[i];
                    }

                    nitrateNew = nitrate + dNitrateDt * dt;
                    nitrate = nitrateNew;

                    // Increment time
                    time = time + dt;

                    // Print to screen and put diagnostics into an array
                    nCount = nCount + 1;
                    if (nCount == nStepOut)
                    {
                        // Print values to screen
                        string screenOut1 = "";
                        screenOut1 += "Time: " + Math.Round(time, 0)
                            + "; Nitrate: " + nitrate.ToString("0.00E0")
                            + "; Biomass grp0:" + biomassN[0].ToString("0.00E0")
                            + "; Biomass grp1:" + biomassN[1].ToString("0.00E0")
                            + "; Biomass grp2:" + biomassN[2].ToString("0.00E0");
                            //+ "; Biomass grp3:" + biomassN[3].ToString("0.00E0")
                            //+ "; Biomass grp4:" + biomassN[4].ToString("0.00E0");
                            //+ "; Biomass grp5:" + biomassN[5].ToString("0.00E0");

                        //Console.WriteLine(screenOut1);

                        // nCount is the counter to decide whether to write output - reset here
                        nCount = 0;

                        // nOut is the position in the output array
                        for (var i = 0; i < biomassN.Length; i++)
                        {
                            biomassN_Out[i, nOut] = biomassN[i];
                        }

                        // Increment the output counter
                        nOut = nOut + 1;

                    }
                }
                // End main time loop
                Console.WriteLine("end off loop");

                // Save final biomasses
                for (int i = 0; i < nLevels; i++)
                {
                    biomassNout[i] = biomassN_Out[i, (nStepOutMax - 1)];
                }
                Console.WriteLine("Plankton nitrate biomass after 30 days:");
                Console.WriteLine("[{0}]", string.Join(", ", biomassNout));

                // Calculate carbon biomass
                double[] nitrateToCarbon = new double[nLevels];

                for (int i = 0; i < nLevels; i++)
                {
                    // Calculate micromoles carbon per liter (from Redfield proportions)
                    nitrateToCarbon[i] = biomassNout[i] * (106.0 / 16.0);
                    // Convert to grams carbon per cubic meter
                    biomassC[i] = nitrateToCarbon[i] * 12.0107 * 1.0E-6 * 1.0e3;
                    // Calculate grams per square km (assuming m3 = m2)
                    biomassC[i] *= _MsqToKmSqConversion;
                    // Calculate grams per cell
                    biomassC[i] *= cellEnvironment["Cell Area"][0];
                // convert to per month
                biomassC[i] *= Utilities.ConvertTimeUnits("month", "day");

                }
                Console.WriteLine("Plankton carbon biomass after 30 days:");
                Console.WriteLine("[{0}]", string.Join(", ", biomassC));

                return biomassC;

            }
        
    }
}

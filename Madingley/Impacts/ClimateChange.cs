using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// Adjusts cell climate parameters to simulate the impacts of climate change
    /// </summary>
    class ClimateChange
    {
        /// <summary>
        /// Constructor for the climate change class
        /// </summary>
        public ClimateChange()
        {
        }

        public void ApplyTemperatureScenario(SortedList<string, double[]> cellEnvironment,
            Tuple<string, double, double> temperatureScenario, uint currentTimestep, uint currentMonth, uint burninSteps,
            uint impactSteps, Boolean impactCell)
        {
            if (impactCell)
            {
                if (temperatureScenario.Item1 == "no")
                {
                }
                else if (temperatureScenario.Item1 == "constant")
                {
                    // Check to see whether this time step is at the transition between burn-in and impact
                    if ((currentTimestep >= (burninSteps + 1)) && (currentTimestep <= (burninSteps + 12)))
                    {
                        // Increment temperature
                        cellEnvironment["Temperature"][currentMonth] += temperatureScenario.Item2;
                    }
                    else
                    {
                    }
                }
                else if (temperatureScenario.Item1 == "temporary")
                {
                    // Check to see whether this time step is at the transition between burn-in and impact
                    // or at the transition between impact and post-impat
                    if ((currentTimestep >= (burninSteps + 1)) && (currentTimestep <= (burninSteps + 12)))
                    {
                        // Increment temperature
                        cellEnvironment["Temperature"][currentMonth] += temperatureScenario.Item2;
                    }
                    else if ((currentTimestep >= (burninSteps + impactSteps + 1)) && (currentTimestep <= (burninSteps + impactSteps + 12)))
                    {
                        // Revert temperature to original value
                        cellEnvironment["Temperature"][currentMonth] -= temperatureScenario.Item2;
                    }
                    else
                    {
                    }
                }
                else if (temperatureScenario.Item1 == "escalating")
                {
                    // If this is the first time step, add items to the cell environment to store the original cell temperature
                    if (currentTimestep == 0)
                    {
                        cellEnvironment.Add("Original Temperature", new double[12]);
                        for (int m = 0; m < 12; m++)
                        {
                            cellEnvironment["Original Temperature"][m] = cellEnvironment["Temperature"][m];
                        }
                    }
                    // If the spin-up period has been completed, then increment cell temperature
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if (currentTimestep > burninSteps)
                    {
                        cellEnvironment["Temperature"][currentMonth] = Math.Min((cellEnvironment["Original Temperature"][currentMonth] + 5.0),
                            cellEnvironment["Temperature"][currentMonth] + ((((currentTimestep - burninSteps) / 12) + 1) *
                            temperatureScenario.Item2));
                        // cellEnvironment["Temperature"][currentMonth] += gridCellStocks[actingStock].TotalBiomass *
                        //     (Math.Min(5.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));
                    }
                }
                else if (temperatureScenario.Item1 == "pb")
                {
                    // If the spin-up period (100 years) has been completed, then increment cell temperature
                    // according to the number of time-steps that have elapsed since the spin-up ended, for a period of 100 years until the 
                    // maximum thermal increase is reached
                    if ((currentMonth == 0) && (currentTimestep >= (burninSteps)) && (currentTimestep < (burninSteps+impactSteps)))
                    {
                        //if ((currentTimestep > (1199)) && (currentTimestep < 2400))
                        //{
                        //cellEnvironment["Temperature"][currentMonth] += (temperatureScenario.Item2 / 1200) * (currentTimestep - 1199);


                        for (int ii = 0; ii < 12; ii++)
                        {
                            cellEnvironment["Temperature"][ii] = cellEnvironment["BaselineTemperature"][ii] + ((temperatureScenario.Item2 / impactSteps) * (currentTimestep+ii - burninSteps));
                        }

                        cellEnvironment["AnnualTemperature"][0] = cellEnvironment["Temperature"].Average();
                    }

                }
                else
                {
                    Debug.Fail("There is no method for the climate change (temperature) scenario specified");
                }

            }
        }
    }
}

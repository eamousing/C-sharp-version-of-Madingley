using System;
using System.Collections.Generic;

namespace Madingley
{
    //Instantiate one of these for each grid cell
    public class ApplyFishingCatches
    {
            /// <summary>
            /// Instance of the class to perform general functions
            /// </summary>
            private UtilityFunctions Utilities;

        // Slope of the selectivity function (Carozza et al. 2017)
        private double cZero = 16.7787;

        // Parameter to convert sigmoidal function in length to mass (Froese et al. 2013)
        private double delta2 = 3.0;

        public ApplyFishingCatches()
        {
                // Initialise the utility functions
                Utilities = new UtilityFunctions();
        }

        // Apply fishing catches to cohorts within a specfic functional gorup
        public double ApplyCatches(GridCellCohortHandler gridCellCohorts, FunctionalGroupDefinitions madingleyCohortDefinitions, SortedList<string, double[]> cellEnvironment, string globalModelTimeStepUnit, int latIndex, int lonIndex)
        {

            double TotalCatch;
            double CatchDeficit = 0.0;

            // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
            double SaupTimeUnitConversion = Utilities.ConvertTimeUnits("year", globalModelTimeStepUnit);
            Console.WriteLine("Time change factor = :" + SaupTimeUnitConversion);
            Console.ReadKey();

            // Loop through the targeted functional groups
            // Small omnivorous ectotherms
            int[] SmallOmniEcto = madingleyCohortDefinitions.GetFunctionalGroupIndex("group description", "ecto omni mall it", false);

            if (SmallOmniEcto != null)
            {
                // Get total catch
                TotalCatch = cellEnvironment["SmallEctoCatch"][0];
                Console.WriteLine("Original total catch:" + TotalCatch);
                TotalCatch = ConvertSAUPUnitsToMadingley(TotalCatch, SaupTimeUnitConversion);
                Console.WriteLine("Revised total catch: " + TotalCatch);
                Console.ReadKey();

                // Run the fishing
                RunFishing(SmallOmniEcto, TotalCatch, madingleyCohortDefinitions, gridCellCohorts, 1.0, ref CatchDeficit);

                Console.WriteLine("Catch deficit: " + CatchDeficit);
                Console.ReadKey();

            } else
            {
                Console.WriteLine("Small omnivorous ectotherms group not found for ApplyFishingCatches class. Note that group description should be: ecto omni small it");
                Console.ReadKey();
            }

            // Medium omnivorous ectotherms
            int[] MedOmniEcto = madingleyCohortDefinitions.GetFunctionalGroupIndex("group description", "ecto omni medium it", false);

            if (MedOmniEcto != null)
            {
                // Get total catch
                TotalCatch = ConvertSAUPUnitsToMadingley(cellEnvironment["MedEctoCatch"][0], SaupTimeUnitConversion);

                // Run the fishing
                RunFishing(MedOmniEcto, TotalCatch, madingleyCohortDefinitions, gridCellCohorts, 0.5, ref CatchDeficit);

            }
            else
            {
                Console.WriteLine("Medium omnivorous ectotherms group not found for ApplyFishingCatches class. Note that group description should be: ecto omni medium it");
                Console.ReadKey();
            }

            // Large omnivorous ectotherms
            int[] LargeOmniEcto = madingleyCohortDefinitions.GetFunctionalGroupIndex("group description", "ecto omni large it", false);

            if (LargeOmniEcto != null)
            {
                // Get total catch
                TotalCatch = ConvertSAUPUnitsToMadingley(cellEnvironment["LgEctoCatch"][0], SaupTimeUnitConversion);

                // Run the fishing
                RunFishing(LargeOmniEcto, TotalCatch, madingleyCohortDefinitions, gridCellCohorts, 0.25, ref CatchDeficit);

            }
            else
            {
                Console.WriteLine("Large omnivorous ectotherms group not found for ApplyFishingCatches class. Note that group description should be: ecto omni large it");
                Console.ReadKey();
            }

            // Large omnivorous ectotherms
            int[] CarnEcto = madingleyCohortDefinitions.GetFunctionalGroupIndex("group description", "ecto carn it", false);

            if (CarnEcto != null)
            {
                // Get total catch
                TotalCatch = ConvertSAUPUnitsToMadingley(cellEnvironment["CarnEctoCatch"][0], SaupTimeUnitConversion);

                // Run the fishing
                RunFishing(CarnEcto, TotalCatch, madingleyCohortDefinitions, gridCellCohorts, 0.1, ref CatchDeficit);

            }
            else
            {
                Console.WriteLine("Carnivorous ectotherms group not found for ApplyFishingCatches class. Note that group description should be: ecto carn it");
                Console.ReadKey();
            }

            Console.WriteLine("Total catch deficit this time step: " + CatchDeficit);

            return CatchDeficit;
        }

        // Calculate the mass that can be caught from each cohort based on its biomass and the selectivity function
        private double CalculatePotentialMassCaught(GridCellCohortHandler gridCellCohorts, int functionalGroup, int cohortNumber, double maxFGBodyMass, double eMj)
        {
            // Get the body mass of individuals in this cohort
            double PreyBodyMass = gridCellCohorts[functionalGroup][cohortNumber].IndividualBodyMass;
            double PreyReproductiveBodyMass = gridCellCohorts[functionalGroup][cohortNumber].IndividualReproductivePotentialMass;

            if ((gridCellCohorts[functionalGroup][cohortNumber].CohortAbundance > 0) && (PreyBodyMass > 0))
            {
                // Calculate the fraction of individuals caught and their total biomass based on the selectivity given target body mass
                return gridCellCohorts[functionalGroup][cohortNumber].CohortAbundance * CalculateSelectivity(PreyBodyMass, eMj, maxFGBodyMass) * (PreyBodyMass + PreyReproductiveBodyMass);
            }
            else
                return 0.0;
        }

        // Apply potential mass caught
        private void ApplyCatch(GridCellCohortHandler gridCellCohorts, int functionalGroup, int cohortNumber, double maxFGBodyMass, double eMj, double CatchScaleFactor)
        {
            // Get the body mass of individuals in this cohort
            double PreyBodyMass = gridCellCohorts[functionalGroup][cohortNumber].IndividualBodyMass;
            double PreyReproductiveBodyMass = gridCellCohorts[functionalGroup][cohortNumber].IndividualReproductivePotentialMass;

            if ((gridCellCohorts[functionalGroup][cohortNumber].CohortAbundance > 0) && (PreyBodyMass > 0))
            {
                // Calculate the fraction of individuals caught and their total biomass based on the selectivity given target body mass
                gridCellCohorts[functionalGroup][cohortNumber].CohortAbundance -= gridCellCohorts[functionalGroup][cohortNumber].CohortAbundance * CalculateSelectivity(PreyBodyMass,  eMj, maxFGBodyMass) * CatchScaleFactor;
            }
        }

        private double CalculateCatchScaling(double totalModelCatch, double totalEmpiricalCatch, ref double catchDeficit)
        {
            double CatchScaleFactor;

            // Now adjust catch as appropriate between model catch and empirical catch
            if (totalModelCatch >= totalEmpiricalCatch)
            {
                CatchScaleFactor = totalEmpiricalCatch / totalModelCatch;
            }
            else
            {
                CatchScaleFactor = 1.0;
                catchDeficit += totalEmpiricalCatch - totalModelCatch;
            }

            return CatchScaleFactor;
        }

        private void RunFishing(int[] functionalGroupsToApplyCatch, double totalCatchToApply, FunctionalGroupDefinitions madingleyCohortDefinitions, GridCellCohortHandler gridCellCohorts, double dMj, ref double catchDeficit)
        {
            double CatchScaleFactor;
            double TotalModelCatch = 0.0;

            // Loop through to work out how much of each cohort should be removed based on the catch function. Note that the total can then either exceed or fall short of actual SAUP catch
            foreach (int FunctionalGroup in functionalGroupsToApplyCatch)
            {
                // Get the max body mass of this FG

                double maxFGBodyMass = madingleyCohortDefinitions.GetBiologicalPropertyOneFunctionalGroup("maximum mass", FunctionalGroup);

                // Loop over cohorts within the functional group and calculate potential catch
                for (int i = 0; i < gridCellCohorts[i].Count; i++)
                {
                    TotalModelCatch += CalculatePotentialMassCaught(gridCellCohorts, FunctionalGroup, i, maxFGBodyMass, dMj);
                }

                // Now adjust catch as appropriate between model catch and empirical catch
                CatchScaleFactor = CalculateCatchScaling(TotalModelCatch, totalCatchToApply, ref catchDeficit);

                // Loop over cohorts within the functional group and apply catch
                for (int i = 0; i < gridCellCohorts[i].Count; i++)
                {
                    ApplyCatch(gridCellCohorts, FunctionalGroup, i, maxFGBodyMass, dMj, CatchScaleFactor);
                }

                Console.WriteLine("Total model catch: " + TotalModelCatch);
                Console.ReadKey();
            }
        }
            
        // Convert SAUP catch units (tonnes / year / grid cell) to Madingley units (g / month (typically) / grid cell
        private double ConvertSAUPUnitsToMadingley(double SAUPCatch, double timeUnitConversion)
        {
            double MadingleyCatch;

                // Convert from yearly catches to monthly. Assume even distribution of catch over a year
                MadingleyCatch = SAUPCatch * timeUnitConversion;

            // Convert from tonnes to grams
                MadingleyCatch = MadingleyCatch * 1000000;

            return MadingleyCatch;
        }

        // Calculate the selectivity curve for a particular body mass 
        private double CalculateSelectivity(double cohortBodyMass, double dParameter, double maxBiomassFG)
        {

            double thresholdMass = dParameter * maxBiomassFG * 0.6198;

            // Claculate as per the general approach of Carozza et al. 2017
            return Math.Pow((1 + Math.Pow(cohortBodyMass / thresholdMass, -cZero/delta2)), -1);

        }


    }
}

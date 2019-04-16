using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.IO;

namespace Madingley
{
    /// <summary>
    /// Class for seeding grid cells with cohorts and stocks
    /// </summary>
    public class SeedGridCells
    {
        /// <summary>
        /// Instance of random number generator to take a time-dependent seed
        /// </summary>
       private NonStaticSimpleRNG RandomNumberGenerator = new NonStaticSimpleRNG();


        // Constructor
        public SeedGridCells()
        {
            ;
        }
        
        /// <summary>
        /// Seed grid cell with cohorts, as specified in the model input files
        /// </summary>
        /// <param name="functionalGroups">The functional group definitions for cohorts in the grid cell</param>
        /// <param name="cellEnvironment">The environment in the grid cell</param>
        /// <param name="globalDiagnostics">A list of global diagnostic variables</param>
        /// <param name="nextCohortID">YThe unique ID to assign to the next cohort produced</param>
        /// <param name="tracking">boolean to indicate if cohorts are to be tracked in this model</param>
        /// <param name="totalCellTerrestrialCohorts">The total number of cohorts to be seeded in each terrestrial grid cell</param>
        /// <param name="totalCellMarineCohorts">The total number of cohorts to be seeded in each marine grid cell</param>
        /// <param name="DrawRandomly">Whether the model is set to use random draws</param>
        /// <param name="ZeroAbundance">Set this parameter to 'true' if you want to seed the cohorts with zero abundance</param>
        public void SeedInitialGridCellCohorts(GridCellCohortHandler gridCellCohorts, ref FunctionalGroupDefinitions functionalGroups, ref SortedList<string, double[]>
            cellEnvironment, SortedList<string, double> globalDiagnostics, Int64 nextCohortID, Boolean tracking, double totalCellTerrestrialCohorts,
            double totalCellMarineCohorts, Boolean DrawRandomly, Boolean ZeroAbundance)
        {
            // Set the seed for the random number generator from the system time
            RandomNumberGenerator.SetSeedFromSystemTime();

            // Write out initial cohort information
           // StreamWriter tempsw = new StreamWriter("C://users//derekt//desktop//adult_juvenile_masses.txt");
           // tempsw.WriteLine("functional group\tadult mass\tjuvenilemass\tbiomass\tabundance\toptimal pbr");

            // Define local variables
            double CohortJuvenileMass;
            double CohortAdultMassRatio;
            double CohortAdultMass;
            double ExpectedLnAdultMassRatio;
            int[] FunctionalGroupsToUse;
            double NumCohortsThisCell;

            // Get the minimum and maximum possible body masses for organisms in each functional group
            double[] MassMinima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("minimum mass");
            double[] MassMaxima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("maximum mass");
            string[] NutritionSource = functionalGroups.GetTraitValuesAllFunctionalGroups("nutrition source");

            double[] ProportionTimeActive = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("proportion suitable time active");

            Int64 CohortIDIncrementer = nextCohortID;

            // Check which realm the cell is in
            if (cellEnvironment["Realm"][0] == 1.0)
            {
                // Get the indices of all terrestrial functional groups 
                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
                NumCohortsThisCell = totalCellTerrestrialCohorts;
            }
            else
            {
                // Get the indices of all marine functional groups
                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
                NumCohortsThisCell = totalCellMarineCohorts;
            }
            Debug.Assert(cellEnvironment["Realm"][0] > 0.0, "Missing realm for grid cell");

            if (NumCohortsThisCell > 0)
            {
                // Loop over all functional groups in the model
                for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
                {

                    // Create a new list to hold the cohorts in the grid cell
                    gridCellCohorts[FunctionalGroup] = new List<Cohort>();

                    // If it is a functional group that corresponds to the current realm, then seed cohorts
                    if (FunctionalGroupsToUse.Contains(FunctionalGroup))
                    {
                        // Loop over the initial number of cohorts
                        double NumberOfCohortsInThisFunctionalGroup = 1.0;
                        if (!ZeroAbundance)
                        {
                            NumberOfCohortsInThisFunctionalGroup = functionalGroups.GetBiologicalPropertyOneFunctionalGroup("initial number of gridcellcohorts", FunctionalGroup);
                        }


                        double MinFGMass = functionalGroups.GetBiologicalPropertyOneFunctionalGroup("minimum mass", FunctionalGroup);
                        double MaxFGMass = functionalGroups.GetBiologicalPropertyOneFunctionalGroup("maximum mass", FunctionalGroup);

                        int counter = 0;
                        int CounterIncrementRange = (int)Math.Ceiling((int)NumberOfCohortsInThisFunctionalGroup / 10.0);
                        double BasicJuveMinMass = MinFGMass;

                        for (int jj = 0; jj < NumberOfCohortsInThisFunctionalGroup; jj++)
                        {
                            // Check whether the model is set to randomly draw the body masses of new cohorts
                            if (DrawRandomly)
                            {
                            } else
                            {
                                // Use the same seed for the random number generator every time
                                RandomNumberGenerator.SetSeed((uint)(jj + 1), (uint)((jj + 1) * 3));
                            }
                            // Seeding works differently for small things
                            //if (MaxFGMass <= 1.01)
                            //{
                            //    CohortJuvenileMass = BasicJuveMinMass + BasicJuveMinMass * counter * 4.0/10.0;
                            //    counter++;

                            //    if (counter == CounterIncrementRange)
                            //    {
                            //        counter = 0;
                            //        BasicJuveMinMass = BasicJuveMinMass * 4.0;
                            //        if (BasicJuveMinMass > MaxFGMass)
                            //        {
                            //            BasicJuveMinMass = CohortJuvenileMass;
                            //            BasicJuveMinMass = BasicJuveMinMass + BasicJuveMinMass * 2.0 / 10.0;
                            //        }
                            //    }

                            //    // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                            //    // within the bounds of the minimum and maximum body masses for the functional group
                            //    CohortAdultMass = Math.Pow(10, (RandomNumberGenerator.GetUniform() * (Math.Log10(MaxFGMass) - Math.Log10(50 * CohortJuvenileMass)) + Math.Log10(50 * CohortJuvenileMass)));
                                
                            //}
                            //else
                            {
                                // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                                // within the bounds of the minimum and maximum body masses for the functional group
                                CohortAdultMass = Math.Pow(10, (RandomNumberGenerator.GetUniform() * (Math.Log10(MassMaxima[FunctionalGroup]) - Math.Log10(50 * MassMinima[FunctionalGroup])) + Math.Log10(50 * MassMinima[FunctionalGroup])));

                                // Calculate body mass ratios
                                CalculateBodyMassRatios(out ExpectedLnAdultMassRatio, out CohortAdultMassRatio, out CohortJuvenileMass, ref cellEnvironment, CohortAdultMass, MassMinima[FunctionalGroup], RandomNumberGenerator);
                            }

                            // Get optimal prey body size
                            double OptimalPreyBodySizeRatio = CalculateOptimalPreyBodySizeRatio(ref cellEnvironment, functionalGroups.GetTraitNames("Diet", FunctionalGroup) == "allspecial", RandomNumberGenerator);


                            // Calculate cohort abundance
                            double NewBiomass;
                            double NewAbund = CalculateCohortAbundance(ref cellEnvironment, ZeroAbundance, NumCohortsThisCell, CohortJuvenileMass, out NewBiomass, 1.0);

                            // Calculate trophic index of cohort
                            double TrophicIndex = CalculateTrophicIndex(NutritionSource[FunctionalGroup]);

                            // Write out properties of the selected cohort
                            //tempsw.WriteLine(FunctionalGroup.ToString() + '\t' + CohortAdultMass.ToString() + '\t' +
                            //     CohortJuvenileMass.ToString() + '\t' + NewBiomass.ToString() + '\t' +
                            //     NewAbund.ToString() + '\t' + OptimalPreyBodySizeRatio.ToString());

                            // Initialise the new cohort with the relevant properties
                            Cohort NewCohort = new Cohort((byte)FunctionalGroup, CohortJuvenileMass, CohortAdultMass, CohortJuvenileMass, NewAbund,
                            OptimalPreyBodySizeRatio, (ushort)0, ProportionTimeActive[FunctionalGroup], ref CohortIDIncrementer, TrophicIndex, tracking);

                            // Add the new cohort to the list of grid cell cohorts
                            gridCellCohorts[FunctionalGroup].Add(NewCohort);

                            // Incrememt the variable tracking the total number of cohorts in the model
                            globalDiagnostics["NumberOfCohortsInModel"]++;
                            
                        }

                    }
                }

            }
            else
            {
                // Loop over all functional groups in the model
                for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
                {
                    // Create a new list to hold the cohorts in the grid cell
                    gridCellCohorts[FunctionalGroup] = new List<Cohort>();
                }
            }

           // tempsw.Dispose();
        }

        /// <summary>
        /// Seed grid cell with zooplankton
        /// </summary>
        /// <param name="functionalGroups">The functional group definitions for cohorts in the grid cell</param>
        /// <param name="cellEnvironment">The environment in the grid cell</param>
        /// <param name="globalDiagnostics">A list of global diagnostic variables</param>
        /// <param name="nextCohortID">YThe unique ID to assign to the next cohort produced</param>
        /// <param name="tracking">boolean to indicate if cohorts are to be tracked in this model</param>
        /// <param name="totalCellTerrestrialCohorts">The total number of cohorts to be seeded in each terrestrial grid cell</param>
        /// <param name="totalCellMarineCohorts">The total number of cohorts to be seeded in each marine grid cell</param>
        /// <param name="DrawRandomly">Whether the model is set to use random draws</param>
        /// <param name="ZeroAbundance">Set this parameter to 'true' if you want to seed the cohorts with zero abundance</param>
        public void SeedGridCellCohortsPerTimeStep(GridCellCohortHandler gridCellCohorts, ref FunctionalGroupDefinitions functionalGroups, ref SortedList<string, double[]>
            cellEnvironment, SortedList<string, double> globalDiagnostics, ref Int64 partial, Boolean tracking, double totalCellTerrestrialCohorts,
            double totalCellMarineCohorts, Boolean DrawRandomly, Boolean ZeroAbundance, uint CellNum)
        {
            // Set the seed for the random number generator from the system time
            RandomNumberGenerator.SetSeedFromSystemTime();

            // Write out initial cohort information
            //StreamWriter tempsw = new StreamWriter("C://users//derekt//desktop//ts_cohorts.txt");
            //tempsw.WriteLine("functional group\tadult mass\tjuvenilemass\tbiomass\tabundance\toptimal pbr");

            // Define local variables
            double CohortJuvenileMass;
            double CohortAdultMassRatio;
            double CohortAdultMass;
            double ExpectedLnAdultMassRatio;
            int[] FunctionalGroupsToUse;
            double NumCohortsThisCell;

            // Get the minimum and maximum possible body masses for organisms in each functional group
            double[] MassMinima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("minimum mass");
            double[] MassMaxima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("maximum mass");
            string[] NutritionSource = functionalGroups.GetTraitValuesAllFunctionalGroups("nutrition source");

            double[] ProportionTimeActive = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("proportion suitable time active");

            // Check which realm the cell is in
            if (cellEnvironment["Realm"][0] == 1.0)
            {
                return;
            }
            else
            {
                // Get the indices of all marine functional groups
                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
                NumCohortsThisCell = 100 * FunctionalGroupsToUse.Count();
            }
            Debug.Assert(cellEnvironment["Realm"][0] > 0.0, "Missing realm for grid cell");

            // Loop over all functional groups in the model
            for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
            {
                // If it is a functional group that corresponds to the current realm, then seed cohorts
                if (FunctionalGroupsToUse.Contains(FunctionalGroup))
                {
                    // Loop over the initial number of cohorts
                    double NumberOfCohortsInThisFunctionalGroup = 1.0;
                    if (!ZeroAbundance)
                    {
                        NumberOfCohortsInThisFunctionalGroup = 0.2*functionalGroups.GetBiologicalPropertyOneFunctionalGroup("initial number of gridcellcohorts", FunctionalGroup);
                    }
                    for (int jj = 0; jj < NumberOfCohortsInThisFunctionalGroup; jj++)
                    {
                        // Check whether the model is set to randomly draw the body masses of new cohorts
                        if (DrawRandomly)
                        {
                        }
                        else
                        {
                            // Use the same seed for the random number generator every time
                            RandomNumberGenerator.SetSeed((uint)(jj + 1) + CellNum, (uint)((jj + 1) * 3) + CellNum);
                        }
                        // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                        // within the bounds of the minimum and maximum body masses for the functional group
                        CohortAdultMass = Math.Pow(10, (RandomNumberGenerator.GetUniform() * (Math.Log10(MassMaxima[FunctionalGroup]) - Math.Log10(50 * MassMinima[FunctionalGroup])) + Math.Log10(50 * MassMinima[FunctionalGroup])));

                        // Calculate body mass ratios
                        CalculateBodyMassRatios(out ExpectedLnAdultMassRatio, out CohortAdultMassRatio, out CohortJuvenileMass, ref cellEnvironment, CohortAdultMass, MassMinima[FunctionalGroup], RandomNumberGenerator);

                        // Only seeding new zooplankton
                        if (CohortJuvenileMass > 1.0)
                        {

                        }
                        else
                        {
                            // Get optimal prey body size
                            double OptimalPreyBodySizeRatio = CalculateOptimalPreyBodySizeRatio(ref cellEnvironment, functionalGroups.GetTraitNames("Diet", FunctionalGroup) == "allspecial", RandomNumberGenerator);

                            // Calculate cohort abundance
                            double NewBiomass;
                            double NewAbund = CalculateCohortAbundance(ref cellEnvironment, ZeroAbundance, NumCohortsThisCell, CohortJuvenileMass, out NewBiomass, 0.01);

                            // Calculate trophic index of cohort
                            double TrophicIndex = CalculateTrophicIndex(NutritionSource[FunctionalGroup]);

                            // Write out properties of the selected cohort
                            //tempsw.WriteLine(FunctionalGroup.ToString() + '\t' + CohortAdultMass.ToString() + '\t' +
                            //     CohortJuvenileMass.ToString() + '\t' + NewBiomass.ToString() + '\t' +
                            //     NewAbund.ToString() + '\t' + OptimalPreyBodySizeRatio.ToString());

                            // Initialise the new cohort with the relevant properties
                            Cohort NewCohort = new Cohort((byte)FunctionalGroup, CohortJuvenileMass, CohortAdultMass, CohortJuvenileMass, NewAbund,
                            OptimalPreyBodySizeRatio, (ushort)0, ProportionTimeActive[FunctionalGroup], ref partial, TrophicIndex, tracking);

                            //Console.WriteLine(NewCohort.CohortID[0]);

                            // Add the new cohort to the list of grid cell cohorts
                            gridCellCohorts[FunctionalGroup].Add(NewCohort);

                            // Increment the variable tracking the total number of cohorts in the model
                            globalDiagnostics["NumberOfCohortsInModel"]++;
                        }
                    }

                }
            }

            //tempsw.Dispose();
        }

        private double CalculateOptimalPreyBodySizeRatio(ref SortedList<string, double[]> cellEnvironment, 
            Boolean isBaleenWhale, NonStaticSimpleRNG RandomNumberGenerator)
        {
            // Terrestrial and marine organisms have different optimal prey/predator body mass ratios
            if (cellEnvironment["Realm"][0] == 1.0)
                // Optimal prey body size 10%
                return (Math.Max(0.01, RandomNumberGenerator.GetNormal(0.1, 0.02)));
            else
            {
                if (isBaleenWhale)
                {
                    // Note that for this group
                    // it is actually (despite the name) not an optimal prey body size ratio, but an actual body size.
                    // This is because it is invariant as the predator (filter-feeding baleen whale) grows.
                    // See also the predation classes.
                    return (Math.Max(1, RandomNumberGenerator.GetNormal(2, 2)));
                }
                else
                {
                    // Optimal prey body size or marine organisms is 10%
                    return (Math.Max(0.01, RandomNumberGenerator.GetNormal(0.1, 0.02)));
                }

            }

        }

        private void CalculateBodyMassRatios(out double ExpectedLnAdultMassRatio, out double CohortAdultMassRatio, out double CohortJuvenileMass, 
            ref SortedList<string, double[]> cellEnvironment, double CohortAdultMass, double massMinimaFG, NonStaticSimpleRNG RandomNumberGenerator)
        {
            //Variable for altering the juvenile to adult mass ratio for marine cells when handling certain functional groups eg baleen whales
            double Scaling = 0.0;

            // Draw from a log-normal distribution with mean 10.0 and standard deviation 5.0, then add one to obtain 
            // the ratio of adult to juvenile body mass, and then calculate juvenile mass based on this ratio and within the
            // bounds of the minimum and maximum body masses for this functional group
            if (cellEnvironment["Realm"][0] == 1.0)
            {
                do
                {
                    ExpectedLnAdultMassRatio = 2.24 + 0.13 * Math.Log(CohortAdultMass);
                    CohortAdultMassRatio = 1.0 + RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < massMinimaFG);
            }
            // In the marine realm, have a greater difference between the adult and juvenile body masses, on average
            else
            {
                uint Counter = 0;
                Scaling = 0.2;
                // Use the scaling to deal with baleen whales not having such a great difference
                do
                {

                    ExpectedLnAdultMassRatio = 2.5 + Scaling * Math.Log(CohortAdultMass);
                    CohortAdultMassRatio = 1.0 + 10 * RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                    Counter++;
                    if (Counter > 10)
                    {
                        Scaling -= 0.01;
                        Counter = 0;
                    }
                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < massMinimaFG);
            }


        }


        private double CalculateTrophicIndex(string nutritionSource)
        {
                            switch (nutritionSource)
                            {
                                case "herbivore":
                                    return(2.0);
                                    break;
                                case "omnivore":
                                    return(2.5);
                                    break;
                                case "carnivore":
                                    return(3);
                                    break;
                                default:
                                    Debug.Fail("Unexpected nutrition source trait value when assigning trophic index");
                                    return(0.0);
                                    break;
                            }

}

        private double CalculateCohortAbundance(ref SortedList<string, double[]> cellEnvironment, Boolean zeroAbundance, double numCohortsThisCell, 
            double cohortJuvenileMass, out double NewBiomass, double multiplier)
        {

            double TotalNewBiomass = 0.0;
            // 3000*(0.6^log(mass)) gives individual cohort biomass density in g ha-1
            // * 100 to give g km-2
            // * cell area to give g grid cell
            //*3300/NumCohortsThisCell scales total initial biomass in the cell to some approximately reasonable mass
            //double NewBiomass = (3300 / numCohortsThisCell) * 100 * 3000 *
            //Math.Pow(0.6, (Math.Log10(cohortJuvenileMass))) * (cellEnvironment["Cell Area"][0]);
            NewBiomass = (3300 / numCohortsThisCell) * 100 * 3000 *
Math.Pow(0.6, (Math.Log10(cohortJuvenileMass))) * (cellEnvironment["Cell Area"][0]) * multiplier;
            TotalNewBiomass += NewBiomass;
            double NewAbund = 0.0;
            if (!zeroAbundance)
            {
                NewAbund = NewBiomass / cohortJuvenileMass;
            }
            NewAbund = NewAbund / 100.0;
            return NewAbund;
        }


    /// <summary>
    /// Seed grid cell with stocks, as specified in the model input files
    /// </summary>
    /// <param name="functionalGroups">A reference to the stock functional group handler</param>
    /// <param name="cellEnvironment">The environment in the grid cell</param>
    /// <param name="globalDiagnostics">A list of global diagnostic variables for the model grid</param>
    public void SeedInitialGridCellStocks(ref FunctionalGroupDefinitions functionalGroups, ref SortedList<string, double[]>
            cellEnvironment, SortedList<string, double> globalDiagnostics, 
        GridCellStockHandler gridCellStocks)
        {
            // Set the seed for the random number generator from the system time
            RandomNumberGenerator.SetSeedFromSystemTime();

            Stock NewStock;

            // Define local variables
            int[] FunctionalGroupsToUse;

            // Get the individual body masses for organisms in each stock functional group
            double[] IndividualMass = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("individual mass");

            // Check which realm the cell is in
            if (cellEnvironment["Realm"][0] == 1.0 && cellEnvironment["Precipitation"][0] != cellEnvironment["Missing Value"][0] && cellEnvironment["Temperature"][0] != cellEnvironment["Missing Value"][0])
            {
                // Get the indices of all terrestrial functional groups 
                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
            }
            else if (cellEnvironment["Realm"][0] == 2.0 && cellEnvironment["NPP"][0] != cellEnvironment["Missing Value"][0])
            {
                // Get the indices of all marine functional groups
                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
            }
            else
            {
                // For cells without a realm designation, no functional groups will be used
                FunctionalGroupsToUse = new int[0];
            }

            // Loop over all functional groups in the model
            for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
            {
                // Create a new list to hold the stocks in the grid cell
                gridCellStocks[FunctionalGroup] = new List<Stock>();

                // If it is a functional group that corresponds to the current realm, then seed the stock
                if (FunctionalGroupsToUse.Contains(FunctionalGroup))
                {
                    if (cellEnvironment["Realm"][0] == 1.0)
                    {
                        // An instance of the terrestrial carbon model class
                        RevisedTerrestrialPlantModel PlantModel = new RevisedTerrestrialPlantModel();


                        // Calculate predicted leaf mass at equilibrium for this stock
                        double LeafMass = PlantModel.CalculateEquilibriumLeafMass(cellEnvironment, functionalGroups.GetTraitNames("leaf strategy", FunctionalGroup) == "deciduous");


                        double InitialFractionalArea = 0.0;
                        if (functionalGroups.GetTraitNames("impact state", FunctionalGroup).Equals("primary"))
                        {
                            InitialFractionalArea = cellEnvironment["Fprimary"][0];
                        }
                        else if (functionalGroups.GetTraitNames("impact state", FunctionalGroup).Equals("secondary"))
                        {
                            InitialFractionalArea = cellEnvironment["Fsecondary"][0];
                        }
                        else
                        {
                            //All HANPPlc comes from impacted lands
                            InitialFractionalArea = 1 - (cellEnvironment["Fprimary"][0] + cellEnvironment["Fsecondary"][0]);
                        }

                        // Initialise the new stock with the relevant properties
                        NewStock = new Stock((byte)FunctionalGroup, IndividualMass[FunctionalGroup], LeafMass,InitialFractionalArea, functionalGroups.GetTraitNames("Stock name", FunctionalGroup));

                        // Add the new stock to the list of grid cell stocks
                        gridCellStocks[FunctionalGroup].Add(NewStock);

                        // Increment the variable tracking the total number of stocks in the model
                        globalDiagnostics["NumberOfStocksInModel"]++;


                    }
                    else if (FunctionalGroupsToUse.Contains(FunctionalGroup))
                    {

                        double TotalPhytoBiomass = 1E8;

                        // Initialise the new stock with the relevant properties
                        NewStock = new Stock((byte)FunctionalGroup, IndividualMass[FunctionalGroup], TotalPhytoBiomass, 1.0, functionalGroups.GetTraitNames("Stock name", FunctionalGroup));
                        // Get current temperature in grid cell
                        double temperature = cellEnvironment["Temperature"][0];

                        // Calculate the coefficients for phytoplankton abundance-size curve
                        double no3 = cellEnvironment["no3"][0];

                        AutotrophProcessor a = new AutotrophProcessor();

                        double[] BinEdges = new double[100];
                        double[] BinCentres = new double[BinEdges.Length-1];
                        int b = 0;
                        for (double i = -9; i <= 0; i++)
                        {
                            //BinEdges[b] = Math.Log10(0.5 * Math.Pow(10.0, i));
                            for (double s = 1; s <= 10; s+=1.0 )
                                BinEdges[b++] = ((s/10) * Math.Pow(10.0, i));
                        }
                        b = 0;
                        for (int i = 0; i < BinEdges.Length-1; i++ )
                        {
                           
                                BinCentres[b++] = (BinEdges[i] + BinEdges[i+1])/2.0;
                        }

                        for (int i = 0; i < BinEdges.Length; i++)
                        {
                            BinEdges[i] = Math.Log10(BinEdges[i]);
                        }

                        for (int i = 0; i < BinCentres.Length; i++)
                        {
                            BinCentres[i] = Math.Log10(BinCentres[i]);
                        }

                        //b = 0;
                        //for (int i = -9; i <= -1; i++, b++)
                        //{
                        //    BinCentres[b] = i;
                        //}

                        double[] NSFCentres = a.GetPhytoDistributionEnvironment(temperature, no3, BinEdges, BinCentres);
                        
                        for (b = 0; b < NSFCentres.Length; b++) NSFCentres[b] = TotalPhytoBiomass * NSFCentres[b];

                        NewStock.SetSizeDistribution(BinCentres, BinEdges, NSFCentres);

                        // Add the new stock to the list of grid cell stocks
                        gridCellStocks[FunctionalGroup].Add(NewStock);

                        // Increment the variable tracking the total number of stocks in the model
                        globalDiagnostics["NumberOfStocksInModel"]++;

                    }
                    else
                    {
                    }

                }

            }
        }


        }
    }

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Madingley
{
    public class CohortTracker : Tracker
    {
        private StreamWriter CohortFlowsWriter;

        private TextWriter SyncedCohortFlowsWriter;

        private float[] Latitudes;
        private float[] Longitudes;

        private string OutputPath;
        private string FileName;
        private string FileSuffix;

        // Lists to hold information to record. Everything is done at the individual level,
        // i.e. PredationEaten is the amount eaten per individual
        private List<string[]> CohortIdentifierStrings;
        private List<uint[]> CohortIdentifiers;
        private List<float> PredationEaten;
        private List<float> HerbivoryEaten;
        private List<float[]> FixedCohortProperties;
        private List<float> MetabolicCosts;
        private List<float> GrowthRates;

        private int MaxNumberFunctionalGroups;

        /// <summary>
        /// Set up tracker for outputting properties of the eating process between functional groups
        /// </summary>
        /// 
        public CohortTracker(float[] latitudes, float[] longitudes, 
            FunctionalGroupDefinitions fGDefinitions, FunctionalGroupDefinitions stockFGDefinitions,
            string outputFileSuffix, string outputPath, string fileName)
        {
            Latitudes = latitudes;
            Longitudes = longitudes;

            FileName = fileName;
            OutputPath = outputPath;
            FileSuffix = outputFileSuffix;

            // Assign a number and name to each functional group (FG as defined by output, not by model)
            AssignFunctionalGroups(fGDefinitions, stockFGDefinitions);

            // Initialize array to hold mass flows among functional groups
            MaxNumberFunctionalGroups = Math.Max(NumberMarineFGsForTracking, NumberTerrestrialFGsForTracking);

            CohortIdentifierStrings = new List<string[]>();
            CohortIdentifiers = new List<uint[]>();
            FixedCohortProperties = new List<float[]>();
            PredationEaten = new List<float>();
            HerbivoryEaten = new List<float>();
            MetabolicCosts = new List<float>();
            GrowthRates = new List<float>();
        }

        /// <summary>
        /// Open the tracking file for writing to
        /// </summary>
        /// 
        public override void OpenTrackerFile()
        {
            CohortFlowsWriter = new StreamWriter(OutputPath + FileName + FileSuffix + ".txt");
            SyncedCohortFlowsWriter = TextWriter.Synchronized(CohortFlowsWriter);

            // Identifier properties: Lat index, lon index, cohort ID
            // Identified string properties: FG
            // Double flows: Predation, herbivory
            // Fixed double properties: Mass, abundance, metabolic cost, growth

            SyncedCohortFlowsWriter.WriteLine
                ("Latitude\tLongitude\ttime_step\tfunctional_group\tcohort_ID\tbody_mass_g\tabundance\tpredation_eaten_perind_g\therbivory_eaten_perind_g\tmetabolic_cost_perind_g\tgrowth_perind_g");
        }

        public void RecordGeneralCohortInformation(uint latIndex, uint lonIndex, Cohort cohort, 
            MadingleyModelInitialisation madingleyInitialisation, FunctionalGroupDefinitions stockFunctionalGroupDefinitions,
            string cohortOrStockName, double cohortOrStockBodyMass, Boolean Marine)
            // todo(erik): not all arguments are used.
        {
            CohortIdentifiers.Add(new uint[3] { latIndex, lonIndex, cohort.CohortID[0] });
            if(Marine)
                CohortIdentifierStrings.Add(new string[1] { MarineFGsForTracking.Keys.ToArray()[DetermineFunctionalGroup(cohortOrStockName, Marine)] });
            else
                CohortIdentifierStrings.Add(new string[1] { TerrestrialFGsForTracking.Keys.ToArray()[DetermineFunctionalGroup(cohortOrStockName, Marine)] });

            FixedCohortProperties.Add(new float[] { (float)cohort.IndividualBodyMass, (float)cohort.CohortAbundance });
        }

        public void RecordMetabolicCost(double metabolicCostPerIndividual)
        {
            MetabolicCosts.Add((float)metabolicCostPerIndividual);
        }

        public void RecordGrowth(double growthPerIndividual)
        {
            GrowthRates.Add((float)growthPerIndividual);
        }

        /// <summary>
        /// Record an eating (predation or herbivory) event
        /// </summary>
        /// 
        public void RecordPredationFlow(uint latIndex, uint lonIndex, Cohort cohort, double massEaten, Boolean herbivory, Boolean marineCell)
        {
            if(!herbivory)
            {
                PredationEaten.Add((float)massEaten);
            }
            else
            {
                HerbivoryEaten.Add((float)massEaten);
            }
        }

        /// <summary>
        /// Write flows of matter among functional groups to the output file at the end of the time step
        /// </summary>
        /// 
        public override void WriteToTrackerFile(uint currentTimeStep, ModelGrid madingleyModelGrid, uint numLats, uint numLons, 
            MadingleyModelInitialisation initialisation, bool MarineCell)
        {
            int NumCohortsToWrite = CohortIdentifiers.Count;

            for(int i = 0; i < NumCohortsToWrite; i++)
            {
                string newline = 
                    Convert.ToString(madingleyModelGrid.GetCellLatitude(CohortIdentifiers.ElementAt(i)[0])) + '\t' +
                    Convert.ToString(madingleyModelGrid.GetCellLongitude(CohortIdentifiers.ElementAt(i)[1])) + '\t' +
                    Convert.ToString(currentTimeStep) + '\t' +
                    CohortIdentifierStrings.ElementAt(i)[0] + '\t' +
                    CohortIdentifiers.ElementAt(i)[2] + '\t' + // 1 or 2?
                    FixedCohortProperties.ElementAt(i)[0] + '\t' +
                    FixedCohortProperties.ElementAt(i)[1] + '\t' +
                    PredationEaten.ElementAt(i) + '\t' +
                    HerbivoryEaten.ElementAt(i) + '\t' +
                    MetabolicCosts.ElementAt(i) + '\t' +
                    GrowthRates.ElementAt(i);

                SyncedCohortFlowsWriter.WriteLine(newline);
            }

            // Reset lists
            CohortIdentifiers = new List<uint[]>();
            CohortIdentifierStrings = new List<string[]>();
            FixedCohortProperties = new List<float[]>();
            PredationEaten = new List<float>();
            HerbivoryEaten = new List<float>();
            MetabolicCosts = new List<float>();
            GrowthRates = new List<float>();
        }

        /// <summary>
        /// Close the file that has been written to
        /// </summary>
        /// 
        public override void CloseTrackerFile()
        {
            CohortFlowsWriter.Dispose();
            SyncedCohortFlowsWriter.Dispose();
        }
    }
}
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Madingley
{
    public class HighResFGTracker : Tracker
    {
        private StreamWriter FGFlowsWriter;

        private TextWriter SyncedFGFlowsWriter;
        
        private string OutputPath;
        private string FileName;
        private string FileSuffix;

        //predator mass, prey mass, initial abundance, abundance eaten
        private List<Tuple<double, double, double, double>> FeedingInteractionsMatrix;
        //From FG number, to FG number, 
        private List<Tuple<int, int>> FunctionalGroupsMatrix;

        private int MaxNumberFunctionalGroups;

        /// <summary>
        /// Set up tracker for outputting properties of the eating process between functional groups
        /// </summary>
        /// <param name="latitudes">Array of latitudes</param>
        /// <param name="longitudes">Array of longitudes</param>
        /// <param name="fGDefinitions">Cohort functional group definitions</param>
        /// <param name="stockDefinitions">Stock function group definitions</param>
        /// <param name="outputFileSuffix">Output file suffix</param>
        /// <param name="outputPath">Output directory path</param>
        /// <param name="fileName">Output file name</param>
        public HighResFGTracker(FunctionalGroupDefinitions fGDefinitions, FunctionalGroupDefinitions stockDefinitions, string outputFileSuffix,
            string outputPath, string fileName)
        {
            FileName = fileName;
            OutputPath = outputPath;
            FileSuffix = outputFileSuffix;

            // Assign a number and name to each functional group (FG as defined by output, not by model)
            AssignFunctionalGroups(fGDefinitions, stockDefinitions);

            // Initialise array to hold mass flows among functional groups
            MaxNumberFunctionalGroups = Math.Max(NumberMarineFGsForTracking, NumberTerrestrialFGsForTracking);

            // Predator body mass, prey body mass, % prey eaten
            FeedingInteractionsMatrix = new List<Tuple<double, double, double, double>>();
            FunctionalGroupsMatrix = new List<Tuple<int, int>>();
        }

        /// <summary>
        /// Open the tracking file for writing
        /// </summary>
        /// 
        public override void OpenTrackerFile()
        {
            FGFlowsWriter = new StreamWriter(OutputPath + FileName + FileSuffix + ".txt");
            SyncedFGFlowsWriter = TextWriter.Synchronized(FGFlowsWriter);
            SyncedFGFlowsWriter.WriteLine("Time_step\tFromIndex\tToIndex\tPredatorMass\tPreyMass\tInitialAbundance\tAbundanceEaten");
        }

        /// <summary>
        /// Record a flow of biomass between two functional groups (as specified by the the tracker)
        /// </summary>
        /// <param name="latIndex"></param>
        /// <param name="lonIndex"></param>
        /// <param name="predatorCohortOrStockName">Predator cohort or stock name</param>
        /// <param name="preyCohortOrStockName">Prey cohort or stock name</param>
        /// <param name="massEaten">Biomass eaten which flows from one functional group to another</param>
        /// <param name="marineCell">Whether this is a marine cell</param>
        public void RecordFGFlow(string predatorCohortOrStockName, string preyCohortOrStockName, 
            double predatorBodyMass, double preyBodyMass, double InitialAbundance, double AbundanceEaten, Boolean marineCell)
        {
            int fromIndex = 0;
            int toIndex = 0;

            // Get the functional group that the mass is flowing to
            toIndex = DetermineFunctionalGroup(predatorCohortOrStockName, marineCell);

            // Get the functional group that the mass is flowing from
            fromIndex = DetermineFunctionalGroup(preyCohortOrStockName, marineCell);

            FunctionalGroupsMatrix.Add(new Tuple<int, int>(fromIndex, toIndex));

            //predator mass, prey mass, mass eaten, percent cohort eaten
            FeedingInteractionsMatrix.Add(new Tuple<double, double, double, double>(predatorBodyMass, preyBodyMass, InitialAbundance, AbundanceEaten));
    }

        /// <summary>
        /// Write flows of matter among functional groups to the output file at the end of the time step
        /// </summary>
        /// 
        public override void WriteToTrackerFile(uint currentTimeStep, ModelGrid madingleyModelGrid, uint numLats, uint numLons, MadingleyModelInitialisation initialisation, Boolean MarineCell)
        {
            //predator mass, prey mass, mass eaten, percent cohort eaten
            //private List<Tuple<double, double, double, double>> FeedingInteractionsMatrix;
            
            //From FG number, to FG number, 
            //private List<Tuple<int, int>> FunctionalGroupsMatrix;

            for (int i = 0; i < FunctionalGroupsMatrix.Count(); i++)
            {
                if (MarineCell)
                {
                    SyncedFGFlowsWriter.WriteLine(Convert.ToString(currentTimeStep) +
                        '\t' + MarineFGsForTracking.Keys.ToArray()[FunctionalGroupsMatrix[i].Item1] + '\t' + 
                        MarineFGsForTracking.Keys.ToArray()[FunctionalGroupsMatrix[i].Item2] + '\t' +
                        Convert.ToString(FeedingInteractionsMatrix[i].Item1) + '\t' + Convert.ToString(FeedingInteractionsMatrix[i].Item2) + '\t' +
                        Convert.ToString(FeedingInteractionsMatrix[i].Item3) + '\t' + Convert.ToString(FeedingInteractionsMatrix[i].Item4));
                } else
                {
                    SyncedFGFlowsWriter.WriteLine(Convert.ToString(currentTimeStep) +
                          '\t' + TerrestrialFGsForTracking.Keys.ToArray()[FunctionalGroupsMatrix[i].Item1] + '\t' + TerrestrialFGsForTracking.Keys.ToArray()[FunctionalGroupsMatrix[i].Item2] + '\t' +
                          Convert.ToString(FeedingInteractionsMatrix[i].Item1) + '\t' + Convert.ToString(FeedingInteractionsMatrix[i].Item2) + '\t' +
                          Convert.ToString(FeedingInteractionsMatrix[i].Item3) + '\t' + Convert.ToString(FeedingInteractionsMatrix[i].Item4));
                }
            }

            // Reset arrays to hold mass flows among trophic levels
            FeedingInteractionsMatrix.Clear();
            FunctionalGroupsMatrix.Clear();
        }

        /// <summary>
        /// Close the file that has been written to
        /// </summary>
        /// 
        public override void CloseTrackerFile()
        {
            FGFlowsWriter.Dispose();
        }
    }
}

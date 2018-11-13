using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;



// using RDotNet;
// using RserveCli;

namespace Madingley
{
    /// <summary>
    /// Calculates ecosystem-level metrics
    /// </summary>
    public class EcosytemMetrics
    {
        /// <summary>
        /// The trophic level values associated with each of the trophic level bins
        /// </summary>
        private float[] _TrophicIndexBinValues;
        /// <summary>
        /// Get and set the trophic level values associated with each of the trophic level bins
        /// </summary>
        public float[] TrophicIndexBinValues
        {
            get { return _TrophicIndexBinValues; }
            set { _TrophicIndexBinValues = value; }
        }

        /// <summary>
        /// The number of trophic level bins to use in calculating ecosystem-level metrics
        /// </summary>
        private int _NumberTrophicBins;

        /// <summary>
        /// Get and set the number of trophic level bins to use in calculating ecosystem-level metrics
        /// </summary>
        public int NumberTrophicBins
        {
            get { return _NumberTrophicBins; }
            set { _NumberTrophicBins = value; }
        }

        /// <summary>
        /// Instance of the connection to R
        /// </summary>
        private Process _RServeProcess;


        //Define upp and lower limits for trophic index
        private double MaxTI = 40.0;
        private double MinTI = 1.0;
       
        /// <summary>
        /// Constructor for the ecosystem metrics class: sets up trophic level bins
        /// </summary>
        public EcosytemMetrics()
        {

            float TrophicIndexBinWidth = 0.2f;
            float LowestTophicIndex = 1.0f;
            float HighestTrophicIndex = 5.0f;
            NumberTrophicBins = (int)((HighestTrophicIndex - LowestTophicIndex) / TrophicIndexBinWidth);
            _TrophicIndexBinValues = new float[NumberTrophicBins];

            for (int i = 0; i < TrophicIndexBinValues.Length; i++)
            {
                _TrophicIndexBinValues[i] = LowestTophicIndex + (TrophicIndexBinWidth * i);
            }

            /*Console.WriteLine("Opening R connection");

            
            int port = int.Parse((6).ToString() + (2*scenarioIndex).ToString() + (simulationNumber).ToString() + cellIndex.ToString());

            ProcessStartInfo StartInfo = new ProcessStartInfo();
            StartInfo.CreateNoWindow = true;
            StartInfo.UseShellExecute = false;
            StartInfo.FileName = "C:/Users/mikeha/Documents/R/win-library/2.14/Rserve/libs/x64/RServe.exe";
            StartInfo.WindowStyle = ProcessWindowStyle.Hidden;
            StartInfo.Arguments = "--RS-port " + port.ToString();
            
            //Start an instance of RServe to talk with R
            _RServeProcess = Process.Start(StartInfo);

            R = new RConnection(new System.Net.IPAddress(new byte[] { 127, 0, 0, 1 }),port);
            //RConn.VoidEval("install.packages(/"FD/")");
            R.VoidEval("library(FD)");
            */
        }

        /// <summary>
        /// Closes the connection to R (currently disabled)
        /// </summary>
        public void CloseRserve()
        {
            //_RServeProcess.Kill();
        }

        /// <summary>
        /// Calculates the mean trophic level of all individual organisms in a grid cell
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of cell indices in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell in the list of cells to run</param>
        /// <returns>The mean trophic level of individuals in the grid cell</returns>
        public double CalculateMeanTrophicLevelCell(ModelGrid ecosystemModelGrid,List<uint[]> cellIndices, int cellIndex)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            double BiomassWeightedTI = 0.0;
            double TotalBiomass = 0.0;
            double CohortBiomass = 0.0;

            foreach (var CohortList in CellCohorts)
            {
                foreach (Cohort c in CohortList)
                {
                    CohortBiomass = (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    BiomassWeightedTI += CohortBiomass * c.TrophicIndex;
                    TotalBiomass += CohortBiomass;
                }
            }

            return BiomassWeightedTI/TotalBiomass;
        }


        /// <summary>
        /// Return the distribution of biomasses among trophic level bins
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of cell indices to be run in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell in the list of cells to be run</param>
        /// <returns>The distribution of biomasses among trophic level bins</returns>
        public double[] CalculateTrophicDistribution(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            double[] TrophicIndexBinMasses = new double[NumberTrophicBins];
            int BinIndex;


            foreach (var CohortList in CellCohorts)
            {
                foreach (Cohort c in CohortList)
                {
                    BinIndex = _TrophicIndexBinValues.ToList().IndexOf(_TrophicIndexBinValues.Last(x => x < c.TrophicIndex));
                    TrophicIndexBinMasses[BinIndex] += (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                }
            }

            return TrophicIndexBinMasses;
        }


        public double[] CalculateFunctionalRichness(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortDefinitions, 
            List<uint[]> cellIndices, int cellIndex, string trait)
        {

            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            double MinCurrentTraitValue = double.MaxValue;
            double MaxCurrentTraitValue = double.MinValue;
            double MinModelTraitValue = 0.0;
            double MaxModelTraitValue = 0.0;

            switch (trait.ToLower())
            {
                case "biomass":

                    foreach (var CohortList in CellCohorts)
                    {
                        foreach (var cohort in CohortList)
                        {

                            if (cohort.IndividualBodyMass < MinCurrentTraitValue) MinCurrentTraitValue = cohort.IndividualBodyMass;

                            if (cohort.IndividualBodyMass > MaxCurrentTraitValue) MaxCurrentTraitValue = cohort.IndividualBodyMass;

                        }
                    }


                    //Define upper and lower limits for body mass
                    MinModelTraitValue = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("minimum mass").Min();
                    MaxModelTraitValue = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("maximum mass").Max();
                    break;
                case "trophic index":
                    foreach (var CohortList in CellCohorts)
                    {
                        foreach (var cohort in CohortList)
                        {

                            if (cohort.TrophicIndex < MinCurrentTraitValue) MinCurrentTraitValue = cohort.TrophicIndex;

                            if (cohort.TrophicIndex > MaxCurrentTraitValue) MaxCurrentTraitValue = cohort.TrophicIndex;

                        }
                    }


                    //Define upper and lower limits for body mass
                    MinModelTraitValue = MinTI;
                    MaxModelTraitValue = MaxTI;

                    break;
                default:
                    Debug.Fail("Trait not recognised in calculation of ecosystem metrics: " + trait);
                    break;
            }

            Debug.Assert((MaxModelTraitValue - MinModelTraitValue) > 0.0, "Division by zero or negative model trait values in calculation of functional richness");

            double[] NewArray = {(MaxCurrentTraitValue-MinCurrentTraitValue)/(MaxModelTraitValue-MinModelTraitValue),MinCurrentTraitValue,MaxCurrentTraitValue};

            return NewArray;
        }

        /// <summary>
        /// Calculate trophic evenness using the Rao Index
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of indices of cells to be run in the current simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to be run</param>
        /// <returns>Trophic evenness</returns>
        public double CalculateFunctionalEvennessRao(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortDefinitions,
            List<uint[]> cellIndices, int cellIndex, string trait)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);

            double[] EvennessValues = new double[2];

            double[,] Distances = new double[CellCohorts.GetNumberOfCohorts(), CellCohorts.GetNumberOfCohorts()];

            double[] FunctionalTrait = new double[CellCohorts.GetNumberOfCohorts()];
            double MaxModelTraitValue=0;
            double MinModelTraitValue=0;

            // Construct a vector of cohort biomass (in case we want to weight by them)
            double[] CohortTotalBiomasses = new double[CellCohorts.GetNumberOfCohorts()];


            int CohortNumberCounter = 0;
            switch (trait.ToLower())
            {
                case "biomass":
                    for (int fg = 0; fg < CellCohorts.Count; fg++)
                    {
                        foreach (Cohort c in CellCohorts[fg])
                        {

                            FunctionalTrait[CohortNumberCounter] = c.IndividualBodyMass;
                            CohortTotalBiomasses[CohortNumberCounter] = (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;

                            CohortNumberCounter++;
                        }
                    }

                    //Define upper and lower limits for body mass
                    MinModelTraitValue = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("minimum mass").Min();
                    MaxModelTraitValue = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("maximum mass").Max();
                    break;
                case "trophic index":
                    for (int fg = 0; fg < CellCohorts.Count; fg++)
                    {
                        foreach (Cohort c in CellCohorts[fg])
                        {

                            FunctionalTrait[CohortNumberCounter] = c.IndividualBodyMass;
                            CohortTotalBiomasses[CohortNumberCounter] = (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;

                            CohortNumberCounter++;
                        }
                    }
                    MinModelTraitValue = MinTI;
                    MaxModelTraitValue = MaxTI;
                    break;
            }


            Distances = CalculateDistanceMatrix(FunctionalTrait, MaxModelTraitValue, MinModelTraitValue);

            return RaoEntropy(Distances, CohortTotalBiomasses);

        }

        /// <summary>
        /// Calculate the ecosystem metabolism for an individual grid cell. This is the amount of respiration per gram of biomass in the grid cell; i.e. the units are g C / g / time-step. A value of 
        /// one would indicate that metabolic cost is exactly equivalent to total respiring biomass.
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of indices of cells to be run in the current simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to be run</param>
        /// <returns></returns>
        public double CalculateBiomassWeightedSystemMetabolism(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex)
        {
            double MetabolicLoss = 0.0;
            double TotalBiomass = 0.0;

            // Retrieve the metabolic loss this timestep
            bool varExists;
            MetabolicLoss = ecosystemModelGrid.GetEnviroLayer("Respiratory CO2 Pool Per Timestep", 0, cellIndices[cellIndex][0], cellIndices[cellIndex][1], out varExists);
            TotalBiomass = ecosystemModelGrid.GetEnviroLayer("Respiring Biomass Pool Per Timestep", 0, cellIndices[cellIndex][0], cellIndices[cellIndex][1], out varExists);


            return MetabolicLoss / TotalBiomass;
        }
        /// <summary>
        /// Calculates the arithmetic community weighted mean body mass
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of indices of cells to be run in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to be run</param>
        /// <returns>arithmetic community weighted mean body mass</returns>
        public double CalculateArithmeticCommunityMeanBodyMass(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex)
        {

            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            double CumulativeAbundance = 0.0;
            double CumulativeBiomass = 0.0;

            //Retrieve the biomass
            foreach (var CohortList in CellCohorts)
            {
                foreach (Cohort c in CohortList)
                {
                    CumulativeBiomass += (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    CumulativeAbundance += c.CohortAbundance;
                }
            }

            double CWAMBM = (CumulativeBiomass / CumulativeAbundance);

            return (CWAMBM);

        }

        /// <summary>
        /// Calculates the geometric community weighted mean body mass
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of indices of cells to be run in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to be run</param>
        /// <returns>geometric community weighted mean body mass</returns>
        public double CalculateGeometricCommunityMeanBodyMass(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex)
        {

            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            double CumulativeAbundance = 0.0;
            double CumulativeLogBiomass = 0.0;
            
            //Retrieve the biomass
            foreach (var CohortList in CellCohorts)
            {
                foreach (Cohort c in CohortList)
                {
                    CumulativeLogBiomass += Math.Log(c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    CumulativeAbundance += c.CohortAbundance;
                }
            }

            double CWGMBM = Math.Exp(CumulativeLogBiomass / CumulativeAbundance);

            return (CWGMBM);

        }

        public double CalculateGeometricCommunityMeanBodyMassbyFG(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex, int functionalGroup)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);

            double CumulativeAbundance = 0.0;
            double CumulativeLogBiomass = 0.0;

                foreach (Cohort c in CellCohorts[functionalGroup])
                {
                    CumulativeLogBiomass += Math.Log(c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    CumulativeAbundance += c.CohortAbundance;
                }

            return Math.Exp(CumulativeLogBiomass / CumulativeAbundance);
        }

        /// <summary>
        /// Calculates trophic evenness using the FRO Index of Mouillot et al.
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cellIndices">The list of indices of cells to be run in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to be run</param>
        /// <returns>Trophic evenness</returns>
        /// <remarks>From Mouillot et al (2005) Functional regularity: a neglected aspect of functional diversity, Oecologia</remarks>
        public double CalculateTrophicEvennessFRO(ModelGrid ecosystemModelGrid, List<uint[]> cellIndices, int cellIndex)
        {

            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);
            List<double[]> TrophicIndexBiomassDistribution = new List<double[]>();
            double[] TIBiomass;
            double[] EW;

            foreach (var CohortList in CellCohorts)
            {
                foreach (Cohort c in CohortList)
                {
                    TIBiomass = new double[2];
                    TIBiomass[0] = c.TrophicIndex;
                    TIBiomass[1] = (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    TrophicIndexBiomassDistribution.Add(TIBiomass);
                }
            }

            TrophicIndexBiomassDistribution = TrophicIndexBiomassDistribution.OrderBy(x => x[0]).ToList();


            //Use the Mouillot Evenness index - Functional Regularity Index or FRO
            //From Mouillot et al (2005) Functional regularity: a neglected aspect of functional diversity, Oecologia

            EW = new double[TrophicIndexBiomassDistribution.Count];
            double TotalEW = 0.0 ;

            for (int ii = 0; ii < TrophicIndexBiomassDistribution.Count-1; ii++)
            {
                EW[ii] = (TrophicIndexBiomassDistribution[ii + 1][0] - TrophicIndexBiomassDistribution[ii][0]) / (TrophicIndexBiomassDistribution[ii + 1][1] + TrophicIndexBiomassDistribution[ii][1]);
                TotalEW += EW[ii];
            }

            double FRO = 0.0;

            for (int ii = 0; ii < TrophicIndexBiomassDistribution.Count - 1; ii++)
            {
                FRO += Math.Min(EW[ii]/TotalEW,1.0/(TrophicIndexBiomassDistribution.Count-1));
            }

            return FRO;
        }

        /// <summary>
        /// Calculates functional diversity of cohorts in a grid cell as functional richness and functional diveregence (using the Rao Index)
        /// </summary>
        /// <param name="ecosystemModelGrid">The model grid</param>
        /// <param name="cohortDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="cellIndices">The list of cell indices in the current model simulation</param>
        /// <param name="cellIndex">The index of the current cell within the list of cells to run</param>
        /// <returns>A pair of values representing the functional richness and functional divergence (functional richness currently disabled!)</returns>
        public double[] CalculateFunctionalDiversity(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortDefinitions, 
            List<uint[]> cellIndices, int cellIndex)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);

            //Variable to hold the functional richness value for the current cohorts
            double FunctionalRichness;
            //Variable to hold the functional divergence value for the current cohorts
            double RaoFunctionalDivergence = 0.0;
            double[,] Distances= new double[CellCohorts.GetNumberOfCohorts(), CellCohorts.GetNumberOfCohorts()];

            List<string> AllTraitNames = cohortDefinitions.GetAllTraitNames().ToList();

            AllTraitNames.Remove("realm");
            AllTraitNames.Remove("heterotroph/autotroph");
            AllTraitNames.Remove("diet");
            string[] TraitNames = AllTraitNames.ToArray();


            //Define upper and lower limits for body mass
            double MinMass = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("minimum mass").Min();
            double MaxMass = cohortDefinitions.GetBiologicalPropertyAllFunctionalGroups("maximum mass").Max();
            //Define upp and lower limits for trophic index
            double MaxTI = 40.0;
            double MinTI = 1.0;

            // Construct an array of functional trait values for each cohort
            // Rows are specific cohorts
            // Columns are the functional traits (these include different types:
            //      quantative: current mass, trophic index
            //      nominal: diet, reproductive strategy, mobility, metabolism
            Tuple<double[], string[]>[] CohortFunctionalTraits = new Tuple<double[], string[]>[CellCohorts.GetNumberOfCohorts()];
            double[] IndividualBodyMasses = new double[CellCohorts.GetNumberOfCohorts()];
            double[] TrophicIndex = new double[CellCohorts.GetNumberOfCohorts()];
            string[][] CohortNominalTraitValues= new string[TraitNames.Length][];

            for (int i = 0; i < TraitNames.Length; i++)
			{
			    CohortNominalTraitValues[i] = new string[CellCohorts.GetNumberOfCohorts()];
			}

            // Construct a vector of cohort biomass (in case we want to weight by them)
            double[] CohortTotalBiomasses = new double[CellCohorts.GetNumberOfCohorts()];

            
            string[] TraitValues = new string[TraitNames.Length];
            double[] QuantitativeTraitValues= new double[2];
            int CohortNumberCounter = 0;
            for (int fg = 0; fg < CellCohorts.Count; fg++)
			{
                foreach (Cohort c in CellCohorts[fg])
                {
                    TraitValues = cohortDefinitions.GetTraitValues(TraitNames, fg);
                    for (int ii = 0; ii < TraitValues.Length; ii++)
                    {
			            CohortNominalTraitValues[ii][CohortNumberCounter] = TraitValues[ii];
                    }


                    IndividualBodyMasses[CohortNumberCounter] = c.IndividualBodyMass;
                    TrophicIndex[CohortNumberCounter] = c.TrophicIndex;
 
                    QuantitativeTraitValues[0] = c.IndividualBodyMass;
                    QuantitativeTraitValues[1] = c.TrophicIndex;

                    CohortFunctionalTraits[CohortNumberCounter] = new Tuple<double[], string[]>(QuantitativeTraitValues, TraitValues);
                    
                    CohortTotalBiomasses[CohortNumberCounter] = (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                    
                    CohortNumberCounter++;
                }
            }
            
            List<double[,]> DistanceList = new List<double[,]>();

            DistanceList.Add(CalculateDistanceMatrix(IndividualBodyMasses, MaxMass, MinMass));
            DistanceList.Add(CalculateDistanceMatrix(TrophicIndex, MaxTI, MinTI));
            foreach (string[] t in CohortNominalTraitValues)
            {
                DistanceList.Add(CalculateDistanceMatrix(t));
            }

            Distances = CalculateAggregateDistance(DistanceList);

            RaoFunctionalDivergence = RaoEntropy(Distances, CohortTotalBiomasses);

            return new double[] {0.0,RaoFunctionalDivergence};
            

        }

        /// <summary>
        /// Functional richness in the 3D space: log10 Adult body mass, trophic level, metabolic pathway
        /// </summary>
        /// <param name="ecosystemModelGrid"></param>
        /// <param name="cohortDefinitions"></param>
        /// <param name="cellIndices"></param>
        /// <param name="cellIndex"></param>
        /// <returns></returns>
        public double FunctionalRichness(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortDefinitions,
            List<uint[]> cellIndices, int cellIndex)
        {
            //Get the cohorts for the specified cell
            GridCellCohortHandler CellCohorts = ecosystemModelGrid.GetGridCellCohorts(cellIndices[cellIndex][0], cellIndices[cellIndex][1]);

            //Find the min and max adult body mass for
            // Carnivores, omnivores and herbivores, & endotherms or ectotherms
            
            //Find FG inds of endotherms
            string[] Traits = new string[] { "nutrition source", "endo/ectotherm" };
            string[] TraitVals;

            TraitVals = new string[]{"carnivore","endotherm"};
            double[] carn_end_min_max = this.FindMinMaxAdultMass(CellCohorts,cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals,true));

            TraitVals = new string[] { "carnivore", "ectotherm" };
            double[] carn_ect_min_max= this.FindMinMaxAdultMass(CellCohorts, cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals, true));

            TraitVals = new string[] { "omnivore", "endotherm" };
            double[] omn_end_min_max = this.FindMinMaxAdultMass(CellCohorts, cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals, true));

            TraitVals = new string[] { "omnivore", "ectotherm" };
            double[] omn_ect_min_max = this.FindMinMaxAdultMass(CellCohorts, cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals, true));

            TraitVals = new string[] { "herbivore", "endotherm" };
            double[] herb_end_min_max = this.FindMinMaxAdultMass(CellCohorts, cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals, true));

            TraitVals = new string[] { "herbivore", "ectotherm" };
            double[] herb_ect_min_max = this.FindMinMaxAdultMass(CellCohorts, cohortDefinitions.GetFunctionalGroupIndex(Traits, TraitVals, true));

            //Construct points for these twelve points
            // <Trophic Level, Metabolism, log10 Adult mass>
            Point3D[] points = new Point3D[12];

            // Ects = 1, Ends = 2
            // Herb = 1, Omn = 2, Carn = 3
            int i = 0;

            points[i++] = new Point3D(1.0, 1.0, Math.Log10(herb_ect_min_max[0]));//0
            points[i++] = new Point3D(1.0, 1.0, Math.Log10(herb_ect_min_max[1]));//1
            points[i++] = new Point3D(1.0, 2.0, Math.Log10(herb_end_min_max[0]));//3
            points[i++] = new Point3D(1.0, 2.0, Math.Log10(herb_end_min_max[1]));//4

            points[i++] = new Point3D(2.0, 1.0, Math.Log10(omn_ect_min_max[0]));//5
            points[i++] = new Point3D(2.0, 1.0, Math.Log10(omn_ect_min_max[1]));//6
            points[i++] = new Point3D(2.0, 2.0, Math.Log10(omn_end_min_max[0]));//7
            points[i++] = new Point3D(2.0, 2.0, Math.Log10(omn_end_min_max[1]));//8

            points[i++] = new Point3D(3.0, 1.0, Math.Log10(carn_ect_min_max[0]));//9
            points[i++] = new Point3D(3.0, 1.0, Math.Log10(carn_ect_min_max[1]));//10
            points[i++] = new Point3D(3.0, 2.0, Math.Log10(carn_end_min_max[0]));//11
            points[i++] = new Point3D(3.0, 2.0, Math.Log10(carn_end_min_max[1]));//12

            //Should be 20 simplices
            Triangle[] Mesh = new Triangle[20];
            i = 0;
            Mesh[i++] = new Triangle(points[6], points[2], points[0]);//6,  2,  0
            Mesh[i++] = new Triangle(points[4], points[6], points[0]);// 4,  6,  0
            Mesh[i++] = new Triangle(points[10], points[4], points[8]);//10,  4,  8
            Mesh[i++] = new Triangle(points[10], points[4], points[6]);//[10,  4,  6
            Mesh[i++] = new Triangle(points[5], points[1], points[3]);//5,  1,  3
            Mesh[i++] = new Triangle(points[5], points[9], points[3]);// 5,  9,  3
            Mesh[i++] = new Triangle(points[7], points[11], points[9]);// 7, 11,  9
            Mesh[i++] = new Triangle(points[7], points[9], points[3]);//7,  9,  3]
            Mesh[i++] = new Triangle(points[1], points[3], points[2]);//1,  3,  2
            Mesh[i++] = new Triangle(points[1], points[2], points[0]);// 1,  2,  0
            Mesh[i++] = new Triangle(points[5], points[1], points[0]);// 5,  1,  0
            Mesh[i++] = new Triangle(points[5], points[4], points[8]);// 5,  4,  8
            Mesh[i++] = new Triangle(points[5], points[9], points[8]);//  5,  9,  8
            Mesh[i++] = new Triangle(points[5], points[4], points[0]);//  5,  4,  0
            Mesh[i++] = new Triangle(points[11], points[9], points[8]);// [11,  9,  8
            Mesh[i++] = new Triangle(points[11], points[10], points[8]);// [11, 10,  8
            Mesh[i++] = new Triangle(points[7], points[3], points[2]);//  7,  3,  2
            Mesh[i++] = new Triangle(points[7], points[6], points[2]);//  7,  6,  2
            Mesh[i++] = new Triangle(points[7], points[10], points[6]);// 7, 10,  6
            Mesh[i++] = new Triangle(points[7], points[11], points[10]);//  7, 11, 10



            return VolumeOfMeshVertex(Mesh,points[0]);
        }


        public double VolumeOfMeshVertex(Triangle[] triangles,Point3D p0)
        {

            foreach (Triangle t in triangles)
            {
                t.ToTetrahedron(p0);
            }

            var vols = from t in triangles
                       select (Math.Abs(DotProduct(Minus(t.P1,t.P4),CrossProduct(Minus(t.P2,t.P4),Minus(t.P3,t.P4))))/6.0);

            return Math.Abs(vols.Sum());
        }

        
        public double DotProduct(Point3D a, Point3D b)
        {
            return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
        }

        public Point3D CrossProduct(Point3D a, Point3D b)
        {
            return new Point3D(x: (a.Y*b.Z) - (a.Z*b.Y),y: (a.Z*b.X)-(a.X*b.Z), z: (a.X*b.Y)-(a.Y*b.X));
        }

        public Point3D Minus(Point3D a, Point3D b)
        {
            return new Point3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        private double[] FindMinMaxAdultMass(GridCellCohortHandler CellCohorts, int[] fgs)
        {
            double MinAM = double.MaxValue;
            double MaxAM = double.MinValue;

            foreach (var fg in fgs)
            {
                foreach (Cohort c in CellCohorts[fg])
                {
                    if (c.AdultMass < MinAM) MinAM = c.AdultMass;
                    if (c.AdultMass > MaxAM) MaxAM = c.AdultMass;
                }
            }

            return new double[] { MinAM, MaxAM };
        }


        private double[,] CalculateDistanceMatrix(double[] continuousTrait, double traitMaxVal, double traitMinVal)
        {
            double[,] D = new double[continuousTrait.Length, continuousTrait.Length];
            double Range = traitMaxVal - traitMinVal;

            for (int ii = 0; ii < continuousTrait.Length; ii++)
            {
                for (int jj = ii; jj < continuousTrait.Length; jj++)
                {
                    D[ii, jj] = Math.Abs(continuousTrait[ii] - continuousTrait[jj]) / Range;
                    D[jj, ii] = D[ii, jj];
                }
            }

            return D;
        }

        private double[,] CalculateDistanceMatrix(string[] nominalTrait)
        {
            int NumberOfCohorts = nominalTrait.Length;
            double[,] D = new double[NumberOfCohorts, NumberOfCohorts];

            for (int ii = 0; ii < NumberOfCohorts; ii++)
            {
                for (int jj = ii; jj < NumberOfCohorts; jj++)
                {
                    if (nominalTrait[ii] == nominalTrait[jj])
                    {
                        D[ii, jj] = 1.0;
                        D[jj, ii] = 1.0;
                    }
                    else
                    {
                        D[ii, jj] = 0.0;
                        D[jj, ii] = 0.0;
                    }
                }
            }

            return D;
        }


        private double[,] CalculateAggregateDistance(List<double[,]> distanceList)
        {
            int NumberOfTraits = distanceList.Count();
            int NumberOfCohorts = distanceList[0].GetLength(0);
            double[,] D = new double[NumberOfCohorts, NumberOfCohorts];

            foreach (double[,] d in distanceList)
            {
                for (int ii = 0; ii < d.GetLength(0); ii++)
                {
                    for (int jj = ii; jj < d.GetLength(1); jj++)
                    {
                        D[ii, jj] += d[ii, jj] / NumberOfTraits;
                        D[jj, ii] = D[ii, jj];
                    }
                }
            }

            return D;
        }

        private double RaoEntropy(double[,] d, double[] b)
        {
            double TotalB = 0.0;
            double R = 0.0;

            for (int ii = 0; ii < b.Length; ii++)
            {
                TotalB += b[ii];
            }

            for (int ii = 0; ii < d.GetLength(0)-1; ii++)
            {
                for (int jj = 1; jj < d.GetLength(0); jj++)
                {
                    R += d[ii,jj]*b[ii]*b[jj]/(TotalB*TotalB);
                }   
            }

            return R;

        }

    }

    public class Mesh
    {
        public Triangle[] _Mesh;

        public Mesh(Triangle[] mesh)
        {
            _Mesh = mesh;
        }
    }

    public class Triangle
    {
        public Point3D P1;
        public Point3D P2;
        public Point3D P3;
        //Change into tetrahedrons
        public Point3D P4;

        public Triangle(Point3D p1,Point3D p2,Point3D p3)
        {
            P1 = p1;
            P2 = p2;
            P3 = p3;
        }

        public void ToTetrahedron(Point3D p4)
        {
            P4 = p4;
        }

        public Point3D SurfaceNormal()
        {

            //U = p2 - p1
            //V = p3 - p1
            Point3D U = new Point3D(P2.X - P1.X,P2.Y - P1.Y, P2.Z-P1.Z);
            Point3D V = new Point3D(P3.X-P1.X,P3.Y-P1.Y,P3.Z-P3.Y);
            double Nx = (U.Y * V.Z) - (U.Z * V.Y);
            double Ny = (U.Z * V.X) - (U.X * V.Z);
            double Nz = (U.X * V.Y) - (U.Y * V.X);

            return new Point3D(Nx, Ny, Nz);
        }

    }

    public class Point3D
    {
        public double X;

        public double Y;
        public double Z;

        public Point3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }


    }



}

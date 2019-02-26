using System;
using System.Collections.Generic;
using System.Linq;

namespace Madingley
{
    /// <summary>
    /// A revised version of the herbivory process, written November 2011
    /// </summary>
    public partial class RevisedHerbivory: IEatingImplementation
    {
        /// <summary>
        /// Holds the thread-local variables to track numbers of extinctions and productions of cohorts
        /// </summary>
        /// <todoD>Needs a little tidying and checking of access levels</todoD>
        public class ThreadLockedParallelVariables
        {
            /// <summary>
            /// Thread-local variables to track numbers of cohort extinctions and productions
            /// </summary>
            public int Extinctions, Productions;
        }

        /// <summary>
        /// Scalar to convert from the time step units used by this herbivory implementation to global model time step units
        /// </summary>
        private double _DeltaT;
        /// <summary>
        /// Return the scalar to convert from the time step units used by this herbivory implementation to global model time step units
        /// </summary>
        public double DeltaT { get { return _DeltaT; } }

        /// <summary>
        /// The proportion of time that a herbivore cohort spends eating
        /// </summary>
        private double _ProportionOfTimeEating;
        /// <summary>
        /// Get or set the proportion of time that a herbivore cohort devotes to eating behaviours
        /// </summary>
        public double ProportionTimeEating
        { 
            get { return _ProportionOfTimeEating; }
            set { _ProportionOfTimeEating = value; }
        }

        /// <summary>
        /// Jagged array mirroring the grid cell stocks to store the biomasses eaten in herbivory
        /// </summary>
        private double[][][] _BiomassesEaten;
        /// <summary>
        /// Get the jagged array storing the biomasses eaten in herbivory
        /// </summary>
        public double[][][] BiomassesEaten
        { get { return _BiomassesEaten; } }

        /// <summary>
        /// Jagged array mirroring the grid cell stocks to store the potential biomasses eaten (given the rate of encounter) in herbivory
        /// </summary>
        private double[][][] _PotentialBiomassesEaten;
        /// <summary>
        /// Get the jagged array storing the potential biomasses eaten in herbivory
        /// </summary>
        public double[][][] PotentialBiomassesEaten
        { get { return _PotentialBiomassesEaten; } }

        /// <summary>
        /// List of autotroph functional group indices to be eaten in herbivory
        /// </summary>
        private int[] _FunctionalGroupIndicesToEat;
        /// <summary>
        /// Get the list of autotroph functional group indices to eat in herbivory
        /// </summary>
        public int[] FunctionalGroupIndicesToEat
        { get { return _FunctionalGroupIndicesToEat; } }

        /// <summary>
        ///The total biomass eaten by the acting cohort 
        /// </summary>
        private double _TotalBiomassEatenByCohort;
        /// <summary>
        /// Get the total biomass eaten by the acting cohort
        /// </summary>
        public double TotalBiomassEatenByCohort
        {
            get { return _TotalBiomassEatenByCohort; }
        }

        /// <summary>
        /// Cumulative number of time units to handle all of the potential biomass eaten from all autotroph stocks
        /// </summary>
        private double _TimeUnitsToHandlePotentialFoodItems;
        /// <summary>
        /// Get and set the cumulative number of time units to handle all of the potential biomass eaten from all autotroph stocks
        /// </summary>
        public double TimeUnitsToHandlePotentialFoodItems
        {
            get { return _TimeUnitsToHandlePotentialFoodItems; }
            set { _TimeUnitsToHandlePotentialFoodItems = value; }

        }

        /// <summary>
        /// The area (in square km) of the grid cell
        /// </summary>
        private double _CellArea;
        /// <summary>
        /// Get and set the area of the grid cell
        /// </summary>
        public double CellArea 
        { 
            get { return _CellArea; } 
            set { _CellArea = value; } 
        }

        /// <summary>
        /// The area of the current grid cell in hectares
        /// </summary>
        private double _CellAreaHectares;
        /// <summary>
        /// Get and set the area of the current grid cell in hectares
        /// </summary>
        public double CellAreaHectares
        {
            get { return _CellAreaHectares; }
            set { _CellAreaHectares = value; }
        }

        /// <summary>
        /// Individual body mass of herbivores
        /// </summary>
        private double _BodyMassHerbivore;

        // Holds the edible plant mass available
        private double EdibleMass;
        // Holds the scaling to get from exstant autotroph biomass to the edible mass
        private double EdibleScaling;

        /// <summary>
        /// Double to hold the log-optimal prey body size ratio for the acting herbivore cohort
        /// </summary>
        /// 
        private double _HerbivoreLogOptimalPreyBodySizeRatio;

        private double _AttackRateTemperatureScalar;
        private double _HandlingTimeTemperatureScalar;

        /// <summary>
        /// String to hold current stock type
        /// </summary>
        /// 
        private string _PhytoStockType;

        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        private UtilityFunctions Utilities;

        double InitialAbundance;
        double AbundanceEaten;

        double _MaxdistanceOptimalPreyPredRatio = 0.7;


        /// <summary>
        /// Constructor for herbivory: assigns all parameter values
        /// </summary>
        /// <param name="cellArea">The area (in square km) of the grid cell</param>
        /// <param name="globalModelTimeStepUnit">The time step unit used in the model</param>
        public RevisedHerbivory(double cellArea, string globalModelTimeStepUnit)
        {
            InitialiseParametersHerbivory();

            // Initialise the utility functions
            Utilities = new UtilityFunctions();

            // Calculate the scalar to convert from the time step units used by this implementation of herbivory to the global model time step units
            _DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, _TimeUnitImplementation);
                        
            // Store the specified cell area in this instance of this herbivory implementation
            _CellArea = cellArea;
            _CellAreaHectares = cellArea * 100;


            
        }

        /// <summary>
        /// Initialises herbivory implementation each time step
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the current grid cell</param>
        /// <param name="gridCellStocks">The stocks in the current grid cell</param>
        /// <param name="madingleyCohortDefinitions">The definitions for cohorts in the model</param>
        /// <param name="madingleyStockDefinitions">The definitions for stocks in the model</param>
        /// <remarks>This only works if: a) herbivory is initialised in every grid cell; and b) if parallelisation is done by latitudinal strips
        /// It is critical to run this every time step</remarks>
        public void InitializeEatingPerTimeStep(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions)
        {
            // Get the functional group indices of all autotroph stocks
            _FunctionalGroupIndicesToEat = madingleyStockDefinitions.GetFunctionalGroupIndex("Heterotroph/Autotroph", "Autotroph", false);          
        }

        /// <summary>
        /// Calculate the potential biomass that could be gained through herbivory on each grid cell autotroph stock
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the grid cell</param>
        /// <param name="gridCellStocks">The stocks in the grid cell</param>
        /// <param name="actingCohort">The acting cohort</param>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        /// <param name="madingleyCohortDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="madingleyStockDefinitions">The functional group definitions for stocks  in the model</param>
        public void GetEatingPotentialTerrestrial(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, 
            int[] actingCohort, SortedList<string, double[]> cellEnvironment, FunctionalGroupDefinitions madingleyCohortDefinitions, 
            FunctionalGroupDefinitions madingleyStockDefinitions, uint currentMonth, uint currentTimeStep)
        {
            // Set the total biomass eaten by the acting cohort to zero
            _TotalBiomassEatenByCohort = 0.0;

            // Get the individual body mass of the acting cohort
            _BodyMassHerbivore = gridCellCohorts[actingCohort].IndividualBodyMass;
            
            // Set the total number of units to handle all potential biomass eaten to zero
            _TimeUnitsToHandlePotentialFoodItems = 0.0;

            // Initialise the jagged arrays to hold the potential and actual biomass eaten in each of the grid cell autotroph stocks
            _BiomassesEaten = new double[gridCellStocks.Count][][];
            _PotentialBiomassesEaten = new double[gridCellStocks.Count][][];

            _AttackRateTemperatureScalar = 
                (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", actingCohort[0]) == "ectotherm") ?
                Math.Exp(AttackRateActivationEnergy*(cellEnvironment["Temperature"][currentMonth]+273.15-_ReferenceTemperature)/
                (_BoltzmannConstant*_ReferenceTemperature*(cellEnvironment["Temperature"][currentMonth]+273.15))): 1.0;


            _HandlingTimeTemperatureScalar =
                (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", actingCohort[0]) == "ectotherm") ?
                Math.Exp(HandlingTimeActivationEnergy * (cellEnvironment["Temperature"][currentMonth] + 273.15 - _ReferenceTemperature) /
                (_BoltzmannConstant * _ReferenceTemperature * (cellEnvironment["Temperature"][currentMonth] + 273.15))) : 1.0;


            // Loop over rows in the jagged arrays and initialise each vector
            for (int i = 0; i < gridCellStocks.Count; i++)
            {
                _BiomassesEaten[i] = new double[gridCellStocks[i].Count][];
                _PotentialBiomassesEaten[i] = new double[gridCellStocks[i].Count][];
                for (int j = 0; j < gridCellStocks[i].Count; j++)
                {
                    _BiomassesEaten[i][j] = new double[1];
                    _PotentialBiomassesEaten[i][j] = new double[1];
                }
            }

            

            // Loop over functional groups that can be eaten
            foreach (int FunctionalGroup in _FunctionalGroupIndicesToEat)
            {
                // Loop over stocks within the functional group
                for (int i = 0; i < gridCellStocks[FunctionalGroup].Count; i++)
                {
                    // Get the mass from this stock that is available for eating (assumes only 10% is edible)
                    EdibleMass = gridCellStocks[FunctionalGroup][i].TotalBiomass* 0.1;

                    // Calculate the potential biomass eaten from this stock by the acting cohort
                    _PotentialBiomassesEaten[FunctionalGroup][i][1] = CalculatePotentialBiomassEatenTerrestrial(EdibleMass, _BodyMassHerbivore) * _AttackRateTemperatureScalar;

                    // Add the time required to handle the potential biomass eaten from this stock to the cumulative total for all stocks
                    _TimeUnitsToHandlePotentialFoodItems += _PotentialBiomassesEaten[FunctionalGroup][i][1] *
                        CalculateHandlingTimeTerrestrial(_BodyMassHerbivore)*_HandlingTimeTemperatureScalar;
                    
                }
            }

        }

        /// <summary>
        /// Calculate the potential biomass that could be gained through herbivory on each grid cell autotroph stock
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the grid cell</param>
        /// <param name="gridCellStocks">The stocks in the grid cell</param>
        /// <param name="actingCohort">The acting cohort</param>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        /// <param name="madingleyCohortDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="madingleyStockDefinitions">The functional group definitions for stocks  in the model</param>
        public void GetEatingPotentialMarine(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, 
            int[] actingCohort, SortedList<string, double[]> cellEnvironment, FunctionalGroupDefinitions madingleyCohortDefinitions, 
            FunctionalGroupDefinitions madingleyStockDefinitions, uint currentMonth, uint currentTimeStep)
        {
            // Set the total biomass eaten by the acting cohort to zero
            _TotalBiomassEatenByCohort = 0.0;

            // Get current temperature in grid cell
            double temperature = cellEnvironment["Temperature"][currentMonth];

            // Calculate the coefficients for phytoplankton abundance-size curve
            double no3 = cellEnvironment["no3"][currentMonth];
            double[] slopeCoefficient = new double[2];

            //Tuple<double, double> regCoefs = GetPhytoAbundSizeSlope(temperature, no3);

            // Get the individual body mass of the acting cohort
            _BodyMassHerbivore = gridCellCohorts[actingCohort].IndividualBodyMass;

            // Set the total number of units to handle all potential biomass eaten to zero
            _TimeUnitsToHandlePotentialFoodItems = 0.0;

            // Initialise the jagged arrays to hold the potential and actual biomass eaten in each of the grid cell autotroph stocks
            _BiomassesEaten = new double[gridCellStocks.Count][][];
            _PotentialBiomassesEaten = new double[gridCellStocks.Count][][];

            _AttackRateTemperatureScalar =
                (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", actingCohort[0]) == "ectotherm") ?
                Math.Exp(AttackRateActivationEnergy * (cellEnvironment["Temperature"][currentMonth] + 273.15 - _ReferenceTemperature) /
                (_BoltzmannConstant * _ReferenceTemperature * (cellEnvironment["Temperature"][currentMonth] + 273.15))) : 1.0;


            _HandlingTimeTemperatureScalar =
                (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", actingCohort[0]) == "ectotherm") ?
                Math.Exp(HandlingTimeActivationEnergy * (cellEnvironment["Temperature"][currentMonth] + 273.15 - _ReferenceTemperature) /
                (_BoltzmannConstant * _ReferenceTemperature * (cellEnvironment["Temperature"][currentMonth] + 273.15))) : 1.0;



            // Loop over rows in the jagged arrays and initialise each vector
            for (int i = 0; i < gridCellStocks.Count; i++)
            {
                _BiomassesEaten[i] = new double[gridCellStocks[i].Count][];
                _PotentialBiomassesEaten[i] = new double[gridCellStocks[i].Count][];
                for (int j = 0; j < gridCellStocks[i].Count; j++)
                {
                    _BiomassesEaten[i][j] = new double[gridCellStocks[i][j].SizeStructure.Length];
                    _PotentialBiomassesEaten[i][j] = new double[gridCellStocks[i][j].SizeStructure.Length];
                }
            }

            _HerbivoreLogOptimalPreyBodySizeRatio = gridCellCohorts[actingCohort[0]][actingCohort[1]].LogOptimalPreyBodySizeRatio;
            //
            double LogHerbivoreOptimatePreySize = Math.Log10(Math.Exp(_HerbivoreLogOptimalPreyBodySizeRatio) *_BodyMassHerbivore);
            double UpperHerbivorePreySizeWindow = Math.Log10(Math.Exp(_HerbivoreLogOptimalPreyBodySizeRatio) *_BodyMassHerbivore * Math.Exp(_MaxdistanceOptimalPreyPredRatio));
            double LowerHerbivorePreySizeWindow = Math.Log10(Math.Exp(_HerbivoreLogOptimalPreyBodySizeRatio) *_BodyMassHerbivore * Math.Exp(-_MaxdistanceOptimalPreyPredRatio));

            // Calculate the inidividual herbivory rate per unit autotroph mass-density per hectare
            SpecificHerbivoreKillRateConstant = CalculateIndividualHerbivoryRatePerHectareMarine(_BodyMassHerbivore);

            // Loop over functional groups that can be eaten
            foreach (int FunctionalGroup in _FunctionalGroupIndicesToEat)
            {
                // Loop over stocks within the functional group
                for (int i = 0; i < gridCellStocks[FunctionalGroup].Count; i++)
                {
                    //Calculate the stock mass bins that this herbivore can eat from - assume that the herbivore has a window of biomass magnitude around the optimum
                    var inds = Enumerable.Range(0, gridCellStocks[FunctionalGroup][i].SizeStructure.Length).
                        Where(b => (gridCellStocks[FunctionalGroup][i].SizeBinEdges[b] >= LowerHerbivorePreySizeWindow ) &&
                            (gridCellStocks[FunctionalGroup][i].SizeBinCentres[b] <= UpperHerbivorePreySizeWindow));

                    //Find the bin in which the optimal prey size fits
                    //var inds = Enumerable.Range(0, gridCellStocks[FunctionalGroup][i].SizeBinEdges.Length-1).
                    //    Where(b => (gridCellStocks[FunctionalGroup][i].SizeBinEdges[b] < LogHerbivoreOptimatePreySize) &&
                    //        (gridCellStocks[FunctionalGroup][i].SizeBinEdges[b+1] > LogHerbivoreOptimatePreySize));

                    //Now loop over these stock size bin inds
                    foreach (var b in inds)
                    {
                        EdibleMass = gridCellStocks[FunctionalGroup][i].SizeStructure[b];

                        _PhytoStockType = madingleyStockDefinitions.GetTraitNames("stock name", FunctionalGroup);
                        // Calculate the potential biomass eaten from this stock by the acting cohort
                        _PotentialBiomassesEaten[FunctionalGroup][i][b] = CalculatePotentialBiomassEatenMarine(EdibleMass, _BodyMassHerbivore,Math.Pow(10,gridCellStocks[FunctionalGroup][i].SizeBinCentres[b]),
                            _HerbivoreLogOptimalPreyBodySizeRatio, _PhytoStockType)*_AttackRateTemperatureScalar;

                        // Add the time required to handle the potential biomass eaten from this stock to the cumulative total for all stocks
                        _TimeUnitsToHandlePotentialFoodItems += _PotentialBiomassesEaten[FunctionalGroup][i][b] *
                            CalculateHandlingTimeMarine(_BodyMassHerbivore, Math.Pow(10,gridCellStocks[FunctionalGroup][i].SizeBinCentres[b])) * _HandlingTimeTemperatureScalar;
                    }
/*                    // Get the mass from this stock that is available for eating (assumes all marine autotrophic organisms are edible)
                    //EdibleMass = gridCellStocks[FunctionalGroup][i].TotalBiomass * 0.1;
                    EdibleMass = gridCellStocks[FunctionalGroup][i].TotalBiomass;

                    _PhytoStockType = madingleyStockDefinitions.GetTraitNames("stock name", FunctionalGroup);
                    // Calculate the potential biomass eaten from this stock by the acting cohort
                    _PotentialBiomassesEaten[FunctionalGroup][i] = CalculatePotentialBiomassEatenMarine(EdibleMass, _BodyMassHerbivore,
                        _HerbivoreLogOptimalPreyBodySizeRatio, _PhytoStockType, temperature, regCoefs,
                        gridCellStocks[FunctionalGroup][i].SizeBinCentres,gridCellStocks[FunctionalGroup][i].SizeStructure)*_AttackRateTemperatureScalar;

                    // Add the time required to handle the potential biomass eaten from this stock to the cumulative total for all stocks
                    _TimeUnitsToHandlePotentialFoodItems += _PotentialBiomassesEaten[FunctionalGroup][i] *
                        CalculateHandlingTimeMarine(_BodyMassHerbivore)*_HandlingTimeTemperatureScalar;
*/
                }
            }

        }

        /// <summary>
        /// Calculate the actual amount eaten in herbivory, apply the changes to the eaten autotroph stocks, and update deltas for the herbivore cohort
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in this grid cell</param>
        /// <param name="gridCellStocks">The stocks in this grid cell</param>
        /// <param name="actingCohort">The acting cohort</param>
        /// <param name="cellEnvironment">The environmental conditions in this grid cell</param>
        /// <param name="deltas">The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell</param>
        /// <param name="madingleyCohortDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="madingleyStockDefinitions">The functional group definitions for stocks in the model</param>
        /// <param name="trackProcesses">An instance of ProcessTracker to hold diagnostics for herbivory</param>
        /// <param name="currentTimestep">The current model time step</param>
        /// <param name="specificLocations">Whether the model is being run for specific locations</param>
        /// <param name="outputDetail">The level of output detail being used in this model run</param>
        /// <param name="initialisation">The Madingley Model initialisation</param>
        public void RunEating(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, int[] actingCohort, SortedList<string, double[]>
            cellEnvironment, Dictionary<string, Dictionary<string, double>> deltas, FunctionalGroupDefinitions madingleyCohortDefinitions,
            FunctionalGroupDefinitions madingleyStockDefinitions, ProcessTracker trackProcesses, FunctionalGroupTracker functionalTracker, HighResFGTracker highResFGTracker,
            CohortTracker cohortTracker, uint currentTimestep, Boolean specificLocations, string outputDetail, MadingleyModelInitialisation initialisation)
        {

            EdibleScaling = 1.0;
            if (cellEnvironment["Realm"][0] == 1.0) EdibleScaling = 0.1;

            // Loop over autotroph functional groups that can be eaten
            foreach (int FunctionalGroup in _FunctionalGroupIndicesToEat)
            {
                // Loop over stocks within the functional groups
                for (int i = 0; i < gridCellStocks[FunctionalGroup].Count; i++)
                {
                    for (int b = 0; b < _BiomassesEaten[FunctionalGroup][i].Length; b++)
                    {
                        if (cellEnvironment["Realm"][0] == 2.0)
                        {
                            // Get the mass from this stock that is available for eating (assumes only 10% is edible in the terrestrial realm)
                            EdibleMass = gridCellStocks[FunctionalGroup][i].SizeStructure[b] ;
                        }
                        else
                        {
                            EdibleMass = gridCellStocks[FunctionalGroup][i].TotalBiomass * EdibleScaling;
                        }

                        // Calculate the biomass actually eaten from this stock by the acting cohort
                        _BiomassesEaten[FunctionalGroup][i][b] = CalculateBiomassesEaten(_PotentialBiomassesEaten[FunctionalGroup][i][b],
                            _TimeUnitsToHandlePotentialFoodItems, gridCellCohorts[actingCohort].CohortAbundance, EdibleMass);

                        AbundanceEaten = _BiomassesEaten[FunctionalGroup][i][b];
                        InitialAbundance = EdibleMass;

                        gridCellCohorts[actingCohort].TrophicIndex += _BiomassesEaten[FunctionalGroup][i][b];

                        InitialAbundance = gridCellStocks[FunctionalGroup][i].SizeStructure[b];
                        AbundanceEaten = _BiomassesEaten[FunctionalGroup][i][b];

                        if (cellEnvironment["Realm"][0] == 2.0)
                        {
                            // Get the mass from this stock that is available for eating (assumes only 10% is edible in the terrestrial realm)
                            gridCellStocks[FunctionalGroup][i].SizeStructure[b] -= _BiomassesEaten[FunctionalGroup][i][b];
                            gridCellStocks[FunctionalGroup][i].TotalBiomass -= _BiomassesEaten[FunctionalGroup][i][b];
                        }
                        else
                        {
                            // Remove the biomass eaten from the autotroph stock
                            gridCellStocks[FunctionalGroup][i].TotalBiomass -= _BiomassesEaten[FunctionalGroup][i][b];
                        }

                        // If track processes has been specified, then track flow between primary producers and herbivores
                        if ((trackProcesses.TrackProcesses) && (initialisation.HighResSlowTrackingOn == true) && (currentTimestep >= initialisation.TimeStepToStartProcessTrackers))
                        {
                            // Track flows between functional groups
                            highResFGTracker.RecordFGFlow(
                                madingleyCohortDefinitions.GetTraitNames("group description", gridCellCohorts[actingCohort].FunctionalGroupIndex),
                                madingleyStockDefinitions.GetTraitNames("stock name", FunctionalGroup), gridCellCohorts[actingCohort].IndividualBodyMass,
                                Math.Pow(10,gridCellStocks[FunctionalGroup][i].SizeBinCentres[b]),
                                InitialAbundance, AbundanceEaten, cellEnvironment["Realm"][0] == 2.0);

                        }

                    }

                    // If track processes has been specified, then track flow between primary producers and herbivores
                    if((trackProcesses.TrackProcesses) && (currentTimestep >= initialisation.TimeStepToStartProcessTrackers))
                    {
                        // Track flows between functional groups
                        functionalTracker.RecordFGFlow((uint)cellEnvironment["LatIndex"][0], (uint)cellEnvironment["LonIndex"][0],
                            madingleyCohortDefinitions.GetTraitNames("group description", gridCellCohorts[actingCohort].FunctionalGroupIndex),
                            madingleyStockDefinitions.GetTraitNames("stock name", FunctionalGroup),
                            _BiomassesEaten[FunctionalGroup][i].Sum(), cellEnvironment["Realm"][0] == 2.0);
                    }

                    // If the model is being run for specific locations and if track processes has been specified, then track the mass flow between
                    // primary producer and herbivore (this is the old tracker)
                    if (specificLocations && trackProcesses.TrackProcesses && (currentTimestep >= initialisation.TimeStepToStartProcessTrackers))
                    {
                        trackProcesses.RecordHerbivoryMassFlow(currentTimestep, _BodyMassHerbivore, _BiomassesEaten[FunctionalGroup][i].Sum());
                    }

                    // If track processes has been specified and the output detail level is set to high and the model is being run for specific locations,
                    // then track the flow of mass between trophic levels
                    if (trackProcesses.TrackProcesses && (outputDetail == "high") && specificLocations)
                    {
                        trackProcesses.TrackHerbivoryTrophicFlow((uint)cellEnvironment["LatIndex"][0], (uint)cellEnvironment["LonIndex"][0],
                            gridCellCohorts[actingCohort].FunctionalGroupIndex, madingleyCohortDefinitions, _BiomassesEaten[FunctionalGroup][i].Sum(), _BodyMassHerbivore, initialisation, cellEnvironment["Realm"][0] == 2.0);

                    }


                    // Check that the biomass eaten is not a negative value
                    // Commented out for purposes of speed
                    //Debug.Assert(_BiomassesEaten[FunctionalGroup][i] >= 0,
                    //    "Herbivory negative for this herbivore cohort" + actingCohort);
                    
                    // Add the biomass eaten and assimilated by an individual to the delta biomass for the acting cohort
                    deltas["biomass"]["herbivory"] += _BiomassesEaten[FunctionalGroup][i].Sum() * AssimilationEfficiency / gridCellCohorts[actingCohort].CohortAbundance;

                    // Move the biomass eaten but not assimilated by an individual into the organic matter pool
                    deltas["organicpool"]["herbivory"] += _BiomassesEaten[FunctionalGroup][i].Sum() * (1 - AssimilationEfficiency);
                

                }

                // Check that the delta biomass from eating for the acting cohort is not negative
                // Commented out for the purposes of speed
                //Debug.Assert(deltas["biomass"]["herbivory"] >= 0, "Delta biomass from herbviory is negative");
                //Debug.Assert(!double.IsNaN(deltas["biomass"]["herbivory"]), "Delta biomass from herbviory is NaN");

                // Calculate the total biomass eaten by the acting (herbivore) cohort
                _TotalBiomassEatenByCohort = deltas["biomass"]["herbivory"] * gridCellCohorts[actingCohort].CohortAbundance;

                

            }
        }


    }


}

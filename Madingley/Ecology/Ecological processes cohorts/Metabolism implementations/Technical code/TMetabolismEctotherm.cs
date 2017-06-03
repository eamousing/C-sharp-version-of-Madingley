﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Madingley
{
    /// <summary>
    /// A formulation of the metabolism process for Ectothermic organisms
    /// </summary>
    /// <remarks>Functional form and parameters taken from fitted relationship in Brown's (2004) Metabolic Theory of Ecology.
    /// Currently mass assigned to reproductive potential is not metabolised
    /// Assumes that ectothermic organisms have a body temperature equal to the ambient temperature,
    /// therefore metabolising at that ambient temperature</remarks>
    public partial class MetabolismEctotherm : IMetabolismImplementation
    {

        /// <summary>
        /// Scalar to convert from the time units used by this metabolism implementation to the global model time step units
        /// </summary>
        private double _DeltaT;
        /// <summary>
        /// Get the scalar to convert from the time units used by this metabolism implementation to the global model time step units
        /// </summary>
        public double DeltaT  { get  {  return _DeltaT;  } }

        /// <summary>
        /// Constant to convert temperature in degrees Celsius to temperature in Kelvin
        /// </summary>
        private double _TemperatureUnitsConvert;

        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        private UtilityFunctions Utilities;

        /// <summary>
        /// Whether the proportion of time that the cohort is active has been recalculated this time step
        /// </summary>
        private Boolean _ProportionTimeActiveCalculatedThisTimestep;
        /// <summary>
        /// Get whether the proportion of time that the cohort is active has been recalculated this time step
        /// </summary>
        public Boolean ProportionTimeActiveCalculatedThisTimestep
        {
            get { return _ProportionTimeActiveCalculatedThisTimestep; }
        }

        /// <summary>
        /// Constructor for metabolism: assigns all parameter values
        /// </summary>
        public MetabolismEctotherm(string globalModelTimeStepUnit)
        {
            
            // Initialise ecological parameters for metabolism
            InitialiseMetabolismParameters();

            // Initialise the utility functions
            Utilities = new UtilityFunctions();

            // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
            _DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, _TimeUnitImplementation);

            _ProportionTimeActiveCalculatedThisTimestep = false;
        }

        /// <summary>
        /// Run metabolism for the acting cohort
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the current grid cell</param>
        /// <param name="gridCellStocks">The stocks in the current grid cell</param>
        /// <param name="actingCohort">The position of the acting cohort in the jagged array of grid cell cohorts</param>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        /// <param name="deltas">The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell</param>
        /// <param name="madingleyCohortDefinitions">The definitions for cohort functional groups in the model</param>
        /// <param name="madingleyStockDefinitions">The definitions for the stock functional groups in the model</param>
        /// <param name="currentTimestep">The current model time step</param>
        /// <param name="currentMonth">The current model month</param>
        public void RunMetabolism(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, 
            int[] actingCohort, SortedList<string, double[]> cellEnvironment, Dictionary<string, Dictionary<string, double>> 
            deltas, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, 
            uint currentTimestep, uint currentMonth)
        {

            // Add in the total biomass which is respiring to the appropriate delta
            deltas["biomass"]["respiring biomass"] += gridCellCohorts[actingCohort].IndividualBodyMass * gridCellCohorts[actingCohort].CohortAbundance;

            // Calculate metabolic loss for an individual and add the value to the delta biomass for metabolism
            deltas["biomass"]["metabolism"] = -CalculateIndividualMetabolicRate(gridCellCohorts[actingCohort].IndividualBodyMass,
                cellEnvironment["Temperature"][currentMonth] + _TemperatureUnitsConvert, gridCellCohorts[actingCohort].ProportionTimeActive) * _DeltaT;


            // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
            deltas["biomass"]["metabolism"] = Math.Max(deltas["biomass"]["metabolism"], -(gridCellCohorts[actingCohort].IndividualBodyMass + deltas["biomass"]["predation"] + deltas["biomass"]["herbivory"]));

            // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
            deltas["respiratoryCO2pool"]["metabolism"] = -deltas["biomass"]["metabolism"] * gridCellCohorts[actingCohort].CohortAbundance;

        }

    }
}

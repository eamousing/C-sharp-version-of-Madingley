using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

using System.IO;

namespace Madingley
{
    /// <summary>
    /// A class to specify, initialise and run ecological processes pertaining to stocks
    /// </summary>
    class EcologyStock
    {
        /// <summary>
        /// An instance of the nutrient plankton model
        /// </summary>
        NutrientsPlanktonModel NPZModel;

        /// <summary>
        /// An instance of the Autotroph Processor for this model
        /// </summary>
        //AutotrophProcessor MarineNPPtoAutotrophStock;

        /// <summary>
        /// An instance of the plant model class
        /// </summary>
        RevisedTerrestrialPlantModel DynamicPlantModel;

        /// <summary>
        /// An instance of the class for human appropriation of NPP
        /// </summary>
        HumanAutotrophMatterAppropriation HANPP;



        public void InitializeEcology()
        {
            // Initiate the NPZ model
            NPZModel = new NutrientsPlanktonModel();

            //Initialize the autotrophprocessor
            //MarineNPPtoAutotrophStock = new AutotrophProcessor();

            // Initialise the plant model
            DynamicPlantModel = new RevisedTerrestrialPlantModel();

            // Initialise the human NPP appropriation class
            HANPP = new HumanAutotrophMatterAppropriation();

        }

        /// <summary>
        /// Run the NPZ model in a particular grid cell
        /// </summary>
        /// <param name="cellEnvironment"></param>
        /// <param name="currentMonth"></param>
        /// <param name="stockDefinitions"></param>
        /// <param name="gridCellStocks"></param>
        /// <returns></returns>
        public double[] RunWithinCellNPZModel(SortedList<string, double[]> cellEnvironment, uint currentMonth, GridCellStockHandler gridCellStocks)
        {
            // Run the NPZ model
            return NPZModel.RunNPZModel(cellEnvironment, currentMonth);
        }

        /// <summary>
        /// Run ecological processes that operate on stocks within a single grid cell
        /// </summary>
        ///<param name="gridCellStocks">The stocks in the current grid cell</param>
        ///<param name="actingStock">The acting stock</param>
        ///<param name="cellEnvironment">The stocks in the current grid cell</param>
        ///<param name="environmentalDataUnits">List of units associated with the environmental variables</param>
        ///<param name="humanNPPScenario">The human appropriation of NPP scenario to apply</param>
        ///<param name="madingleyStockDefinitions">The functional group definitions for stocks in the model</param>
        ///<param name="currentTimeStep">The current model time step</param>
        ///<param name="burninSteps">The number of time steps to spin the model up for before applying human impacts</param>
        ///<param name="impactSteps">The number of time steps to apply human impacts for</param>
        ///<param name="globalModelTimeStepUnit">The time step unit used in the model</param>
        ///<param name="trackProcesses">Whether to track properties of ecological processes</param>
        ///<param name="tracker">An instance of the ecological process tracker</param>
        ///<param name="globalTracker">An instance of the global process tracker</param>
        ///<param name="currentMonth">The current model month</param>
        ///<param name="outputDetail">The level of detail to use in outputs</param>
        ///<param name="specificLocations">Whether to run the model for specific locations</param>
        ///<param name="impactCell">Whether this cell should have human impacts applied</param>
        public void RunWithinCellEcology(GridCellStockHandler gridCellStocks, int[] actingStock, SortedList<string, double[]> cellEnvironment,
            SortedList<string, string> environmentalDataUnits, Tuple<string, double, double> humanNPPScenario,
            FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, 
            uint currentTimeStep, uint burninSteps, uint impactSteps,uint recoverySteps, uint instantStep, uint numInstantSteps, string globalModelTimeStepUnit, Boolean trackProcesses, 
            ProcessTracker tracker, FunctionalGroupTracker functionalTracker, GlobalProcessTracker globalTracker, uint currentMonth, 
            string outputDetail, bool specificLocations, Boolean impactCell, double currentPlanktonStock)
        {
            if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "marine")
            {
                // Update biomasses of marine stocks after running the NPZ model
                UpdateMarineStocks(currentPlanktonStock, actingStock, gridCellStocks);
                
                // Run the autotroph processor
                //MarineNPPtoAutotrophStock.ConvertNPPToAutotroph(madingleyCohortDefinitions, madingleyStockDefinitions, cellEnvironment, gridCellStocks, actingStock, environmentalDataUnits["LandNPP"], 
                //    environmentalDataUnits["OceanNPP"], currentTimeStep,globalModelTimeStepUnit,tracker, functionalTracker, globalTracker ,outputDetail,specificLocations,currentMonth);

                
            }
            else if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "terrestrial")
            {

                // Run the dynamic plant model to update the leaf stock for this time step
                double WetMatterNPP = DynamicPlantModel.UpdateLeafStock(cellEnvironment, gridCellStocks, actingStock, currentTimeStep, madingleyStockDefinitions.
                    GetTraitNames("leaf strategy", actingStock[0]).Equals("deciduous"), globalModelTimeStepUnit, tracker, globalTracker, currentMonth,
                    outputDetail, specificLocations);
                        
                double fhanpp = HANPP.RemoveHumanAppropriatedMatter(WetMatterNPP, cellEnvironment, humanNPPScenario, gridCellStocks, actingStock, 
                    currentTimeStep,burninSteps,impactSteps,recoverySteps,instantStep, numInstantSteps,impactCell, globalModelTimeStepUnit);

                // Apply human appropriation of NPP
                gridCellStocks[actingStock].TotalBiomass += WetMatterNPP * (1.0 - fhanpp);
                if (globalTracker.TrackProcesses)
                {
                    globalTracker.RecordHANPP((uint)cellEnvironment["LatIndex"][0], (uint)cellEnvironment["LonIndex"][0], (uint)actingStock[0],
                        fhanpp);
                }

                if (gridCellStocks[actingStock].TotalBiomass < 0.0) gridCellStocks[actingStock].TotalBiomass = 0.0;

            }
            else
            {
                Debug.Fail("Stock must be classified as belonging to either the marine or terrestrial realm");
            }
        }

        /// <summary>
        /// Updates marine stock biomasses
        /// </summary>
        /// <param name="StockPlanktonBiomass"></param>
        /// <param name="actingStock"></param>
        /// <param name="gridCellStocks"></param>
        public void UpdateMarineStocks(double StockPlanktonBiomass, int[] actingStock, GridCellStockHandler gridCellStocks)
        {
            gridCellStocks[actingStock].TotalBiomass += StockPlanktonBiomass;
        }
    }
}

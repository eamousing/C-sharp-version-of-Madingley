using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// A revised version of the herbivory process, written November 2011
    /// </summary>
    public partial class RevisedHerbivory: IEatingImplementation
    {
        /// <summary>
        /// The time unit associated with this herbivory implementation and its parameters
        /// </summary>
        private string _TimeUnitImplementation;
        /// <summary>
        /// Get the time unit associated with this herbivory implementation and its parameters
        /// </summary>
        public string TimeUnitImplementation { get { return _TimeUnitImplementation; } }

        /// <summary>
        /// The assimilation efficiency of eaten autotroph mass into herbivore body mass
        /// </summary>
        private double _AssimilationEfficiency;
        /// <summary>
        /// Get and set the assimilation efficiency of eaten autotroph mass into herbivore body mass
        /// </summary>
        public double AssimilationEfficiency
        {
            get { return _AssimilationEfficiency; }
            set { _AssimilationEfficiency = value; }
        }

        /// <summary>
        /// The scalar of the relationship between handling time and the function of herbivore mass for the terrestrial realm
        /// </summary>
        private double _HandlingTimeScalarTerrestrial;
        /// <summary>
        /// Get the scalar of the relationship between handling time and the function of herbivore mass for the terrestrial realm
        /// </summary>
        public double HandlingTimeScalarTerrestrial
        { get { return _HandlingTimeScalarTerrestrial; } }

        /// <summary>
        /// The scalar of the relationship between handling time and the function of herbivore mass for the marine realm
        /// </summary>
        private double _HandlingTimeScalarMarine;
        /// <summary>
        /// Get the scalar of the relationship between handling time and the function of herbivore mass for the marine realm
        /// </summary>
        public  double HandlingTimeScalarMarine
        { get { return _HandlingTimeScalarMarine; } }
        

        /// <summary>
        /// The exponent applied to herbivore mass in the handling time relationship for the terrestrial realm
        /// </summary>
        private double _HandlingTimeExponentTerrestrial;
        /// <summary>
        /// Get the exponent applied to herbivore mass in the handling time relationship for the terrestrial realm
        /// </summary>
        public double HandlingTimeExponentTerrestrial
        { get { return _HandlingTimeExponentTerrestrial; } }

        /// <summary>
        /// The exponent applied to herbivore mass in the handling time relationship for the marine realm
        /// </summary>
        private double _HandlingTimeExponentMarine;
        /// <summary>
        /// The exponent applied to herbivore mass in the handling time relationship for the marine realm
        /// </summary>
        public double HandlingTimeExponentMarine
        { get { return _HandlingTimeExponentMarine; } }
        

        /// <summary>
        /// Reference mass of plant matter for calculating handling times
        /// </summary>
        public double _ReferenceMass;
        /// <summary>
        /// Get the reference mass of plant matter for calculating handling times
        /// </summary>
        public double ReferenceMass {  get { return _ReferenceMass; } }
        
        /// <summary>
        /// The maximum herbivory rate for a herbivore of 1 g in the terrestrial realm
        /// </summary>
        private double _HerbivoryRateConstantTerrestrial;
        /// <summary>
        /// Get the maximum herbivory rate for a herbivore of 1 g
        /// </summary>
        public double HerbivoryRateConstantTerrestrial
        { get { return _HerbivoryRateConstantTerrestrial; } }

        /// <summary>
        /// The maximum herbivory rate for a herbivore of 1 g in the marine realm
        /// </summary>
        private double _HerbivoryRateConstantMarine;
        /// <summary>
        /// Get the maximum herbivory rate for a herbivore of 1 g
        /// </summary>
        public double HerbivoryRateConstantMartine
        { get { return _HerbivoryRateConstantMarine; } }

        /// <summary>
        /// The exponent to apply to body mass in the relationship between body mass and herbivory rate
        /// </summary>
        private double _HerbivoryRateMassExponent;
        /// <summary>
        /// Get and set the exponent to apply to body mass in the relationship between body mass and herbivory rate
        /// </summary>
        public double HerbivoryRateMassExponent
        {   get { return _HerbivoryRateMassExponent; }      }

        /// <summary>
        /// The exponent applied to prey density when calculating attack rates for organisms in the terrestrial realm
        /// </summary>
        private double _AttackRateExponentTerrestrial;
        /// <summary>
        /// Get and set the exponent applied to prey density when calculating attack rates for organisms in the terrestrial realm
        /// </summary>
        public double AttackRateExponentTerrestrial
        {    get { return _AttackRateExponentTerrestrial; }  }

        /// <summary>
        /// The exponent applied to prey density when calculating attack rates for organisms in the marine realm
        /// </summary>
        private double _AttackRateExponentMarine;

        /// <summary>
        /// Get and set the exponent applied to prey density when calculating attack rates for organisms in the marine realm
        /// </summary>
        public double AttackRateExponentMarine
        {  get { return _AttackRateExponentMarine; }  }

        // Variable to hold the instantaneous fraction of the autotroph stock biomass that is eaten
        double InstantFractionEaten;


        /// <summary>
        /// Set and get the standard deviation in attack rates around the optimal prey to predator mass ratio
        /// </summary>
        ///
        private double _FeedingPreferenceStandardDeviation;
        public double FeedingPreferenceStandardDeviation { get { return _FeedingPreferenceStandardDeviation; } }

        private double RelativeFeedingPreference;

        /// <summary>
        /// Pre-calculate the maximum kill rate for a specific predator of 1 g on prey of an optimal size
        /// </summary>
        private double _SpecificHerbivoreKillRateConstant;

        /// <summary>
        /// Get the pre-calculated maximum kill rate for a specific predator of 1 g on prey of an optimal size
        /// </summary>
        public double SpecificHerbivoreKillRateConstant
        {
            get { return _SpecificHerbivoreKillRateConstant; }
            set { _SpecificHerbivoreKillRateConstant = value; }
        }

        private double _AttackRateActivationEnergy;
        public double AttackRateActivationEnergy{ get { return _AttackRateActivationEnergy; } }

        private double _HandlingTimeActivationEnergy;
        public double HandlingTimeActivationEnergy { get { return _HandlingTimeActivationEnergy; } }

        private double _ReferenceTemperature;
        public double ReferenceTemperature { get { return _ReferenceTemperature; } }

        private double _BoltzmannConstant;
        public double BoltzmannConstant { get { return _BoltzmannConstant; } }
        


        /// <summary>
        /// An instance of the simple random number generator class
        /// </summary>
        /// 
        private NonStaticSimpleRNG RandomNumberGenerator = new NonStaticSimpleRNG();

        public void InitialiseParametersHerbivory()
        {
            _TimeUnitImplementation =
                EcologicalParameters.TimeUnits[(int)EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.TimeUnitImplementation"]];
            _HandlingTimeScalarTerrestrial = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Terrestrial.HandlingTimeScalar"];
            _HandlingTimeScalarMarine = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Marine.HandlingTimeScalar"];
            _HandlingTimeExponentTerrestrial = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Terrestrial.HandlingTimeExponent"];
            _HandlingTimeExponentMarine = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Marine.HandlingTimeExponent"];
            _HerbivoryRateConstantTerrestrial = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Terrestrial.RateConstant"];
            _HerbivoryRateConstantMarine = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Marine.RateConstant"];
            _HerbivoryRateMassExponent = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.RateMassExponent"];
            _AttackRateExponentTerrestrial = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Terrestrial.AttackRateExponent"];
            _AttackRateExponentMarine = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.Marine.AttackRateExponent"];
            _ReferenceMass = EcologicalParameters.Parameters["Herbivory.RevisedHerbivory.HandlingTimeReferenceMass"];
            // todo: make an entry in the parameter file for herbivory
            _FeedingPreferenceStandardDeviation = EcologicalParameters.Parameters["Predation.RevisedPredation.FeedingPreferenceStandardDeviation"];
            _AttackRateActivationEnergy = EcologicalParameters.Parameters["Feeding.TemperatureDependence.AttackRateActivationEnergy"];
            _HandlingTimeActivationEnergy = EcologicalParameters.Parameters["Feeding.TemperatureDependence.HandlingTimeActivationEnergy"];
            _ReferenceTemperature = EcologicalParameters.Parameters["Feeding.TemperatureDependence.ReferenceTemperature"];
            _BoltzmannConstant = EcologicalParameters.Parameters["BoltzmannConstant"];
        }


        /// <summary>
        /// Write out the values of the parameters to an output file
        /// </summary>
        /// <param name="sw">A streamwriter object to write the parameter values to</param>
        public void WriteOutParameterValues(StreamWriter sw)
        {
            // Write out parameters
            sw.WriteLine("Herbivory\tTimeUnitImplementation\t" + Convert.ToString(_TimeUnitImplementation));
            sw.WriteLine("Herbivory\tReferenceMass_g\t" + Convert.ToString(_ReferenceMass));
            sw.WriteLine("Herbivory\tHandlingTimeScalarTerrestrial\t" + Convert.ToString(_HandlingTimeScalarTerrestrial));
            sw.WriteLine("Herbivory\tHandlingTimeScalarMarine\t" + Convert.ToString(_HandlingTimeScalarMarine));
            sw.WriteLine("Herbivory\tHandlingTimeExponentTerrestrial\t" + Convert.ToString(_HandlingTimeExponentTerrestrial));
            sw.WriteLine("Herbivory\tHandlingTimeExponentMarine\t" + Convert.ToString(_HandlingTimeExponentMarine));
            sw.WriteLine("Herbivory\tHerbivoryRateConstantTerrestrial\t" + Convert.ToString(_HerbivoryRateConstantTerrestrial));
            sw.WriteLine("Herbivory\tHerbivoryRateConstantMarine\t" + Convert.ToString(_HerbivoryRateConstantMarine));

            sw.WriteLine("Herbivory\t_AttackRateExponentTerrestrial\t" + Convert.ToString(HerbivoryRateMassExponent));
            sw.WriteLine("Herbivory\t_AttackRateExponentMarine\t" + Convert.ToString(HerbivoryRateMassExponent));
            sw.WriteLine("Herbivory\tHerbivoryRateMassExponent\t" + Convert.ToString(HerbivoryRateMassExponent));
            sw.WriteLine("Herbivory\tAttackRateActivationEnergy\t" + Convert.ToString(_AttackRateActivationEnergy));
            sw.WriteLine("Herbivory\tHandlingTimeActivationEnergy\t" + Convert.ToString(_HandlingTimeActivationEnergy));
            sw.WriteLine("Herbivory\tReferenceTemperature\t" + Convert.ToString(_ReferenceTemperature));
            sw.WriteLine("Herbivory\tBoltzmannConstant\t" + Convert.ToString(_BoltzmannConstant));
        }



        /// <summary>
        /// Calculates the potential biomass of an autotroph stock eaten by a herbivore cohort (terrestrial)
        /// </summary>
        /// <param name="autotrophBiomass">The total biomass of the autotroph stock</param>
        /// <param name="herbivoreIndividualMass">The individual body mass of the acting (herbivore) cohort</param>
        /// <returns>The potential biomass eaten by the herbivore cohort</returns>
        private double CalculatePotentialBiomassEatenTerrestrial(double autotrophBiomass, double herbivoreIndividualMass)
        {
            // Calculate the inidividual herbivory rate per unit autotroph mass-density per hectare
            double IndividualHerbivoryRate = CalculateIndividualHerbivoryRatePerHectareTerrestrial(herbivoreIndividualMass);

            // Calculate autotroph biomass density per hectare
            double AutotrophBiomassDensity = autotrophBiomass / _CellAreaHectares;
            
            // Calculate the expected autotroph biomass eaten
            return IndividualHerbivoryRate * Math.Pow(AutotrophBiomassDensity, _AttackRateExponentTerrestrial);
        }

        /// <summary>
        /// Calculates the potential biomass of an autotroph stock eaten by a herbivore cohort (marine)
        /// </summary>
        /// <param name="autotrophBiomass">The total biomass of the autotroph stock</param>
        /// <param name="herbivoreIndividualMass">The individual body mass of the acting (herbivore) cohort</param>
        /// <returns>The potential biomass eaten by the herbivore cohort</returns>
        private double CalculatePotentialBiomassEatenMarine(double autotrophBiomass, double herbivoreIndividualMass,double preyIndividualMass,
            double logOptimalPreyPredatorMassRatio, string cohortPhytoType)//, double temperature, Tuple<double, double> regCoefs, double[] BinCentres, double[] BiomassSizeStructure
        {



            double Alphaij = CalculateIndividualKillingRatePerHectare(preyIndividualMass, herbivoreIndividualMass, logOptimalPreyPredatorMassRatio);

            // Calculate optimum prey size
            //double optPreyPredatorMassRatio = Math.Exp(logOptimalPreyPredatorMassRatio);
            //double optPreySize = optPreyPredatorMassRatio * herbivoreIndividualMass;

            // Calculate proportion of NPP available to the cohort
/*            double propAvailableNPP = 0.0; 
            propAvailableNPP = regCoefs.Item2 * Math.Log10(optPreySize) + regCoefs.Item1;
            propAvailableNPP = Math.Pow(10, propAvailableNPP);
            
            // Adjust total NPP by the proportion available
            double autotrophBiomassBin = autotrophBiomass * propAvailableNPP;
 */

            

            // Calculate autotroph biomass density per hectare
            double AutotrophBiomassDensity = autotrophBiomass / _CellAreaHectares;//autotrophBiomassBin / _CellAreaHectares;

            // Calculate the expected autotroph biomass eaten
            return Alphaij * Math.Pow(AutotrophBiomassDensity, _AttackRateExponentMarine);
        }


        // Original
        /// <summary>
        /// Calculates the killing rate of an individual predator per unit prey density per hectare per time unit 
        /// </summary>
        /// <param name="preyIndividualMass">The body mass of individuals in the prey cohort</param>
        /// <param name="preyMassBinNumber">The mass bin index of the prey</param>
        /// <param name="preyFunctionalGroup">The functional group index of the prey</param>
        /// <param name="predatorIndividualMass">The body mass of individuals in the predator cohort</param>
        /// <param name="logOptimalPreyPredatorMassRatio">The log ratio of optimal prey body mass to predator body mass</param>
        /// <returns>The killing rate of an individual predator per unit prey density per hectare per time unit</returns>
        private double CalculateIndividualKillingRatePerHectare(double preyIndividualMass,
            double herbivoreIndividualMass, double logOptimalPreyPredatorMassRatio)
        {
            //int PreyBinNumber;

            // Calculate the relative feeding preference from a log-normal distribution with mean equal to the optimal 
            // prey to predator ratio and standard deviation as specified
            RelativeFeedingPreference = Math.Exp(-(Math.Pow
                    (((Math.Log(preyIndividualMass / herbivoreIndividualMass) - logOptimalPreyPredatorMassRatio) /
                    _FeedingPreferenceStandardDeviation), 2)));

            // Calculate the individual killing rate
            return _SpecificHerbivoreKillRateConstant * RelativeFeedingPreference;// *BinnedPreyDensities[preyFunctionalGroup, preyMassBinNumber];
        } 

        /// <summary>
        /// Calculate the herbivory rate of an individual herbivore per unit autotroph mass-density per hectare (terrestrial)
        /// </summary>
        /// <param name="herbivoreIndividualMass">Herbivore individual body mass</param>
        /// <returns>The herbivory rate of an individual herbivore per unit autotroph mass-density per hectare</returns>
        private double CalculateIndividualHerbivoryRatePerHectareTerrestrial(double herbivoreIndividualMass)
        {
            // Calculate the individual herbivory rate
            return _HerbivoryRateConstantTerrestrial * Math.Pow(herbivoreIndividualMass, (_HerbivoryRateMassExponent));

        }

        /// <summary>
        /// Calculate the herbivory rate of an individual herbivore per unit autotroph mass-density per hectare (marine)
        /// </summary>
        /// <param name="herbivoreIndividualMass">Herbivore individual body mass</param>
        /// <returns>The herbivory rate of an individual herbivore per unit autotroph mass-density per hectare</returns>
        private double CalculateIndividualHerbivoryRatePerHectareMarine(double herbivoreIndividualMass)
        {
            return _HerbivoryRateConstantMarine * Math.Pow(herbivoreIndividualMass, (_HerbivoryRateMassExponent));
        }

        /// <summary>
        /// Calculate the time taken for a herbivore in the marine realm to handle unit mass (1 g) of autotroph mass
        /// </summary>
        /// <param name="herbivoreIndividualMass">The body mass of an individual herbivore</param>
        /// <returns>The time taken for a herbivore to handle unit mass (1 g) of autotroph mass</returns>
        private double CalculateHandlingTimeMarine(double herbivoreIndividualMass, double preyIndividualMass)
        {
            return preyIndividualMass*_HandlingTimeScalarMarine * Math.Pow((_ReferenceMass / herbivoreIndividualMass), _HandlingTimeExponentMarine);
        }


        /// <summary>
        /// Calculate the time taken for a herbivore in the terrestrial realm to handle unit mass (1 g) of autotroph mass
        /// </summary>
        /// <param name="herbivoreIndividualMass">The body mass of an individual herbivore</param>
        /// <returns>The time taken for a herbivore to handle unit mass (1 g) of autotroph mass</returns>
        private double CalculateHandlingTimeTerrestrial(double herbivoreIndividualMass)
        {
               return _HandlingTimeScalarTerrestrial * Math.Pow((_ReferenceMass / herbivoreIndividualMass), _HandlingTimeExponentTerrestrial);
        }

        /// <summary>
        /// Calculate the actual biomass eaten by a herbivore cohort from an autotroph stock
        /// </summary>
        /// <param name="potentialBiomassEaten">The potential biomass eaten by the herbivore cohort from the autotroph stock given the encounter rate</param>
        /// <param name="totalHandlingTime">The total time that would be taken to handle all encountered autotroph biomass in all autotroph stocks</param>
        /// <param name="herbivoreAbundance">The number of individuals in the acting herbivore cohort</param>
        /// <param name="autotrophBiomass">The total biomass in the autotroph stock</param>
        /// <returns>The biomass eaten by the herbivore cohort from the autotroph stock</returns>
        private double CalculateBiomassesEaten(double potentialBiomassEaten, double totalHandlingTime, double herbivoreAbundance, double autotrophBiomass)
        {
            // Check whether there is any biomass in the autotroph stock
            if (autotrophBiomass > 0.0)
            {
                // Calculate the instantaneous fraction of the autotroph stock eaten
                InstantFractionEaten = herbivoreAbundance * ((potentialBiomassEaten / (1 + totalHandlingTime)) / autotrophBiomass);
            }
            else
            {
                // Set the instantaneous fraction of the autotroph stock eaten to zero
                InstantFractionEaten = 0.0;
            }
            
            // Return the total  biomass of the autotroph stock eaten
            return autotrophBiomass * (1 - Math.Exp(-InstantFractionEaten * _DeltaT * _ProportionOfTimeEating));
        }

    }
}

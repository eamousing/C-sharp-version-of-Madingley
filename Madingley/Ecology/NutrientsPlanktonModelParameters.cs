using System;
using System.Collections.Generic;
using System.IO;

namespace Madingley
{
    /// <summary>
    /// Class to read in parameters necessary to run the Nutrient-Plankton model
    /// </summary>
    public static class NutrientPlanktonModelParameters
    {
        /// <summary>
        /// Initialize a directory to hold the parameter values
        /// </summary>
        public static Dictionary<string, double> NPZParameters;

        /// <summary>
        /// Method for reading in parameter values and output a copy to the output directory
        /// </summary>
        /// <param name="parametersFile">Name of the file containing the parameters</param>
        /// <param name="outputPath">Path of the output directory</param>
        public static void ReadNPZModelParameters(string parametersFile, string outputPath)
        {
            // Copy the parameter values to the output directory
            System.IO.File.Copy("input/Model setup/Ecological definition files/" + parametersFile, outputPath + parametersFile, true);

            // Read the parameter values into a directory
            NPZParameters = new Dictionary<string, double>();
            StreamReader r_npz = new StreamReader("input/Model setup/Ecological definition files/" + parametersFile);
            string l;
            char[] comma = ",".ToCharArray();

            string[] f;

            l = r_npz.ReadLine();
            while(!r_npz.EndOfStream)
            {
                l = r_npz.ReadLine();
                // Split field by commas
                f = l.Split(comma);
                // Lists of the fields
                NPZParameters.Add(f[0], Convert.ToDouble(f[1]));
            }
        }
    }
}

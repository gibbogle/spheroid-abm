#include <qstring.h>
#include "params.h"

Params::Params()
{
    PARAM_SET params[] = {

/*
{"BC_AVIDITY_MEDIAN", 1.0, 0.1, 10.0,
"BCR avidity median parameter",
"BCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
(BCR stimulation rate is proportional to the product of BC avidity and antigen load.)"},

{"BC_AVIDITY_SHAPE", 1.1, 1.01, 3.0,
"BCR avidity shape parameter",
"BCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
The shape value must be greater than 1, and values close to 1 give distributions that are close to normal."},
*/

{"NX", 100, 100, 300,
"Lattice size",
"Dimension of the lattice (number of sites in X,Y and Z directions).  Typically 5*BLOB_RADIUS is OK."},

{"INITIAL_COUNT", 3000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"DIVIDE_TIME_MEDIAN", 18, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"NDAYS", 5.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[days]"},

{"NT_CONC", 10, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver."},

{"NMM3", 120000, 0, 0,
"Cells/cubic mm",
"Number of cells per cubic mm of non-necrotic tumour."},

{"FLUID_FRACTION", 0.5, 0, 0,
"Fluid fraction",
"Fraction of non-necrotic tumour that is extracellular fluid."},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.05, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"DEATH_THRESHOLD", 0.0001, 0, 0,
"O2 death threshold",
"Concentration of O2 at which a cell begins to experience anoxia leading to cell death.\n\
[mM]"},

{"THRESHOLD_FACTOR", 1.5, 0, 0,
"O2 death threshold factor",
"Multiplying factor for death threshold of O2 concentration."},

{"T_HYPOXIC_LIMIT", 3.0, 0, 0,
"Hypoxic time limit",
"Length of time under hypoxia (O2 < death threshold) after which a cell dies.\n\
[h]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 2, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "OXYGEN diffusion coeff",
 "OXYGEN diffusion coefficient"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "OXYGEN boundary concentration",
 "OXYGEN boundary concentration"},

{"OXYGEN_CONSUMPTION", 2.3e-16, 0, 0,
 "OXYGEN consumption rate",
 "OXYGEN consumption rate"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_DIFF_COEF", 6.0e-7, 0, 0,
 "GLUCOSE diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_BDRY_CONC", 9.0, 0, 0,
 "GLUCOSE boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSUMPTION", 3.8e-17, 0, 0,
 "GLUCOSE consumption rate",
 "GLUCOSE consumption rate"},

{"USE_DRUG_A", 0, 0, 1,
"Use Drug A?",
"Drug A is simulated"},

{"DRUG_A_DIFF_COEF", 6.0e-7, 0, 0,
 "DRUG_A diffusion coeff",
 "DRUG_A diffusion coefficient"},

{"DRUG_A_BDRY_CONC", 0.01, 0, 0,
 "DRUG_A boundary concentration",
 "DRUG_A boundary concentration"},

{"DRUG_A_CONSUMPTION", 0, 0, 0,
 "DRUG_A consumption rate",
 "DRUG_A consumption rate"},

{"USE_DRUG_B", 0, 0, 1,
"Use Drug B?",
"Drug B is simulated"},

{"DRUG_B_DIFF_COEF", 6.0e-7, 0, 0,
 "DRUG_B diffusion coeff",
 "DRUG_B diffusion coefficient"},

{"DRUG_B_BDRY_CONC", 9.0, 0, 0,
 "DRUG_B boundary concentration",
 "DRUG_B boundary concentration"},

{"DRUG_B_CONSUMPTION", 0, 0, 0,
 "DRUG_B consumption rate",
 "DRUG_B consumption rate"},

{"DRUG_A_DECAY", 0, 0, 0,
 "DRUG_A decay",
 "DRUG_A boundary conc decays with the specified half-life"},

{"DRUG_A_HALFLIFE", 2.0, 0, 0,
 "DRUG_A half-life",
 "DRUG_A half-life"},

{"DRUG_B_DECAY", 0, 0, 0,
 "DRUG_B decay",
 "DRUG_B boundary conc decays with the specified half-life"},

{"DRUG_B_HALFLIFE", 2.0, 0, 0,
 "DRUG_B half-life",
 "DRUG_B half-life"},

{"SN30K_KMET0", 1.54, 0, 0,
 "Kmet0",
 "SN30000 max value of 1st order rate constant for metabolism under zero oxygen"},

{"SN30K_C1", 0.0, 0, 0,
 "C1",
 "SN30000 C1 in function for oxygen-dependence of rate of metabolism"},

{"SN30K_C2", 1.0, 0, 0,
 "C2",
 "SN30000 C2 in function for oxygen-dependence of rate of metabolism"},

{"SN30K_KO2", 1.14, 0, 0,
 "KO2",
 "SN30000 KO2 in function for oxygen-dependence of rate of metabolism"},

{"SN30K_GAMMA", 1.0, 0, 0,
 "gamma",
 "SN30000 gamma"},

{"SN30K_KLESION", 0.001, 0, 0,
 "Klesion",
 "SN30000 Klesion is the parameter that converts total metabolite into lesion level"},

{"KILL_O2_CONC", 0.0, 0, 0,
 "O2 conc",
 "SN30000 constant O2 concentration in kill experiment"},

{"KILL_DRUG_CONC", 0.01, 0, 0,
 "Drug conc",
 "SN30000 constant drug concentration in kill experiment"},

{"KILL_DURATION", 60, 0, 0,
 "Duration",
 "SN30000 duration of kill experiment"},

{"KILL_FRACTION", 0.9, 0, 0,
 "Kill fraction",
 "SN30000 fraction of cells killed in the experiment"},

{"INPUT_FILE", 0, 0, 0,
"spheroid_fixed.inpdata",
"The auxiliary input file contains data that (almost!) never changes"}

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}

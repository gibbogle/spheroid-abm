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

{"FLUID_FRACTION", 0.25, 0, 0,
"Fluid fraction",
"Fraction of non-necrotic tumour that is extracellular fluid."},

{"MEDIUM_VOLUME", 1.0, 0, 0,
"Medium volume",
"Volume of the medium in which the spheroid is growing."},

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

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Diffusion coeff",
 "Constituent diffusion coefficient"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSUMPTION", 2.3e-16, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

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

{"USE_SN30K", 0, 0, 1,
"Use SN30000?",
"SN30000 is simulated"},

{"SN30K_BDRY_CONC", 0.01, 0, 0,
 "SN30000 boundary concentration",
 "SN30000 concentration in the medium"},

{"SN30K_DECAY", 0, 0, 0,
 "SN30000 decay",
 "SN30000 boundary conc decays with the specified half-life"},

{"SN30K_HALFLIFE", 2.0, 0, 0,
 "SN30000 half-life",
 "SN30000 half-life"},

{"SN30K_METABOLITE", 0, 0, 0,
 "SN30000 metabolite",
 "SN30000 simulate metabolite"},

{"SN30K_DIFF_COEF", 6.0e-7, 0, 0,
 "SN30000 diffusion coeff",
 "SN30000 diffusion coefficient"},

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

{"SN30K_KILL_O2_CONC", 0.0, 0, 0,
 "O2 conc",
 "SN30000 constant O2 concentration in kill experiment"},

{"SN30K_KILL_DRUG_CONC", 0.01, 0, 0,
 "Drug conc",
 "SN30000 constant drug concentration in kill experiment"},

{"SN30K_KILL_DURATION", 60, 0, 0,
 "Duration",
 "SN30000 duration of kill experiment"},

{"SN30K_KILL_FRACTION", 0.9, 0, 0,
 "Kill fraction",
 "SN30000 fraction of cells killed in the experiment"},

{"USE_DRUG_B", 0, 0, 1,
"Use Drug B?",
"Drug B is simulated"},

{"DRUG_B_BDRY_CONC", 9.0, 0, 0,
 "DRUG_B boundary concentration",
 "DRUG_B boundary concentration"},

{"DRUG_B_DECAY", 0, 0, 0,
 "DRUG_B decay",
 "DRUG_B boundary conc decays with the specified half-life"},

{"DRUG_B_HALFLIFE", 2.0, 0, 0,
 "DRUG_B half-life",
 "DRUG_B half-life"},

{"DRUG_B_METABOLITE", 0, 0, 0,
 "DRUG_B metabolite",
 "DRUG_B simulate metabolite"},

{"USE_TREATMENT_FILE", 0, 0, 1,
"Use treatment file?",
"Treatment programme is specified in the treatment file"},

{"TREATMENT_FILE", 0, 0, 0,
"treatment.data",
"The treatment file contains data describing the drug and radiation dosing schedule"}

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

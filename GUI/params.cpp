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

{"INITIAL_COUNT", 15000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"DIVIDE_TIME_MEDIAN", 24, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 1, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 7.0, 0.0, 30.0,
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

{"NMM3", 500000, 0, 0,
"Cells/cubic mm",
"Number of cells per cubic mm of non-necrotic tumour."},

{"FLUID_FRACTION", 0.5, 0, 0,
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

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_FACTOR", 1.5, 0, 0,
"Anoxia threshold factor",
"A cell begins to experience anoxia leading to cell death at the O2 concentration given by this multiplying factor times the Michaelis-Menten threshold value."},

{"ANOXIA_TAG_TIME", 3.0, 0, 0,
"Anoxic time limit",
"Length of time under hypoxia (O2 < anoxic threshold) after which a cell is tagged to die of anoxia.\n\
[h]"},

{"ANOXIA_DEATH_TIME", 3.0, 0, 0,
"Anoxic death delay time",
"Time taken for a cell to die after it is tagged to die of anoxia.\n\
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

{"OXYGEN_CELL_DIFF", 200, 0, 0,
 "Membrane diff coef",
 "Cell membrane diffusion coefficient Kd"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSUMPTION", 2.3e-16, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_DIFF_COEF", 6.0e-7, 0, 0,
 "GLUCOSE diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_CELL_DIFF", 20, 0, 0,
 "Membrane diff coef",
 "Cell membrane diffusion coefficient Kd"},

{"GLUCOSE_BDRY_CONC", 9.0, 0, 0,
 "GLUCOSE boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSUMPTION", 3.8e-17, 0, 0,
 "GLUCOSE consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"USE_DRUG_A", 0, 0, 1,
"Use SN30000?",
"SN30000 is simulated"},

{"DRUG_A_NAME", 0, 0, 0,
"SN30000",
"Name of drug A"},

{"DRUG_A_BDRY_CONC", 0.0, 0, 0,
 "Drug A boundary concentration",
 "Drug A concentration in the medium"},

{"DRUG_A_DECAY", 0, 0, 0,
 "Drug A decay",
 "Drug A conc decays with the specified half-life"},

{"DRUG_A_SIMULATE_METABOLITE", 0, 0, 0,
 "Drug A metabolite",
 "Drug A simulate metabolite"},

{"DRUG_A_METABOLITE_DECAY", 0, 0, 0,
 "Drug A metabolite decay",
 "Drug A metabolite conc decays with the specified half-life"},

{"SN30K_DIFF_COEF", 6.0e-7, 0, 0,
 "Diffusion coeff",
 "SN30000 diffusion coefficient (also used for metabolite)"},

{"SN30K_CELL_DIFF", 5, 0, 0,
 "Membrane diff coef",
 "Cell membrane diffusion coefficient Kd (also used for metabolite)"},

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

{"SN30K_HALFLIFE", 2.0, 0, 0,
 "Half-life",
 "SN30000 half-life (hours)"},

{"SN30K_METABOLITE_HALFLIFE", 2.0, 0, 0,
 "Metabolite half-life",
 "SN30000 metabolite half-life (hours)"},

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

{"DRUG_B_NAME", 0, 0, 0,
"",
"Name of drug B"},

{"DRUG_B_BDRY_CONC", 0.0, 0, 0,
 "DRUG_B boundary concentration",
 "DRUG_B boundary concentration"},

{"DRUG_B_DECAY", 0, 0, 0,
 "DRUG_B decay",
 "DRUG_B boundary conc decays with the specified half-life"},

{"DRUG_B_SIMULATE_METABOLITE", 0, 0, 0,
 "DRUG_B metabolite",
 "DRUG_B simulate metabolite"},

{"DRUG_B_METABOLITE_DECAY", 0, 0, 0,
 "Drug B metabolite decay",
 "Drug B metabolite conc decays with the specified half-life"},

{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"SPCRAD", 500.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

{"USE_EXTRA", 0, 0, 1,
"Use extra conc",
"Use extracellular O2 concentration to determine cell death"},

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

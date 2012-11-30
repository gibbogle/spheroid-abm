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

{"BC_COGNATE_FRACTION", 0.0001, 0, 0,
"B cell cognate fraction",
"The fraction of B cells that are cognate, i.e. recognize and respond to the antigen."},

{"BC_STIM_RATE_CONSTANT", 1, 0.0, 100.0,
"BCR stimulation rate constant",
"Rate constant Ks for BCR stimulation, where:\n\
rate of BCR stimulation = Ks*(BCR avidity)*(antigen load)\n\
[molecules/min]"},
	
{"BC_STIM_HALFLIFE", 24.0, 0.0, 100.0,
"BCR stimulation halflife",
"Integrated BCR stimulation decays with a specified halflife. \n\
[hours]"},

{"DIVIDE1_MEDIAN", 6.0, 0.0, 100.0,
"1st division time median parameter",
"The time taken for the first B cell division, after full activation, has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE1_SHAPE", 1.1, 1.01, 3.0,
"1st division time shape parameter",
"The time taken for the first B cell division, after full activation, has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE2_MEDIAN", 5.0, 0, 100.0,
"Later division time median parameter",
"The time taken for later B cell divisions has a lognormal distribution, described by the median and shape parameters.\n\
[hours]"},

{"DIVIDE2_SHAPE", 1.1, 1.01, 3.0,
"Later division time shape parameter",
"The time taken for later B cell divisions has a lognormal distribution, described by the median and shape parameters."},

{"MOTILITY_BETA", 0.06, 0.0, 0.0,
"Motility speed parameter",
"B cell motility is described by speed and persistence parameters, each in the range 0 - 1. Median B cell speed is roughly proportional to MOTILITY_BETA."},

{"MOTILITY_RHO", 0.76, 0.0, 0.0,
"Motility persistence parameter",
"B cell motility is described by speed and persistence parameters, each in the range 0 - 1. MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next."},

{"GCCPROB_GENLT3",0.0,0,0,
"Probability for generation < 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation < 3."},

{"GCCPROB_GENEQ3",0.1,0,0,
"Probability for generation = 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation = 3."},

{"GCCPROB_GENGT3",0.3,0,0,
"Probability for generation > 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation > 3."},

{"PLASMAPROB_GENLT3",0,0,0,
"Probability for generation < 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation < 3."},

{"PLASMAPROB_GENEQ3",0,0,0,
"Probability for generation = 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation = 3."},

{"PLASMAPROB_GENGT3",0.4,0,0,
"Probability for generation > 3",
"Probabilities of differentiation into GC cells or Plasma cells on cell division yielding generation > 3."},

*/
{"NX", 100, 100, 300,
"Lattice size",
"Dimension of the lattice (number of sites in X,Y and Z directions).  Typically 5*BLOB_RADIUS is OK."},

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

/*
{"BLOB_RADIUS", 10.0, 5.0, 50.0,
"Initial blob size",
"The major radius of the initial ellipsoidal blob of B cells, as number of sites.  (18.5 -> 12k sites)"},

{"BC_FRACTION", 0.6, 0.4, 0.8,
"B cell fraction",
"Fraction of the follicular volume occupied by B cells."},

{"FLUID_FRACTION", 0.1, 0.05, 0.2,
"Fluid fraction",
"Fraction of the paracortical volume occupied by fluid."},

{"USE_TRAFFIC", 1, 0, 1,
"B cell trafficking?",
"B cell trafficking is simulated (ingress and egress)"},

{"USE_TCELLS", 0, 0, 1,
"Simulate T cells?",
"CD4 T cells are explicitly simulated (i.e. they occupy space))"},

{"USE_CROWDING_CORRECTION", 0, 0, 1,
"Apply crowding correction?",
"Apply crowding correction based on gradient of cell density."},

{"USE_EXIT_CHEMOTAXIS", 0, 0, 1,
"B cell exit chemotaxis?",
"S1P-modulated B cell chemotaxis towards exit portals is simulated"},

{"COMPUTED_OUTFLOW", 0, 0, 1,
"Compute B cell outflow limit?",
"The upper bound on B cell outflow is computed together with inflow.  The alternative is to permit (probabilistic) egress of any cell at a portal."},

{"RESIDENCE_TIME", 15.0, 10.0, 30.0,
"B cell residence time",
"B cell residence time.\n\
[hours]"},

{"INFLAMM_DAYS1", 3.5, 0, 0,
"Inflammation plateau duration",
"Period over which the level of inflammation signal from the periphery is constant.\n\
[days]"},

{"INFLAMM_DAYS2", 4.5, 0, 0,
"Inflammation cessation time",
"Time at which the level of inflammation signal from the periphery goes to zero.\n\
[days]"},

{"INFLAMM_LEVEL", 0.0, 0.0, 1.0,
"Inflammation level",
"The plateau inflammation signal level."},

{"USE_S1PR1", 1, 0, 1,
"Use S1PR1?",
"S1P chemotactic attraction is simulated"},

{"USE_S1PR2", 1, 0, 1,
"Use S1PR2?",
"S1P chemotactic repulsion is simulated"},

{"S1P_BDRY_0",0, 0, 1,
"Use S1P bdry secretion?",
"Use S1P bdry secretion rate (otherwise use bdry concentration)"},

{"S1P_BDRY_RATE", 1, 0, 0,
 "S1P boundary rate",
 "Influx rate of S1P at the follicle surface"},

{"S1P_BDRY_CONC", 100, 0, 0,
 "S1P boundary conc",
 "Concentration of S1P at the follicle surface"},

{"S1P_DIFF_COEFF", 1.0, 0, 0,
 "S1P diffusion coeff",
 "S1P diffusion coefficient"},

{"S1P_HALFLIFE", 10, 0, 0,
 "S1P halflife",
 "S1P halflife (hours)"},

{"S1P_STRENGTH_POS", 0.3, 0, 0,
 "S1P pos strength",
 "Relative strength of S1P chemotactic attraction (S1PR1)"},

{"S1P_STRENGTH_NEG", 1, 0, 0,
 "S1P neg strength",
 "Relative strength of S1P chemotactic repulsion (S1PR2)"},

{"S1PR1_1", 1.0, 0, 0},
{"S1PR1_2", 0.2, 0, 0},
{"S1PR1_3", 0.2, 0, 0},
{"S1PR1_4", 0.0, 0, 0},
{"S1PR1_5", 0.2, 0, 0},

{"S1PR1_sat_threshold", 0.9, 0, 0},
{"S1PR1_refractory_time", 5, 0, 0},

{"S1PR2_1", 0.0, 0, 0},
{"S1PR2_2", 0.0, 0, 0},
{"S1PR2_3", 0.0, 0, 0},
{"S1PR2_4", 2.0, 0, 0},
{"S1PR2_5", 0.0, 0, 0},

{"S1PR2_sat_threshold", 0.9, 0, 0},
{"S1PR2_refractory_time", 5, 0, 0},

{"USE_CCR7", 1, 0, 1,
"Use CCR7?",
"CCL21 chemotaxis is simulated"},

{"CCL21_BDRY_0", 0, 0, 1,
"Use CCL21 bdry secretion?",
"Use CCL21 bdry secretion rate (otherwise use bdry concentration)"},

{"CCL21_BDRY_RATE", 1, 0, 0,
 "CCL21 boundary rate",
 "Influx rate of CCL21 at the follicle surface"},

{"CCL21_BDRY_CONC", 100, 0, 0,
 "CCL21 boundary conc",
 "Concentration of CCL21 at the follicle surface"},

{"CCL21_DIFF_COEFF", 1.0, 0, 0,
 "CCL21 diffusion coeff",
 "CCL21 diffusion coefficient"},

{"CCL21_HALFLIFE", 10, 0, 0,
 "CCL21 halflife",
 "CCL21 halflife (hours)"},

{"CCL21_STRENGTH", 0.2, 0, 0,
 "CCL21 strength",
 "Relative strength of CCL21 chemotactic influence"},

{"CCR7_1", 1.0, 0, 0},
{"CCR7_2", 2.0, 0, 0},
{"CCR7_3", 0.5, 0, 0},
{"CCR7_4", 0.5, 0, 0},
{"CCR7_5", 0.5, 0, 0},

{"CCR7_sat_threshold", 0.9, 0, 0},
{"CCR7_refractory_time", 5, 0, 0},

{"USE_EBI2", 1, 0, 1,
"Use EBI2?",
"Oxysterol chemotaxis is simulated"},

{"OXY_BDRY_0", 0, 0, 1,
"Use OXY bdry secretion?",
"Use OXY bdry secretion rate"},

{"OXY_BDRY_RATE", 1, 0, 0,
 "OXY boundary rate",
 "Influx rate of Oxysterol at the follicle surface"},

{"OXY_BDRY_CONC", 100, 0, 0,
 "OXY boundary conc",
 "Concentration of Oxysterol at the follicle surface"},

{"OXY_DIFF_COEFF", 1.0, 0, 0,
 "OXY diffusion coeff",
 "Oxysterol diffusion coefficient"},

{"OXY_HALFLIFE", 5, 0, 0,
 "OXY halflife",
 "Oxysterol halflife (hours)"},

{"OXY_STRENGTH", 0.1, 0, 0,
 "OXY strength",
 "Relative strength of Oxysterol chemotactic influence"},

{"EBI2_1", 1.0, 0, 0},
{"EBI2_2", 2.0, 0, 0},
{"EBI2_3", 2.0, 0, 0},
{"EBI2_4", 0.2, 0, 0},
{"EBI2_5", 1.0, 0, 0},

{"EBI2_sat_threshold", 0.9, 0, 0},
{"EBI2_refractory_time", 5, 0, 0},

{"USE_CXCR5", 1, 0, 1,
"Use CXCR5?",
"CXCL13 chemotaxis is simulated"},

{"CXCL13_BDRY_0", 0, 0, 1,
"Use FDC CXCL13 secretion?",
"Use FDC CXCL13 secretion rate"},

{"CXCL13_BDRY_RATE", 1, 0, 0,
 "CXCL13 boundary rate",
 "Secretion rate of CXCL13 at the FDC surface"},

{"CXCL13_BDRY_CONC", 100, 0, 0,
 "CXCL13 boundary conc",
 "Concentration of CXCL13 at the FDC surface"},

{"CXCL13_DIFF_COEFF", 1.0, 0, 0,
 "CXCL13 diffusion coeff",
 "CXCL13 diffusion coefficient"},

{"CXCL13_HALFLIFE", 10, 0, 0,
 "CXCL13 halflife",
 "CXCL13 halflife (hours)"},

{"CXCL13_STRENGTH", 0.2, 0, 0,
 "CXCL13 strength",
 "Relative strength of CXCL13 chemotactic influence"},

{"CXCR5_1", 1.0, 0, 0},
{"CXCR5_2", 1.0, 0, 0},
{"CXCR5_3", 1.0, 0, 0},
{"CXCR5_4", 1.5, 0, 0},
{"CXCR5_5", 0.0, 0, 0},

{"CXCR5_sat_threshold", 0.9, 0, 0},
{"CXCR5_refractory_time", 5, 0, 0},

{"NFDC", 100, 0, 0,
"Number of FDCs",
"Initial number of FDCs in the follicle.  Each FDC occupies 7 lattice sites."},

{"NMRC", 100, 0, 0,
"Number of MRCs",
"Initial number of MRCs in the follicle.  Each MRC occupies 7 lattice sites."},

{"BASE_EXIT_PROB", 0.007, 0.0, 0.0,
"Base B cell exit probability",
"A cell located at a boundary site on the lower exit surface has a probability of egress in a time step.\n\
The specified probability applies to residence time Tres = 15 hr, and will be scaled appropriately for a different Tres value."},

*/

{"DIVIDE_TIME_MEDIAN", 18, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"NDAYS", 1.0, 0.0, 30.0,
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

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

/*
{"SHOW_NONCOGNATE", 0, 0, 1,
"Display non-cognate cells?",
"Display a representative fraction of the non-cognate B cells"},

{"DISPLAY_FRACTION", 0, 0, 1,
"Fraction to display",
"Fraction of non-cognate B cells to display."},
*/

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

    {"DRUG_A_BDRY_CONC", 9.0, 0, 0,
     "DRUG_A boundary concentration",
     "DRUG_A boundary concentration"},

    {"DRUG_A_CONSUMPTION", 3.8e-17, 0, 0,
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

    {"DRUG_B_CONSUMPTION", 3.8e-17, 0, 0,
     "DRUG_B consumption rate",
     "DRUG_B consumption rate"},

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

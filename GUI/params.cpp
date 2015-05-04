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

{"GUI_VERSION_NAME", 0, 0, 0,
 "GUI0.00",
 "GUI program version number."},

{"DLL_VERSION_NAME", 0, 0, 0,
 "DLL0.00",
 "DLL version number."},

{"NX", 150, 0, 0,
"Lattice size",
"Dimension of the lattice (number of sites in X,Y and Z directions).  Typically 5*BLOB_RADIUS is OK."},

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"DIVIDE_TIME_MEDIAN", 24, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 0, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 10.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[mins]"},

{"NT_CONC", 1, 0, 0,
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

{"UNSTIRRED_LAYER", 0.001, 0, 0,
"Unstirred layer width",
"Thickness of the unstirred layer around the spheroid."},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.1, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_THRESHOLD", 0.15, 0, 0,
"Anoxia threshold",
"A cell begins to experience anoxia leading to cell death at the O2 concentration given by this threshold value."},

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

{"NCELLTYPES", 2, 0, 0,
"Number of cell types",
"Maximum number of cell types in the spheroid.  The initial percentage of each type must be specified"},

{"CELLPERCENT_1", 100, 0, 100,
"Percentage of cell type 1",
"Percentage of cell type 1"},

//{"CELLDISPLAY_1", 1, 0, 1,
//"Display cell type 1",
//"Display cell type 1"},

{"CELLPERCENT_2", 0, 0, 100,
"Percentage of cell type 2",
"Percentage of cell type 2"},

//{"CELLDISPLAY_2", 1, 0, 1,
//"Display cell type 2",
//"Display cell type 2"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SHOW_PROGENY", 0, 0, 0,
 "Show descendents of cell #",
 "All the descendents of cell with the specified ID are highlighted.  (0 = no selection)"},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Spheroid diffusion coeff",
 "Constituent diffusion coefficient in the spheroid"},

{"OXYGEN_MEDIUM_DIFF", 2.5e-5, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"OXYGEN_CELL_DIFF", 600, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion constant Kd"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"OXYGEN_CONSUMPTION", 6.25e-17, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"OXYGEN_HILL_N", 1, 1, 2,
 "Hill function N",
 "Oxygen uptake rate Hill function N"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUCOSE_CELL_DIFF", 1200, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kd"},

{"GLUCOSE_BDRY_CONC", 9.0, 0, 0,
 "Boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUCOSE_CONSUMPTION", 3.8e-17, 0, 0,
 "Max consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUCOSE_HILL_N", 1, 1, 2,
 "Hill function N",
 "Glucose uptake rate Hill function N"},

{"USE_TRACER", 0, 0, 1,
"Use Tracer?",
"Tracer is simulated"},

{"TRACER_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TRACER diffusion coefficient"},

{"TRACER_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"TRACER_CELL_DIFF", 20, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kd"},

{"TRACER_BDRY_CONC", 1.0, 0, 0,
 "Boundary concentration",
 "TRACER boundary concentration"},

{"TRACER_CONSTANT", 1, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"TRACER_CONSUMPTION", 0, 0, 0,
 "Consumption rate",
 "TRACER consumption rate"},

{"TRACER_MM_KM", 0, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"TRACER_HILL_N", 0, 0, 2,
 "Hill function N",
 "Tracer uptake rate Hill function N"},

//==========================
// TPZ-type drug parameters
//==========================

{"USE_TPZ_DRUG", 0, 0, 1,
"Use TPZ-type?",
"TPZ-type drug is simulated"},

{"TPZ_DRUG_NAME", 0, 0, 0,
"SN30000",
"Name of TPZ-type drug"},

{"TPZ_DRUG_BDRY_CONC", 0.0, 0, 0,
 "Boundary concentration",
 "Drug concentration in the medium"},

{"TPZ_DRUG_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"TPZ_DRUG_SIMULATE_METABOLITE", 0, 0, 0,
 "Simulate metabolites",
 "Simulate drug metabolites"},

{"TPZ_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TPZ-type drug diffusion coefficient in the spheroid"},

{"TPZ_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "TPZ-type drug diffusion coefficient in the medium"},

{"TPZ_CELL_DIFF_IN", 5, 0, 0,
 "Membrane diff (in)",
 "TPZ-type drug cell membrane inwards permeability constant Kd"},

{"TPZ_CELL_DIFF_OUT", 5, 0, 0,
 "Membrane diff (out)",
 "TPZ-type drug cell membrane outwards permeability constant Kd"},

{"TPZ_HALFLIFE", 2.0, 0, 0,
 "Half-life",
 "TPZ-type drug half-life (hours)"},

{"TPZ_DIFF_COEF_MET1", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TPZ-type drug metabolite #1 diffusion coefficient in the spheroid"},

{"TPZ_MEDIUM_DIFF_MET1", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "TPZ-type drug metabolite #1 diffusion coefficient in the medium"},

{"TPZ_CELL_DIFF_IN_MET1", 5, 0, 0,
 "Membrane diff (in)",
 "TPZ-type drug metabolite #1 cell membrane inwards permeability constant Kd"},

{"TPZ_CELL_DIFF_OUT_MET1", 5, 0, 0,
 "Membrane diff (out)",
 "TPZ-type drug metabolite #1 cell membrane outwards permeability constant Kd"},

{"TPZ_HALFLIFE_MET1", 2.0, 0, 0,
 "Half-life",
 "TPZ-type drug metabolite #1 half-life (hours)"},

{"TPZ_DIFF_COEF_MET2", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TPZ-type drug metabolite #2 diffusion coefficient in the spheroid"},

{"TPZ_MEDIUM_DIFF_MET2", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "TPZ-type drug metabolite #2 diffusion coefficient in the medium"},

{"TPZ_CELL_DIFF_IN_MET2", 5, 0, 0,
 "Membrane diff (in)",
 "TPZ-type drug metabolite #2 cell membrane inwards permeability constant Kd"},

{"TPZ_CELL_DIFF_OUT_MET2", 5, 0, 0,
 "Membrane diff (out)",
 "TPZ-type drug metabolite #2 cell membrane outwards permeability constant Kd"},

{"TPZ_HALFLIFE_MET2", 2.0, 0, 0,
 "Half-life",
 "TPZ-type drug metabolite #2 half-life (hours)"},

//------------------------------------------------------
// Cell type 1 metabolism and kill experiment parameters
//------------------------------------------------------

// Parent
{"TPZ_KMET0_CELL1", 1.54, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL1", 1.0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL1", 1.14, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL1", 0, 0, 0,
 "Vmax",
    "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism: Kmet0 -> Kmet0 + Vmax/(Km + C)"},

{"TPZ_KM_CELL1", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL1", 0.001, 0, 0,
 "Klesion    ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 1 TPZ drug kill experiment parameters

    {"TPZ_KILL_MODEL_CELL1", 1, 0, 0,
     "Kill model",
     "Model of TPZ-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"TPZ_KILL_O2_CONC_CELL1", 0.0, 0, 0,
     "O2 conc",
     "TPZ-type drug constant O2 concentration in kill experiment"},

    {"TPZ_KILL_DRUG_CONC_CELL1", 0.01, 0, 0,
     "Drug conc",
     "TPZ-type drug constant drug concentration in kill experiment"},

    {"TPZ_KILL_DURATION_CELL1", 60, 0, 0,
     "Duration",
     "TPZ-type drug duration of kill experiment"},

    {"TPZ_KILL_FRACTION_CELL1", 0.9, 0, 0,
     "Kill fraction",
     "TPZ-type drug fraction of cells killed in the experiment"},

// Metabolite 1
{"TPZ_KMET0_CELL1_MET1", 0.5, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL1_MET1", 0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL1_MET1", 1, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL1_MET1", 0, 0, 0,
 "Vmax",
 "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KM_CELL1_MET1", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL1_MET1", 0.001, 0, 0,
 "Klesion  ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},

// Metabolite 2
{"TPZ_KMET0_CELL1_MET2", 0.5, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL1_MET2", 0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL1_MET2", 1, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL1_MET2", 0, 0, 0,
 "Vmax",
 "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KM_CELL1_MET2", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL1_MET2", 0.001, 0, 0,
 "Klesion  ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},

//------------------------------------------------------
// Cell type 2 metabolism and kill experiment parameters
//------------------------------------------------------

// Parent
{"TPZ_KMET0_CELL2", 1.54, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL2", 1.0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL2", 1.14, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL2", 0, 0, 0,
 "Vmax",
 "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KM_CELL2", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL2", 0.001, 0, 0,
 "Klesion    ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 2 TPZ drug kill experiment parameters

    {"TPZ_KILL_MODEL_CELL2", 1, 0, 0,
     "Kill model",
     "Model of TPZ-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"TPZ_KILL_O2_CONC_CELL2", 0.0, 0, 0,
     "O2 conc",
     "TPZ-type drug constant O2 concentration in kill experiment"},

    {"TPZ_KILL_DRUG_CONC_CELL2", 0.01, 0, 0,
     "Drug conc",
     "TPZ-type drug constant drug concentration in kill experiment"},

    {"TPZ_KILL_DURATION_CELL2", 60, 0, 0,
     "Duration",
     "TPZ-type drug duration of kill experiment"},

    {"TPZ_KILL_FRACTION_CELL2", 0.9, 0, 0,
     "Kill fraction",
     "TPZ-type drug fraction of cells killed in the experiment"},

// Metabolite 1
{"TPZ_KMET0_CELL2_MET1", 0.5, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL2_MET1", 0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL2_MET1", 1, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL2_MET1", 0, 0, 0,
 "Vmax",
 "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KM_CELL2_MET1", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL2_MET1", 0.001, 0, 0,
 "Klesion  ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},

// Metabolite 2
{"TPZ_KMET0_CELL2_MET2", 0.5, 0, 0,
 "Kmet0",
 "TPZ-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"TPZ_C2_CELL2_MET2", 0, 0, 0,
 "C2",
 "TPZ-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KO2_CELL2_MET2", 1, 0, 0,
 "KO2",
 "TPZ-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"TPZ_VMAX_CELL2_MET2", 0, 0, 0,
 "Vmax",
 "TPZ-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KM_CELL2_MET2", 1, 0, 0,
 "Km",
 "TPZ-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"TPZ_KLESION_CELL2_MET2", 0.001, 0, 0,
 "Klesion  ",
 "TPZ-type drug Klesion is the parameter that converts total metabolite into lesion level"},


//==========================
// DNB-type drug parameters
//==========================

{"USE_DNB_DRUG", 0, 0, 1,
"Use DNB-type drug?",
"DNB-type drug is simulated"},

{"DNB_DRUG_NAME", 0, 0, 0,
"PR104A",
"Name of DNB-type drug"},

{"DNB_DRUG_BDRY_CONC", 0.0, 0, 0,
 "Drug boundary concentration",
 "Drug concentration in the medium"},

{"DNB_DRUG_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"DNB_DRUG_SIMULATE_METABOLITE", 0, 0, 0,
 "Simulate metabolites",
 "Simulate drug metabolites"},

{"DNB_DIFF_COEF", 4.42e-7, 0, 0,
 "Spheroid diffusion coeff",
 "DNB-type drug diffusion coefficient in the spheroid"},

{"DNB_MEDIUM_DIFF", 5.0e-6, 0, 0,
 "Medium diffusion coeff",
 "DNB-type drug diffusion coefficient in the medium"},

{"DNB_CELL_DIFF_IN", 2, 0, 0,
 "Membrane diff (in)",
 "DNB-type drug cell membrane inwards permeability constant Kd"},

{"DNB_CELL_DIFF_OUT", 2, 0, 0,
 "Membrane diff (out)",
 "DNB-type drug cell membrane outwards permeability constant Kd"},

{"DNB_HALFLIFE", 10.0, 0, 0,
 "Half-life",
 "DNB-type drug half-life (hours)"},

{"DNB_DIFF_COEF_MET1", 6.24e-7, 0, 0,
 "Spheroid diffusion coeff",
 "DNB-type drug metabolite #1 diffusion coefficient in the spheroid"},

{"DNB_MEDIUM_DIFF_MET1", 5.0e-6, 0, 0,
 "Medium diffusion coeff",
 "DNB-type drug metabolite #1 diffusion coefficient in the medium"},

{"DNB_CELL_DIFF_IN_MET1", 1, 0, 0,
 "Membrane diff (in)",
 "DNB-type drug metabolite #1 cell membrane inwards permeability constant Kd"},

{"DNB_CELL_DIFF_OUT_MET1", 1, 0, 0,
 "Membrane diff (out)",
 "DNB-type drug metabolite #1 cell membrane outwards permeability constant Kd"},

{"DNB_HALFLIFE_MET1", 2.0, 0, 0,
 "Half-life",
 "DNB-type drug metabolite #1 half-life (hours)"},

{"DNB_DIFF_COEF_MET2", 4.29e-7, 0, 0,
 "Spheroid diffusion coeff",
 "DNB-type drug metabolite #2 diffusion coefficient in the spheroid"},

{"DNB_MEDIUM_DIFF_MET2", 5.0e-6, 0, 0,
 "Medium diffusion coeff",
 "DNB-type drug metabolite #2 diffusion coefficient in the medium"},

{"DNB_CELL_DIFF_IN_MET2", 1, 0, 0,
 "Membrane diff (in)",
 "DNB-type drug metabolite #2 cell membrane inwards permeability constant Kd"},

{"DNB_CELL_DIFF_OUT_MET2", 1, 0, 0,
 "Membrane diff (out)",
 "DNB-type drug metabolite #2 cell membrane outwards permeability constant Kd"},

{"DNB_HALFLIFE_MET2", 0.1, 0, 0,
 "Half-life",
 "DNB-type drug metabolite #2 half-life (hours)"},

//------------------------------------------------------
// Cell type 1 metabolism and kill experiment parameters
//------------------------------------------------------

// Parent
{"DNB_KMET0_CELL1", 0.95, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL1", 1.0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL1", 0.126, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL1", 0, 0, 0,
 "Vmax",
    "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism: Kmet0 -> Kmet0 + Vmax/(Km + C)"},

{"DNB_KM_CELL1", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL1", 0.001, 0, 0,
 "Klesion    ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 1 DNB drug kill experiment parameters

    {"DNB_KILL_MODEL_CELL1", 0, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 4 = K x Ci, 5 = K x Ci^2"},

    {"DNB_KILL_O2_CONC_CELL1", 0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL1", 0, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL1", 0, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL1", 0, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

// Metabolite 1
{"DNB_KMET0_CELL1_MET1", 0.313, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL1_MET1", 0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL1_MET1", 1, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL1_MET1", 0, 0, 0,
 "Vmax",
 "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"DNB_KM_CELL1_MET1", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL1_MET1", 0.001, 0, 0,
 "Klesion  ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 1 DNB metabolite 1 kill experiment parameters

    {"DNB_KILL_MODEL_CELL1_MET1", 4, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"DNB_KILL_O2_CONC_CELL1_MET1", 0.0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL1_MET1", 0.01, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL1_MET1", 60, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL1_MET1", 0.9, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

// Metabolite 2
{"DNB_KMET0_CELL1_MET2", 0.72, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL1_MET2", 0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL1_MET2", 1, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL1_MET2", 0, 0, 0,
 "Vmax",
 "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"DNB_KM_CELL1_MET2", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL1_MET2", 0.001, 0, 0,
 "Klesion  ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},


    {"DNB_KILL_MODEL_CELL1_MET2", 4, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"DNB_KILL_O2_CONC_CELL1_MET2", 0.0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL1_MET2", 0.01, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL1_MET2", 60, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL1_MET2", 0.9, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

//------------------------------------------------------
// Cell type 2 metabolism and kill experiment parameters
//------------------------------------------------------

// Parent
{"DNB_KMET0_CELL2", 0.95, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL2", 1.0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL2", 0.126, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL2", 0, 0, 0,
 "Vmax",
 "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"DNB_KM_CELL2", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL2", 0.001, 0, 0,
 "Klesion    ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 2 DNB drug kill experiment parameters

    {"DNB_KILL_MODEL_CELL2", 0, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"DNB_KILL_O2_CONC_CELL2", 0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL2", 0, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL2", 0, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL2", 0, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

// Metabolite 1
{"DNB_KMET0_CELL2_MET1", 0.313, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL2_MET1", 0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL2_MET1", 1, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL2_MET1", 0, 0, 0,
 "Vmax",
 "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"DNB_KM_CELL2_MET1", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL2_MET1", 0.001, 0, 0,
 "Klesion  ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 2 DNB metabolite 1 kill experiment parameters

    {"DNB_KILL_MODEL_CELL2_MET1", 4, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"DNB_KILL_O2_CONC_CELL2_MET1", 0.0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL2_MET1", 0.01, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL2_MET1", 60, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL2_MET1", 0.9, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

// Metabolite 2
{"DNB_KMET0_CELL2_MET2", 0.72, 0, 0,
 "Kmet0",
 "DNB-type drug max value of 1st order rate constant for metabolism under zero oxygen"},

{"DNB_C2_CELL2_MET2", 0, 0, 0,
 "C2",
 "DNB-type drug C2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_KO2_CELL2_MET2", 1, 0, 0,
 "KO2",
 "DNB-type drug KO2 in function for oxygen-dependence of rate of metabolism"},

{"DNB_VMAX_CELL2_MET2", 0, 0, 0,
 "Vmax",
 "DNB-type drug Vmax in function for oxygen-dependence of rate of metabolism"},

{"DNB_KM_CELL2_MET2", 1, 0, 0,
 "Km",
 "DNB-type drug Km in function for oxygen-dependence of rate of metabolism"},

{"DNB_KLESION_CELL2_MET2", 0.001, 0, 0,
 "Klesion  ",
 "DNB-type drug Klesion is the parameter that converts total metabolite into lesion level"},

    // Cell type 2 DNB metabolite 2 kill experiment parameters

    {"DNB_KILL_MODEL_CELL2_MET2", 4, 0, 0,
     "Kill model",
     "Model of DNB-type drug killing: 1 = K x metabolism, 2 = K x Ci x metabolism, 3 = K x metabolism^2"},

    {"DNB_KILL_O2_CONC_CELL2_MET2", 0.0, 0, 0,
     "O2 conc",
     "DNB-type drug constant O2 concentration in kill experiment"},

    {"DNB_KILL_DRUG_CONC_CELL2_MET2", 0.01, 0, 0,
     "Drug conc",
     "DNB-type drug constant drug concentration in kill experiment"},

    {"DNB_KILL_DURATION_CELL2_MET2", 60, 0, 0,
     "Duration",
     "DNB-type drug duration of kill experiment"},

    {"DNB_KILL_FRACTION_CELL2_MET2", 0.9, 0, 0,
     "Kill fraction",
     "DNB-type drug fraction of cells killed in the experiment"},

//==========================
// Radiotherapy parameters
//==========================

{"RADIATION_ALPHA_H", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},


{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"GROWTH_FRACTION_1", 0.25, 0, 0,
"Growth fraction threshold 1",
"Growth fraction threshold 1"},

{"GROWTH_FRACTION_2", 0.1, 0, 0,
"Growth fraction threshold 2",
"Growth fraction threshold 2"},

{"GROWTH_FRACTION_3", 0.01, 0, 0,
"Growth fraction threshold 3",
"Growth fraction threshold 3"},

{"SPCRAD", 200.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

{"USE_EXTRA", 0, 0, 1,
"Use extra conc",
"Use extracellular O2 concentration to determine cell death"},

{"USE_RELAX", 1, 0, 1,
"Use O2 relaxation solver",
"Use over- and under-relaxation to solve reaction-diffusion for oxygen"},

{"USE_PAR_RELAX", 1, 0, 1,
"Use parallel O2 relaxation solver",
"Use over- and under-relaxation to solve reaction-diffusion for oxygen, with parallelized over-relaxation"},

//{"USE_RADIATION", 1, 0, 1,
//"Use radiation?",
//"Treatment with radiation"},

//{"USE_TREATMENT_FILE", 0, 0, 1,
//"Use treatment file?",
//"Treatment programme is specified in the treatment file"},

//{"TREATMENT_FILE_NAME", 0, 0, 0,
//"treatment.data",
//"The treatment file contains data describing the drug and radiation dosing schedule"},

{"USE_DROP", 0, 0, 1,
"Account for drop deformation",
"Account for drop deformation when it is released to sit at the bottom of the well"},

{"NDROP", 1000, 0, 0,
"Dropping cell count",
"Number of cells in the spheroid when it is dropped."},

{"DROP_ALPHA", 0.4, 0, 0,
"Contact_diameter/diameter",
"Drop parameter alpha = initial (surface contact diameter)/(blob diameter).  Must be < 1."},

{"DROP_BETA", 0.6, 0, 0,
"Height/diameter",
"Drop parameter beta = initial (blob height)/(blob diameter).  Must be < 1."},

    {"SAVE_PROFILE_DATA",0,0,1,
     "Save profile data",
     "Save data for profile plots as a specified interval"},

    {"SAVE_PROFILE_DATA_FILE_NAME",0,0,0,
     "profile_data",
     "Base file name for saving profile data"},

    {"SAVE_PROFILE_DATA_INTERVAL",0,0,0,
     "Interval",
     "Time interval for saving profile data"},

    {"SAVE_PROFILE_DATA_NUMBER",1,0,0,
     "Number",
     "Number of times to save profile data"},


// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted

{"DUMMY_HYPOXIA_THRESHOLD", 0, 0, 0,
"Hypoxia threshold",
"Select the intracellular O2 level below which the cell is counted as hypoxic"},

{"DUMMY_GROWTH_FRACTION", 0, 0, 0,
"Growth fraction",
"Select the threshold fraction of average growth rate (i.e. with no nutrient limits) used to count slow-growing cells"},

// Time-series plots
    {"nlive",                     1, 0,1,"","Number of live cells"},
    {"nanoxiadead",               1, 0,1,"","Total number of cells that have been killed by anoxia"},
    {"ndrugAdead",                1, 0,1,"","Total number of cells that have been killed by drugA"},
    {"ndrugBdead",                0, 0,1,"","Total number of cells that have been killed by drugB"},
    {"nradiationdead",            1, 0,1,"","Total number of cells that have been killed by radiation"},
    {"nanoxiatagged",             1, 0,1,"","Current number of cells tagged to die by anoxia"},
    {"ndrugAtagged",              0, 0,1,"","Current number of cells tagged to die by drugA"},
    {"ndrugBtagged",              0, 0,1,"","Current number of cells tagged to die by drugB"},
    {"nradiationtagged",          1, 0,1,"","Current number of cells tagged to die by radiation"},
    {"diameter",                  1, 0,1,"","Spheroid diameter (um)"},
    {"volume",                    0, 0,1,"","Spheroid volume (mm3)"},
    {"hypoxicfraction",           1, 0,1,"","Fraction of cells with oxygen level below the specified threshold for hypoxia"},
    {"growthfraction",            1, 0,1,"","Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits"},
    {"necroticfraction",          1, 0,1,"","Percentage of the spheroid that is necrotic = (number of vacant sites)/(number of sites taken up by the spheroid)"},
    {"platingefficiency",         0, 0,1,"","Percentage of live cells that are viable"},
    {"mediumoxygen",              0, 0,1,"","Average concentration of oxygen in the medium (far-field)"},
    {"mediumglucose",             0, 0,1,"","Average concentration of glucose in the medium (far-field)"},
    {"mediumTPZdrug",             0, 0,1,"","Average concentration of TPZ drug in the medium (far-field)"},
    {"mediumDNBdrug",             0, 0,1,"","Average concentration of DNB drug in the medium (far-field)"},

// Profile plots
    {"MULTI",                     1, 0,1,"","Selected constituent on a line through the blob centre"},
    {"CFSE",                      0, 0,1,"","Extracellular CFSE concentration on a line through the blob centre"},
    {"Oxygen",                    0, 0,1,"","Extracellular oxygen concentration on a line through the blob centre"},
    {"Glucose",                   0, 0,1,"","Extracellular glucose concentration on a line through the blob centre"},
    {"Tracer",                    0, 0,1,"","Extracellular tracer concentration on a line through the blob centre"},
    {"TPZdrug",                   0, 0,1,"","Extracellular TPZ drug concentration on a line through the blob centre"},
    {"TPZmetab1",                 0, 0,1,"","Extracellular TPZ metabolite 1 concentration on a line through the blob centre"},
    {"TPZmetab2",                 0, 0,1,"","Extracellular TPZ metabolite 2 concentration on a line through the blob centre"},
    {"DNBdrug",                   0, 0,1,"","Extracellular DNB drug concentration on a line through the blob centre"},
    {"DNBmetab1",                 0, 0,1,"","Extracellular DNB metabolite 1 concentration on a line through the blob centre"},
    {"DNBmetab2",                 0, 0,1,"","Extracellular DNB metabolite 2 concentration on a line through the blob centre"},
    {"growthrate",                1, 0,1,"","Cell growth rate on a line through the blob centre"},
    {"cellvolume",                1, 0,1,"","Cell volume fraction on a line through the blob centre"},
    {"O2byvolume",                0, 0,1,"","Cell volume fraction on a line through the blob centre"},
// Distribution plots
//    {"Oxygen",                    0, 0,1,"","Probability distribution of extracellular oxygen concentration"},
//    {"cellvolume",                0, 0,1,"","Probability distribution of cell volume fraction"}

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

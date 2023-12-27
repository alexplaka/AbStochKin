URLS=[
"abstochkin/index.html",
"abstochkin/het_calcs.html",
"abstochkin/process.html",
"abstochkin/graphing.html",
"abstochkin/agentstatedata.html",
"abstochkin/utils.html",
"abstochkin/base.html",
"abstochkin/de_calcs.html",
"abstochkin/simulation.html"
];
INDEX=[
{
"ref":"abstochkin",
"url":0,
"doc":"Agent-based (or Particle-based) Stochastic Kinetics."
},
{
"ref":"abstochkin.het_calcs",
"url":1,
"doc":"Some functions for calculating metrics of population heterogeneity."
},
{
"ref":"abstochkin.het_calcs.get_het_processes",
"url":1,
"doc":"Filter the heterogeneous processes from a given list of processes. A process is heterogeneous if any of its parameters are defined as such.",
"func":1
},
{
"ref":"abstochkin.het_calcs.richness",
"url":1,
"doc":"Calculate the species richness, or how many subspecies a species population comprises.",
"func":1
},
{
"ref":"abstochkin.het_calcs.idx_het",
"url":1,
"doc":"Calculate the Index of Heterogeneity (\u03c8), defined as the probability that two randomly chosen agents (without replacement) from a species population belong to different subspecies. - A homogeneous population returns 0. - A heterogeneous population with two distinct subspecies of equal fractional abundance (\u03c7=0.5)  approaches 0.5 as the population size increases. - A fully heterogeneous population returns 1.",
"func":1
},
{
"ref":"abstochkin.het_calcs.info_het",
"url":1,
"doc":"Information-theoretic measure of population heterogeneity. - A homogeneous population returns 0. - A heterogeneous population with two distinct subspecies of equal fractional abundance (\u03c7=0.5) returns ln(2). Note that this is true regardless of the population size. - For a fully heterogeneous population, the measure increases with population size and has no upper limit.",
"func":1
},
{
"ref":"abstochkin.process",
"url":2,
"doc":"Define a process of the form Reactants -> Products."
},
{
"ref":"abstochkin.process.Process",
"url":2,
"doc":"Define a unidirectional process: Reactants -> Products, where the Reactants and Products are specified using standard chemical notation. That is, stoichiometric coefficients (integers) and species names are specified. For example: 2A + B -> C. Attributes      reactants : dict The reactants of a given process are specified with key-value pairs describing each species name and its stoichiometric coefficient, respectively. products : dict The products of a given process are specified with key-value pairs describing each species name and its stoichiometric coefficient, respectively. k : float, int, list of floats, tuple of floats The  microscopic rate constant(s) for the given process. The data type of  k determines the \"structure\" of the population as follows: - A homogeneous population: if  k is a single value (float or int), then the population is assumed to be homogeneous with all agents of the reactant species having kinetics defined by this value. - A heterogeneous population with a distinct number of subspecies (each with a corresponding  k value): if  k is a list of floats, then the population is assumed to be heterogeneous with a number of subspecies equal to the length of the list. - A heterogeneous population with normally-distributed  k values: If  k is a tuple whose length is 2, then the population is assumed to be heterogeneous with a normally distributed  k value. The two entries in the tuple represent the mean and standard deviation (in that order) of the desired normal distribution. order : int The order of the process (or the molecularity of an elementary process). It is the sum of the stoichiometric coefficients of the reactants. species : set of strings A set of all species in a process. reacts_ : list of strings A list containing all the reactants in a process. prods_ : list of strings A list containing all the products in a process. Methods    - from_string Class method for creating a Process object from a string."
},
{
"ref":"abstochkin.process.Process.from_string",
"url":2,
"doc":"Create a process from a string. Parameters      proc_str : str A string describing the process in standard chemical notation (e.g., 'A + B -> C') k : float or int or list of floats or 2-tuple of floats The rate constant for the given process. If  k is a float or int, then the process is homogeneous. If  k is a list, then the population of the reactants constsists of distinct subspecies or subinteractions depending on the order. If  k is a 2-tuple, then the constant is normally-distributed with a mean and standard deviation specified in the tuple's elements. sep : str, default: '->' Specifies the characters that distinguish the reactants from the products. The default is '->'. The code also treats   > as a default, if it's present in  proc_str . Notes   - - Species names should not contain spaces, dashes, and should start with a non-numeric character. - Zeroth order processes should be specified by an empty space or 'None'. Examples     >>> Process.from_string(\"2A + B -> X\", 0.3) >>> Process.from_string(\" -> Y\", 0.1)  for a 0th order (birth) process. >>> Process.from_string(\"Protein_X -> None\", 0.15)  for a 1st order degradation process.",
"func":1
},
{
"ref":"abstochkin.process.ReversibleProcess",
"url":2,
"doc":"Define a reversible process. Attributes      k_rev : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the reverse process. is_heterogeneous_rev : bool Denotes if the parameter  k_rev exhibits heterogeneity (distinct subspecies/interactions or normally-distributed). Notes   - A ReversibleProcess object gets split into two Process objects (forward and reverse process) when the algorithm runs."
},
{
"ref":"abstochkin.process.ReversibleProcess.from_string",
"url":2,
"doc":"Create a reversible process from a string. Parameters      proc_str : str A string describing the process in standard chemical notation (e.g., 'A + B  C') k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the forward process. k_rev : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the reverse process. sep : str, default: ' ' Specifies the characters that distinguish the reactants from the products. The default is ' '. The code also treats    as a default, if it's present in  proc_str . Notes   - - Species names should not contain spaces, dashes, and should start with a non-numeric character. Examples     >>> ReversibleProcess.from_string(\"2A + B  X\", 0.3, k_rev=0.2)",
"func":1
},
{
"ref":"abstochkin.process.MichaelisMentenProcess",
"url":2,
"doc":"Define a process that obeys Michaelis-Menten kinetics. Attributes      catalyst : str Name of the species acting as a catalyst for this process. Km : float or int or list of floats or 2-tuple of floats  Microscopic Michaelis constant. Corresponds to the number of  catalyst agents that would produce half-maximal activity. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . is_heterogeneous_Km : bool Denotes if the parameter  Km exhibits heterogeneity (distinct subspecies/interactions or normally-distributed)."
},
{
"ref":"abstochkin.process.MichaelisMentenProcess.from_string",
"url":2,
"doc":"Create a Michaelis-Menten process from a string. Parameters      proc_str : str A string describing the process in standard chemical notation (e.g., 'A + B -> C') k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the given process. If  k is a float or int, then the process is homogeneous. If  k is a list, then the population of the reactants constsists of distinct subspecies or subinteractions depending on the order. If  k is a 2-tuple, then the constant is normally-distributed with a mean and standard deviation specified in the tuple's elements. catalyst : str Name of species acting as a catalyst. Km : float or int or list of floats or 2-tuple of floats  Microscopic Michaelis constant for the process. Heterogeneity in this parameter is determined by the type of  Km , using the same rules as for parameter  k . sep : str, default: '->' Specifies the characters that distinguish the reactants from the products. The default is '->'. The code also treats   > as a default, if it's present in  proc_str . Notes   - - Species names should not contain spaces, dashes, and should start with a non-numeric character. - Zeroth order processes should be specified by an empty space or 'None'. Examples     >>> MichaelisMentenProcess.from_string(\"A -> X\", k=0.3, catalyst='E', Km=10) >>> MichaelisMentenProcess.from_string(\"A -> X\", k=0.3, catalyst='alpha', Km=(10, 1 ",
"func":1
},
{
"ref":"abstochkin.process.RegulatedProcess",
"url":2,
"doc":"Define a process that is regulated. This class allows a Process to be defined in terms of how it is regulated. If there is only one regulating species, then the parameters have the same type as would be expected for a homogeneous/heterogeneous process. If there are multiple regulating species, then all parameters are a list of their expected type, with the length of the list being equal to the number of regulating species. Attributes      k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the given process. It is the  basal rate constant in the case of activation (or the minimum  k value) and the maximum rate constant in the case of repression. regulating_species : str or list of str Name of the regulating species. Multiple species can be specified as comma-separated in a string or a list of strings with the species names. alpha : float or int or list[float or int] Parameter denoting the degree of activation/repression. - 0  1: activation alpha is a multiplier: in the case of activation, the maximum rate constant value will be  alpha  k . In the case of repression, the minimum rate constant value will be  alpha  k . K50 : float or int or list of floats or 2-tuple of floats or list[float or int or list of floats or 2-tuple of floats]  Microscopic constant that corresponds to the number of  regulating_species agents that would produce half-maximal activation/repression. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . nH : float or int or list[float or int] Hill coefficient for the given process. Indicates the degree of cooperativity in the regulatory interaction. is_heterogeneous_K50 : bool or list of bool Denotes if the parameter  K50 exhibits heterogeneity (distinct subspecies/interactions or normally-distributed). regulation_type : str or list of str The type of regulation for this process based on the value of alpha: 'activation' or 'repression' or 'no regulation'. Notes   - Allowing a 0th order process to be regulated. However, heterogeneity in  k and  K50 (or any other parameters) is not allowed for such a process."
},
{
"ref":"abstochkin.process.RegulatedProcess.from_string",
"url":2,
"doc":"Create a regulated process from a string. Parameters      proc_str : str A string describing the process in standard chemical notation (e.g., 'A + B -> C') k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the given process. It is the  basal rate constant in the case of activation (or the minimum  k value) and the maximum rate constant in the case of repression. If  k is a float or int, then the process is homogeneous. If  k is a list, then the population of the reactants constsists of distinct subspecies or subinteractions depending on the order. If  k is a 2-tuple, then the constant is normally-distributed with a mean and standard deviation specified in the tuple's elements. Note that  k cannot be zero for this form of regulation. regulating_species : str or list of str Name of the regulating species. alpha : float or int or list[float or int] Parameter denoting the degree of activation/repression. - 0  1: activation alpha is a multiplier: in the case of activation, the maximum rate constant value will be  alpha  k . In the case of repression, the minimum rate constant value will be  alpha  k . K50 : float or int or list of floats or 2-tuple of floats or list of each of the previous types  Microscopic constant that corresponds to the number of  regulating_species agents that would produce half-maximal activation/repression. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . nH : float or int or list[float or int] Hill coefficient for the given process. Indicates the degree of cooperativity in the regulatory interaction. sep : str, default: '->' Specifies the characters that distinguish the reactants from the products. The default is '->'. The code also treats   > as a default, if it's present in  proc_str . Notes   - - Species names should not contain spaces, dashes, and should start with a non-numeric character. - Zeroth order processes should be specified by an empty space or 'None'. Examples     >>> RegulatedProcess.from_string(\"A -> X\", k=0.2, regulating_species='X', alpha=2, K50=10, nH=1) >>> RegulatedProcess.from_string(\"A -> X\", k=0.3, regulating_species='X', alpha=0.5, K50=[10, 15], nH=2) >>> RegulatedProcess.from_string(\"A + B -> X\", k=0.5, regulating_species='B, X', alpha=[2, 0], K50=[(15, 5), [10, 15 , nH=[1, 2])",
"func":1
},
{
"ref":"abstochkin.process.RegulatedMichaelisMentenProcess",
"url":2,
"doc":"Define a process that is regulated and obeys Michaelis-Menten kinetics. This class allows a Michaelis-Menten Process to be defined in terms of how it is regulated. If there is only one regulating species, then the parameters have the same type as would be expected for a homogeneous/heterogeneous process. If there are multiple regulating species, then all parameters are a list of their expected type, with the length of the list being equal to the number of regulating species. Attributes      k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the given process. It is the  basal rate constant in the case of activation (or the minimum  k value) and the maximum rate constant in the case of repression. regulating_species : str or list of str Name of the regulating species. Multiple species can be specified as comma-separated in a string or a list of strings with the species names. alpha : float or int or list[float or int] Parameter denoting the degree of activation/repression. - 0  1: activation alpha is a multiplier: in the case of activation, the maximum rate constant value will be  alpha  k . In the case of repression, the minimum rate constant value will be  alpha  k . K50 : float or int or list of floats or 2-tuple of floats or list[float or int or list of floats or 2-tuple of floats]  Microscopic constant that corresponds to the number of  regulating_species agents that would produce half-maximal activation/repression. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . nH : float or int or list[float or int] Hill coefficient for the given process. Indicates the degree of cooperativity in the regulatory interaction. is_heterogeneous_K50 : bool or list of bool Denotes if the parameter  K50 exhibits heterogeneity (distinct subspecies/interactions or normally-distributed). regulation_type : str or list of str The type of regulation for this process based on the value of alpha: 'activation' or 'repression' or 'no regulation'. catalyst : str Name of the species acting as a catalyst for this process. Km : float or int or list of floats or 2-tuple of floats  Microscopic Michaelis constant. Corresponds to the number of  catalyst agents that would produce half-maximal activity. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . is_heterogeneous_Km : bool Denotes if the parameter  Km exhibits heterogeneity (distinct subspecies/interactions or normally-distributed). Notes   - Currently only implemented for 1st order processes. 0th order processes cannot obey Michaelis-Menten kinetics and 2nd order Michaelis-Menten processes are not implemented yet."
},
{
"ref":"abstochkin.process.RegulatedMichaelisMentenProcess.from_string",
"url":2,
"doc":"Create a regulated Michaelis-Menten process from a string. Parameters      proc_str : str A string describing the process in standard chemical notation (e.g., 'A + B -> C') k : float or int or list of floats or 2-tuple of floats The  microscopic rate constant for the given process. It is the  basal rate constant in the case of activation (or the minimum  k value) and the maximum rate constant in the case of repression. If  k is a float or int, then the process is homogeneous. If  k is a list, then the population of the reactants constsists of distinct subspecies or subinteractions depending on the order. If  k is a 2-tuple, then the constant is normally-distributed with a mean and standard deviation specified in the tuple's elements. Note that  k cannot be zero for this form of regulation. regulating_species : str or list of str Name of the regulating species. alpha : float or int or list[float or int] Parameter denoting the degree of activation/repression. - 0  1: activation alpha is a multiplier: in the case of activation, the maximum rate constant value will be  alpha  k . In the case of repression, the minimum rate constant value will be  alpha  k . K50 : float or int or list of floats or 2-tuple of floats or list of each of the previous types  Microscopic constant that corresponds to the number of  regulating_species agents that would produce half-maximal activation/repression. Heterogeneity in this parameter is determined by the type of  K50 , using the same rules as for parameter  k . nH : float or int or list[float or int] Hill coefficient for the given process. Indicates the degree of cooperativity in the regulatory interaction. catalyst : str Name of species acting as a catalyst. Km : float or int or list of floats or 2-tuple of floats  Microscopic Michaelis constant for the process. Heterogeneity in this parameter is determined by the type of  Km , using the same rules as for parameter  k . sep : str, default: '->' Specifies the characters that distinguish the reactants from the products. The default is '->'. The code also treats   > as a default, if it's present in  proc_str . Notes   - - Species names should not contain spaces, dashes, and should start with a non-numeric character. - Zeroth order processes should be specified by an empty space or 'None'. Examples     >>> RegulatedMichaelisMentenProcess.from_string(\"A -> X\", k=0.2, regulating_species='X', alpha=2, K50=10, nH=1, catalyst='E', Km=15) >>> RegulatedMichaelisMentenProcess.from_string(\"A -> X\", k=0.3, regulating_species='A', alpha=0.5, K50=[10, 15], nH=2, catalyst='C', Km=5)",
"func":1
},
{
"ref":"abstochkin.process.NullSpeciesNameError",
"url":2,
"doc":"Error when the species name is an empty string."
},
{
"ref":"abstochkin.process.update_all_species",
"url":2,
"doc":"Categorize all species in a list of processes. Extract all species from a list of processes. Then categorize each of them as a reactant or product and list the process(es) it takes part in. Parameters      procs : tuple A tuple of objects of type  Process or its subclasses. Returns    - tuple all_species : set of strings A set of all species present in the processes. procs_by_reactant : dict A dictionary whose keys are the species that are reactants in one or more processes. The value for each key is a list of processes. procs_by_product : dict A dictionary whose keys are the species that are products in one or more processes. The value for each key is a list of processes.",
"func":1
},
{
"ref":"abstochkin.graphing",
"url":3,
"doc":"Graphing for AbStochKin simulations."
},
{
"ref":"abstochkin.graphing.Graph",
"url":3,
"doc":"Graphing class for displaying the results of AbStochKin simulations. Notes   - To successfully use the LaTeX engine for rendering text on Linux, run the following command in a terminal:  sudo apt install cm-super ."
},
{
"ref":"abstochkin.graphing.Graph.setup_spines_ticks",
"url":3,
"doc":"Make only the left and bottom spines/axes visible on the graph and place major ticks on them. Also set the minor ticks.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.plot_ODEs",
"url":3,
"doc":"Plot the deterministic trajectories of all species obtained by obtaining the solution to a system of ODEs. Parameters      de_data : DEcalcs object Data structure containing all the data related to solving the system of ODEs. num_pts : int, default: 1000, optional Number of points used to calculate DE curves at. Used to approximate a smooth/continuous curve. species : sequence of strings, default: (), optional An iterable sequence of strings specifying the species names to plot. If no species are specified (the default), then all species trajectories are plotted. ax_loc : tuple, optional If the figure is made up of subplots, specify the location of the axis to draw the data at. Ex: for two subplots, the possible values of  ax_loc are (0, ) and (1, ). That's because the  self.ax object is a 1-D array. For figures with multiple rows and columns of subplots, a 2-D tuple is needed.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.plot_trajectories",
"url":3,
"doc":"Graph simulation time trajectories.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.plot_avg_std",
"url":3,
"doc":"Graph simulation average trajectories and 1-standard-deviation envelopes.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.plot_eta",
"url":3,
"doc":"Graph the coefficient of variation.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.plot_het_metrics",
"url":3,
"doc":"Graph species- and process-specific metrics of population heterogeneity.",
"func":1
},
{
"ref":"abstochkin.graphing.Graph.savefig",
"url":3,
"doc":"Save the figure as a file.",
"func":1
},
{
"ref":"abstochkin.agentstatedata",
"url":4,
"doc":"Class for storing the state of all agents of a certain species during an AbStochKin simulation."
},
{
"ref":"abstochkin.agentstatedata.AgentStateData",
"url":4,
"doc":"Class for storing the state of all agents of a certain species during an AbStochKin simulation. Attributes      p_init : int The initial population size of the species whose data is represented in an  AgentStateData object. max_agents : int The maximum number of agents for the species whose data is represented in an  AgentStateData object. reps : int The number of times the AbStochKin algorithm will repeat a simulation. This will be the length of the  asv list. fill_state: int asv_ini : numpy.ndarray Agent-State Vector (asv) is a species-specific 2-row vector to monitor agent state according to Markov's property. This array is the initial asv, i.e., at  t=0 . The array shape is  (2, max_agents) asv : list of numpy.ndarray A list of length  reps with copies of  asv_ini . Each simulation run uses its corresponding entry in  asv to monitor the state of all agents."
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.p_init",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.max_agents",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.reps",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.fill_state",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.asv_ini",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.asv",
"url":4,
"doc":""
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.apply_markov_property",
"url":4,
"doc":"The future state of the system depends only on its current state. This method is called at the end of each time step in an AbStochKin simulation. Therefore, the new agent-state vector becomes the current state.",
"func":1
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.cleanup_asv",
"url":4,
"doc":"Empty the contents of the array  asv .",
"func":1
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.get_vals_o1",
"url":4,
"doc":"Get random values in [0,1) at a given time step for agents of a given state. Agents of other states have a value of zero. Get probability values at a given time step for agents of a given state. Agents of other states have a transition probability of zero. Notes   - Note that only elements of the  asv that have the same state in the previous and current time steps are considered. This is to ensure that agents that have already transitioned to a different state in the current time step are not reconsidered for a possible transition.",
"func":1
},
{
"ref":"abstochkin.agentstatedata.AgentStateData.get_vals_o2",
"url":4,
"doc":"Get random values in [0,1) at a given time step for interactions between agents of a given state. Agents of other states have a value of zero. Get probability values at a given time step for interactions between agents of a given state. Interactions of agents in other states have a transition probability of zero. Notes   - Note that only elements of the  asv that have the same state in the previous and current time steps are considered. This is to ensure that agents that have already transitioned to a different state in the current time step are not reconsidered for a possible transition.",
"func":1
},
{
"ref":"abstochkin.utils",
"url":5,
"doc":"Some utility functions on generating random number streams, measuring a function's runtime, calculating a statistic measuring the goodness of fit when comparing time series data, and performing unit conversion of kinetic parameters."
},
{
"ref":"abstochkin.utils.rng_streams",
"url":5,
"doc":"Generate independent streams of random numbers spawned from the same initial seed. Parameters      n : int number of generators/streams to generate. random_state : int, optional initial seed from which n new seeds are spawned. Returns    - list[PCG64DXSM Generator objects] List of  n generators of independent streams of random numbers Notes   - See https: numpy.org/doc/stable/reference/random/parallel.html for more info. On PCG64DXSM: - https: numpy.org/doc/stable/reference/random/upgrading-pcg64.html upgrading-pcg64 - https: numpy.org/doc/stable/reference/random/bit_generators/pcg64dxsm.html Examples     >>> rng = rng_streams(5)  make 5 random number generators >>> a = rng[1].integers(1, 10, 100)",
"func":1
},
{
"ref":"abstochkin.utils.measure_runtime",
"url":5,
"doc":"Decorator for measuring the duration of a function's execution.",
"func":1
},
{
"ref":"abstochkin.utils.r_squared",
"url":5,
"doc":"Compute the coefficient of determination, \\(R^2\\). In the case of comparing the average AbStochKin-simulated species trajectory to its deterministic trajectory. Since the latter is only meaningful for a homogeneous population, \\(R^2\\) should be close to  1 for a simulated homogeneous process. For a heterogeneous process, it can be interpreted as how close the simulated trajectory is to the deterministic trajectory of a  homogeneous process. In this case, \\(R^2\\) would not be expected to be close to  1 and the importance of looking at this metric is questionable. Parameters      actual : numpy.array Actual data obtained through a simulation. theoretical : numpy.array Theoretical data to compare the actual data to. Returns    - float The coefficient of determination, \\(R^2\\).",
"func":1
},
{
"ref":"abstochkin.utils.macro_to_micro",
"url":5,
"doc":"Convert a kinetic parameter value from macroscopic to microscopic form. The ABK algorithm uses microscopic kinetic constants, thus necessitating the conversion of any molar quantities to their microscopic counterpart. For a kinetic parameter, the microscopic form is interpreted as the number of transition events per second (or whatever the time unit may be). For a molar quantity, its microscopic form is the number of particles in the given volume. Parameters      macro_val : float or int The value of the parameter to be converted, expressed in terms of molar quantities. volume : float or int The volume, in liters, in which the process that the given parameter value is a descriptor of. order : int, default: 0 The order of the process whose kinetic parameter is to be converted. The default value of 0 is for parameters (such as Km or K50) whose units are molarity. Returns     float A kinetic parameter is returned in units of reciprocal seconds. A molar quantity is returned as the number of particles in the given volume. Notes   - - A kinetic parameter for a 1st order process will remain unchanged because its units are already reciprocal seconds. Reference     - Plakantonakis, Alex. \u201cAgent-based Kinetics: A Nonspatial Stochastic Method for Simulating the Dynamics of Heterogeneous Populations.\u201d OSF Preprints, 26 July 2019. Web. Section 2.1.",
"func":1
},
{
"ref":"abstochkin.base",
"url":6,
"doc":"Base class, AbStochKin, for initializing and storing all data for performing stochastic simulations using the Agent-based Kinetics method. A simulation project can be initialized and run as follows: Example    - >>> from abstochkin import AbStochKin >>> sim = AbStochKin() >>> sim.add_process_from_str('A -> ', 0.2)  degradation process >>> sim.simulate(p0={'A': 100}, t_max=20) >>>  All data for the above simulation is stored in  sim.sims[0] . >>> >>>  Now set up a new simulation without actually running it. >>> sim.simulate(p0={'A': 10}, t_max=10, n=50, run=False) >>>  All data for the new simulation is stored in  sim.sims[1] . >>>  The simulation can then be manually run using methods >>>  documented in the class  Simulation ."
},
{
"ref":"abstochkin.base.AbStochKin",
"url":6,
"doc":"Base class for Agent-based Kinetics (AbStochKin) simulator. Attributes      time_unit : str, default : sec, optional A string of the time unit to be used for describing the kinetics of the given processes. processes : list A list of the processes that the AbStochKin object has. het_processes : list A list of the processes where population heterogeneity in one of the parameters is to be modeled. This list is a subset of the  processes attribute. sims : list A list of all simulations performed for the given set of processes. Each member of the list is an object of the  Simulation class and contains all data for that simulation."
},
{
"ref":"abstochkin.base.AbStochKin.add_processes_from_file",
"url":6,
"doc":"Add a batch of processes from a text file.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.extract_process_from_str",
"url":6,
"doc":"Extract a process and all of its specified parameters from a string. This functions parses a string specifying all values and parameters needed to define a process. It then creates a Process object based on the extracted data.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.add_process_from_str",
"url":6,
"doc":"Add a process by specifying a string: 'reactants -> products'. Additional arguments determine if a specialized process (such as a reversible, regulated, or Michaelis-Menten process) is to be defined.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.add_process",
"url":6,
"doc":"Add a process by using a dictionary for the reactants and products. Additional arguments determine if a specialized process (such as a reversible, regulated, or Michaelis-Menten process) is to be defined.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.del_process_from_str",
"url":6,
"doc":"Delete a process by specifying a string: 'reactants -> products'.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.del_process",
"url":6,
"doc":"Delete a process by using a dictionary for the reactants and products.",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.simulate",
"url":6,
"doc":"Start an AbStochKin simulation by creating an instance of the class  Simulation . The resulting object is appended to the list in the class attribute  AbStochKin.sims . Parameters      p0 : dict[str: int] Dictionary specifying the initial population sizes of all species in the given processes. t_max : float or int Numerical value of the end of simulated time in the units specified in the class attribute  AbStochKin.time_unit . dt : float, default: 0.1, optional The duration of the time interval that the simulation's algorithm considers. The current implementation only supports a fixed time step interval whose value is  dt . n : int, default: 100, optional The number of repetitions of the simulation to be performed. random_seed : int, default: 19, optional A number used to seed the random number generator. solve_odes : bool, default: True, optional Specify whether to numerically solve the system of ODEs defined from the given set of processes. ode_method : str, default: RK45, optional Available ODE methods: RK45, RK23, DOP853, Radau, BDF, LSODA. run : bool, default: True, optional Specify whether to run an AbStochKin simulation. show_plots : bool, default: True, optional Specify whether to graph the results of the AbStochKin simulation. multithreading : bool, default: True, optional Specify whether to parallelize the simulation using multithreading. If  False , the ensemble of simulations is run sequentially. max_agents_by_species : None or dict, default: dict Specification of the maximum number of agents that each species should have when running the simulation. If  None , that a default approach will be taken by the class  Simulation and the number for each species will be automatically determined (see method  Simulation._setup_runtime_data() for details). The entries in the dictionary should be  species name (string): number (int) . max_agents_multiplier : float or int, default: 2 This parameter is used to calculate the maximum number of agents of each species that the simulation engine allocates memory for. This be determined by multiplying the maximum value of the ODE time trajectory for this species by the multiplier value specified here. _return_simulation : bool Determines if the  self.simulate method returns a  Simulation object or appends it to the list  self.sims . Returning a  Simulation object is needed when calling the method  simulate_series_in_parallel .",
"func":1
},
{
"ref":"abstochkin.base.AbStochKin.simulate_series_in_parallel",
"url":6,
"doc":"Perform a series of simulations in parallel by initializing separate processes. Each process runs a simulation and appends a  Simulation object in the list  self.sims . Parameters      series_kwargs : list of dict A list containing dictionaries of the keyword arguments for performing each simulation in the series. The number of elements in the list is the number of simulations that will be run. max_workers : int, default: None The maximum number of processes to be used for performing the given series of simulations. If None, then as many worker processes will be created as the machine has processors. Examples     - Run a series of simulations by varying the initial population size of A. >>> from abstochkin import AbStochKin >>> sim = AbStochKin() >>> sim.add_process_from_str(\"A -> B\", 0.3, catalyst='E', Km=10) >>> series_kwargs = [{\"p0\": {'A': i, 'B': 0, 'E': 10}, \"t_max\": 10} for i in range(40, 51)] >>> sim.simulate_series_in_parallel(series_kwargs)",
"func":1
},
{
"ref":"abstochkin.de_calcs",
"url":7,
"doc":"Perform deterministic calculations on a set of processes. Construct the ordinary differential equations (ODEs) describing the system and obtain a numerical solution."
},
{
"ref":"abstochkin.de_calcs.DEcalcs",
"url":7,
"doc":"Perform deterministic calculations given the processes specified in an AbStochKin simulation object. Attributes      p0 : dict[str: int] Dictionary specifying the initial population sizes of all species in the given processes. t_min : float or int Numerical value of the start of simulated time in the units specified in the class attribute  time_unit . t_max : float or int Numerical value of the end of simulated time in the units specified in the class attribute  time_unit . processes : list of Process objects A list of all the processes that define the system of ODEs. Each process is a  Process object. ode_method: str Method for the ODE solver to use. time_unit: str A string of the time unit to be used for describing the kinetics of the given processes."
},
{
"ref":"abstochkin.de_calcs.DEcalcs.setup_ODEs",
"url":7,
"doc":"Set up the system of ODEs to be used for computing the deterministic trajectory of all the species in the given processes. The equations consist of sympy objects. Parameters      agent_based : bool, default: True If True, set up the agent-based (or microscopic) form of ODEs. For instance, for the process  2X -> Y , the ODE for species  X would include an  X(X - 1) term instead of  X^2 (the canonical form). If False, the canonical form of the ODEs is constructed. Notes   - The rate constant,  k , for a given process is taken to be the mean of  k , unless  k was defined to be normally-distributed, in which case  k is a 2-tuple and  k[0] is the specified mean. Implementation note: Building up the ODE expressions by separately iterating over the processes for products and reactants. This is to properly handle 0th order processes for product species. For example, for the 0th order process '  > X' with rate constant k1, the ODE is dX/dt = k1.",
"func":1
},
{
"ref":"abstochkin.de_calcs.DEcalcs.get_term_multiplier",
"url":7,
"doc":"Generate the multiplicative term (or terms) needed for generating the correct algebraic expressions for specialized processes (such as Michaelis-Menten and regulated processes).",
"func":1
},
{
"ref":"abstochkin.de_calcs.DEcalcs.solve_ODEs",
"url":7,
"doc":"Solve system of ordinary differential equations (ODEs). Notes   - Using the solver  scipy.integrate.solve_ivp , whose method can be one of the following: - RK45 : Explicit Runge-Kutta method of order 5(4). - RK23 : Explicit Runge-Kutta method of order 3(2). - DOP853 : Explicit Runge-Kutta method of order 8. - Radau : Implicit Runge-Kutta method of the Radau IIA family of order 5 - BDF : Implicit multistep variable-order (1 to 5) method based on a backward differentiation formula for the derivative approximation. - LSODA : Adams/BDF method with automatic stiffness detection and switching. Documentation       - https: docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html",
"func":1
},
{
"ref":"abstochkin.de_calcs.DEcalcs.get_fixed_pts",
"url":7,
"doc":"Not currently implemented.",
"func":1
},
{
"ref":"abstochkin.simulation",
"url":8,
"doc":"Perform an Agent-based Kinetics simulation. This module contains the code for the class  Simulation , which, along with the  SimulationMethodsMixin class, does everything that is needed to run an Agent-based Kinetics simulation and store its results. The class  AgentStateData is used by a  Simulation object to store and handle some of the necessary runtime data."
},
{
"ref":"abstochkin.simulation.Simulation",
"url":8,
"doc":"Run an Agent-based Kinetics simulation. Attributes      p0 : dict[str: int] Dictionary specifying the initial population sizes of all species in the given processes. t_max : float or int Numerical value of the end of simulated time in the units specified in the class attribute  AbStochKin.time_unit . dt : float The duration of the time interval that the simulation's algorithm considers. The current implementation only supports a fixed time step interval whose value is  dt . n : int The number of repetitions of the simulation to be performed. random_state : float or int A number used to seed the random number generator. use_multithreading : bool Specify whether to parallelize the simulation using multithreading. If  False , the ensemble of simulations is run sequentially. max_agents : dict Specification of the maximum number of agents that each species should have when running the simulation. An empty dictionary signifies that a default approach will be taken and the number for each species will be automatically determined (see method  _setup_runtime_data() for details). The entries in the dictionary should be  species name (string): number (int) . max_agents_multiplier : float or int Set the number of possible agents for each species to be the maximum value of the ODE solution for the species times  max_agents_multiplier . time_unit : str A string of the time unit to be used for describing the kinetics of the given processes. The parameters below are not class attributes, but are part of a  Simulation object's initialization to trigger specific actions to be automatically performed. Note that these actions can also be performed manually by calling the appropriate methods once a class object has been instantiated. Other Parameters         do_solve_ODEs : bool If  True , attempt to numerically solve the system of ODEs defined from the given set of processes. If  False , do not attempt to solve the ODEs and do not run the simulation. ODE_method : str Method to use when attempting to solve the system of ODEs (if  do_solve_ODEs is  True ). Available ODE methods: RK45, RK23, DOP853, Radau, BDF, LSODA. do_run : bool Specify whether to run the AbStochKin simulation. If  False , then a  Simulation object is created but the simulation is not run. A user can then manually run it by calling the class method  run_simulation() . show_graphs : bool Specify whether to show graphs of the results."
},
{
"ref":"abstochkin.simulation.Simulation.run_simulation",
"url":8,
"doc":"Run an ensemble of simulations and compute statistics of simulation data.",
"func":1
},
{
"ref":"abstochkin.simulation.Simulation.graph_results",
"url":8,
"doc":"Make graphs of the results. Parameters      species_to_show :  None or list of string(s), default:  None If  None , data for all species are plotted. graphs_to_show :  None or string or list of string(s), default:  None If  None , all graphs are shown. If a string is given then the graph that matches the string is shown. A list of strings shows all the graphs specified in the list. Graph specifications: 'avg' : Plot the average trajectories and their one-standard deviation envelopes. The ODE trajectories are also shown. 'traj' : Plot the species trajectories of individual simulations. 'ode' : Plot the deterministic species trajectories, obtained by numerically solving the ODEs. 'eta' : Plot the coefficient of variation (CoV) and the CoV assuming that all processes a species participates in obey Poisson statistics. 'het' : Plot species- and process-specific metrics of heterogeneity ( k and  \u03c8 ) and their one-standard deviation envelopes.",
"func":1
}
]
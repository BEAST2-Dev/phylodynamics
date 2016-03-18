package beast.phylodynamics;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.coalescent.TreeIntervals;

import java.io.File;
import java.util.*;

import static java.lang.System.loadLibrary;

/**
 @author Louis du Plessis

 Model: Tanja Stadler, Gabriel Leventhal
 ExpoTree library: Gabriel Leventhal

 * Usage:
        - The external expoTree library needs to be compiled separately.  This requires a Fortran, C and C++ compiler.
          (and some more libraries).  So far it has only been tested on Mac OS X.
          (will be changed soon to a jar in the libs directory - check code for details)

 * Issues:
        - Get a null pointer exception when I use TreeIntervalsInput and then try to get the nodes from the tree
          recovered from TreeIntervalsInput.get().TreeInput

 * To do / Questions:
        - SHIFTS SHIFTS SHIFTS!!
        - Find out what all the different node types are and where to use them (ask Gabriel)
        - Nr of susceptibles at the beginning (ask Gabriel)  n
        - How are shifts implemented? (ask Gabriel - and Denise)
            Seems that in BDSIR the shifts for different variables can take place at different times, eg. shift times
            lambda are not necessarily the shifting times for mu, etc.
            For now all parameter shifts in BDSIS need to be at the same times and each parameter needs the same nr of
            shifts.
            Is it necessary that the number of shifts are always specified a priori?  and the shift times?  Couldn't
            these be optimized as well?
        - Is it necessary to always recalculate the branching times?  Intervals?
        - Make fat calculation node?  implement store/restore...  (Ask Denise)
        - Properly implement requiresRecalculation
        - Check canHandleTipDates (for now copied BDSKY)

 java -classpath Beast2-read-only/:Phylodynamics/:../lib/* -Djava.library.path=Phylodynamics/lib/ beast.app.BeastMCMC ../../phylodynamics/exples/testBDSISEstimation.xml

 */



@Description("Calculate the exact likelihood of sampling a birth-death tree under an SIS model while accounting for "+
             "saturation effects (density dependence)")
@Citation("Gabriel E. Leventhal, Huldrych F. Guenthard, Sebastian Bonhoeffer and Tanja Stadler. Using an epidemiological " +
          "model for phylogenetic inference reveals density-dependence in HIV transmission. Molecular Biology and " +
          "Evolution. (2013). ")
public class ExactBDSIS extends TreeDistribution {

    protected static final int SAMPLE_EVENT = 0;
    protected static final int BRANCH_EVENT = 1;
    protected static final int SWITCH_EVENT = 3;
    protected static final int SURVIVAL     = 99;


    /******************************************************************************************************************/
    /* Input parameters for the epidemic (SIS model) */
    public Input<RealParameter> popSizeParameterInput =
            new Input<RealParameter>("popSize", "Population size, N in the SIS model", Input.Validate.REQUIRED);  // K
    public Input<RealParameter> rhoParameterInput =
            new Input<RealParameter>("rho", "Proportion of tips sampled at present time (default = 0)"); // rho

    // Set 1  (beta, mu, psi)
    public Input<RealParameter> infectionRateParameterInput =
            new Input<RealParameter>("infectionRate", "Infection rate parameter, beta in SIS model (lambda = beta*N)");
    public Input<RealParameter> deathRateParameterInput =
            new Input<RealParameter>("deathRate", "Death rate parameter, mu in SIS model");
    public Input<RealParameter> samplingRateParameterInput =
            new Input<RealParameter>("samplingRate", "Sampling rate per individual, psi in SIS model");

    // Set 2 (R0, delta, s)
    public Input<RealParameter> R0ParameterInput =
            new Input<RealParameter>("R0", "The basic reproductive number (R0 = beta*S0)/(mu+psi)", Input.Validate.XOR, infectionRateParameterInput);
    public Input<RealParameter> becomeUninfectiousRateParameterInput =
            new Input<RealParameter>("becomeUninfectiousRate", "Rate at which individuals become uninfectious, delta (throuch recovery or sampling)", Input.Validate.XOR, deathRateParameterInput);
    public Input<RealParameter> samplingProportionParameterInput =
            new Input<RealParameter>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRateParameterInput);

    /******************************************************************************************************************/
    /* Input parameters for phylogeny */
    /* TreeInput XOR TreeIntervalsInput inherited from ancestor */
    public Input<RealParameter> originParameterInput =
            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", Input.Validate.REQUIRED);

    public Input<RealParameter> eventTimesParameterInput =
            new Input<RealParameter>("eventTimes", "Times of events in the tree");

    /*
        Branching event -
        Sampling event  -
        Shift in        -
    */
    public Input<IntegerParameter> eventTypesParameterInput =
            new Input<IntegerParameter>("eventTypes", "Types of nodes in the tree");

    public Input<IntegerParameter> extantAtRootParameterInput =
            new Input<IntegerParameter>("extantAtRoot", "Number of extant lineages at the root");  // optional

    /******************************************************************************************************************/
    /* Random input parameters */
    public Input<BooleanParameter> survivalProbParameter =
            new Input<BooleanParameter>("survival", "Condition on the likelihood of observing the tree");  // False by default

    public Input<BooleanParameter> rescaleParameter =
            new Input<BooleanParameter>("rescale", "Rescale during calculation (improves accuracy)"); // True by default

    /******************************************************************************************************************/
    /* Shadow parameters */
    RealParameter    popSizeParameter, rhoParameter,
                     infectionRateParameter, deathRateParameter, samplingRateParameter,
                     R0Parameter, becomeUninfectiousRateParameter, samplingProportionParameter,
                     originParameter, eventTimesParameter;
    IntegerParameter eventTypesParameter, extantAtRootParameter;

    /* Model parameters */
    protected boolean       transform;           // Transform input parameters
    protected TreeInterface tree       = null;
    protected TreeIntervals intervals  = null;
    protected double [] K              = null,   // Population sizes
                        beta           = null,   // Infection rates
                        mu             = null,   // Death rates
                        psi            = null,   // Sampling rates through time
                        eventTimes     = null,   // Array of all event times on the tree
                        eventTimesTemp = null;   // Used for easy switching
    protected double    rho            = 0,      // Proportion of tips sampled at present
                        origin         = 0;      // Origin of the epidemic
    protected int []    eventTypes     = null,   // Different types of nodes in the tree ()
                        eventTypesTemp = null,   // Used for easy switching
                        survTypes      = null;   // Node types for calculating survival probability
    protected int       extantAtRoot   = 0,      // Number of extant lineages at the root
                        extantAtPresent= 0,      // Number of extant lineages at the present time
                        maxExtant      = 0,      // Maxmimum number of extant lineages
                        totalIntervals = 0,      // Number of SIS parameter intervals (parameter shifts = intervals - 1)
                        totalEvents    = 0;      // Number of events (tree and parameter shifts)
    protected boolean   survival       = false,  // Calculate survival probability
                        rescaling      = true;   // Use rescaling in calculations
    protected int       verbose        = 0;      // Not an input parameter, just for debugging (different levels)


    /******************************************************************************************************************/
    /*                                                                                                                */
    /* Initialization of native library                                                                               */
        /* Load expoTree library - relying on class name for now, else use input options                              */
        /*                                                                                                            */
        /* For now I'm using a dirty hack to load the native library and bundle it with the code.  As it is it only   */
        /* works in Mac OS X 10.x                                                                                     */
        /*                                                                                                            */
        /* The plan is to bundle different versions of the library for different platforms/architectures into a jar   */
        /* which will be distributed with the Phylodynamics library                                                   */
        /* I will include the project to build the jar library with Gabriel' expoTree library:                        */
        /*   https://bitbucket.org/gaberoo/expotree/src                                                               */
        /*                                                                                                            */
        /* Something along the lines of this:                                                                         */
        /*   http://stackoverflow.com/questions/2937406/how-to-bundle-a-native-library-and-a-jni-library-inside-a-jar */
        /*   http://frommyplayground.com/how-to-load-native-jni-library-from-jar/                                     */

    private static final String MODULE_NAME  = "Phylodynamics";
    private static final String LIBRARY_NAME = "expoTree";

    static {
        // Input option to add to library path (shouldn't be necessary if classpath set up right)
        // -Djava.library.path=/Users/louis/Documents/Projects/BirthDeath/Beast2Dev/beast2-read-only/build/phylodynamics/lib/

        String [] classpath = System.getProperty("java.class.path").split(":");
        String os = System.getProperty("os.name");
        boolean loaded = false;

        for (String p : classpath) {
            if (p.indexOf(MODULE_NAME) > 0 && (new File(p+"/lib/lib"+LIBRARY_NAME+".jnilib").exists() ||
                                               new File(p+"/lib/lib"+LIBRARY_NAME+".so").exists()        )) {
                if (os.equals("Mac OS X"))
                    System.load(p+"/lib/lib"+LIBRARY_NAME+".jnilib");
                else if (os.equals("Linux"))
                    System.load(p+"/lib/lib"+LIBRARY_NAME+".so");
                else
                    throw new RuntimeException("Operating system not supported (No native library for ExpoTree compiled yet...)");
                loaded = true;
                break;
            }
        }
        if (!loaded) {
            System.err.println("Library expoTree not found in "+MODULE_NAME+" lib directory...");
            System.err.println("Trying java.library.path (set command line argument or copy to run directory");
            loadLibrary(LIBRARY_NAME); // libexpoTree.jnilib or libexpoTree.so
        }
    }

    /* Native C method */
    private native double [] DirtyExpoTree(double[] K, double[] beta, double[] mu, double [] psi, double rho,
                                           double[] branchtimes, int[] nodetypes,
                                           boolean survival, int extant, boolean rescaling, int verbose);

    /******************************************************************************************************************/


    public void initAndValidate() {

        super.initAndValidate();

        // Assign shadow parameters
        popSizeParameter                = popSizeParameterInput.get();
        rhoParameter                    = rhoParameterInput.get();
        infectionRateParameter          = infectionRateParameterInput.get();
        deathRateParameter              = deathRateParameterInput.get();
        samplingRateParameter           = samplingRateParameterInput.get();
        R0Parameter                     = R0ParameterInput.get();
        becomeUninfectiousRateParameter = becomeUninfectiousRateParameterInput.get();
        samplingProportionParameter     = samplingProportionParameterInput.get();
        originParameter                 = originParameterInput.get();
        eventTimesParameter             = eventTimesParameterInput.get();
        eventTypesParameter             = eventTypesParameterInput.get();
        extantAtRootParameter           = extantAtRootParameterInput.get();

        // Get the tree or if eventTimes are specified ignore the tree
        if (eventTimesParameter == null) {
            tree = treeInput.get();
            if (tree == null) {
                throw new RuntimeException("Need to specify the tree or the eventTimes and eventTypes!");
            }
            totalEvents = tree.getLeafNodeCount()+tree.getInternalNodeCount();
        } else {
            totalEvents = eventTimesParameter.getDimension();

            if (eventTypesParameter == null) {
                throw new RuntimeException("If eventTimes are specified then eventTypes also need to be specified!");
            }

            System.err.println("Using eventTimes and eventTypes inputs - Tree input will be ignored! (not recommended)");
        }

        // Check if parameters need to be transformed or not
        if (infectionRateParameter != null &&
            deathRateParameter     != null &&
            samplingRateParameter  != null) {

            transform = false;
        } else
        if (R0Parameter                     != null &&
            becomeUninfectiousRateParameter != null &&
            samplingProportionParameter     != null) {

            transform = true;
        } else {
            throw new RuntimeException("Either specify infectionRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }

        // Check dimensions
        totalIntervals = popSizeParameter.getDimension();
        totalEvents   += totalIntervals;
        if ((!transform && (infectionRateParameter.getDimension()          != totalIntervals ||
                            deathRateParameter.getDimension()              != totalIntervals ||
                            samplingRateParameter.getDimension()           != totalIntervals))  ||
             (transform && (R0Parameter.getDimension()                     != totalIntervals ||
                            becomeUninfectiousRateParameter.getDimension() != totalIntervals ||
                            samplingProportionParameter.getDimension()     != totalIntervals))) {
            throw new RuntimeException("SIS model parameter dimensions not consistent");
        }

        if (eventTypesParameter != null && eventTypesParameter.getDimension() != totalEvents) {
            throw new RuntimeException("Wrong number of node types (expected "+totalEvents+")");
        }


        // Read constant parameter values (only need to be assigned once)
        if (survivalProbParameter.get() != null) {
            survival = survivalProbParameter.get().getValue();
        }

        if (rescaleParameter.get() != null) {
            rescaling = rescaleParameter.get().getValue();
        }

        if (extantAtRootParameter != null) {
            extantAtRoot = extantAtRootParameter.getValue();
        }

        // Assign memory (only necessary to do this once unless nr of intervals/events change)
        K           = new double[totalIntervals];
        beta        = new double[totalIntervals];
        mu          = new double[totalIntervals];
        psi         = new double[totalIntervals];

        eventTimes     = new double[totalEvents];
        eventTypes     = new int[totalEvents];
        eventTimesTemp = new double[totalEvents];
        eventTypesTemp = new int[totalEvents];
        if (survival) survTypes = new int[totalEvents];

        // Check that initial values are within absolute limits
        checkAbsoluteBounds();
    }



    /**
     * Check if all parameters are within bounds
     * These bounds may not be exceeded regardless of the limits set by the user!
     * (only checked at initialization)
     *
     * Set default limits for parameters if limits not specified already.
     * This should keep parameters within limits when being scaled
     * However still need to check manually that initial values are within limits
     *
     * rho    in [0,1]
     * K      in [max co-existing lineages, inf)
     *
     * beta   in (0,inf)
     * mu     in [0,inf)
     * psi    in [0,inf)
     *
     * R0     in (0,inf)
     * delta  in [0,inf)
     * s      in [0,1]
     *
     * extant in [0,inf)   (constant)
     *
     * @throws Exception
     */
    public void checkAbsoluteBounds() {

        // Check if the limits set by the user are sensible
        if (rhoParameter != null) {
            rhoParameter.setBounds(Math.max(0., rhoParameter.getLower()),
                                   Math.min(1., rhoParameter.getUpper()));
        }

        // TODO: Need to change this for shifts (maxExtant becomes an array)
        extantAtPresent = extantAtRoot;
        for (int i = totalEvents-1; i >= 0; --i) {
            switch (eventTypes[i]) {
                case 1:
                    ++extantAtPresent;
                    break;
                case 0:
                case 2:
                case 4:
                    --extantAtPresent;
                    break;
                default:
                    break;
            }
            if (maxExtant < extantAtPresent) maxExtant = extantAtPresent;
        }
        popSizeParameter.setLower(Math.max(maxExtant, popSizeParameter.getLower()));

        if (!transform) {
            infectionRateParameter.setLower(Math.max(0., infectionRateParameter.getLower()));
            deathRateParameter.setLower(Math.max(0., deathRateParameter.getLower()));
            samplingRateParameter.setLower(Math.max(0., samplingRateParameter.getLower()));
        } else {
            R0Parameter.setLower(Math.max(0.,R0Parameter.getLower()));
            becomeUninfectiousRateParameter.setLower(Math.max(0.,becomeUninfectiousRateParameter.getLower()));
            samplingProportionParameter.setLower(Math.max(0.,samplingProportionParameter.getLower()));
        }
        if (extantAtRootParameter != null)
            extantAtRootParameter.setLower(Math.max(0,extantAtRootParameter.getLower()));


        // Copy parameter values to primitives and check that no absolute bounds are exceeded
        if (updateRatesAndTimes(tree) >= 0) {

            if (rho < 0 || rho > 1) {
                throw new RuntimeException("Illegal parameter value (rho = "+rho+")");
            }

            for (int i = 0; i < totalIntervals; i++) {
                if (K[i] < 0 || beta[i] <= 0 ||  mu[i] < 0.0 || psi[i] < 0.0) {
                    throw new RuntimeException("Illegal parameter values\nK\tbeta\tmu\tpsi\n"+K[i]+"\t"+beta[i]+"\t"+mu[i]+"\t"+psi[i]);
                }
            }

            if (extantAtRoot < 0) {
                throw new RuntimeException("Illegal parameter value (extantAtRoot = "+extantAtRoot+")");
            }
        } else
            throw new RuntimeException("Illegal parameter value (origin > tree height)");
    }


    protected double updateRatesAndTimes(TreeInterface tree) {

        // TODO: Initial number of susceptibles
        int s0;

        // TODO: Default value for Rho?
        // TODO: Can rho be optimized or is it constant?
        // Convert SIS parameters to primitives
        rho      = (rhoParameter != null) ? rhoParameter.getValue() : 0;
        origin   = originParameter.getValue();
        if (origin < tree.getRoot().getHeight())
            return -1;

        if (transform) {
            for (int i = 0; i < totalIntervals; i++) {
                K[i]    = popSizeParameter.getValue(i);
                s0      = (int) Math.floor(K[i])-1;
                beta[i] = R0Parameter.getValue(i)*becomeUninfectiousRateParameter.getValue(i)*K[i]/s0;
                psi[i]  = samplingProportionParameter.getValue(i)*becomeUninfectiousRateParameter.getValue(i);
                mu[i]   = becomeUninfectiousRateParameter.getValue(i) - psi[i];
            }
        } else {
            for (int i = 0; i < totalIntervals; i++) {
                K[i]    = popSizeParameter.getValue(i);
                beta[i] = infectionRateParameter.getValue(i);
                mu[i]   = deathRateParameter.getValue(i);
                psi[i]  = samplingRateParameter.getValue(i);
            }
        }

        // Collect event times and types (using the tree is preferred)
        if (eventTimesParameter != null) {
            for (int i = 0; i < totalEvents; i++) {
                eventTimes[i] = eventTimesParameter.getValue(i);
                eventTypes[i] = eventTypesParameter.getValue(i);
                if (survival) survTypes[i] = (eventTypes[i] == 3) ? 3 : 99;
            }
            eventTimes[totalEvents-1] = origin;
            eventTypes[totalEvents-1] = BRANCH_EVENT;
            if (survival) survTypes[totalEvents-1] = 99;
        } else {
            collectEvents(tree);
        }

        if (verbose > 0) printModelState();

        return 0;
    }



    @Override
    public double calculateLogP() {
            final TreeInterface tree = treeInput.get();
            logP = calculateTreeLogLikelihood(tree);
            return logP;
    } // calculateLogP



    double calculateTreeLogLikelihood(TreeInterface tree) {

        int lkidx = 1+extantAtRoot;

        if (updateRatesAndTimes(tree) >= 0) {

            double [] loglk = DirtyExpoTree(K, beta, mu, psi, rho,
                                            eventTimes, eventTypes,
                                            false, extantAtRoot, rescaling, verbose);

            if (survival) {
                double [] survlk = DirtyExpoTree(K, beta, mu, psi, rho,
                                               eventTimes, survTypes,
                                               survival, extantAtRoot, rescaling, verbose);
                return loglk[lkidx] - survlk[1];
            } else
                return loglk[lkidx];

        } else
            return Double.NEGATIVE_INFINITY;
    }

    /**
     * Gather all of the event times and types from the tree and the parameter inputs
     *
     * TODO: Since the parameter shift times wont change and the tree has a dirty() function this should only be done if
     * the tree is dirty
     *
     * @param t
     */
    protected void collectEvents(TreeInterface t) {

        Comparator<Integer> eventTimeComparator = new Comparator<Integer>() {
            public int compare(Integer i, Integer j) {
                return Double.compare(eventTimesTemp[i], eventTimesTemp[j]);
            }
        };

        Integer [] indices = new Integer[totalEvents];
        int i = 0;

        // TODO: Get parameter shift events here

        // Read branching and sampling events from tree and add origin at the end
        for (Node n : t.getNodesAsArray()) {
            indices[i]   = i;
            eventTypesTemp[i] = (n.isLeaf()) ? SAMPLE_EVENT : BRANCH_EVENT;
            eventTimesTemp[i] = n.getHeight();
            i++;
        }
        indices[i]   = i;
        eventTypesTemp[i] = BRANCH_EVENT;
        eventTimesTemp[i] = origin;

        // Get times order and sort
        Arrays.sort(indices, eventTimeComparator);
        for (i = 0; i < totalEvents; i++) {
            eventTimes[i] = eventTimesTemp[indices[i]];
            eventTypes[i] = eventTypesTemp[indices[i]];
            if (survival) survTypes[i] = (eventTypes[i] == 3) ? 3 : 99;
        }

        //System.out.println("Times:");
        //for (double t : eventTimes) System.out.print(t+"\t");

        //System.out.println("\nTypes: ");
        //for (int t : eventTypes) System.out.print(t+"\t");

    }


    @Override
    /**
     * TODO: Dirty for now, force recalculation need to be made more elegant later
     */
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    /**
     * TODO: For now copied BDSKY, need to check later
     */
    public boolean canHandleTipDates() {
        return (rhoParameterInput.get() == null);
    }

    private void printModelState() {
        System.out.println("Model state");

        System.out.printf("\tintervals:\t\t%8d\n",totalIntervals);
        System.out.print("\tK:\t\t\t");
        for (int i = 0; i < K.length; i++) System.out.printf("\t%8.5f", K[i]);
        System.out.print("\n\tbeta:\t\t");
        for (int i = 0; i < beta.length; i++) System.out.printf("\t%8.5f", beta[i]);
        System.out.print("\n\tmu:\t\t\t");
        for (int i = 0; i < mu.length; i++) System.out.printf("\t%8.5f",mu[i]);
        System.out.print("\n\tpsi:\t\t");
        for (int i = 0; i < psi.length; i++) System.out.printf("\t%8.5f",psi[i]);
        System.out.printf("\n\trho:\t\t\t%8.3f",rho);

        System.out.printf("\n\tevents:\t\t\t%8d\n",totalEvents);
        System.out.print("\teventTimes:   ");
        for (int i = 0; i < eventTimes.length; i++) System.out.printf("\t%8.5f",eventTimes[i]);
        System.out.print("\n\teventTypes:   ");
        for (int i = 0; i < eventTypes.length; i++) System.out.printf("\t%8d",eventTypes[i]);
        if (survival) {
            System.out.print("\n\tsurvTypes:   ");
            for (int i = 0; i < survTypes.length; i++) System.out.printf("\t%8d",survTypes[i]);
        }
        System.out.printf("\n\ttextantAtRoot:\t%8d\n",extantAtRoot);
        System.out.printf("\textantAtPresent:%8d\n",extantAtPresent);
        System.out.printf("\tmaxExtant:\t\t%8d\n",maxExtant);

        //System.out.println("\tTree: "+tree.toString());

        System.out.println("\tsurvival:\t\t\t" +(survival ? "True" : "False"));
        System.out.println("\trescaling:\t\t\t"+(rescaling ? "True" : "False"));
        System.out.println("\tverbose:\t\t\t"+(verbose > 0 ? verbose : "False"));
    }

}

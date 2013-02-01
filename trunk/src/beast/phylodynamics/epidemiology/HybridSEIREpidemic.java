package beast.phylodynamics.epidemiology;

import java.util.List;
import java.util.Arrays;
import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Loggable;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Valuable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.parameter.BooleanParameter;
import beast.util.Randomizer;
import beast.evolution.tree.Tree;
import cern.jet.random.engine.MersenneTwister;

/**
 * @author dkuh004
 *         Date: Aug 22, 2011
 *         Time: 3:39:19 PM
 */
public class HybridSEIREpidemic extends CalculationNode implements Loggable {



    public Input<Valuable> S0 =
            new Input<Valuable>("S0", "The numbers of susceptible individuals", Validate.REQUIRED);
    public Input<IntegerParameter> E0 =
            new Input<IntegerParameter>("E0", "The numbers of exposed individuals");

    public Input<Tree> m_tree =
            new Input<Tree>("tree", "The phylogenetic tree on which to condition the trajectory",  Input.Validate.REQUIRED);
    public Input<RealParameter> orig_root =
            new Input<RealParameter>("orig_root", "The origin of infection x0", Input.Validate.REQUIRED);

    public Input<Integer> Nt = new Input<Integer>("Nt", "Number of timesteps in trajectory simulation", 10001);
    public Input<Integer> Nsamples = new Input<Integer>("Nsamples", "Number of samples = dimension of dS", Validate.REQUIRED);
    public Input<Double> alpha = new Input<Double>("alpha", "threshold parameter for the transition to the SSA", Validate.REQUIRED);

    public Input<RealParameter> exposeRate = new Input<RealParameter>("expose", "The rate at which individuals become exposed");
    public Input<Valuable> birthRateScalar = new Input<Valuable>("birth", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time", Validate.REQUIRED);
    public Input<Valuable> deathRateScalar = new Input<Valuable>("death", "The deathRate vector with deathRates between times", Validate.REQUIRED);
    public Input<Valuable> samplingRate = new Input<Valuable>("sampling", "The sampling rate per individual", Validate.REQUIRED);      // psi
    // the interval times for sampling rate
    public Input<RealParameter> samplingRateChangeTimesInput =
            new Input<RealParameter>("samplingRateChangeTimes", "The times t_i specifying when sampling rate or sampling proportion changes occur", (RealParameter) null);
    public Input<BooleanParameter> reverseTimeArrays =
            new Input<BooleanParameter>("reverseTimeArrays", "True if the time arrays are given in backwards time (from the present back to root). Order: 1) birth 2) death 3) sampling 4) rho. Default false." +
                    "Careful, rate array must still be given in FORWARD time (root to tips).");


    public Input<Boolean> useExposedBoolean = new Input<Boolean>("useExposed", "Boolean, system is SIR if false, beast.epidemiology.SEIR if true", false);

    public Input<String> simulationType = new Input<String>("simulationType", "which method to use for simulation: hybrid or SAL?", "hybrid", new String[]{"SAL", "hybrid"});

    // Trajectory parameters:
    Double expose;
    Double infect;
    Double[] recover;
    int e0;
    int i0 = 1;
    int r0 = 0;
    Boolean useExposed;

    public Integer[] dS;
    public Integer[] dE;
    public Integer[] dR;
    double tau;
    public double[] times;

    Integer[] store_dS;
    Integer[] store_dE;
    Integer[] store_dR;
    double store_tau;
    int breachCount;
    int noBreachCount;


    SEIR_simulator hybridTauleapSEIR;

    public void initAndValidate() throws Exception {

        useExposed = useExposedBoolean.get();
        if (useExposed && (E0.get()==null || exposeRate.get() == null) ) throw new Exception("HybridSEIREpidemic: When useExposed=true, exposeRate and E0 must be specified.");
        if (useExposed){
            expose = exposeRate.get().getValue();
            e0 = E0.get().getValue();
        }
        else{
            expose=0.;
            e0 = 0;
        }

        infect = birthRateScalar.get().getArrayValue();

        // allow the sampling rate to change over time
        recover = new Double[samplingRate.get().getDimension()];
        for (int i=0; i<recover.length; i++){

            recover[i] = deathRateScalar.get().getArrayValue(i) +  samplingRate.get().getArrayValue(i);
        }

        times=new double[recover.length-1];
        if (samplingRateChangeTimesInput.get()!= null) updateSamplingChangeTimes();

        breachCount = 0;
        noBreachCount = 0;

        initializeTrajectory(Nsamples.get(), alpha.get());
    }

    public void initializeTrajectory(int Nsamples, double alpha){

        SEIRState x0 = new SEIRState((int) S0.get().getArrayValue() - e0 - i0 - r0, e0, i0, r0, 0.0);

        if (simulationType.get().equals("hybrid"))
            hybridTauleapSEIR = new HybridTauleapSEIR(x0, expose, infect, recover[0], false, alpha);
        else
            hybridTauleapSEIR = new SALTauleapSEIR(x0, expose, infect, recover, false, alpha, new MersenneTwister((int)Randomizer.getSeed()));

        hybridTauleapSEIR.setState(x0);
//        List<SEIRState> trajectory = null;

        dS = new Integer[Nsamples-1];
        Arrays.fill(dS, 0);
        dE = new Integer[Nsamples-1];
        dR = new Integer[Nsamples-1];


    }

    void updateSamplingChangeTimes(){

        Double[] samplingChangeTimes = samplingRateChangeTimesInput.get().getValues();

        if (reverseTimeArrays.get()!=null && reverseTimeArrays.get().getValue(2)){
            //times are reversed
            for (int i=0; i<times.length; i++){
                times[i] = samplingChangeTimes[times.length-i-1];
            }
        } else {
            for (int i=0; i<times.length; i++){
                times[i] = samplingChangeTimes[i];
            }

        }

    }


    /**
     * @param S0 the initial number of susceptibles
     * @param Nsamples the number of samples to be taken from the time series
     * @param ntaxa number of taxa in phylogeny
     * @param Nt  Number of time-steps
     * @param useExposed
     * @param T
     * @param maxLoop
     * @param times
     * @return if the generation of a time series was succesful
     */
     public Boolean generateTrajectory(int S0, int Nsamples, int ntaxa, int Nt, Boolean useExposed, double T, int maxLoop, double[] times) {

        SEIRState x0 = new SEIRState(S0 - e0 - i0 - r0, e0, i0, r0, 0.0);

        List<SEIRState> trajectory = null;

        int loopCount = 0;

        while (loopCount < maxLoop && trajectory == null ) {

            try{
                loopCount++;
                hybridTauleapSEIR.setState(x0);
                trajectory = null;
                trajectory = hybridTauleapSEIR.genTrajectory(T, Nt, Nsamples, ntaxa, true, times);

            } catch (Exception e){
                /*System.out.println(e.getMessage());*/}
        }



        if (trajectory == null )  {
            return false;
        }

        dS = new Integer[Nsamples];
        dE = new Integer[Nsamples];
        dR = new Integer[Nsamples];

        SEIRState former = trajectory.get(0);

        for (int i = 0; i < Nsamples; i++) {
            SEIRState current = trajectory.get(i+1);

            dS[i] = former.S - current.S;

            if (useExposed)
                dE[i] = current.E - former.E;

            dR[i] = current.R - former.R;

            former = current;
        }

        return true;

    }

    @Override
    public void store() {
        super.store();
        if (dS != null){
            store_dS = dS.clone();
            store_dE = dE.clone();
            store_dR = dR.clone();
            store_tau = tau;
        }
    }


    /** Restore internal calculations
     *
     * This is called when a proposal is rejected
     **/
    @Override
    public void restore() {
        super.restore();
        if (store_dS != null){

            final Integer[] tmpS;
            final Integer[] tmpE;
            final Integer[] tmpR;

            tmpS =  store_dS;
            tmpE =  store_dE;
            tmpR =  store_dR;

            store_dS = dS;
            store_dE = dE;
            store_dR = dR;

            dS = tmpS;
            dE = tmpE;
            dR = tmpR;

            final double tmp = store_tau;
            store_tau = tau;
            tau = tmp;

        }
    }


    //
    public Boolean init(Double birth,Double expose, Double[] death, Double[] psi, double T, int ntaxa, int maxLoop, int intervals, double[] times){

        return refresh((int) S0.get().getArrayValue(), useExposed ? E0.get().getValue(): 0, birth, useExposed?expose:0., death, psi, T, ntaxa, maxLoop, intervals, times);

    }


    public Boolean refresh(int S0, int E0, Double birth,Double exposed, Double[] death, Double[] psi, double T, int ntaxa, int maxLoop, int intervals, double[] times){

        tau = T / intervals;

        infect = birth;

        // allow the sampling rate to change over time
        recover = new Double[psi.length];

        for (int i=0; i<recover.length; i++){

            recover[i] = death[i] +  psi[i];
        }

        if (samplingRateChangeTimesInput.get()!= null) updateSamplingChangeTimes();

        if (useExposed){
            expose = exposed;
            e0 = E0;
        }

        hybridTauleapSEIR.setRates(expose, infect, recover, alpha.get());

        return generateTrajectory(S0, intervals, ntaxa, Nt.get(), useExposed, T, maxLoop, times);

    }

    public Integer[] get_dS(){
        return dS;
    }

    public Integer[] get_dE(){
        return dE;
    }

    public Integer[] get_dR(){
        return dR;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    public void init(PrintStream printStream) throws Exception {
        printStream.print("breachCount\tnoBreachCount\tnoBreachProportion\t");
    }

    public void log(int i, PrintStream printStream) {
        printStream.print(breachCount+"\t"+noBreachCount+"\t"+(noBreachCount/((breachCount+noBreachCount)*1.))+"\t");
    }

    public void close(PrintStream printStream) {
        // nothing to do
    }
}

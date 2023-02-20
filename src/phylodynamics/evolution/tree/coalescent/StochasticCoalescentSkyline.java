package phylodynamics.evolution.tree.coalescent;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.util.Binomial;

/**
 * User: Denise Kuehnert
 * Date: 10.11.13
 * Time: 16:49
 */

@Description("Coalescent skyline plot coupled with forward simulated SIR epidemiological model (aka the coalescent BDSIR)")
public class StochasticCoalescentSkyline extends TreeDistribution {

    public Input<RealParameter> susceptiblePopulationInput = new Input<RealParameter>("susceptiblePopulation", "susceptible population size at times eventTimes",
            Input.Validate.REQUIRED);

    public Input<RealParameter> infectedPopulationInput = new Input<RealParameter>("infectedPopulation", "infected population size at times eventTimes",
            Input.Validate.REQUIRED);

    public Input<RealParameter> transmissionRateInput = new Input<RealParameter>("transmissionRate", "the transmission rate β",
            Input.Validate.REQUIRED);

    public Input<RealParameter> eventTimeInput = new Input<RealParameter>("eventTimes", "event times from first infection to present (i.e. in forward time)",
            Input.Validate.REQUIRED);

    TreeIntervals intervals;
    double transmission;
    Double[] susceptiblePopulation;
    Double[] infectedPopulation;
    Double[] eventTimes;
    Double[] S;
    Double[] I;
    Boolean[] isCoalescent;
    double[] coalescentTimes;
    double[] times;

    int dimension;

    boolean m_bIsPrepared = false;


    public void initAndValidate() {
        if (treeInput.get() != null) {
            throw new RuntimeException("only tree intervals (not tree) should not be specified");
        }
        intervals = treeIntervalsInput.get();

        if (infectedPopulationInput.get().getDimension()!=eventTimeInput.get().getDimension()) {
            throw new RuntimeException("infectedPopulation size input should have same dimension as eventTime input");
        }

        prepare();
    }

    public void prepare() {

        transmission = transmissionRateInput.get().getValue();

        coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);
        eventTimes =  eventTimeInput.get().getValues();
        infectedPopulation = infectedPopulationInput.get().getValues();
        susceptiblePopulation = susceptiblePopulationInput.get().getValues();

        dimension = eventTimes.length + coalescentTimes.length;

        S = new Double[dimension];
        I = new Double[dimension];
        times = new double[dimension];
        isCoalescent = new Boolean[dimension];

        int eventIndex = eventTimes.length-1;
        int coalescentIndex = 0;

        double T = Math.max(eventTimes[eventIndex],intervals.getTotalDuration());

        for (int i=0; i<dimension && coalescentIndex<coalescentTimes.length; i++) {

            if (eventIndex>=0 && (T-eventTimes[eventIndex])<coalescentTimes[coalescentIndex]){

                times[i] = T-eventTimes[eventIndex];
                I[i] = infectedPopulation[eventIndex];
                S[i] = susceptiblePopulation[eventIndex];
                isCoalescent[i] = false;
                eventIndex--;

            }
            else {
                times[i] = coalescentTimes[coalescentIndex];
                S[i] = susceptiblePopulation[Math.min(eventIndex+1,eventTimes.length-1)];
                I[i] = infectedPopulation[Math.min(eventIndex+1,eventTimes.length-1)];

                isCoalescent[i] = true;

                if(eventIndex>=0 && eventTimes[eventIndex]==coalescentTimes[coalescentIndex]) {
                    I[i] = infectedPopulation[eventIndex];
                    eventIndex--;
                    dimension--;
                }
                coalescentIndex++;

            }
        }

        m_bIsPrepared = true;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    protected boolean requiresRecalculation() {
        m_bIsPrepared = false;
        return true;
    }

    @Override
    public void store() {
        m_bIsPrepared = false;
        super.store();
    }

    @Override
    public void restore() {
        m_bIsPrepared = false;
        super.restore();
    }


    /**
     * Calculates the log likelihood of this set of coalescent intervals
     */
    @Override
    public double calculateLogP() {
        if (!m_bIsPrepared) {
            prepare();
        }

        logP = 0.0;

        double currentTime = 0.0;
        double time;

        for (int j = 1 ; j < times.length && currentTime<intervals.getTotalDuration(); j++) {

            time = times[j] - currentTime;

            logP += calculateIntervalLikelihood((2 * transmission * (S[j-1]+S[j])/2.) / (I[j-1]+I[j])/2., time,
                    lineageCountAtTime(currentTime,intervals.treeInput.get()), isCoalescent[j]);

            currentTime += time;
        }
        return logP;
    }

    public static double calculateIntervalLikelihood(double coalescentRate, double width,
                                                     int lineageCount, Boolean isCoalescent) {

        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * width * coalescentRate; //exponential waiting time

        if (isCoalescent!=null && isCoalescent)  {

            like += Math.log(coalescentRate);
        }

        return like;
    }


    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, TreeInterface tree) {

        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;

        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() > time) count -= 1;
        }
        return count;
    }

}


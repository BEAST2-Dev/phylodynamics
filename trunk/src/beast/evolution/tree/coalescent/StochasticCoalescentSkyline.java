package beast.evolution.tree.coalescent;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.math.Binomial;

/**
 * User: Denise
 * Date: 10.11.13
 * Time: 16:49
 */
public class StochasticCoalescentSkyline extends TreeDistribution {

    public Input<RealParameter> infectedPopulationInput = new Input<RealParameter>("infectedPopulation", "infected population size at times eventTimes",
            Input.Validate.REQUIRED);

    public Input<RealParameter> eventTimeInput = new Input<RealParameter>("eventTimes", "event times from first infection to present (i.e. in forward time)",
            Input.Validate.REQUIRED);

    TreeIntervals intervals;
    Double[] infectedPopulation;
    Double[] eventTimes;
    Double[] I;
    Boolean[] isCoalescent;
    double[] coalescentTimes;
    double[] times;

    int dimension;

    boolean m_bIsPrepared = false;

    // This pseudo-constructor is only used for junit tests
    public StochasticCoalescentSkyline() {
    }

    public void initAndValidate() throws Exception {
        if (treeInput.get() != null) {
            throw new Exception("only tree intervals (not tree) should not be specified");
        }
        intervals = treeIntervalsInput.get();

        if (infectedPopulationInput.get().getDimension()!=eventTimeInput.get().getDimension()) {
            throw new Exception("infectedPopulation size input should have same dimension as eventTime input");
        }

        prepare();
    }

    public void prepare() {

        coalescentTimes = intervals.getCoalescentTimes(coalescentTimes);
        eventTimes =  eventTimeInput.get().getValues();
        infectedPopulation = infectedPopulationInput.get().getValues();

        dimension = eventTimes.length + coalescentTimes.length;

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
                isCoalescent[i] = false;
                eventIndex--;

            }
            else {
                times[i] = coalescentTimes[coalescentIndex];
                I[i] = infectedPopulation[Math.min(eventIndex+1,eventTimes.length-1)];

                isCoalescent[i] = intervals.getIntervalType(coalescentIndex)==IntervalType.COALESCENT;

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
    public double calculateLogP() throws Exception {
        if (!m_bIsPrepared) {
            prepare();
        }

        logP = 0.0;

        double currentTime = 0.0;
        double time;

        for (int j = 1 ; j < times.length && currentTime<=intervals.getTotalDuration(); j++) {

            time = times[j];

            logP += calculateIntervalLikelihood(I[j-1], time, currentTime,
                    lineageCountAtTime(time,intervals.treeInput.get()), isCoalescent[j]); // todo: change I[j-1] to (beta S I)/(I^2)!!! (like in Volz)

            currentTime += time;
        }
        return logP;
    }

    public static double calculateIntervalLikelihood(double popSize, double width,
                                                     double timeOfPrevCoal, int lineageCount, Boolean isCoalescent) {

        final double timeOfThisCoal = width + timeOfPrevCoal;

        final double intervalArea = (timeOfThisCoal - timeOfPrevCoal) / popSize;

        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * intervalArea;

        if (isCoalescent)  {
            final double demographic = Math.log(popSize);//demogFunction.getLogDemographic(timeOfThisCoal);
            like += -demographic;
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


package beast.phylodynamics.epidemiology;


import java.util.*;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;

import beast.evolution.tree.coalescent.PopulationFunction;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.evolution.tree.coalescent.IntervalList;
import beast.evolution.tree.coalescent.IntervalType;

/**
 * @author Alexei Drummond
 * @author Alex Popinga
 */

@Description("Calculates the probability of a beast.tree conditional on a population size function. " +
        "Note that this does not take the number of possible tree interval/tree topology combinations " +
        "in account, in other words, the constant required for making this a proper distribution that integrates " +
        "to unity is not calculated (partly, because we don't know how for sequentially sampled data).")
public class StochasticCoalescent extends TreeDistribution {

    public Input<SIRPopulationFunction> popSizeInput = new Input<SIRPopulationFunction>("populationModel", "A population size model", Validate.REQUIRED);
    public Input<Integer> minimumNumberOfTrajectories = new Input<Integer>("minTraj", "The minimum number of trajectories over which to average the coalescent probability, i.e. a value of 1 means one trajectory is used, (defaults 1).", 1);
    public Input<Integer> minimumNumberOfSuccesses = new Input<Integer>("minTrajSuccess",
            "minimum number of trajectories that must span the entire tree. (i.e. default is 1)", 1);
    public Input<Integer> maxTries = new Input<Integer>("maxTries",
            "maximum number of trajectory simulations attempted. (i.e. default is 10000). Reject if reached.", 1000);


    TreeIntervals treeIntervals;

    public int minTraj = 1;
    public int minTrajSuccess = 1;

    private int lastEnsembleSize = 0;

    public StochasticCoalescent() {}

    /**
     *
     * @param treeIntervals the tree intervals to calculate coalescent probability for
     * @param populationModel the population function
     * @param minEnsembleSize the minimum number of trajectories to simulation from the population function
     * @param minSuccessfulTraj the minimum number of "successful" trajectories, i.e that have population sizes that are strictly positive over the whole tree.
     * @throws Exception
     */
    public StochasticCoalescent(TreeIntervals treeIntervals, SIRPopulationFunction populationModel, int minEnsembleSize, int minSuccessfulTraj) throws Exception {
        initByName("treeIntervals", treeIntervals, "populationModel", populationModel, "minTraj", minEnsembleSize, "minTrajSuccess", minSuccessfulTraj);
    }

    @Override
    public void initAndValidate() throws Exception {
        treeIntervals = treeIntervalsInput.get();
        if (treeIntervals == null) {
            throw new Exception("Expected treeIntervals to be specified");
        }

        minTraj = minimumNumberOfTrajectories.get();
        minTrajSuccess = minimumNumberOfSuccesses.get();
        calculateLogP();
    }

    /**
     * @return the size of the ensemble used to calculate the last probability density (i.e. last call to calculateLogP).
     */
    public int getLastEnsembleSize() {
        return lastEnsembleSize;
    }

    /**
     * do the actual calculation *
     * 
     * @return log of HR
     * @throws java.lang.Exception
     */
    @Override
    public double calculateLogP() throws Exception {

        SIRPopulationFunction popFunction = popSizeInput.get();

        // if stochastic, average log-likelihoods
        if (popFunction.volzSIR.get() instanceof StochasticSIR) {

            ArrayList<Double> logps = new ArrayList<Double>();

            int numTraj = 0;
            int maxTriesBeforeReject = maxTries.get();

            while (numTraj < minTraj || logps.size() < minTrajSuccess) {
                boolean fail = popFunction.simulateTrajectory();

                double logp = Double.NEGATIVE_INFINITY;
                if (!fail) {
                    logp = calculateLogLikelihood(treeIntervals, popFunction);
                }

                if (!Double.isNaN(logp) && !Double.isInfinite(logp)) {
                    logps.add(logp);
                }
                numTraj += 1;

                if (numTraj >= maxTriesBeforeReject) {
                    // skip the party early, no luck here
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }
            }
            lastEnsembleSize = numTraj;

            if (logps.isEmpty()) return Double.NEGATIVE_INFINITY;

            // find maximum likelihood, make it the 'new zero'
            double ML = Collections.max(logps);

            Collections.sort(logps);
            double sum = 0.0;

            // shift array elements according to newML and exponentiate
            for (int j = 0; j < logps.size(); j++) {
                logps.set(j, Math.exp(logps.get(j) - ML));

                double logp = logps.get(j);
                sum += logp;
            }

            // average probability
            double logAve = Math.log(sum/lastEnsembleSize);

            logP = logAve + ML;

        } else {
            logP = calculateLogLikelihood(treeIntervals, popFunction);
        }

        if (Double.isInfinite(logP)) {
            logP = Double.NEGATIVE_INFINITY;
        }

        return logP;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO this should eventually sample a coalescent tree conditional on population size function
        throw new UnsupportedOperationException("This should eventually sample a coalescent tree conditional on population size function.");
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(treeIntervalsInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return popSizeInput.get().getParameterIds();
    }


    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     *
     * @param intervals       the intervals whose likelihood is computed
     * @param popSizeFunction the population size function
     * @return the log likelihood of the intervals given the population size function
     */
    public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction) {
        return calculateLogLikelihood(intervals, popSizeFunction, 0.0);
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a population size function.
     *
     * @param intervals       the intervals whose likelihood is computed
     * @param popSizeFunction the population size function
     * @param threshold       the minimum allowable coalescent interval size; negative infinity will be returned if
     *                        any non-zero intervals are smaller than this
     * @return the log likelihood of the intervals given the population size function
     */
    public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {

        double logL = 0.0;

        double startTime = 0.0;
        final int n = intervals.getIntervalCount();
        for (int i = 0; i < n; i++) {

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            final double intervalArea = popSizeFunction.getIntegral(startTime, finishTime);
            if (intervalArea == 0 && duration != 0) {
                return Double.NEGATIVE_INFINITY;
            }
            final int lineageCount = intervals.getLineageCount(i);

            final double kChoose2 = Binomial.choose2(lineageCount);
            // common part
            logL += -kChoose2 * intervalArea;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                final double demographicAtCoalPoint = popSizeFunction.getPopSize(finishTime);

                // if value at end is many orders of magnitude different than mean over interval reject the interval
                // This is protection against cases where ridiculous infinitesimal
                // population size at the end of a linear interval drive coalescent values to infinity.

                if (duration == 0.0 || demographicAtCoalPoint * (intervalArea / duration) >= threshold) {
                    //                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
                    logL -= Math.log(demographicAtCoalPoint);
                } else {
                    // remove this at some stage
                    //  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
                    return Double.NEGATIVE_INFINITY;
                }
            }
            startTime = finishTime;
        }

        return logL;
    }

    @Override
    protected boolean requiresRecalculation() {
        return ((CalculationNode) popSizeInput.get()).isDirtyCalculation() || super.requiresRecalculation();
    }

    @Override
    public boolean isStochastic() {
        return popSizeInput.get().volzSIR.get() instanceof StochasticSIR;
    }
}

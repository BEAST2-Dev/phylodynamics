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

    public Input<PopulationFunction> popSizeInput = new Input<PopulationFunction>("populationModel", "A population size model", Validate.REQUIRED);
    public Input<Integer> REP = new Input<Integer>("REP", "Average probabilities, i.e. a value of 1 means no averaging is done, (defaults 1).", 1);
    public Input<Boolean> allowNaN = new Input<Boolean>("allowNaN", "Indicate whether to allow NaN likelihoods (defaults false).", false);

    TreeIntervals treeIntervals;

    public int REPS = REP.get();
    public boolean allowNaNs = allowNaN.get();

    public StochasticCoalescent() {}

    public StochasticCoalescent(TreeIntervals treeIntervals, PopulationFunction populationModel) throws Exception {
        initByName("treeIntervals", treeIntervals, "populationModel", populationModel);
    }

    @Override
    public void initAndValidate() throws Exception {
        treeIntervals = treeIntervalsInput.get();
        if (treeIntervals == null) {
            throw new Exception("Expected treeIntervals to be specified");
        }
        calculateLogP();
    }


    /**
     * do the actual calculation *
     */
    @Override
    public double calculateLogP() throws Exception {

        PopulationFunction popFunction = popSizeInput.get();

        // if stochastic, average log-likelihoods
        if (popFunction instanceof StochasticSIRPopulationFunction) {

            StochasticSIRPopulationFunction scSIR = (StochasticSIRPopulationFunction)popFunction;

            ArrayList<Double> logps = new ArrayList<Double>();

            int failCount = 0;
            while (logps.size() < REPS) {
                scSIR.simulateStochasticTrajectory();

                double logp = calculateLogLikelihood(treeIntervals, popSizeInput.get());

                if (allowNaNs || (!Double.isNaN(logp) && !Double.isInfinite(logp))) {
                    logps.add(logp);
                } else {
                    failCount += 1;
                }
            }

            // find maximum likelihood, make it the 'new zero'
            double ML = Collections.max(logps);
            double newML = 0.0 - ML;

            Collections.sort(logps);
            double ave = 0.0;

            // shift array elements according to newML and exponentiate all the things
            for (int j = 0; j < logps.size(); j++) {
                logps.set(j, Math.exp(logps.get(j) + newML));

                double logp = logps.get(j);
                ave += logp;
            }

            // average probability
            double logAve = Math.log(ave/logps.size());
            logP = logAve + ML;

            if (failCount > 0) {
                // do something special to correct for the proportion of failures
                // someone with more brains than AJD should check this

                double probSurvivalOfTrajectory = (double)REPS/(double)(failCount+REPS);

                logP += Math.log(probSurvivalOfTrajectory);
            }


        } else {
            logP = calculateLogLikelihood(treeIntervals, popSizeInput.get());
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
}

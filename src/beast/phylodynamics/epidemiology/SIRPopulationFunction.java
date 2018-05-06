package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.util.Arrays;
import java.util.List;


/**
 * @author Timothy Vaughan
 * @author Alex Popinga
 */

@Description("Population function for performing likelihood calculations based on Volz 2012.")
public class SIRPopulationFunction extends PopulationFunction.Abstract {


    public Input<CoalescentSIR> volzSIR = new Input<CoalescentSIR>("volzSIR", "Volz parameters for coalescent SIR");

    public Input<Boolean> oldMethodInput = new Input<Boolean>(
            "oldMethod",
            "Use old (slow) method to evaluate intensity.  Default false.", false);

    public SIRPopulationFunction() {
    }

    public SIRPopulationFunction(CoalescentSIR ssir) throws Exception {
        initByName("volzSIR", ssir);
    }

    /**
     * Simulate a new trajectory, irrespective of dirtiness.
     * 
     * @return true if the simulation failed
     */
    public boolean simulateTrajectory() {
        return volzSIR.get().simulateTrajectory();
    }

    // Implementation of abstract methods
    @Override
    public List<String> getParameterIds() {

        String[] parameterIds = new String[]{
                volzSIR.get().R0.get().getID(),
                volzSIR.get().originParameter.get().getID(),
                volzSIR.get().n_S_Parameter.get().getID(),
                volzSIR.get().gammaParameter.get().getID(),
        };

        return Arrays.asList(parameterIds);
    }

    /**
     * "Effective" population size for calculating coalescent likelihood
     * under Volz's model.  This is only used if numerical integration
     * of the intensity is employed.
     *
     * @param t
     * @return NI(t)/(2*beta*NS(t))
     */
    @Override
    public double getPopSize(double t) {

        double popSize = volzSIR.get().getPopSize(t);
        if (!volzSIR.get().reject)
            return popSize;
        else
            return Double.POSITIVE_INFINITY; // Causes -infinity coalescent prob.
    }


    @Override
    public double getIntegral(double start, double finish) {
        
        double res;
        
        // Approximation for the case of very small intervals:
        if (finish - start < 1e-15)
            res = (finish - start) / getPopSize(0.5 * (start + finish));
        else {
            if (oldMethodInput.get())
                res = getNumericalIntegral(start, finish);
            else
                res = super.getIntegral(start, finish);
        }
        
        if (!volzSIR.get().reject)
            return res;
        else
            return Double.POSITIVE_INFINITY; // causes -infinity coalescent prob.
    }

    /**
     * Uses pre-calculated intensity values to estimate the intensity
     * at time t. For times within the domain of the calculated values,
     * linear interpolation is used.  For times outside this range, the
     * initial or final effectivePopSizeTraj value is used to extrapolate.
     *
     * @param t
     * @return
     */

    @Override
    public double getIntensity(double t) {
        return volzSIR.get().getIntensity(t);
    }

    @Override
    public double getInverseIntensity(final double intensity) {
        return volzSIR.get().getInverseIntensity(intensity);
    }

}
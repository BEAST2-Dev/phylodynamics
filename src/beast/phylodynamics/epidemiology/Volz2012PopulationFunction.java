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
public class Volz2012PopulationFunction extends PopulationFunction.Abstract {

    public Input<VolzSIR> volzSIRInput = new Input<VolzSIR>("volzSIR", "Volz parameters", Input.Validate.REQUIRED);

    public Input<Boolean> oldMethodInput = new Input<Boolean>(
            "oldMethod",
            "Use old (slow) method to evaluate intensity.  Default false.", false);



    // Implementation of abstract methods
    public List<String> getParameterIds() {

        String[] parameterIds = new String[]{
                volzSIRInput.get().betaParameter.get().getID(),
                volzSIRInput.get().originParameter.get().getID(),
                volzSIRInput.get().n_S_Parameter.get().getID(),
                volzSIRInput.get().gammaParameter.get().getID(),

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
          return volzSIRInput.get().getPopSize(t);
    }

    @Override
    public double getIntegral(double start, double finish) {

        // Approximation for the case of very small intervals:

        if (finish - start < 1e-15)
            return (finish-start)/getPopSize(0.5*(start+finish));

        if (oldMethodInput.get())
            return getNumericalIntegral(start, finish);
        else
            return super.getIntegral(start, finish);
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
        if (volzSIRInput.get().reject == true) {
            return Double.NEGATIVE_INFINITY;
        }
        else return volzSIRInput.get().getIntensity(t);
    }

    @Override
    public double getInverseIntensity(final double intensity) {
        return volzSIRInput.get().getInverseIntensity(intensity);
    }

}
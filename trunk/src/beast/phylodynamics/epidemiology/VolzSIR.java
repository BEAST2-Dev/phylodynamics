package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BrentSolver;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Timothy Vaughan
 * @author Alexei Drummond
 */
@Description("Population function based on Volz 1999 `coalescent' likelihood.")
public class VolzSIR extends PopulationFunction.Abstract implements Loggable {
    public Input<RealParameter> n_S_Parameter = new Input<RealParameter>("n_S0",
            "the number of susceptibles at time of origin (defaults to 1000). ");
    public Input<RealParameter> betaParameter = new Input<RealParameter>("beta",
            "the mass action rate of infection.");
    public Input<RealParameter> gammaParameter = new Input<RealParameter>("gamma",
            "the per-infected rate of recovery.");
    public Input<RealParameter> originParameter = new Input<RealParameter>("origin",
            "the time before the root that the first infection occurred.");

    //public Input<Double> integrationStepInput = new Input<Double>("integrationStep",
    //        "length of integration time step.", Validate.REQUIRED);
    public Input<Integer> integrationStepCount = new Input<Integer>("integrationStepCount",
            "number of integration time steps to use.", Validate.REQUIRED);
    public Input<Double> finishingThresholdInput = new Input<Double>("finishingThreshold",
            "Integration will finish when infected pop drops below this.", 1.0);

    public Input<Integer> statesToLogInput = new Input<Integer>("statesToLog",
            "Number of states to log. (Default 100.)", 100);

    public Input<Boolean> oldMethodInput = new Input<Boolean>(
            "oldMethod",
            "Use old (slow) method to evaluate intensity.  Default false.", false);

    //
    // Private stuff
    //

    private boolean dirty;
    private List<Double> effectivePopSizeTraj, intensityTraj;
    private List<Double> NStraj, NItraj;
    private double tIntensityTrajStart;

    private double dt;
    private int Nt;

    //
    // Public stuff
    //

    @Override
    public void initAndValidate() throws Exception {
        if (betaParameter.get() != null) {
            betaParameter.get().setBounds(
                    Math.max(0.0, betaParameter.get().getLower()),
                    betaParameter.get().getUpper());
        }
        if (gammaParameter.get() != null) {
            gammaParameter.get().setBounds(
                    Math.max(0.0, gammaParameter.get().getLower()),
                    gammaParameter.get().getUpper());
        }
        if (originParameter.get() != null) {
            originParameter.get().setBounds(
                    Math.max(0.0, originParameter.get().getLower()),
                    originParameter.get().getUpper());
        }
        if (n_S_Parameter.get() != null) {
            n_S_Parameter.get().setBounds(
                    Math.max(1.0, n_S_Parameter.get().getLower()),
                    n_S_Parameter.get().getUpper());
        }

        effectivePopSizeTraj = new ArrayList<Double>();
        intensityTraj = new ArrayList<Double>();
        NStraj = new ArrayList<Double>();
        NItraj = new ArrayList<Double>();

        Nt = integrationStepCount.get();
        //dt = integrationStepInput.get();

        dirty = true;
        update();

    }

    /**
     * Update deterministic trajectory.
     */
    private void update() {
        if (!dirty)
            return;

        // Short-hand for model parameters:
        double beta = betaParameter.get().getValue();
        double gamma = gammaParameter.get().getValue();
        double NS0 = n_S_Parameter.get().getValue();

        double threshNI = finishingThresholdInput.get();


        if (gamma < beta * NS0) {
            // Estimate length of epidemic:
            double elEstimate = Math.log(NS0) * (1.0 / (beta * NS0 - gamma) + 1.0 / gamma);

            // Use estimate to choose step size:
            dt = elEstimate / Nt;
        } else
            // Insurance against outrageous parameter combinations.
            dt = 0.1;

        // Clear old trajectory
        NStraj.clear();
        NItraj.clear();
        effectivePopSizeTraj.clear();

        // Set up initial conditions:
        double NS = NS0;
        double NI = 1.0;
        double effectivePopSize = NI / (2.0 * beta * NS);

        NStraj.add(NS);
        NItraj.add(NI);
        effectivePopSizeTraj.add(effectivePopSize);

        // Integrate trajectory:
        double dNSdt, dNIdt = 0.0;
        do {

            double NSmid = NS;
            double NImid = NI;
            for (int iter = 0; iter < 3; iter++) {
                dNSdt = -beta * NSmid * NImid;
                dNIdt = beta * NSmid * NImid - gamma * NImid;

                NSmid = NS + 0.5 * dt * dNSdt;
                NImid = NI + 0.5 * dt * dNIdt;
            }
            NS = 2.0 * NSmid - NS;
            NI = 2.0 * NImid - NI;

            effectivePopSize = NI / (2.0 * beta * NS);

            NStraj.add(NS);
            NItraj.add(NI);
            effectivePopSizeTraj.add(effectivePopSize);
        } while (dNIdt > 0 || NI > threshNI);

        // Switch effective pop size to reverse time:
        Collections.reverse(effectivePopSizeTraj);

        // Estimate intensity on integration lattice:
        intensityTraj.clear();

        double intensity = 0.0;
        intensityTraj.add(intensity);
        for (int i = 0; i < effectivePopSizeTraj.size(); i++) {
            intensity += dt / effectivePopSizeTraj.get(i);
            intensityTraj.add(intensity);
        }

        // Start of integral is 0.5*dt from end of forward-time integration.
        tIntensityTrajStart = originParameter.get().getValue()
                - dt * (effectivePopSizeTraj.size() - 1) - 0.5 * dt;

        dirty = false;
    }


    // Implementation of abstract methods

    @Override
    public List<String> getParameterIds() {

        String[] parameterIds = new String[]{
                betaParameter.get().getID(),
                originParameter.get().getID()

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
        update();

        // Choose which index into integration lattice to use:        
        int tidx = (int) Math.floor((t - tIntensityTrajStart) / dt);

        // Use initial or final state of trajectory if t outside the bounds of the
        // simulation.  This is a CLUDEGE to deal with trees which don't fit
        // the trajectories at all.
        if (tidx >= effectivePopSizeTraj.size())
            return 1e-10 * effectivePopSizeTraj.get(effectivePopSizeTraj.size() - 1);
        else if (tidx < 0)
            return effectivePopSizeTraj.get(0);

        return effectivePopSizeTraj.get(tidx);
    }

    @Override
    public double getIntegral(double start, double finish) {

        if (start == finish)
            return 0.0;

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
        update();

        if (t < tIntensityTrajStart) {
            return -(tIntensityTrajStart - t) / effectivePopSizeTraj.get(0);
        } else {
            if (t > originParameter.get().getValue() + 0.5 * dt) {
                return intensityTraj.get(intensityTraj.size() - 1)
                        + (t - originParameter.get().getValue() + 0.5 * dt) * 1e10;//intensityTraj.get(intensityTraj.size()-1)
//                        + (t-(originParameter.get().getValue()+0.5*dt))
//                        /effectivePopSizeTraj.get(effectivePopSizeTraj.size()-1);
            } else {
                int idx = (int) Math.floor((t - tIntensityTrajStart) / dt);
                double alpha = (t - tIntensityTrajStart - dt * idx) / dt;
                return intensityTraj.get(idx) * (1.0 - alpha)
                        + intensityTraj.get(idx + 1) * alpha;
            }
        }
    }

    @Override
    public double getInverseIntensity(final double intensity) {

        // given a value for the intensity, find a time

        UnivariateRealFunction intensityFunction = new UnivariateRealFunction() {
            public double value(double x) throws FunctionEvaluationException {
                return getIntensity(x) - intensity;
            }
        };

        BrentSolver solver = new BrentSolver();
        try {
            return solver.solve(intensityFunction, 0, n_S_Parameter.get().getValue() / betaParameter.get().getValue());
        } catch (MaxIterationsExceededException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException("Max iterations (" + e.getMaxIterations() + ") exceeded:" + e.getMessage());
        } catch (FunctionEvaluationException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException(e.getMessage());
        }
    }

    /*
     * CalculationNode interface
     */

    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    public void restore() {
        dirty = true;
        super.restore();
    }
    
    /*
     * Loggable interface
     */

    @Override
    public void init(PrintStream out) throws Exception {
        out.print("dt\t");
        for (int i = 0; i < statesToLogInput.get(); i++) {
            out.format("S%d\t", i);
            out.format("I%d\t", i);
            out.format("t%d\t", i);
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {

        out.format("%g\t", dt);

        double tend = NStraj.size() * dt;

        double delta = tend / (statesToLogInput.get() - 1);

        for (int i = 0; i < statesToLogInput.get(); i++) {

            double t = delta * i;
            int tidx = (int) Math.round(t / dt);

            if (tidx >= NStraj.size())
                tidx = NStraj.size() - 1;

            out.format("%g\t", NStraj.get(tidx));
            out.format("%g\t", NItraj.get(tidx));
            out.format("%g\t", t);
        }
    }

    @Override
    public void close(PrintStream out) {
    }
}
package beast.phylodynamics.epidemiology;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NumberIsTooSmallException;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;


/**
 * @author Timothy Vaughan
 * @author Alexei Drummond
 * @author Alex Popinga
 */
public abstract class VolzSIR extends CalculationNode implements Loggable {

    public Input<RealParameter> n_S_Parameter = new Input<RealParameter>("n_S0",
            "the number of susceptibles at time of origin (defaults to 1000).", Input.Validate.REQUIRED);
    public Input<RealParameter> betaParameter = new Input<RealParameter>("beta",
            "the mass action rate of infection.", Input.Validate.REQUIRED);
    public Input<RealParameter> gammaParameter = new Input<RealParameter>("gamma",
            "the per-infected rate of recovery + sampling rate.", Input.Validate.REQUIRED);
    public Input<RealParameter> originParameter = new Input<RealParameter>("origin",
            "the time before the root that the first infection occurred.", Input.Validate.REQUIRED);

    public Input<Integer> storedStateCount = new Input<Integer>("integrationStepCount",
            "number of integration time steps to use (defaults to 1000).", 1000);

    public Input<Double> finishingThresholdInput = new Input<Double>("finishingThreshold",
            "Integration will finish when infected pop drops below this.", 1.0);
    public Input<Double> maxSimLengthInput = new Input<Double>("maxSimLength",
            "Maximum length of simulation. (Default 1000.)", 1000.0);


    public List<Double> NStraj;
    public List<Double> NItraj;
    public boolean reject = false;
    public List<Double> effectivePopSizeTraj;
    public List<Double> intensityTraj;
    public double tIntensityTrajStart;
    public double dt;

    protected boolean dirty;
    protected ContinuousOutputModel integrationResults;

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

        dirty = true;
        update();

    }

    protected abstract boolean update();

    /**
     * Obtain dt for integration results at discrete locations
     *
     * @param beta
     * @param gamma
     * @param NS0
     * @return
     */
    public double simulateTrajectory(final double beta, final double gamma, final double NS0) {

          // Equations of motion:
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {

            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot) throws MaxCountExceededException, DimensionMismatchException {
                double S = y[0];
                double I = y[1];
                ydot[0] = -beta * S * I;
                ydot[1] = beta * S * I - gamma * I;
            }
        };

        try {
            AdaptiveStepsizeIntegrator integrator = new HighamHall54Integrator(1E-10, 100, 0.5, 0.01);

            integrationResults = new ContinuousOutputModel();
            integrator.addStepHandler(integrationResults);

            integrator.addEventHandler(new EventHandler() {

                @Override
                public void init(double t0, double[] y, double t) { };

                @Override
                public double g(double t, double[] y) {
                    return y[1] - finishingThresholdInput.get();
                }

                @Override
                public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
                    if (!increasing)
                        return Action.STOP;
                    else
                        return Action.CONTINUE;
                }

                @Override
                public void resetState(double d, double[] doubles) { };

            }, 1.0, 0.1, 10);

            integrator.addEventHandler(new EventHandler() {

                @Override
                public void init(double t0, double[] y, double t) { };

                @Override
                public double g(double t, double[] y) {
                    return y[1];
                }

                @Override
                public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
                    return EventHandler.Action.STOP;
                }

                @Override
                public void resetState(double d, double[] doubles) { };

            }, 0.1, 0.1, 10);


            double [] y0 = new double[2];
            y0[0] = NS0;
            y0[1] = 1.0;
            double [] y = new double[2];

            // Integrate SIR model ODEs:
            try {
                integrator.integrate(ode, 0, y0, maxSimLengthInput.get(), y);
            } catch (MaxCountExceededException e) {
                reject = true;
            }

        } catch (NumberIsTooSmallException tse) {
            reject=true;
        }

        // Obtain integration results at discrete locations
        dt = integrationResults.getFinalTime()/storedStateCount.get();

        return dt;
    }

    /**
     * "Effective" population size for calculating coalescent likelihood
     * under Volz's model.  This is only used if numerical integration
     * of the intensity is employed.
     *
     * @param t
     * @return NI(t)/(2*beta*NS(t))
     */

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

    /**
     * Uses pre-calculated intensity values to estimate the intensity
     * at time t. For times within the domain of the calculated values,
     * linear interpolation is used.  For times outside this range, the
     * initial or final effectivePopSizeTraj value is used to extrapolate.
     *
     * @param t
     * @return
     */

    public double getIntensity(double t) {
        update();

        if (t < tIntensityTrajStart) {
            return -(tIntensityTrajStart - t) / effectivePopSizeTraj.get(0);
        } else {
            if (t > originParameter.get().getValue() + 0.5 * dt) {
                return
                        intensityTraj.get(intensityTraj.size() - 1)
                                + (t - originParameter.get().getValue() + 0.5 * dt) * 1e10;
            } else {
                // Return linear interpolation between values at surrounding
                // lattice points
                int idx = (int) Math.floor((t - tIntensityTrajStart) / dt);
                double alpha = (t - (tIntensityTrajStart + dt * idx)) / dt;
                return intensityTraj.get(idx) * (1.0 - alpha)
                        + intensityTraj.get(idx + 1) * alpha;
            }
        }
    }

    public double getInverseIntensity(final double intensity) {
        update();

        // given a value for the intensity, find a time
        UnivariateRealFunction intensityFunction = new UnivariateRealFunction() {

            @Override
            public double value(double x) throws FunctionEvaluationException {
                return getIntensity(x) - intensity;
            }
        };

        BrentSolver solver = new BrentSolver();
        try {
            return solver.solve(intensityFunction, 0, n_S_Parameter.get().getValue() / betaParameter.get().getValue());
        } catch (MaxIterationsExceededException e) {
            throw new RuntimeException("Max iterations (" + e.getMaxIterations() + ") exceeded:" + e.getMessage());
        } catch (FunctionEvaluationException e) {
            throw new RuntimeException(e.getMessage());
        }
    }

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
        //out.print("dt\t");

        out.print("R_0\t");

        //for (int i = 0; i < statesToLogInput.get(); i++) {
        //    out.format("S%d\t", i);
        //    out.format("I%d\t", i);
        //    out.format("t%d\t", i);
        //}
    }

    @Override
    public void log(int nSample, PrintStream out) {

        //out.format("%g\t", dt);

        out.format("%g\t", betaParameter.get().getValue()*n_S_Parameter.get().getValue()
                /gammaParameter.get().getValue());

        //double tend = NStraj.size() * dt;
        //double delta = tend / (statesToLogInput.get() - 1);

        //for (int i = 0; i < statesToLogInput.get(); i++) {
        //    double t = delta * i;
        //    int tidx = (int) Math.round(t / dt);
        //
        //    if (tidx >= NStraj.size())
        //        tidx = NStraj.size() - 1;
        //    out.format("%g\t", NStraj.get(tidx));
        //    out.format("%g\t", NItraj.get(tidx));
        //    out.format("%g\t", t);
        // }
    }

    @Override
    public void close(PrintStream out) {
    }

    /**
     * DEBUG: Dump intensities estimated at a given
     * number of times ranging between the extreme values of the provided
     * interval list.
     *
     * @param tree tree used to determine range of times to use
     * @param n number of times to use
     * @param ps PrintStream to send results to
     */

}

package beast.phylodynamics.epidemiology;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.evolution.tree.Tree;
import beast.core.Loggable;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;


/**
 * @author Timothy Vaughan
 * @author Alexei Drummond
 */
@Description("Keeps track of trajectories for effective, infected, and susceptible populations," +
             "used by Volz2009 Tree Distribution and Volz2012PopulationFunction Population Function")

public class VolzSIR extends CalculationNode implements Loggable {

    public Input<RealParameter> n_S_Parameter = new Input<RealParameter>("n_S0",
            "the number of susceptibles at time of origin (defaults to 1000).", Input.Validate.REQUIRED);
    public Input<RealParameter> betaParameter = new Input<RealParameter>("beta",
            "the mass action rate of infection.", Input.Validate.REQUIRED);
    public Input<RealParameter> gammaParameter = new Input<RealParameter>("gamma",
            "the per-infected rate of recovery.", Input.Validate.REQUIRED);
    public Input<RealParameter> originParameter = new Input<RealParameter>("origin",
            "the time before the root that the first infection occurred.", Input.Validate.REQUIRED);

    public Input<Integer> storedStateCount = new Input<Integer>("integrationStepCount",
            "number of integration time steps to use (defaults to 1000).", 1000);
    public Input<Double> finishingThresholdInput = new Input<Double>("finishingThreshold",
            "Integration will finish when infected pop drops below this.", 1.0);
    public Input<Double> maxSimLengthInput = new Input<Double>("maxSimLength",
            "Maximum length of simulation. (Default 1000.)", 1000.0);

    public Input<Integer> statesToLogInput = new Input<Integer>("statesToLog",
            "Number of states to log. (Default 100.)", 100);



    protected boolean dirty;
    public List<Double> NStraj, NItraj;
    private ContinuousOutputModel integrationResults;
    public boolean reject = false;

    // index checkers
    int NScounter = 0;
    int NIcounter = 0;

    double decay_x = -0.1;
    double decay_y = -0.1;

    //
    // Public stuff
    //

    public List<Double> effectivePopSizeTraj, intensityTraj;
    public double tIntensityTrajStart;
    public double dt;

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



    /**
     * Update deterministic trajectory.
     */
    protected boolean update() {

        boolean reject = false;

        if (!dirty)
            return reject;

        // Short-hand for model parameters:
        final double beta = betaParameter.get().getValue();
        final double gamma = gammaParameter.get().getValue();
        double NS0 = n_S_Parameter.get().getValue();

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
                ydot[0] = -beta*S*I;
                ydot[1] = beta*S*I - gamma*I;
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
        //reject = false;

        // Obtain integration results at discrete locations
        dt = integrationResults.getFinalTime()/storedStateCount.get();
        NStraj.clear();
        NItraj.clear();
        effectivePopSizeTraj.clear();
        for (int i=0; i<=storedStateCount.get(); i++) {
            double t = i*dt;
            integrationResults.setInterpolatedTime(t);
            double [] thisy = integrationResults.getInterpolatedState();
            
            if (thisy[0]<0.0 || thisy[1]<0.0)
                break;
            NStraj.add(thisy[0]);
            NItraj.add(thisy[1]);
            effectivePopSizeTraj.add(thisy[1]/(2.0*beta*thisy[0]));
        }

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
        tIntensityTrajStart = originParameter.get().getValue() - dt * (effectivePopSizeTraj.size()-1) - 0.5 * dt;

        dirty = false;

        } catch (NumberIsTooSmallException tse) {
            reject = true;
        }
        //reject = false;

        return reject;
    }

    /**
     * @param t the time
     * @return the number of susceptibles at the given time
     */
    public double getNS(double t) {

        // Choose which index to use (note these are around the other way to effective population size!)
        int tidx = NStraj.size() - (int) Math.floor((t - tIntensityTrajStart) / dt) - 1;

        // kludge
        if (tidx >= NStraj.size()) {
            NScounter++;
               // System.out.print("NScounter: " + NScounter + " ");
            decay_x = decay_x - 0.5;
            return Math.exp(decay_x) * NStraj.get(NStraj.size() - 1);

        } else if (tidx < 0)  {
            return NStraj.get(0);
        } else {
            return NStraj.get(tidx);
        }
    }

    /**
     * @param t the time
     * @return the number of infecteds at the given time
     */
    public double getNI(double t) {

        // Choose which index to use (note these are around the other way to effective population size!)
        int tidx = NItraj.size() - (int) Math.floor((t - tIntensityTrajStart) / dt) - 1;

        // kludge
        if (tidx >= NItraj.size())  {
            NIcounter++;
              //  System.out.print("NIcounter: " + NIcounter + " ");
            decay_y = decay_y - 0.5;
            return Math.exp(decay_y) * NItraj.get(NItraj.size() - 1);
        } else if (tidx < 0)  {
            return NItraj.get(0);
        } else {
            return NItraj.get(tidx);
        }
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

    public void init(PrintStream out) throws Exception {
        //out.print("dt\t");

        out.print("R_0\t");

        //for (int i = 0; i < statesToLogInput.get(); i++) {
        //    out.format("S%d\t", i);
        //    out.format("I%d\t", i);
        //    out.format("t%d\t", i);
        //}
    }

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
    public void dumpIntensities(Tree tree, int n, PrintStream ps) {

        double T = tree.getRoot().getHeight();

        ps.println("t intensity");
        for (int i=0; i<n; i++) {
            double t = T*i/((double)(n-1));
            ps.format("%g %g\n", t, getIntensity(t));
        }

    }

}

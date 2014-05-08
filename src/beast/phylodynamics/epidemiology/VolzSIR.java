package beast.phylodynamics.epidemiology;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BrentSolver;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


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

    public Input<RealParameter> R0 =
            new Input<RealParameter>("R0", "The basic reproduction number", Input.Validate.XOR, betaParameter);

    public Input<RealParameter> gammaParameter = new Input<RealParameter>("gamma",
            "the per-infected rate of recovery + sampling rate.", Input.Validate.REQUIRED);
    public Input<RealParameter> originParameter = new Input<RealParameter>("origin",
            "the time before the root that the first infection occurred.", Input.Validate.REQUIRED);

    public Input<Integer> integrationStepCount = new Input<Integer>("integrationStepCount",
            "number of integration time steps to use for ODE solver or tau-leaping (defaults to 1000).", 1000);

    public Input<Double> finishingThresholdInput = new Input<Double>("finishingThreshold",
            "Integration will finish when infected pop drops below this.", 1.0);
    public Input<Double> maxSimLengthInput = new Input<Double>("maxSimLength",
            "Maximum length of simulation. (Default 1000.)", 1000.0);

    public Input<Boolean> logTrajectoriesInput = new Input<Boolean>("logTrajectories",
            "Log entire trajectories. Default false.", false);

    public List<Double> NStraj;
    public List<Double> NItraj;
    public boolean reject = false;
    public List<Double> effectivePopSizeTraj;
    public List<Double> intensityTraj;
    public double tIntensityTrajStart;

    // the smallest unit of time for computation of the population size function
    protected double dt;

    protected boolean dirty;
    protected ContinuousOutputModel integrationResults;

    @Override
    public void initAndValidate() throws Exception {
        if (betaParameter.get() != null) {
            betaParameter.get().setBounds(
                    Math.max(0.0, betaParameter.get().getLower()),
                    betaParameter.get().getUpper());
        }

        if (R0.get() != null) {
            R0.get().setBounds(
                    Math.max(0.0, R0.get().getLower()),
                    R0.get().getUpper());
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
     * @return beta (birth rate) value, possibly calculating from R0, gamma and S0, thereby catering for both parameterizations.
     */
    public double beta() {
        if (betaParameter.get() != null) {
            return betaParameter.get().getValue();
        } else {
            return R0.get().getValue() * gammaParameter.get().getValue() / n_S_Parameter.get().getValue();
        }
    }

    /**
     * @return R0 (fundamental reproductive number) value, possibly calculating from beta, gamma and S0, thereby catering for both parameterizations.
     */
    private double R0() {
        if (R0.get() != null) {
            return R0.get().getValue();
        } else {
            return betaParameter.get().getValue() * n_S_Parameter.get().getValue() / gammaParameter.get().getValue();
        }
    }

    protected boolean update() {
        
        if (!dirty)
            return reject;

        return simulateTrajectory(beta(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
    }

    /**
     * Simulate a stochastic or deterministic trajectory.
     *
     * @param beta
     * @param gamma
     * @param NS0
     * @return
     */
    public abstract boolean simulateTrajectory(final double beta,
            final double gamma, final double NS0);


    /**
     * Simulate trajectory using values of parameters found in Inputs.
     *
     * @return true if the simulated epidemic continues for the entire
     * interval, false otherwise.
     */
    public boolean simulateTrajectory() {
        return simulateTrajectory(beta(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
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

    /**
     * Return inverse intensity for use in coalescent simulation.
     *
     * @param intensity
     * @return
     */
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
            return solver.solve(intensityFunction, 0, n_S_Parameter.get().getValue() / beta());
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

        if (logTrajectoriesInput.get())
            out.print("trajI\ttrajTime\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {

        //out.format("%g\t", dt);

        out.format("%g\t", R0());

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

        if (logTrajectoriesInput.get()) {
            out.print("\"");
            for (int idx=0; idx<NItraj.size(); idx++) {
                if (idx>0)
                    out.print(",");
                out.print(NItraj.get(idx));
            }
            out.print("\"\t\"");
            for (int idx=0; idx<NItraj.size(); idx++) {
                double t = dt*idx;
                if (idx>0)
                    out.print(",");
                out.print(t);
            }
            out.print("\"\t");
        }
    }

    @Override
    public void close(PrintStream out) {
    }
    
}

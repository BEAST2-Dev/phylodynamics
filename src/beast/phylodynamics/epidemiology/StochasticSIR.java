package beast.phylodynamics.epidemiology;

import beast.core.*;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * @author Alex Popinga
 * @author Alexei Drummond
 * @author Timothy Vaughan
 */
@Description("Keeps track of trajectories for effective, infected, and susceptible populations," +
        "used by Volz2009TreeDistribution and DeterministicSIRPopulationFunction")

public class StochasticSIR extends VolzSIR {

    private ContinuousOutputModel integrationResults;

    public double totalItime = 0;
    double storedTotalItime = 0;

    public StochasticSIR() {}

    public StochasticSIR(RealParameter nSO, RealParameter beta, RealParameter gamma, RealParameter origin) throws Exception {

        initByName("n_S0", nSO, "beta", beta, "gamma", gamma, "origin", origin);
    }

    public boolean simulateStochasticTrajectory() {
        return  simulateStochasticTrajectory(betaParameter.get().getValue(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
    }

    /**
     *
     * @param beta
     * @param gamma
     * @param NS0
     * @return true if simulated stochastic trajectory should force a reject
     */
    private boolean simulateStochasticTrajectory(final double beta, final double gamma, double NS0) {


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
                public void init(double t0, double[] y, double t) {
                }

                ;

                @Override
                public double g(double t, double[] y) {
                    return y[1];
                }

                @Override
                public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
                    return EventHandler.Action.STOP;
                }

                @Override
                public void resetState(double d, double[] doubles) {
                }

                ;

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

            // Obtain integration results at discrete locations
            dt = integrationResults.getFinalTime()/storedStateCount.get();
            NStraj.clear();
            NItraj.clear();
            effectivePopSizeTraj.clear();

            totalItime = 0.0;

            int Nt = storedStateCount.get();

            int Ntraj = 1;
            double T = originParameter.get().getValue();
            double[] times = {T};

            int Nsamples = 10;

            // USE SEIR_simulator to populate the following variables: NStraj, NItraj, effectivePopSizeTraj, intensityTraj

            SEIRState state0 = new SEIRState(NS0, 0, 1, 0, 0);
            double alpha = 10;

            SALTauleapSEIR simulator = new SALTauleapSEIR(state0, 0, beta, new Double[] {gamma}, false, alpha, false);


            // List to hold integrated trajectories:
            List<List<SEIRState>> trajectoryList = new ArrayList<List<SEIRState>>();

            // Allocate and zero critical steps list:
            List<Integer> criticalTrajectories = new ArrayList<Integer>();
            for (int i = 0; i < Nsamples; i++) {
                criticalTrajectories.add(0);
            }

            // Integrate trajectories:
            for (int i = 0; i < Ntraj; i++) {
                simulator.setState(state0);
                trajectoryList.add(simulator.genTrajectory(T, Nt, Nsamples, criticalTrajectories, times));
            }

            for (int i = 0; i < Ntraj; i++) {
                //System.out.println("t\tS\tI\tR");
                List<SEIRState> traj = trajectoryList.get(i);

                for (int j = 0; j < Nsamples; j++) {

                    SEIRState state = traj.get(j);
                    //System.out.println(state.time + "\t" + state.S + "\t" + state.I + "\t" + state.R);

                    //Obtain integration results at discrete locations
                    double t = j*dt;
                    integrationResults.setInterpolatedTime(t);

                    if (state.I < 0.0 || state.S < 0.0) {
                        //System.out.println("panic! T="+((double)j/T/(double)Nsamples));
                        //break;
                        return true;
                    }
                    NStraj.add(state.S);
                    NItraj.add(state.I);
                    effectivePopSizeTraj.add(state.I / (2.0 * beta * state.S));

                    totalItime += state.I * (T/Nsamples);
                }
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
            return true;
        }

        return false;
    }

    public void store() {
        super.store();
        storedTotalItime = totalItime;
        dirty = true;
    }

    /**
     * Update stochastic trajectory.
     */
    @Override
    protected boolean update() {

        if (!dirty) {
            // if not dirty succeed without work
            return false;
        }

        final double beta = betaParameter.get().getValue();
        final double gamma = gammaParameter.get().getValue();
        double NS0 = n_S_Parameter.get().getValue();

        return simulateStochasticTrajectory(beta, gamma, NS0);
    }

    @Override
    public void log(int nSample, PrintStream out) {

        // logs R0
        super.log(nSample,out);

        // now log the integral of infecteds

        out.format("%g\t", totalItime);

    }

    public void init(PrintStream out) throws Exception {

        // inits R0
        super.init(out);

        out.print("totalItime\t");
    }

    public void restore() {
        super.restore();
        totalItime = storedTotalItime;
        dirty = true;
    }

}

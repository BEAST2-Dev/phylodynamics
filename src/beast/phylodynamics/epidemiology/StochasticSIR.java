package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

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

    public Input<Integer> numSamplesFromTrajectory = new Input<Integer>("numSamplesFromTrajectory",
            "number of samples taken from trajectory to use in piecewise-constant coalescent pop-size function " +
                    "(defaults to 100). integrationStepCount should be an integer multiple of this.", 100);
    
    public double totalItime = 0;
    double storedTotalItime = 0;

    public StochasticSIR() {
    }

    public StochasticSIR(RealParameter nSO, RealParameter beta, RealParameter gamma, RealParameter origin, int integrationStepCount, int numSamplesFromTrajectory) throws Exception {
        initByName("n_S0", nSO, "beta", beta, "gamma", gamma, "origin", origin, "integrationStepCount", integrationStepCount, "numSamplesFromTrajectory", numSamplesFromTrajectory);
    }

    public boolean simulateStochasticTrajectory() {
        return simulateStochasticTrajectory(betaParameter.get().getValue(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
    }

    /**
     * Simulate a stochastic trajectory using SALTauleapSEIR
     *
     * @param beta
     * @param gamma
     * @param NS0
     * @return true if simulated stochastic trajectory should force a reject
     */
    private boolean simulateStochasticTrajectory(final double beta, final double gamma, double NS0) {

        NStraj.clear();
        NItraj.clear();
        effectivePopSizeTraj.clear();
        intensityTraj.clear();

        totalItime = 0.0;

        int Nt = integrationStepCount.get();
        int Nsamples = numSamplesFromTrajectory.get();

        int Ntraj = 1;
        double T = originParameter.get().getValue();
        double[] times = {T};

        // USE SEIR_simulator to populate the following variables: NStraj, NItraj, effectivePopSizeTraj, intensityTraj

        SEIRState state0 = new SEIRState(
                NS0,       // susceptibles
                0,         // exposed
                1,         // infected
                0,         // recovered
                0.0        // time
        );
        double alpha = 100;

        SALTauleapSEIR simulator = new SALTauleapSEIR(
                state0,                 // initial state for simulation
                0,                      // exposure rate
                beta,                   // infection rate
                new Double[]{gamma},    // recovery rate
                false,                  // use exposed compartment?
                alpha,                  // threshold for determining critical reactions
                false                   // is deterministic?
        );


        // List to hold integrated trajectories:
        List<List<SEIRState>> trajectoryList = new ArrayList<List<SEIRState>>();

        // Integrate trajectories:
        for (int i = 0; i < Ntraj; i++) {
            simulator.setState(state0);
            trajectoryList.add(simulator.genTrajectory(
                    T,        // total time for simulation of trajectory
                    Nt,       // the number of points used for the tau-leaping lattice (i.e. dt for tau-leaping is T / Nt)
                    Nsamples, // the number of points stored in the resulting trajectory object (i.e. dt for returned trajectory is T / Nsamples)
                    times
            ));
        }

        dt = T / (Nsamples - 1);

        for (int i = 0; i < Ntraj; i++) {
            //System.out.println("t\tS\tI\tR");
            List<SEIRState> traj = trajectoryList.get(i);

            for (int j = 0; j < Nsamples; j++) {

                SEIRState state = traj.get(j);

                if (state.I < 1.0 || state.S < 0.0) {
                    // trajectory failed to reach big T
                    return true;
                }
                NStraj.add(state.S);
                NItraj.add(state.I);
                effectivePopSizeTraj.add((state.I - 1) / (2.0 * beta * state.S));

            }
        }

        // Switch effective pop size to reverse time:
        Collections.reverse(effectivePopSizeTraj);

        // Estimate intensity on integration lattice:

        double intensity = 0.0;
        intensityTraj.add(intensity);

        for (Double popSize : effectivePopSizeTraj) {
            intensity += dt / popSize;
            intensityTraj.add(intensity);
        }

        // Start of integral is 0.5*dt from end of forward-time integration.
        tIntensityTrajStart = originParameter.get().getValue() - dt * (effectivePopSizeTraj.size() - 1) - 0.5 * dt;

        dirty = false;

        return false;
    }

    @Override
    public void store() {
        super.store();
        storedTotalItime = totalItime;
        dirty = true;
    }

    /**
     * Update stochastic trajectory.
     * @return 
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
        super.log(nSample, out);

        // now log the integral of infecteds

        out.format("%g\t", totalItime);

    }

    @Override
    public void init(PrintStream out) throws Exception {

        // inits R0
        super.init(out);

        out.print("totalItime\t");
    }

    @Override
    public void restore() {
        super.restore();
        totalItime = storedTotalItime;
        dirty = true;
    }

}

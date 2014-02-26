package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.parameter.RealParameter;

import java.util.Collections;


/**
 * @author Timothy Vaughan
 * @author Alexei Drummond
 * @author Alex Popinga
 */
@Description("Keeps track of trajectories for effective, infected, and susceptible populations," +
        "used by Volz2009TreeDistribution and Volz2012PopulationFunction")

public class DeterministicSIR extends VolzSIR {

    // index checkers
    int NScounter = 0;
    int NIcounter = 0;

    double decay_x = -0.1;
    double decay_y = -0.1;

    public DeterministicSIR() {
    }

    public DeterministicSIR(RealParameter nSO, RealParameter beta, RealParameter gamma, RealParameter origin) throws Exception {
        initByName("n_S0", nSO, "beta", beta, "gamma", gamma, "origin", origin);
    }

    public boolean simulateDeterministicTrajectory() {
        return simulateDeterministicTrajectory(betaParameter.get().getValue(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
    }

    @Override
    public double simulateTrajectory(double beta, double gamma, double NS0) {
        return super.simulateTrajectory(beta, gamma, NS0);
    }

    /**
     * Simulate a deterministic trajectory using coalescent rate described by Volz (2012)
     *
     * @param beta
     * @param gamma
     * @param NS0
     * @return true if simulated stochastic trajectory should force a reject
     */
    private boolean simulateDeterministicTrajectory(final double beta, final double gamma, double NS0) {

        dt = simulateTrajectory(beta, gamma, NS0);

        NStraj.clear();
        NItraj.clear();
        effectivePopSizeTraj.clear();

        for (int i = 0; i <= integrationStepCount.get(); i++) {
            double t = i * dt;
            integrationResults.setInterpolatedTime(t);
            double[] thisy = integrationResults.getInterpolatedState();

            if (thisy[0] < 0.0 || thisy[1] < 0.0)
                break;
            NStraj.add(thisy[0]);
            NItraj.add(thisy[1]);
            effectivePopSizeTraj.add(thisy[1] / (2.0 * beta * thisy[0]));

            //if (thisy[1] < 1) { System.out.print(" Number of infecteds: " + thisy[1] + " " + "Time: " + integrationResults.getFinalTime());
            //}
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
        tIntensityTrajStart = originParameter.get().getValue() - dt * (effectivePopSizeTraj.size() - 1) - 0.5 * dt;

        dirty = false;

        return false;
    }

    public void store() {
        super.store();
        dirty = true;
    }

    /**
     * Update deterministic trajectory.
     */
    @Override
    protected boolean update() {

        if (!dirty)
            return false;

        return simulateDeterministicTrajectory(betaParameter.get().getValue(),
                gammaParameter.get().getValue(), n_S_Parameter.get().getValue());
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

        } else if (tidx < 0) {
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
        if (tidx >= NItraj.size()) {
            NIcounter++;
            //  System.out.print("NIcounter: " + NIcounter + " ");
            decay_y = decay_y - 0.5;
            return Math.exp(decay_y) * NItraj.get(NItraj.size() - 1);
        } else if (tidx < 0) {
            return NItraj.get(0);
        } else {
            return NItraj.get(tidx);
        }
    }

    public void restore() {
        super.restore();
        dirty = true;
    }

}

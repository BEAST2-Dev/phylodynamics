package beast.phylodynamics.epidemiology;

import beast.core.*;

import java.util.Collections;

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
 * @author Alex Popinga
 */
@Description("Keeps track of trajectories for effective, infected, and susceptible populations," +
             "used by Volz2009TreeDistribution and Volz2012PopulationFunction")

public class DeterministicSIR extends VolzSIR {

    private ContinuousOutputModel integrationResults;

    // index checkers
    int NScounter = 0;
    int NIcounter = 0;

    double decay_x = -0.1;
    double decay_y = -0.1;


    private boolean simulateDeterministicTrajectory(final double beta, final double gamma, double NS0) {

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
            tIntensityTrajStart = originParameter.get().getValue() - dt * (effectivePopSizeTraj.size()-1) - 0.5 * dt;

            dirty = false;

        } catch (NumberIsTooSmallException tse) {
            reject = true;
        }

        return reject;
    }

    /**
     * Update deterministic trajectory.
     */
    @Override
    protected boolean update() {

        boolean reject = false;

        if (!dirty)
            return reject;

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

}

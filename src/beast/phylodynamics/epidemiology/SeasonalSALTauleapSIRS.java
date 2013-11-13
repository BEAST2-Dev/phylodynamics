package beast.phylodynamics.epidemiology;

import beast.phylodynamics.util.Stuff;
import beast.util.Randomizer;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.ArrayList;
import java.util.List;

/**
 * User: Denise
 * Date: 12.11.13
 * Time: 11:33
 */
public class SeasonalSALTauleapSIRS implements SEIR_simulator {


    SEIRState state;
    double exposeRate, loseImmunityRate;
    Double[] infectRate, recoverRate;
    double alpha;
    boolean useExposed;
    boolean isDeterministic;
    int season;
    double[] seasons;

    /**
     * Constructor
     *
     * @param state0	initial state of system
     * @param infect	infection rate
     * @param recover	recovery rate
     * @param alpha	threshold for determining critical reactions
     */

    public SeasonalSALTauleapSIRS(SEIRState state0,
                                  double expose, Double[] infect, Double[] recover, Double loseImmunity,
                                  boolean useExposed,
                                  double alpha, Boolean isDeterministic,
                                  double[] seasons) {
        super();
        this.state = state0.copy();
        this.exposeRate = expose;
        this.infectRate = infect;
        this.recoverRate = recover;
        this.useExposed = useExposed;
        this.isDeterministic = isDeterministic;
        this.alpha = isDeterministic ? 0. : alpha;
        this.loseImmunityRate = loseImmunity;
        this.seasons = seasons;
    }


    public void setRates(double expose, double infect, Double[] recover, double loseImmunity, double alpha) {


    }

    public void setRates(double expose, Double[] infect, Double[] recover, double loseImmunity, double alpha) {

        this.exposeRate = expose;
        this.infectRate = infect;
        this.recoverRate = recover;
        this.loseImmunityRate = loseImmunity;
        this.alpha = isDeterministic ? 0. : alpha;

    }

    /**
     * Set system state
     *
     * @param newState
     */
    public void setState(SEIRState newState) {
        state = newState.copy();
    }

    /**
     * Perform one time-step of fixed length (Does not use exposed compartment.)
     *
     * @param	dt	time-step size
     */
    public boolean step(double dt, int index) {

        double t = 0.0;

        boolean critical_step = false;

        while (true) {

            // Calculate propensities:
            double a_infect = infectRate[season] * state.S * state.I;
            double a_recover = recoverRate[index] * state.I;
            double a_loseImmunity = loseImmunityRate * state.R;

            // Calculate 2nd order corrections (SAL):
            double a2_infect = infectRate[season] * state.S * state.I * (infectRate[season] * (state.S - state.I) - recoverRate[index]);
            double a2_recover = recoverRate[index] * state.I * (infectRate[season] * state.S - recoverRate[index] -loseImmunityRate);
            double a2_loseImmunity = loseImmunityRate * state.R * ( -loseImmunityRate - infectRate[season] * state.I);

             /*
              double a2_infect = 0;
              double a2_recover = 0;
              */
            boolean infectIsCrit = false;
            boolean recoverIsCrit = false;
            boolean loseImmunityIsCrit = false;
            double dtCrit = dt - t;
            int critReactionToFire = 0;

            // Determine which reactions are "critical"
            // and calculate next reaction time:

            double lambda_infect = a_infect * dtCrit + 0.5 * a2_infect * dtCrit * dtCrit;
            if ((alpha > 0) && (state.S < lambda_infect + alpha * Math.sqrt(lambda_infect))) {

                // Infection reaction is critical
                infectIsCrit = true;

                // Determine whether infection will fire:
                double thisdt = Randomizer.nextExponential(a_infect);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 1;
                }
            }

            double lambda_recover = a_recover * dtCrit + 0.5 * a2_recover * dtCrit * dtCrit;
            if ((alpha > 0) && (state.I < lambda_recover + alpha * Math.sqrt(lambda_recover))) {

                // Recovery is critical:
                recoverIsCrit = true;

                // Determine whether recovery will fire:
                double thisdt = Randomizer.nextExponential(a_recover);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 2;
                }
            }


            double lambda_loseImmunity = a_loseImmunity * dtCrit + 0.5 * a2_loseImmunity * dtCrit * dtCrit;
            if (loseImmunityRate>0 && (alpha > 0) && (state.R < lambda_loseImmunity + alpha * Math.sqrt(lambda_loseImmunity))) {

                // losing immunity is critical:
                loseImmunityIsCrit = true;

                // Determine whether recovery will fire:
                double thisdt = Randomizer.nextExponential(a_loseImmunity);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 3;
                }
            }

            // Reaction has been marked as critical, so this is a critical step:
            if (infectIsCrit || recoverIsCrit || loseImmunityIsCrit)
                critical_step = true;

            // Update time:
            t += dtCrit;

            // tau-leap non-critical reactions:
            if (infectIsCrit == false) {
                double q = isDeterministic ? (a_infect * dtCrit)
                        : Randomizer.nextPoisson(a_infect * dtCrit + 0.5 * a2_infect * dtCrit * dtCrit);
                state.S -= q;
                state.I += q;
            }

            if (recoverIsCrit == false) {
                double q = isDeterministic ? (a_recover * dtCrit)
                        : Randomizer.nextPoisson(a_recover * dtCrit + 0.5 * a2_recover * dtCrit * dtCrit);
                state.I -= q;
                state.R += q;
            }

            if (loseImmunityRate>0 && loseImmunityIsCrit == false) {
                double q = isDeterministic ? (a_loseImmunity * dtCrit)
                        : Randomizer.nextPoisson(a_loseImmunity * dtCrit + 0.5 * a2_loseImmunity * dtCrit * dtCrit);
                state.R -= q;
                state.S += q;
            }

            // Zero negative populations:
            // (Beginning to think this is the sensible thing to do.)
            if (state.S < 0)
                state.S = 0;
            if (state.I < 0)
                state.I = 0;
            if (state.R < 0)
                state.R = 0;

            // End step if no critical reaction fires:
            if (critReactionToFire == 0)
                break;

            // Implement one critical reaction:
            switch (critReactionToFire) {
                case 1:
                    // Infection
                    state.S -= 1;
                    state.I += 1;
                    break;
                case 2:
                    // Recovery
                    state.I -= 1;
                    state.R += 1;
                    break;

                case 3:
                    // Losing immunity
                    state.R -= 1;
                    state.S += 1;
                    break;

                default:
                    // No critical reaction
                    break;
            }

        }

        // Check for negative populations:
        if (state.S < 0 || state.I < 0 || state.R < 0) {
            throw new RuntimeException("Error: negative population detected. Rejecting trajectory.");
        }


        // Update state time:
        state.time += dt;

        return critical_step;

    }

    /**
     * Perform one fixed-size time step using exposed compartment
     *
     * @param dt	length of time step
     * @return true or false depending on whether step involved "critical"
     * reactions
     */
    public boolean step_exposed(double dt, int index) {


        if (loseImmunityRate>0) throw new NotImplementedException(); //    SEIRS is not yet implemented!

        double t = 0.0;

        boolean critical_step = false;

        while (true) {

            // Calculate propensities:
            double a_expose = exposeRate * state.S * state.I;
            double a_infect = infectRate[season] * state.E;
            double a_recover = recoverRate[index] * state.I;

            // Calculate 2nd order corrections:
            double a2_expose = -exposeRate * exposeRate * state.S * state.I * state.I
                    + exposeRate * infectRate[season] * state.E * state.S
                    - exposeRate * recoverRate[index] * state.I * state.S;
            double a2_infect = infectRate[season] * exposeRate * state.S * state.I
                    - infectRate[season] * infectRate[season] * state.E;
            double a2_recover = recoverRate[index] * infectRate[season] * state.E
                    - recoverRate[index] * recoverRate[index] * state.I;

            // Determine which reactions are "critical"
            // and calculate next reaction time:
            double dtCrit = dt - t;
            int critReactionToFire = 0;

            boolean exposeIsCrit = false;
            double lambda_expose = a_expose * dtCrit + 0.5 * a2_expose * dtCrit * dtCrit;
            if ((alpha > 0) && (state.S < lambda_expose + alpha * Math.sqrt(lambda_expose))) {

                // Exposure reaction is critical
                exposeIsCrit = true;

                // Determine whether exposure will fire:
                double thisdt = Randomizer.nextExponential(a_expose);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 1;
                }
            }

            boolean infectIsCrit = false;
            double lambda_infect = a_infect * dtCrit + 0.5 * a2_infect * dtCrit * dtCrit;
            if ((alpha > 0) && (state.E < lambda_infect + alpha * Math.sqrt(lambda_infect))) {

                // Infection reaction is critical
                infectIsCrit = true;

                // Determine whether infection will fire:
                double thisdt = Randomizer.nextExponential(a_infect);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 2;
                }
            }

            boolean recoverIsCrit = false;
            double lambda_recover = a_recover * dtCrit + 0.5 * a2_recover * dtCrit * dtCrit;
            if ((alpha > 0) && (state.I < lambda_recover + alpha * Math.sqrt(lambda_recover))) {

                // Recovery is critical:
                recoverIsCrit = true;

                // Determine whether recovery will fire:
                double thisdt = Randomizer.nextExponential(a_recover);
                if (thisdt < dtCrit) {
                    dtCrit = thisdt;
                    critReactionToFire = 3;
                }
            }

            // Reaction has been marked as critical, so this is a critical step:
            if (exposeIsCrit || infectIsCrit || recoverIsCrit)
                critical_step = true;

            // Update time:
            t += dtCrit;

            // tau-leap non-critical reactions:
            if (exposeIsCrit == false) {
                int q = (int)Randomizer.nextPoisson(a_expose * dtCrit + 0.5 * a2_expose * dtCrit * dtCrit);
                state.S -= q;
                state.E += q;
            }

            // tau-leap non-critical reactions:
            if (infectIsCrit == false) {
                int q = (int)Randomizer.nextPoisson(a_infect * dtCrit + 0.5 * a2_infect * dtCrit * dtCrit);
                state.E -= q;
                state.I += q;
            }

            if (recoverIsCrit == false) {
                int q = (int)Randomizer.nextPoisson(a_recover * dtCrit + 0.5 * a2_recover * dtCrit * dtCrit);
                state.I -= q;
                state.R += q;
            }

            // Zero negative populations:
            if (state.S < 0)
                state.S = 0;
            if (state.E < 0)
                state.E = 0;
            if (state.I < 0)
                state.I = 0;
            if (state.R < 0)
                state.R = 0;

            // End step if no critical reaction fires:
            if (critReactionToFire == 0)
                break;

            // implement one critical reaction:
            switch (critReactionToFire) {
                case 1:
                    // Exposure
                    state.S -= 1;
                    state.E += 1;
                    break;

                case 2:
                    // Infection
                    state.E -= 1;
                    state.I += 1;
                    break;

                case 3:
                    // Recovery
                    state.I -= 1;
                    state.R += 1;
                    break;
            }

        }

        // Check for negative populations:
        if (state.S < 0 || state.E < 0 || state.I < 0 || state.R < 0) {
            throw new RuntimeException("Error: negative population detected. Rejecting trajectory.");
        }

        // Update state time:
        state.time += dt;

        return critical_step;

    }


    public List<SEIRState> genTrajectory(double T, int Nt, int Nsamples, int ntaxa, Boolean check, double[] times) {

        // Determine time-step size:
        double dt = T / (Nt - 1);

        // Determine number of time steps per sample:
        int stepsPerSample = (Nt - 1) / (Nsamples - 1);

        // Allocate memory for sampled states:
        List<SEIRState> trajectory = new ArrayList<SEIRState>();

        // Sample first state:
        trajectory.add(state.copy());

        //index of rates (in case they change over time)
        int index = 0;

        for (int tidx = 1; tidx < Nt; tidx++) {

            index = Stuff.index(tidx * dt, times);

            season = Stuff.index(tidx * dt, seasons) % 2 ;

            if (useExposed)
                step_exposed(dt, index);
            else
                step(dt, index);

            // Sample if necessary:
            if (tidx % stepsPerSample == 0) {
                if ((Math.round(state.I) < 1) && (!useExposed || (state.E == 0)))
                    //                if (trajectory.get(trajectory.size()-1).S == state.S ) zeroCount++;
                    //                if (check && zeroCount > Nsamples/2.)
                    throw new RuntimeException("Abort simulation. No infecteds left.");

                trajectory.add(state.copy());

            }
        }

        if ((Math.round(state.I) < 1) && (!useExposed || (state.E == 0)))
            throw new RuntimeException("Abort simulation. No infecteds left.");

        return trajectory;

    }

    /**
     * Main method: for debugging only
     *
     * @param args
     */
    public static void main(String[] args) {


        // Simulation parameters:
        int Ntraj = 1;		// Number of trajectories
        int Nt = 10001;			// Number of timesteps
        int Nsamples = 11;		// Number of samples to record
        double T = 10.;		// Length of time of simulation

        //double alpha = 10;		// Critical reaction parameter
        double alpha = 10.;		// No SSA component

        // Model parameters:
        int s0 = 3000;
        int e0 = 0;
        int i0 = 1;
        int r0 = 0;

        double[] times = {T};
        double expose = 0.0;
        Double[] infect = {9. / s0, 2. / s0};
        Double[] recover = {5.};
        double loseImmunity = 2.;

        double[] seasons = {2,4,6,8,10};

        SEIRState x0 = new SEIRState(s0, e0, i0, r0, 0.0);

        // List to hold integrated trajectories:
        List<List<SEIRState>> trajectoryList = new ArrayList<List<SEIRState>>();

        // Create SALTauleapSEIR instance:
        SeasonalSALTauleapSIRS hybridTauleapSEIR = new SeasonalSALTauleapSIRS(x0, expose, infect, recover, loseImmunity, false, alpha, false, seasons);

        // Allocate and zero critical steps list:
        List<Integer> criticalTrajectories = new ArrayList<Integer>();
        for (int i = 0; i < Nsamples; i++)
            criticalTrajectories.add(0);

        // Integrate trajectories:
        for (int i = 0; i < Ntraj; i++) {

            try{

                hybridTauleapSEIR.setState(x0);
                trajectoryList.add(hybridTauleapSEIR.genTrajectory(T, Nt, Nsamples, 10, true, times));

            }catch(Exception e) { i--;}

        }

        for (int i = 0; i < Ntraj; i++) {
            System.out.println("t\tS\tI\tR");
            List<SEIRState> traj = trajectoryList.get(i);

            for (int j = 0; j < Nsamples; j++) {

                SEIRState state = traj.get(j);
                System.out.println(state.time + "\t" + state.S + "\t" + state.I + "\t" + state.R);


            }

        }

        System.exit(0);
    }
}

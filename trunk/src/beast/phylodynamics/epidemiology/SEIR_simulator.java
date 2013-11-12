package beast.phylodynamics.epidemiology;

import java.util.List;

/**
 * @author dkuh004
 *         Date: Jan 27, 2012
 *         Time: 10:50:03 AM
 */
public interface SEIR_simulator {


    public void setRates(double expose, double infect, Double[] recover, double loseImmunity, double alpha);

    public void setState(SEIRState newState);

    public boolean step(double dt, int index);

    public boolean step_exposed(double dt, int index);

//    public List<SEIRState> genTrajectory(double T, int Nt, int Nsamples, List<Integer> criticalTrajectories, double[] times);

    public List<SEIRState> genTrajectory(double T, int Nt, int Nsamples, int ntaxa, Boolean check, double[] times) ;

}

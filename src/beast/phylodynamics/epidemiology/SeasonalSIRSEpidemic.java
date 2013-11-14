package beast.phylodynamics.epidemiology;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.Arrays;

/**
 * User: Denise
 * Date: 12.11.13
 * Time: 11:22
 */
public class SeasonalSIRSEpidemic extends HybridSEIREpidemic {

    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<RealParameter>("birthRateChangeTimes", "The times t_i specifying when birth rate changes occur", (RealParameter) null);

    public double[] birthTimes;

    @Override
    public void initAndValidate() throws Exception {

        if (birthRateScalar.get().getDimension()!=2) throw new Exception("birthrate should have dimesnion 2 (winter rate & summer rate)");

        birthTimes = new double[birthRateChangeTimesInput.get().getDimension()];

        super.initAndValidate();

    }



    @Override
    public void initializeTrajectory(int Nsamples, double alpha){

        SEIRState x0 = new SEIRState((int) S0.get().getArrayValue() - e0 - i0 - r0, e0, i0, r0, 0.0);

        if (simulationType.get().equals("hybrid"))  throw new NotImplementedException();  // Seasonal SIRS is only implemented with SAL simulator

        hybridTauleapSEIR = new SeasonalSALTauleapSIRS(x0, expose, infect, recover, loseImmunity, false, alpha, isDeterministic.get(), birthTimes);

        hybridTauleapSEIR.setState(x0);
        //        List<SEIRState> trajectory = null;

        dS = new Double[Nsamples-1];
        Arrays.fill(dS, 0.);
        dE = new Double[Nsamples-1];
        dR = new Double[Nsamples-1];


    }

    @Override
    public Boolean initTraj(Double[] birth,Double expose, Double[] death, Double[] psi, Double loseImmunity, double T, int ntaxa, int maxLoop, int intervals, double[] times, Boolean outputInfPop){

        return refresh((int) S0.get().getArrayValue(), useExposed ? E0.get().getValue(): 0, birth, useExposed?expose:0., death, psi, loseImmunity, T, ntaxa, maxLoop, intervals, times, outputInfPop);

    }


    @Override
    public Boolean refresh(int S0, int E0, Double[] birth,Double exposed, Double[] death, Double[] psi, Double loseImmunity, double T, int ntaxa, int maxLoop, int intervals, double[] times, Boolean outputInfPop){

        tau = T / intervals;

        infect = birth;

        // allow the sampling rate to change over time
        recover = new Double[psi.length];

        for (int i=0; i<recover.length; i++){

            recover[i] = death[i] +  psi[i];
        }

        if (samplingRateChangeTimesInput.get()!= null) updateSamplingChangeTimes();

        updateBirthChangeTimes(T);

        ( (SeasonalSALTauleapSIRS) hybridTauleapSEIR).seasons = birthTimes;

        if (useExposed){
            expose = exposed;
            e0 = E0;
        }

        ((SeasonalSALTauleapSIRS)hybridTauleapSEIR).setRates(expose, infect, recover, loseImmunity, alpha.get());

        return generateTrajectory(S0, intervals, ntaxa, Nt.get(), useExposed, T, maxLoop, times, outputInfPop);

    }

    void updateBirthChangeTimes(double T){

        Double[] birthChangeTimes = birthRateChangeTimesInput.get().getValues();

        if (reverseTimeArrays.get()!=null && reverseTimeArrays.get().getValue(0)){
            //birthTimes are reversed

            birthTimes[birthTimes.length-1] = T;
            for (int i=1; i<birthTimes.length; i++){
                birthTimes[i] = T-birthChangeTimes[birthTimes.length-i-1];
            }
        } else {
            for (int i=0; i<birthTimes.length; i++){
                birthTimes[i] = birthChangeTimes[i];
            }

        }

    }


}
